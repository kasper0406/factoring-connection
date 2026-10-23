"""DIAGNOSTIC / HARNESS code for agent `newfamily` (knows factors).  Never imported by operational modules.

* `dedekind_brute(a, c)`: s(a, c) = sum_{k=1}^{c-1} ((k/c)) ((ak/c)) as a Fraction.
* `dedekind_split(a, p, q)`: s(a, pq) = s(a, p) + s(a, q) + U(a, pq), with the unit part
      U(a, N) = sum_{k in (Z/N)^*} ((k/N)) ((ak/N))
  (inclusion-exclusion over k = 0 mod p, k = 0 mod q).  Returns the three pieces (exact).
* Salie sums: closed form check and the factoring reduction (Prop SAL in notes/newfamily_results.md).
* Cubic residue symbols (alpha/(N))_3 for alpha in Z[omega], computed from the factorisation, to test Theorem RA.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Dict, List, Optional, Tuple

import mpmath

from src.newfamily.recip import jacobi


def _saw(x: Fraction) -> Fraction:
    if x.denominator == 1:
        return Fraction(0)
    return x - (x.numerator // x.denominator) - Fraction(1, 2)


def dedekind_brute(a: int, c: int) -> Fraction:
    return sum((_saw(Fraction(k, c)) * _saw(Fraction(a * k, c)) for k in range(1, c)), Fraction(0))


def unit_part_brute(a: int, N: int) -> Fraction:
    return sum((_saw(Fraction(k, N)) * _saw(Fraction(a * k, N)) for k in range(1, N) if gcd(k, N) == 1),
               Fraction(0))


def dedekind_split(a: int, p: int, q: int) -> Tuple[Fraction, Fraction, Fraction]:
    N = p * q
    sp, sq = dedekind_brute(a % p, p), dedekind_brute(a % q, q)
    return sp, sq, dedekind_brute(a, N) - sp - sq


# ----------------------------------------------------------------------------------------------- Salie sums

def salie_exact_terms(a: int, b: int, N: int, p: int, q: int) -> List[int]:
    """Square roots v of ab mod N (via the factors; harness only)."""
    from sympy.ntheory import sqrt_mod
    rp = sqrt_mod(a * b % p, p, all_roots=True) or []
    rq = sqrt_mod(a * b % q, q, all_roots=True) or []
    out = []
    for x in rp:
        for y in rq:
            # CRT
            v = (x * q * pow(q, -1, p) + y * p * pow(p, -1, q)) % N
            out.append(v)
    return sorted(set(out))


class SalieOracle:
    """HARNESS oracle for T(c) := S(c^2, 1; N) / (eps_N sqrt(N)) = sum_{v^2 = c^2} e(2v/N)  (a real number:
    2cos(4 pi c/N) + 2cos(4 pi c e/N), e the nontrivial square root of 1 that is 1 mod p, -1 mod q),
    returned with an absolute error <= delta (deterministic pseudo-noise from `seed`).  Computed with mpmath at
    `prec` bits; the noise is added exactly.  Knows p, q.  Counts queries."""

    def __init__(self, p: int, q: int, delta: float, seed: int = 0, prec: int = 0):
        import random
        self.p, self.q, self.N = p, q, p * q
        self.e = (q * pow(q, -1, p) - p * pow(p, -1, q)) % self.N   # 1 mod p, -1 mod q
        self.delta = mpmath.mpf(delta)
        self.rng = random.Random(seed)
        self.prec = prec or (3 * self.N.bit_length() + 64)
        self.queries = 0

    def __call__(self, c: int) -> mpmath.mpf:
        self.queries += 1
        with mpmath.workprec(self.prec):
            N = self.N
            val = 2 * mpmath.cos(4 * mpmath.pi * (c % N) / N) + 2 * mpmath.cos(4 * mpmath.pi * (c * self.e % N) / N)
            noise = (2 * self.rng.random() - 1) * self.delta
            return val + noise


def salie_factor(N: int, oracle, c: int, delta: float) -> Optional[int]:
    """OPERATIONAL decoder (sees N, the oracle, a public multiplier c): Prop SAL.  Queries T(2^k c) for
    k = 0..K, K = bitlen(N) + 2.  From T subtract the known term 2cos(4 pi 2^k c/N) to get d_k ~ ||x_k|| where
    x_k = 2^{k+1} c e / N mod 1, then reconstruct x_0 backwards (x_k = x_{k+1}/2 or (x_{k+1}+1)/2, chosen by
    ||x_k|| = d_k), round to a multiple of 1/N, and take gcds.  Returns a proper factor or None."""
    K = N.bit_length() + 2
    prec = 3 * N.bit_length() + 64
    with mpmath.workprec(prec):
        d = []
        for k in range(K + 1):
            ck = (c << k) % N
            t = oracle(ck) - 2 * mpmath.cos(4 * mpmath.pi * ck / N)
            t = max(min(t / 2, mpmath.mpf(1)), mpmath.mpf(-1))       # cos(2 pi x_k)
            d.append(mpmath.acos(t) / (2 * mpmath.pi))               # ||x_k|| in [0, 1/2]
        for x in (d[K], 1 - d[K]):
            for k in range(K - 1, -1, -1):
                c0, c1 = x / 2, (x + 1) / 2
                x = c0 if abs(min(c0, 1 - c0) - d[k]) <= abs(min(c1, 1 - c1) - d[k]) else c1
            m = int(mpmath.nint(x * N)) % N      # m = 2 c e mod N (or its negative)
            for mm in (m, -m % N):
                ce = mm * pow(2, -1, N) % N      # c e mod N
                for g in (gcd(ce - c, N), gcd(ce + c, N)):
                    if 1 < g < N:
                        return g
    return None


# ------------------------------------------------------------------------------- cubic residue symbols (RA)

def _cubic_char_Fp(x: int, p: int, r: int) -> int:
    """For p = 1 mod 3, r a root of r^2+r+1 = 0 mod p: x^{(p-1)/3} mod p in {1, r, r^2} -> exponent in {0,1,2}."""
    t = pow(x % p, (p - 1) // 3, p)
    if t == 1:
        return 0
    if t == r % p:
        return 1
    if t == r * r % p:
        return 2
    raise ValueError("not a cube root of unity")


def cubic_symbol_prime(u: int, v: int, p: int) -> int:
    """(alpha/(p))_3 for alpha = u + v*omega, p != 3 prime coprime to N(alpha), as exponent k (value omega^k).
    p = 1 mod 3: (p) = pi pibar; (alpha/pi)(alpha/pibar) with the embeddings omega -> r, r^2 (the choice of r
    only swaps the two factors).  p = 2 mod 3: (p) inert, residue field F_p[omega]: alpha^{(p^2-1)/3}."""
    if p % 3 == 1:
        from sympy.ntheory import sqrt_mod
        s = sqrt_mod(-3 % p, p)
        r = (-1 + s) * pow(2, -1, p) % p
        k1 = _cubic_char_Fp(u + v * r, p, r)
        j2 = _cubic_char_Fp(u + v * r * r, p, r)
        # mod pibar = (p, omega - r^2) the symbol omega^{k2} reduces to r^{2 k2}, so k2 = 2 j2 mod 3
        return (k1 + 2 * j2) % 3
    # inert: arithmetic in F_p[w]/(w^2 + w + 1)
    def mul(x, y):
        a, b = x
        c, d = y
        # (a + b w)(c + d w) = ac + (ad + bc) w + bd w^2, w^2 = -1 - w
        return ((a * c - b * d) % p, (a * d + b * c - b * d) % p)
    e = (p * p - 1) // 3
    res, base = (1, 0), (u % p, v % p)
    while e:
        if e & 1:
            res = mul(res, base)
        base = mul(base, base)
        e >>= 1
    return {(1, 0): 0, (0, 1): 1, ((-1) % p, (-1) % p): 2}[res]


def cubic_symbol_N(u: int, v: int, p: int, q: int) -> int:
    return (cubic_symbol_prime(u, v, p) + cubic_symbol_prime(u, v, q)) % 3
