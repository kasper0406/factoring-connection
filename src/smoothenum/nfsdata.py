"""NFS data for the smoothenum prototypes: polynomial pairs, factor bases, exact norms, exact relation tests.

Exact integer arithmetic throughout (gmpy2 / python ints). The polynomial pair either comes from a CADO-NFS
.poly file (polyselect + ropt output; experiments/raw/smoothenum/c61.poly) or from a seeded base-m expansion.
Nothing here factors N. Relation tests factor *norm cofactors* (<= 2^mfb bits), which is the standard NFS
step; cofactor primality uses gmpy2.is_prime (a probabilistic primality test, noted per AGENT_RULES rule 2).

Conventions (CADO): g(x) = Y1 x + Y0, rational norm G(a, b) = Y1 a + Y0 b; algebraic F(a, b) = sum c_i a^i b^(d-i).
A pair (a, b) has b > 0 and gcd(a, b) = 1.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import gcd, isqrt
from typing import Dict, List, Sequence, Tuple

import gmpy2
import numpy as np
from flint import nmod_poly


@dataclass
class PolyPair:
    N: int
    f: List[int]          # c_0 .. c_d
    Y0: int
    Y1: int
    skew: float = 1.0
    name: str = ""

    @property
    def d(self) -> int:
        return len(self.f) - 1

    def G(self, a: int, b: int) -> int:
        return self.Y1 * a + self.Y0 * b

    def F(self, a: int, b: int) -> int:
        return sum(c * a ** i * b ** (self.d - i) for i, c in enumerate(self.f))

    def check_root(self) -> bool:
        """F(-Y0, Y1) == 0 mod N, i.e. m = -Y0/Y1 is a common root mod N."""
        return self.F(-self.Y0, self.Y1) % self.N == 0


def load_cado_poly(path: str) -> PolyPair:
    d: Dict[str, str] = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        k, v = line.split(":", 1)
        d[k.strip()] = v.strip()
    deg = max(int(k[1:]) for k in d if k.startswith("c") and k[1:].isdigit())
    P = PolyPair(N=int(d["n"]), f=[int(d.get(f"c{i}", "0")) for i in range(deg + 1)], Y0=int(d["Y0"]),
                 Y1=int(d["Y1"]), skew=float(d.get("skew", "1")), name=path.split("/")[-1])
    assert P.check_root()
    return P


def base_m_pair(N: int, d: int) -> PolyPair:
    """Base-m pair: m = floor(N^(1/(d+1))) + 1, f = base-m digits of N (non-negative), g = x - m."""
    m = int(gmpy2.iroot(gmpy2.mpz(N), d + 1)[0]) + 1
    f, n = [], N
    for _ in range(d):
        f.append(n % m)
        n //= m
    f.append(n)
    P = PolyPair(N=N, f=f, Y0=-m, Y1=1, skew=1.0, name=f"basem_d{d}")
    assert sum(c * m ** i for i, c in enumerate(P.f)) == N and P.check_root()
    return P


def primes_upto(B: int) -> np.ndarray:
    s = np.ones(B + 1, dtype=bool)
    s[:2] = False
    for p in range(2, isqrt(B) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0].astype(np.int64)


@dataclass
class FactorBase:
    """Per side: arrays (p, r) with a == r*b (mod p) the root condition; projective roots (p | leading coeff
    in b-direction, i.e. b == 0 mod p) are kept separately and ignored by the sieve (b > 0 rows still get the
    exact test)."""
    p: np.ndarray
    r: np.ndarray
    proj: List[int] = field(default_factory=list)


def rational_fb(P: PolyPair, B: int) -> FactorBase:
    ps, rs, proj = [], [], []
    for p in primes_upto(B).tolist():
        if P.Y1 % p == 0:
            proj.append(p)
            continue
        ps.append(p)
        rs.append((-P.Y0 * pow(P.Y1, -1, p)) % p)
    return FactorBase(np.array(ps, np.int64), np.array(rs, np.int64), proj)


def algebraic_fb(P: PolyPair, B: int) -> FactorBase:
    ps, rs, proj = [], [], []
    for p in primes_upto(B).tolist():
        cs = [c % p for c in P.f]
        if cs[-1] == 0:
            proj.append(p)
        poly = nmod_poly(cs, p)
        if poly.degree() <= 0:
            continue
        for root, _mult in poly.roots():
            ps.append(p)
            rs.append(int(root))
    return FactorBase(np.array(ps, np.int64), np.array(rs, np.int64), proj)


def primorial(B: int) -> gmpy2.mpz:
    """Product of all primes <= B (product tree)."""
    xs = [gmpy2.mpz(p) for p in primes_upto(B).tolist()]
    while len(xs) > 1:
        xs = [xs[i] * xs[i + 1] if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
    return xs[0]


class SideTest:
    """Exact las-style side condition: |n| = (B-smooth part) * c with c <= 2^mfb and every prime of c < 2^lpb.
    c is 1, a prime < 2^lpb, or (if mfb allows) a product of two primes < 2^lpb."""

    def __init__(self, B: int, lpb: int, mfb: int):
        self.B, self.lpb, self.mfb = B, lpb, mfb
        self.P = primorial(B)
        self.L = gmpy2.mpz(1) << lpb
        self.Mf = gmpy2.mpz(1) << mfb
        self.calls = 0

    def cofactor(self, n: int) -> gmpy2.mpz:
        t = gmpy2.mpz(abs(n))
        if t == 0:
            return gmpy2.mpz(0)
        g = gmpy2.gcd(self.P % t, t)
        while g != 1:
            t //= g
            g = gmpy2.gcd(t, g)
        return t

    def ok(self, n: int) -> bool:
        self.calls += 1
        c = self.cofactor(n)
        if c == 0:
            return False
        if c == 1:
            return True
        if c > self.Mf:
            return False
        if gmpy2.is_prime(c):
            return c < self.L
        # composite cofactor <= 2^mfb: both primes must be in (B, 2^lpb); c >= B^2 required
        if c < gmpy2.mpz(self.B) ** 2:
            return False
        f = _split_small(int(c))
        return f is not None and all(q < self.L and gmpy2.is_prime(q) for q in f) and len(f) == 2


def _split_small(c: int):
    """Split a small composite cofactor (<= 2^64 in practice) by Pollard rho (Brent). Returns [p, q] or None."""
    if c % 2 == 0:
        return [2, c // 2]
    for s in range(1, 30):
        x = y = 2
        g = 1
        while g == 1:
            x = (x * x + s) % c
            y = (y * y + s) % c
            y = (y * y + s) % c
            g = gcd(abs(x - y), c)
        if g != c:
            return [g, c // g]
    return None
