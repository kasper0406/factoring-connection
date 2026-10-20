"""Algebraic and Aurifeuillian pieces of r^e - 1 and r^e + 1 (notes/snfsrig_results.md section 3).

Operational-safe: given (r, e, sign) the routine computes Phi_d(r) for the relevant d (integer polynomial
evaluation), refines the pieces by gcds, and for r = 2 applies the Aurifeuillian identity
    2^{4k+2} + 1 = (2^{2k+1} - 2^{k+1} + 1)(2^{2k+1} + 2^{k+1} + 1).
It never factors N; the divisors of the small exponent e (e <= log2 N) are found by trial division.
The pieces are then handed to a general-purpose factoring algorithm (e.g. Theorem U); the complexity consequence is
computed exactly by `cost_exponent_ratio`.
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr
from typing import Dict, List, Tuple


def divisors_small(n: int) -> List[int]:
    """Divisors of a small positive integer n by trial division (n is an exponent, n <= log2 N)."""
    if n < 1:
        raise ValueError
    small, large = [], []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            if i * i != n:
                large.append(n // i)
        i += 1
    return small + large[::-1]


def euler_phi_small(n: int) -> int:
    res, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            res -= res // p
        p += 1
    if m > 1:
        res -= res // m
    return res


def mobius_small(n: int) -> int:
    m, p, k = n, 2, 0
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            k += 1
        p += 1
    if m > 1:
        k += 1
    return -1 if k % 2 else 1


def cyclotomic_value(n: int, r: int) -> int:
    """Phi_n(r) = prod_{d | n} (r^d - 1)^{mu(n/d)}, exact integer arithmetic (r >= 2)."""
    num, den = 1, 1
    for d in divisors_small(n):
        mu = mobius_small(n // d)
        if mu == 1:
            num *= r ** d - 1
        elif mu == -1:
            den *= r ** d - 1
    q, rem = divmod(num, den)
    assert rem == 0
    return q


def algebraic_pieces(r: int, e: int, sign: int) -> Dict[int, int]:
    """sign = -1: r^e - 1 = prod_{d | e} Phi_d(r).  sign = +1: r^e + 1 = prod_{d | 2e, d not | e} Phi_d(r).
    Returns {d: Phi_d(r)}; the product equals r^e + sign."""
    if sign == -1:
        ds = divisors_small(e)
    elif sign == 1:
        ds = [d for d in divisors_small(2 * e) if e % d != 0]
    else:
        raise ValueError
    pieces = {d: cyclotomic_value(d, r) for d in ds}
    prod = 1
    for v in pieces.values():
        prod *= v
    assert prod == r ** e + sign
    return pieces


def aurifeuillian_2(n: int) -> Tuple[int, int] | None:
    """For n = 4k + 2: 2^n + 1 = L * M with L, M = 2^{2k+1} -/+ 2^{k+1} + 1.  None otherwise."""
    if n % 4 != 2:
        return None
    k = (n - 2) // 4
    L = 2 ** (2 * k + 1) - 2 ** (k + 1) + 1
    M = 2 ** (2 * k + 1) + 2 ** (k + 1) + 1
    assert L * M == 2 ** n + 1
    return L, M


def refine(pieces: List[int], splitters: List[int]) -> List[int]:
    """Coprime-ish refinement of a list of pieces by gcds with the splitters (keeps the product)."""
    out = []
    for p in pieces:
        parts = [p]
        for s in splitters:
            new = []
            for q in parts:
                g = math.gcd(q, s)
                if 1 < g < q:
                    new.extend([g, q // g])
                else:
                    new.append(q)
            parts = new
        out.extend(parts)
    return out


def special_pieces(r: int, e: int, sign: int) -> List[int]:
    """All pieces >1 of r^e + sign from cyclotomic factors, refined by the Aurifeuillian pair when r = 2,
    sign = +1 and e = 2 mod 4 (and by the pairs of 2^{e'} + 1 for every e' | e with e/e' odd, e' = 2 mod 4)."""
    pieces = [v for v in algebraic_pieces(r, e, sign).values() if v > 1]
    if r == 2 and sign == 1:
        spl = []
        for ep in divisors_small(e):
            if (e // ep) % 2 == 1 and ep % 4 == 2:
                L, M = aurifeuillian_2(ep)
                spl += [L, M]
        pieces = refine(pieces, spl)
    prod = 1
    for v in pieces:
        prod *= v
    assert prod == r ** e + sign
    return [v for v in pieces if v > 1]


def largest_piece_log_ratio(r: int, e: int, sign: int) -> Fr:
    """Exact-rational upper bracket of log(max piece)/log(N) via bit lengths: (bitlen(max)) / (bitlen(N) - 1)."""
    pcs = special_pieces(r, e, sign)
    N = r ** e + sign
    return Fr(max(p.bit_length() for p in pcs), N.bit_length() - 1)


def cost_exponent_ratio(eps: Fr, alpha: Fr = Fr(2, 5)) -> Tuple[Fr, Fr]:
    """L_{N^eps}[alpha, c] = L_N[alpha, c eps^alpha (1 + o(1))]: returns a rational bracket of eps^alpha for
    alpha = 2/5 (eps^{2/5} = (eps^2)^{1/5}, bracketed by bisection on x^5 vs eps^2)."""
    if alpha != Fr(2, 5):
        raise NotImplementedError
    target = Fr(eps) ** 2
    lo, hi = Fr(0), Fr(1)
    for _ in range(60):
        mid = (lo + hi) / 2
        if mid ** 5 < target:
            lo = mid
        else:
            hi = mid
    return lo, hi


def primorial(k: int) -> int:
    out, n, c = 1, 2, 0
    while c < k:
        if all(n % p for p in range(2, int(n ** 0.5) + 1)):
            out *= n
            c += 1
        n += 1
    return out
