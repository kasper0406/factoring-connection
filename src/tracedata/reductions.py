"""OPERATIONAL exact reductions from trace data to a factor of N (tracedata, phase 25).

Every function receives N and the stated data only (never p, q, group orders or factorisations), uses exact integer
arithmetic and returns a nontrivial factor of N, or None.  No factorisation or primality routine is called.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Iterable, List, Optional, Sequence


def _proper(g: int, N: int) -> Optional[int]:
    g = abs(g)
    return g if 1 < g < N else None


def _int_roots_quadratic(A: int, B: int, C: int) -> List[int]:
    """Integer roots of A z^2 + B z + C (A, B, C integers, not all zero)."""
    if A == 0:
        if B == 0:
            return []
        return [-C // B] if C % B == 0 else []
    D = B * B - 4 * A * C
    if D < 0:
        return []
    r = isqrt(D)
    if r * r != D:
        return []
    out = []
    for num in (-B + r, -B - r):
        if num % (2 * A) == 0:
            out.append(num // (2 * A))
    return out


def _split_from_S(N: int, aN: int, S: int) -> Optional[int]:
    """S = q a_p^2 + p a_q^2 and N a_N^2 = (q a_p^2)(p a_q^2): the two roots, then gcd with N."""
    for z in _int_roots_quadratic(1, -S, N * aN * aN):
        g = _proper(gcd(z, N), N)
        if g:
            return g
    return _proper(gcd(S, N), N)


def bc2_factor(N: int, aN: int, aN2: int) -> Optional[int]:
    """Bach--Charles, weight 2: a_{N^2} = (a_p^2 - p)(a_q^2 - q) = a_N^2 - S + N."""
    return _split_from_S(N, aN, aN * aN + N - aN2)


def bc3_factor(N: int, aN: int, aN3: int) -> Optional[int]:
    """a_{p^3} = a_p^3 - 2 p a_p, so a_{N^3} = a_N (a_N^2 - 2 S + 4 N).  Needs a_N != 0."""
    if aN == 0 or aN3 % aN:
        return None
    num = aN * aN + 4 * N - aN3 // aN
    if num % 2:
        return None
    return _split_from_S(N, aN, num // 2)


def linked_pair_factor(N: int, aN_f: int, aN_g: int, c: int = 1, e: int = 2) -> Optional[int]:
    """Two forms with a_p(g) = c a_p(f)^2 - e p at both primes:
    a_N(g) = c^2 a_N^2 - c e S + e^2 N.  (Sym^2: c = 1, e = 1;  CM weight 3 from weight 2: c = 1, e = 2.)"""
    num = c * c * aN_f * aN_f + e * e * N - aN_g
    if num % (c * e):
        return None
    return _split_from_S(N, aN_f, num // (c * e))


def mix_factor(N: int, traces: Sequence[int], M: int) -> Optional[int]:
    """Unordered {a_p, a_q} plus M = (p+1-a_p)(q+1-a_q): try both orders,
    v p^2 + (N + u v - M) p + u N = 0 with u = 1 - a_p, v = 1 - a_q."""
    a, b = traces
    for (x, y) in ((a, b), (b, a)):
        u, v = 1 - x, 1 - y
        for p in _int_roots_quadratic(v, N + u * v - M, u * N):
            g = _proper(p, N)
            if g and N % g == 0:
                return g
    return None


def s_mod_m_from_torsion(N: int, aN: int, m: int) -> int:
    """If E(Q) has a point of order m, a_r = r + 1 (mod m) at both primes, so a_N = N + 1 + s (mod m)."""
    return (aN - 1 - N) % m


def cm_split_bit(aN: int) -> bool:
    """CM curve, (-D/N) = 1: a_N != 0 iff both primes split (a_r != 0 at split r > D)."""
    return aN != 0
