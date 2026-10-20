"""Aligned-cluster sieve (ACS) variants -- OPERATIONAL generators.

Each generator receives only N and public parameters.  It yields candidate relations
    u = prod_l (m - a_l)  (pieces: small integers or factor-base symbols),
    v = u mod N (centered),
built from a split polynomial P(X) = prod (X - a_l) whose top coefficients are aligned with N.
No factorisation of N, no primality-with-factor routine.  (Smoothness tests on v are done by the
harness in run_experiment.py; v is a residual, not N.)

Variants (notes/ptecons_prereg.md):
  cos      M = 2   (Coppersmith-Odlyzko-Schroeppel linear sieve): z = isqrt(N)+1, (z+r1)(z+r2).
  quart_g  M = 4   generic quartic: m = iroot(kN,4)+1, e_1 matched exactly, e_2 in a window.
  quart_s  M = 4   centrally symmetric quartic ((z^2 - r1^2)(z^2 - r2^2)), s = r1^2 + r2^2 window.
  sext_s   M = 6   centrally symmetric sextic, prod (z^2 - r_i^2): sigma_1 matched, sigma_2 window.
"""
from __future__ import annotations

from math import isqrt
from typing import Iterator, List, Tuple

import gmpy2

from .polys import centered

Rel = Tuple[Tuple[int, ...], int, int]   # (pieces, u, v)


def iroot(x: int, k: int) -> int:
    return int(gmpy2.iroot(x, k)[0])


def cos(N: int, A: int, cap: int) -> Iterator[Rel]:
    z = isqrt(N) + 1
    n = 0
    for r1 in range(-A, A + 1):
        for r2 in range(r1, A + 1):
            u = (z + r1) * (z + r2)
            yield (z + r1, z + r2), u, centered(u, N)
            n += 1
            if n >= cap:
                return


def quart_g(N: int, A: int, K: int, W: int, cap: int) -> Iterator[Rel]:
    """P(m) - kN = R0 - e1 m^3 + e2 m^2 - e3 m + e4 with R0 = m^4 - kN; match e1, window e2."""
    n = 0
    for k in range(1, K + 1):
        m = iroot(k * N, 4) + 1
        R0 = m ** 4 - k * N
        E1 = (R0 + m ** 3 // 2) // m ** 3
        R1 = R0 - E1 * m ** 3
        E2 = -((R1 + m ** 2 // 2) // m ** 2)
        # e1 = sum a = E1 ; p2 = sum a^2 = E1^2 - 2 e2 in [E1^2 - 2(E2+W), E1^2 - 2(E2-W)]
        lo_p2, hi_p2 = E1 * E1 - 2 * (E2 + W), E1 * E1 - 2 * (E2 - W)
        for a1 in range(-A, A + 1):
            for a2 in range(a1, A + 1):
                S = E1 - a1 - a2          # a3 + a4 = S, a2 <= a3 <= a4 <= A
                rem_lo, rem_hi = lo_p2 - a1 * a1 - a2 * a2, hi_p2 - a1 * a1 - a2 * a2
                if rem_hi < 0:
                    continue
                # a3 = (S - d)/2, a4 = (S + d)/2 with d >= 0, d = S mod 2, q = (S^2 + d^2)/2
                dlo2, dhi2 = 2 * rem_lo - S * S, 2 * rem_hi - S * S
                if dhi2 < 0:
                    continue
                dlo = isqrt(dlo2 - 1) + 1 if dlo2 > 0 else 0
                dhi = isqrt(dhi2)
                for d in range(dlo + ((dlo - S) & 1), dhi + 1, 2):
                    a3, a4 = (S - d) // 2, (S + d) // 2
                    if a3 < a2 or a4 > A:
                        continue
                    q = a3 * a3 + a4 * a4
                    if rem_lo <= q <= rem_hi:
                        roots = (a1, a2, a3, a4)
                        u = 1
                        for a in roots:
                            u *= (m - a)
                        yield tuple(m - a for a in roots), u, centered(u, N)
                        n += 1
                        if n >= cap:
                            return


def _two_squares(s: int, A: int) -> List[Tuple[int, int]]:
    out = []
    r1 = 0
    while 2 * r1 * r1 <= s and r1 <= A:
        t = s - r1 * r1
        r2 = isqrt(t)
        if r2 * r2 == t and r2 <= A:
            out.append((r1, r2))
        r1 += 1
    return out


def quart_s(N: int, A: int, K: int, W: int, cap: int) -> Iterator[Rel]:
    n = 0
    for k in range(1, K + 1):
        z = iroot(k * N, 4) + 1
        Y = z * z
        R0 = Y * Y - k * N
        s0 = (R0 + Y // 2) // Y
        for s in range(max(0, s0 - W), s0 + W + 1):
            for r1, r2 in _two_squares(s, A):
                u = (Y - r1 * r1) * (Y - r2 * r2)
                yield (z - r1, z + r1, z - r2, z + r2), u, centered(u, N)
                n += 1
                if n >= cap:
                    return


def sext_s(N: int, A: int, K: int, W1: int, W2: int, cap: int) -> Iterator[Rel]:
    """Q(Y) = Y^3 - s1 Y^2 + s2 Y - s3, s_i = e_i(r^2), Y = z^2, z = iroot(kN,6)+1.
    s1 ranges over [S1, S1 + W1] (S1 = nearest), s2 must be within W2 of its aligned target."""
    n = 0
    sq = [r * r for r in range(A + 1)]
    for k in range(1, K + 1):
        z = iroot(k * N, 6) + 1
        Y = z * z
        R0 = Y ** 3 - k * N
        S1 = (R0 + Y * Y // 2) // (Y * Y)
        for s1 in range(max(0, S1), S1 + W1 + 1):
            R1 = R0 - s1 * Y * Y
            S2 = -((R1 + Y // 2) // Y)
            if S2 + W2 < 0:
                continue
            for r1 in range(A + 1):
                if 3 * sq[r1] > s1:
                    break
                for r2 in range(r1, A + 1):
                    t = s1 - sq[r1] - sq[r2]
                    if t < sq[r2]:
                        break
                    r3 = isqrt(t)
                    if r3 * r3 != t or r3 > A:
                        continue
                    x1, x2, x3 = sq[r1], sq[r2], sq[r3]
                    s2 = x1 * x2 + x1 * x3 + x2 * x3
                    if abs(s2 - S2) > W2:
                        continue
                    u = (Y - x1) * (Y - x2) * (Y - x3)
                    yield (z - r1, z + r1, z - r2, z + r2, z - r3, z + r3), u, centered(u, N)
                    n += 1
                    if n >= cap:
                        return


def quart_rand(N: int, A: int, seed: int, cap: int) -> Iterator[Rel]:
    """Control: random quartic multisets in [-A, A] at m = iroot(N,4)+1, NO alignment."""
    import random
    rng = random.Random(f"ptecons|quart_rand|{N}|{seed}")
    m = iroot(N, 4) + 1
    for _ in range(cap):
        roots = sorted(rng.randint(-A, A) for _ in range(4))
        u = 1
        for a in roots:
            u *= (m - a)
        yield tuple(m - a for a in roots), u, centered(u, N)
