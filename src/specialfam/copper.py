"""Known-residue Coppersmith / Howgrave-Graham decoder with interval subdivision.

OPERATIONAL. Adapted (copied, rule 7) from src/trackE/coppersmith.py: same lattice
[May 2010, proof of Thm 3, delta = 1], same rigorous sufficient condition, flint LLL
(fpylll is not installed). Integer roots of the reduced polynomial use flint's
fmpz_poly.roots() over Z (standard last step; never touches N).

Task: find a divisor b of N with b = r + d*x, x in [x_lo, x_hi], b >= L (L public).
Write x = xc + y, |y| <= Y. With c = r*d^{-1} + xc (mod N), f(y) = y + c is monic and
f(y0) = d^{-1} b == 0 (mod b). Rigorous condition (Howgrave-Graham + LLL bound):
    2^{D(D-1)} det^4 D^{2D} < L^{4 m D},  det = N^{m(m+1)/2} Y^{D(D-1)/2}.
Subdividing [x_lo, x_hi] into chunks of width 2Y+1 costs ceil(len/(2Y+1)) LLL calls.
"""
from __future__ import annotations

import math
from math import comb, gcd
from typing import List, Optional, Tuple

import flint

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import NonUnit


def rigorous_ok(N: int, L: int, Y: int, D: int, m: int) -> bool:
    det = N ** (m * (m + 1) // 2) * Y ** (D * (D - 1) // 2)
    return (1 << (D * (D - 1))) * det ** 4 * D ** (2 * D) < L ** (4 * m * D)


def max_rigorous_Y(N: int, L: int, D: int, m: int) -> int:
    """Largest Y >= 0 satisfying rigorous_ok (0 if none). Exact integer binary search
    (rigorous_ok is monotone in Y); no float overflow for huge N (codex specialfam_1 #1)."""
    if not rigorous_ok(N, L, 1, D, m):
        return 0
    lnN, lnL = N.bit_length(), L.bit_length()
    est = (4 * m * D * lnL - 2 * m * (m + 1) * (lnN - 1)) // (2 * D * (D - 1)) + 2
    hi = 1 << max(1, est)
    while rigorous_ok(N, L, hi, D, m):  # safety: grow until it fails
        hi <<= 1
    lo = 1  # invariant: ok(lo) and not ok(hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rigorous_ok(N, L, mid, D, m):
            lo = mid
        else:
            hi = mid
    return lo


def best_params(N: int, L: int, Dmax: int = 13) -> Tuple[int, int, int]:
    """(D, m, Y) maximising the rigorous Y over odd D <= Dmax and m near beta*D, where
    beta = log L / log N (codex specialfam_1 #1: m ~ beta D, not D/2, for unbalanced beta).
    The theorem (May 2010, Thm 1) lets D grow like log N; experiments cap Dmax = 9."""
    beta = (L.bit_length() - 1) / max(1, N.bit_length())
    best = (1, 0, 0)
    for D in range(3, Dmax + 1, 2):
        cands = set(range(max(1, (D - 1) // 2 - 1), min(D, (D - 1) // 2 + 2)))
        c = round(beta * D)
        cands |= {x for x in (c - 1, c, c + 1) if 1 <= x < D}
        for m in sorted(cands):
            Y = max_rigorous_Y(N, L, D, m)
            if Y > best[2]:
                best = (D, m, Y)
    return best


def _basis(N: int, c: int, Y: int, D: int, m: int) -> List[List[int]]:
    fp = [[comb(k, i) * c ** (k - i) for i in range(k + 1)] for k in range(m + 1)]
    Yp = [Y ** j for j in range(D)]
    rows = []
    for i in range(m):
        Ni = N ** (m - i)
        row = [0] * D
        for j, a in enumerate(fp[i]):
            row[j] = Ni * a * Yp[j]
        rows.append(row)
    for j in range(D - m):
        row = [0] * D
        for i2, a in enumerate(fp[m]):
            row[i2 + j] = a * Yp[i2 + j]
        rows.append(row)
    return rows


def decode_chunk(N: int, d: int, r: int, xc: int, Y: int, D: int, m: int,
                 mt: Metrics, n_vectors: int = 2) -> Optional[int]:
    """One lattice call: proper divisor b = r + d*(xc+y), |y| <= Y, or None."""
    g = gcd(d, N)
    if 1 < g < N:
        return g
    dinv = pow(d, -1, N) if d % N != 1 else 1
    c = (r * dinv + xc) % N
    B = _basis(N, c, max(Y, 1), D, m)
    mt.inc("lll_calls")
    red = flint.fmpz_mat(B).lll().tolist()
    Yp = [max(Y, 1) ** j for j in range(D)]
    red.sort(key=lambda v: sum(int(a) * int(a) for a in v))
    for row in red[:n_vectors]:
        coeffs = []
        for j, a in enumerate(row):
            a = int(a)
            coeffs.append(a // Yp[j])
        poly = flint.fmpz_poly(coeffs)
        if poly.degree() < 1:
            continue
        mt.inc("root_calls")
        for y0, _ in poly.roots():
            y0 = int(y0)
            if abs(y0) > Y:
                continue
            b = r + d * (xc + y0)
            mt.inc("gcd")
            h = gcd(b, N)
            if 1 < h < N:
                return h
    return None


def subdivide_search(N: int, d: int, r: int, P_lo: int, P_hi: int, L: int,
                     params: Tuple[int, int, int], mt: Metrics,
                     max_chunks: Optional[int] = None) -> Tuple[Optional[int], int]:
    """Find a divisor b = r + d x of N in [P_lo, P_hi] by chunked Coppersmith calls.

    params = (D, m, Y) rigorous for lower bound L <= every target b. Returns
    (factor or None, chunks used).
    """
    D, m, Y = params
    x_lo = -(-(P_lo - r) // d)
    x_hi = (P_hi - r) // d
    if x_hi < x_lo:
        return None, 0
    width = 2 * Y + 1
    if Y <= 2:  # tiny lattice range: direct enumeration is cheaper
        chunks = 0
        for x in range(x_lo, x_hi + 1):
            chunks += 1
            h = gcd(r + d * x, N)
            mt.inc("gcd")
            if 1 < h < N:
                return h, chunks
        return None, chunks
    nchunks = -(-(x_hi - x_lo + 1) // width)
    used = 0
    for t in range(nchunks):
        if max_chunks is not None and used >= max_chunks:
            break
        xc = x_lo + t * width + Y
        used += 1
        f = decode_chunk(N, d, r, xc, Y, D, m, mt)
        if f:
            return f, used
    return None, used


def n_chunks(d: int, P_lo: int, P_hi: int, Y: int) -> int:
    span = (P_hi - P_lo) // d + 1
    return -(-span // (2 * Y + 1)) if Y > 2 else span
