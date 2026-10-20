"""Coppersmith / Howgrave-Graham divisor-in-interval block oracle (KNOWN algorithm; E5).

OPERATIONAL (sees only N).  Finds every divisor p of N with x <= p < x + X, provided
X satisfies the exact Howgrave-Graham/LLL determinant inequality below (then success is
guaranteed, not heuristic).  This is the classical deterministic N^{1/4 + o(1)} method
(Coppersmith, J. Cryptology 10 (1997); Howgrave-Graham, CaLC 2001 "Approximate integer common
divisors"); it is included to show that the j = 2 ladder's REAL need -- locating divisors of N in a
block of length H = x N^{-1/4} -- already has a polylog block oracle, by lattice reduction rather
than by floor sums.  It does not compute floor sums.

Lattice: f(t) = x + t, root t0 = p - x in [0, X), p | f(t0).  Basis polynomials
  g_i(t) = N^{m-i} f(t)^i  (0 <= i <= m),   h_k(t) = t^k f(t)^m  (1 <= k <= n-m-1),
evaluated at tX (coefficient vectors), lower triangular, det = N^{m(m+1)/2} X^{n(n-1)/2}.
p^m divides every basis polynomial at t0.  LLL (delta = 0.99) gives b1 with
||b1|| <= 2^{(n-1)/4} det^{1/n}; if  sqrt(n) ||b1|| < x^m <= p^m  then b1(t0) = 0 over Z.
Exact integer form of the sufficient condition used here:
  2^{n(n-1)} * det^4 * n^{2n} < x^{4mn}.
"""
from __future__ import annotations

import math
from math import comb, isqrt
from typing import Dict, List, Optional, Tuple

import flint
import gmpy2


def _cond_ok(N: int, x: int, X: int, m: int, n: int) -> bool:
    det4 = N ** (2 * m * (m + 1)) * X ** (2 * n * (n - 1))
    return (1 << (n * (n - 1))) * det4 * n ** (2 * n) < x ** (4 * m * n)


def _logX(N: int, x: int, m: int, n: int) -> float:
    """Float estimate of log max X (parameter choice only; the exact check is _cond_ok)."""
    lN, lx = math.log(N), math.log(x)
    return (4 * m * n * lx - n * (n - 1) * math.log(2) - 2 * m * (m + 1) * lN
            - 2 * n * math.log(n)) / (2 * n * (n - 1))


def best_params(N: int, x: int, nmax: int = 16) -> Tuple[int, int, int]:
    """(m, n, X): lattice parameters and the largest X passing the exact condition."""
    best = (1, 2, 0)
    bestl = -1e300
    for n in range(2, nmax + 1):
        for m in range(1, n):
            l = _logX(N, x, m, n)
            if l > bestl:
                bestl, best = l, (m, n, 0)
    m, n, _ = best
    if bestl < 0:
        return m, n, 0
    # exact: largest X with _cond_ok, by doubling + bisection around the estimate
    lo, hi = 0, max(2, int(math.exp(min(bestl, 700))) * 2 + 2)
    while _cond_ok(N, x, hi, m, n):
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _cond_ok(N, x, mid, m, n):
            lo = mid
        else:
            hi = mid
    return m, n, lo


def _poly_mul(a: List[int], b: List[int]) -> List[int]:
    r = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for k, v in enumerate(b):
                r[i + k] += u * v
    return r


def lattice_block(N: int, x: int, m: int, n: int, X: int, stats: Optional[Dict] = None
                  ) -> List[int]:
    """All divisors p of N with x <= p < x + X (guaranteed when _cond_ok(N, x, X, m, n))."""
    if stats is not None:
        stats["lll_calls"] = stats.get("lll_calls", 0) + 1
        stats["max_dim"] = max(stats.get("max_dim", 0), n)
    f = [x, 1]
    fp = [[1]]
    for i in range(m):
        fp.append(_poly_mul(fp[-1], f))
    rows = []
    for i in range(m + 1):
        g = [c * N ** (m - i) for c in fp[i]]
        rows.append(g)
    for k in range(1, n - m):
        rows.append([0] * k + fp[m])
    Xp = [X ** i for i in range(n)]
    M = flint.fmpz_mat([[(r[i] if i < len(r) else 0) * Xp[i] for i in range(n)] for r in rows])
    L = M.lll()
    out = set()
    for row_i in range(min(3, n)):
        coeffs = [int(L[row_i, i]) // Xp[i] for i in range(n)]
        if not any(coeffs):
            continue
        poly = flint.fmpz_poly(coeffs)
        for r, _mult in poly.roots():
            t = int(r)
            if 0 <= t < X:
                p = x + t
                if 1 < p < N and N % p == 0:
                    out.add(p)
    return sorted(out)


def lattice_factor(N: int, nmax: int = 16, stats: Optional[Dict] = None) -> Optional[int]:
    """Smallest-found nontrivial factor of N by trial division below x1 and lattice blocks on
    [x1, sqrt N], where x1 is the first x with a block length X(x) >= 2.  None if N is prime."""
    if stats is None:
        stats = {}
    stats.setdefault("blocks", 0)
    stats.setdefault("trial_divisions", 0)
    sN = isqrt(N)
    # find x1 by doubling from 2
    x1 = 2
    while x1 <= sN and best_params(N, x1, nmax)[2] < 2:
        x1 *= 2
    x1 = min(x1, sN + 1)
    stats["x1"] = x1
    for d in range(2, x1):
        stats["trial_divisions"] += 1
        if N % d == 0:
            stats["found_in"] = "trial"
            return d
    x = x1
    params = None
    while x <= sN:
        if params is None or stats["blocks"] % 64 == 0:
            params = best_params(N, x, nmax)  # X(x) is increasing in x; refresh periodically
        m, n, X = params
        if not _cond_ok(N, x, X, m, n):
            params = best_params(N, x, nmax)
            m, n, X = params
        X = max(1, min(X, sN - x + 1))
        stats["blocks"] += 1
        ps = lattice_block(N, x, m, n, X, stats)
        if ps:
            stats["found_in"] = "lattice"
            return ps[0]
        x += X
    return None
