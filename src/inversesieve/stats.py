"""Structure statistics for uncovered sets (exact integer arithmetic; floats only for reported ratios).

Every statistic comes with a null model:
  * 'uniform'  -- the same number of distinct integers drawn uniformly from the same interval
                  (seeded numpy Generator); Monte Carlo.
  * 'local'    -- the independent-residue model of the design (exact pair densities), used for
                  pair correlations.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .design import forbidden_mask


# ---------------------------------------------------------------------------
# covering numbers
# ---------------------------------------------------------------------------
def interval_cover_number(xs: np.ndarray, L: int) -> int:
    """Minimum number of closed integer intervals of length L (L consecutive integers) covering xs.
    Greedy from the left is optimal in one dimension."""
    if len(xs) == 0:
        return 0
    xs = np.sort(xs)
    c, i, n = 0, 0, len(xs)
    while i < n:
        end = xs[i] + L - 1
        i = int(np.searchsorted(xs, end, side="right"))
        c += 1
    return c


def ap_triples(xs: np.ndarray, Lmin: int) -> int:
    """Number of 3-subsets {a<b<c} of xs lying in one residue class modulo some M >= Lmin, i.e.
    gcd(b-a, c-a) >= Lmin.  (Coppersmith containers: p = r mod M with M >= P^{1/2}.)"""
    xs = np.sort(np.asarray(xs, dtype=np.int64))
    n, cnt = len(xs), 0
    for i in range(n - 2):
        d = xs[i + 1:] - xs[i]
        g = np.gcd.outer(d, d)
        cnt += int(np.count_nonzero(np.triu(g >= Lmin, 1)))
    return cnt


def ap_pairs_large_modulus(xs: np.ndarray, Lmin: int, P: int) -> int:
    """Pairs with difference d admitting M | d, M >= Lmin, and window P/M <= ... : always true for
    d >= Lmin; reported only as a sanity count (every pair is AP-coverable)."""
    xs = np.sort(xs)
    d = (xs[None, :] - xs[:, None])[np.triu_indices(len(xs), 1)]
    return int(np.count_nonzero(d >= Lmin))


# ---------------------------------------------------------------------------
# class occupancy (low-degree image detector, Green--Harper / Walsh signature)
# ---------------------------------------------------------------------------
def occupancy(xs: np.ndarray, r: int) -> Tuple[int, float]:
    """(#classes mod r occupied, chi^2 statistic of class counts against uniform)."""
    c = np.bincount(np.asarray(xs) % r, minlength=r).astype(np.int64)
    n = int(c.sum())
    exp = n / r
    chi2 = float(((c - exp) ** 2).sum() / exp) if n else 0.0
    return int(np.count_nonzero(c)), chi2


# ---------------------------------------------------------------------------
# additive energy
# ---------------------------------------------------------------------------
def energy_mod(A: np.ndarray, l: int) -> int:
    """Additive energy E(A) = #{(a,b,c,d) in A^4 : a + b = c + d (mod l)} (exact)."""
    A = np.asarray(A, dtype=np.int64)
    s = ((A[:, None] + A[None, :]) % l).ravel()
    r = np.bincount(s, minlength=l).astype(np.int64)
    return int((r * r).sum())


def energy_int(xs: np.ndarray) -> int:
    """Additive energy of a set of integers (exact)."""
    xs = np.asarray(xs, dtype=np.int64)
    s = (xs[:, None] + xs[None, :]).ravel()
    _, c = np.unique(s, return_counts=True)
    return int((c.astype(np.int64) ** 2).sum())


def mult_energy_mod(A: np.ndarray, l: int) -> int:
    """Multiplicative energy #{ab = cd mod l}, A subset of (Z/l)^* (exact)."""
    A = np.asarray(A, dtype=np.int64)
    A = A[A % l != 0]
    s = ((A[:, None] * A[None, :]) % l).ravel()
    _, c = np.unique(s, return_counts=True)
    return int((c.astype(np.int64) ** 2).sum())


# ---------------------------------------------------------------------------
# gaps and pair correlation
# ---------------------------------------------------------------------------
def pair_counts(xs: np.ndarray, D: int) -> np.ndarray:
    """c[d] = #{(x, x+d) both in xs}, d = 0..D (c[0] = |xs|)."""
    s = set(int(v) for v in xs)
    c = np.zeros(D + 1, dtype=np.int64)
    c[0] = len(s)
    for d in range(1, D + 1):
        c[d] = sum(1 for v in s if v + d in s)
    return c


def local_pair_density(N: int, primes, V, D: int, Rs: dict = None) -> np.ndarray:
    """Independent-residue model: g[d] = prod_l (1 - |R u (R - d)| / l) for d = 0..D (g[0] = single density).
    Floats for reporting (each factor is an exact rational)."""
    from .design import forbidden_set
    g = np.ones(D + 1)
    for l in primes:
        R = forbidden_set(N, l, V) if Rs is None else Rs[l]
        diff = ((R[None, :] - R[:, None]) % l).ravel()       # y - x for (x, y) in R^2
        rep = np.bincount(diff[diff <= D], minlength=D + 1)[:D + 1]   # #{x in R : x + d in R}
        bad = 2 * len(R) - rep                                 # |R u (R - d)|  (d = 0: |R|)
        bad[0] = len(R)
        g *= 1.0 - bad / l
    return g


# ---------------------------------------------------------------------------
# Lehman / Coppersmith containers
# ---------------------------------------------------------------------------
def lehman_centres(N: int, h0: int) -> np.ndarray:
    """Stationary points s_v = round(sqrt(A N / B)) for positive primitive weights of height <= h0
    whose stationary point lies in [sqrt(N/2), sqrt(N)] (computed exactly by isqrt). Sees only N."""
    out = []
    for B in range(1, h0 + 1):
        for A in range(1, h0 + 1):
            if gcd(A, B) != 1:
                continue
            if not (B <= 2 * A and A <= B):     # s_v = sqrt(AN/B) in [sqrt(N/2), sqrt(N)] iff 1/2 <= A/B <= 1
                continue
            out.append(isqrt(A * N // B))
    return np.array(sorted(set(out)), dtype=np.int64)


def in_containers(xs: np.ndarray, centres: np.ndarray, half: int) -> np.ndarray:
    """Boolean mask: x within distance `half` of some centre."""
    xs = np.asarray(xs, dtype=np.int64)
    if centres.size == 0:
        return np.zeros(len(xs), dtype=bool)
    i = np.searchsorted(centres, xs)
    lo = np.abs(xs - centres[np.clip(i - 1, 0, len(centres) - 1)])
    hi = np.abs(centres[np.clip(i, 0, len(centres) - 1)] - xs)
    return np.minimum(lo, hi) <= half


def container_measure(centres: np.ndarray, half: int, a0: int, b0: int) -> int:
    """Number of integers of [a0, b0] within distance `half` of some centre (exact)."""
    if centres.size == 0:
        return 0
    iv = sorted((max(a0, int(c) - half), min(b0, int(c) + half)) for c in centres)
    tot, cur_lo, cur_hi = 0, None, None
    for lo, hi in iv:
        if lo > hi:
            continue
        if cur_hi is None or lo > cur_hi + 1:
            if cur_hi is not None:
                tot += cur_hi - cur_lo + 1
            cur_lo, cur_hi = lo, hi
        else:
            cur_hi = max(cur_hi, hi)
    if cur_hi is not None:
        tot += cur_hi - cur_lo + 1
    return tot


def lstar_score(p: int, q: int, hmax: int) -> float:
    """HARNESS: coverage38's Dirichlet-optimal Lehman score L*(p,q) = min_v max(h^2, |Bp-Aq| sqrt(h/P)),
    over positive primitive weights of height <= hmax (float, report only)."""
    P = (p * q) ** 0.5
    best = float("inf")
    for B in range(1, hmax + 1):
        for A in range(1, hmax + 1):
            if gcd(A, B) != 1:
                continue
            h = max(A, B)
            s = max(h * h, abs(B * p - A * q) * (h / P) ** 0.5)
            best = min(best, s)
    return best


# ---------------------------------------------------------------------------
# null samplers
# ---------------------------------------------------------------------------
def uniform_null(rng: np.random.Generator, a0: int, b0: int, k: int) -> np.ndarray:
    return np.sort(rng.choice(np.arange(a0, b0 + 1, dtype=np.int64), size=k, replace=False))
