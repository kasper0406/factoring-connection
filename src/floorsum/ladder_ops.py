"""Operational degree-j ladder (Theorem L_j) with the exact Farey/hull floor-sum oracle.

OPERATIONAL: `ladder_factor_op(N, j)` sees only N.  It is src/analytic/ladder.py's
`ladder_factor` with the oracle made pluggable (copied, not edited: file ownership), and with the
simulated brute-force oracle replaced by `FloorSumConvex` (src/floorsum/farey.py), which is exact
for every degree j because P_j is convex on every block (P_j'' > 0 for 0 <= h < H <= x/2).

Cost instrumentation: the oracle's counters (`evals`, `windows`, `runs`, `bsearch`,
`euclid_calls`, `queries`) are the model cost; `calls` is the Theorem L_j call count.
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt
from typing import Dict, List, Optional

from src.analytic.ladder import _Block, _candidate_count, _locate, iroot_floor
from src.floorsum.farey import Counters, FloorSumConvex


class FareyOracle:
    """Same interface as src.analytic.ladder.Oracle: floor_range(blk, shift, h1, h2)."""

    def __init__(self):
        self.calls = 0
        self.cnt = Counters()
        self.blocks_built = 0

    def floor_range(self, blk: _Block, shift: Fraction, h1: int, h2: int) -> int:
        self.calls += 1
        if h2 <= h1:
            return 0
        cache = getattr(blk, "_fs_cache", None)
        if cache is None:
            cache = {}
            blk._fs_cache = cache
        if shift not in cache:
            # P(h) + u/v = (v num(h) + u den) / (v den)
            u, v = shift.numerator, shift.denominator
            cache[shift] = FloorSumConvex(_ladder_coef(blk, u, v), blk.den * v, 0, blk.H, self.cnt)
            self.blocks_built += 1
        return cache[shift].range_sum(h1, h2)


def _ladder_coef(blk: _Block, u: int, v: int) -> List[int]:
    """Integer coefficients (ascending powers of h) of v*num(h) + u*den.

    In src/analytic/ladder.py, _Block.coef[i] = (-1)^i N x^{j-i} is the coefficient of h^i."""
    c = [v * a for a in blk.coef]
    c[0] += u * blk.den
    return c


def ladder_factor_op(N: int, j: int, c: Fraction = Fraction(1, 2),
                     stats: Optional[Dict] = None) -> Optional[int]:
    """Smallest-found nontrivial factor of N via the degree-j ladder with the exact Farey oracle,
    or None if N has no divisor in [2, sqrt N] (N prime or 1)."""
    if stats is None:
        stats = {}
    for k in ("blocks", "candidates", "trial_divisions", "divisor_tests", "located"):
        stats.setdefault(k, 0)
    orc = FareyOracle()
    sN = isqrt(N)
    B = iroot_floor(N, j + 2)
    x0 = min(sN + 1, -((-B * c.denominator) // c.numerator))
    found = None
    for d in range(2, min(x0, sN + 1)):
        stats["trial_divisions"] += 1
        if N % d == 0:
            found = d
            stats["found_in"] = "trial"
            break
    x = x0
    while found is None and x <= sN:
        H = max(1, (c.numerator * x) // (c.denominator * B))
        H = min(H, sN - x + 1)
        blk = _Block(N, x, H, j)
        stats["blocks"] += 1
        cands: List[int] = []
        tot = _candidate_count(orc, blk, 0, H)
        _locate(orc, blk, 0, H, tot, cands)
        stats["located"] += len(cands)
        if j % 2 == 1 and 0 not in cands:
            cands.append(0)
        if j % 2 == 0 and 0 not in cands:
            cands.append(0)  # R(0) = 0 also for even j; one forced division per block
        stats["candidates"] += len(cands)
        for h in sorted(cands):
            stats["divisor_tests"] += 1
            d = x + h
            if 1 < d < N and N % d == 0:
                found = d
                stats["found_in"] = "block"
                break
        x += H
    stats["oracle_calls"] = orc.calls
    stats["oracle_builds"] = orc.blocks_built
    stats.update({f"op_{k}": v for k, v in orc.cnt.items()})
    return found


def divisor_summatory_ladder(N: int, j: int, c: Fraction = Fraction(1, 2),
                             stats: Optional[Dict] = None) -> int:
    """D(N) = sum_{n <= N} tau(n) = 2 sum_{d <= s} floor(N/d) - s^2 (s = isqrt N), computed with
    the degree-j ladder (Prop D, notes/floorsum_results.md).  OPERATIONAL, exact.

    Per block, sum_{h<H} floor(N/(x+h)) = F(P_j; 0, H) - sum_{h in C} (floor(P_j(h)) - floor(N/(x+h)))
    where C is the located candidate set (it contains every h with floor(P_j) != floor(f)).
    With the Farey oracle the cost is ~ N^{1/3} polylog; with a polylog oracle it would be
    O(N^{1/(j+2)} log^2 N) oracle calls."""
    if stats is None:
        stats = {}
    stats.setdefault("blocks", 0)
    stats.setdefault("located", 0)
    orc = FareyOracle()
    s = isqrt(N)
    B = iroot_floor(N, j + 2)
    x0 = min(s + 1, -((-B * c.denominator) // c.numerator))
    tot = sum(N // d for d in range(1, min(x0, s + 1)))
    x = x0
    while x <= s:
        H = max(1, (c.numerator * x) // (c.denominator * B))
        H = min(H, s - x + 1)
        blk = _Block(N, x, H, j)
        stats["blocks"] += 1
        cands: List[int] = []
        cnt = _candidate_count(orc, blk, 0, H)
        _locate(orc, blk, 0, H, cnt, cands)
        stats["located"] += len(cands)
        F = orc.floor_range(blk, Fraction(0), 0, H)
        corr = 0
        for h in cands:
            corr += blk.num(h) // blk.den - N // (x + h)
        tot += F - corr
        x += H
    stats["oracle_calls"] = orc.calls
    stats.update({f"op_{k}": v for k, v in orc.cnt.items()})
    return 2 * tot - s * s
