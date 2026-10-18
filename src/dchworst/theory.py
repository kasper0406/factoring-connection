"""Closed forms used in notes/dchworst_results.md (model-level, exact where possible).

* worst_case_exponent(c): DCL family with J moduli near Y, weights of height H,
  cost X = H P / Y^2 + J H^2 (+ J), and a coverage GUARANTEE H = sqrt(Y) J^{-c}
  for every pair.  Returns the optimal exponent of N (N = P^2):
        e(c) = min(1/5, (1 - c) / (7/2 - c))       (Theorem W).
* poisson_*: the Poisson heuristic for uncovered pairs.
"""
from __future__ import annotations

from math import exp, log
from typing import Sequence


def worst_case_exponent(c: float) -> float:
    assert 0.0 <= c <= 0.5
    return min(0.2, (1.0 - c) / (3.5 - c))


def worst_case_exponent_numeric(c: float, grid: int = 400) -> float:
    """Brute-force min over y = log_P Y in (0, 1) and theta = log_Y J in [0, 1] of
    max(1 - 1.5y - c theta y, y + theta y (1 - 2c)), halved (exponent of N)."""
    best = 9.0
    for i in range(1, grid):
        y = i / grid
        for k in range(grid + 1):
            th = k / grid
            e = max(1 - 1.5 * y - c * th * y, y + th * y * (1 - 2 * c), th * y)
            best = min(best, e)
    return best / 2


def lambda_design(nV: int, primes: Sequence[int]) -> float:
    """lambda = |V| * sum 1/(l-1): expected number of covering (v, l) events."""
    return nV * sum(1.0 / (l - 1) for l in primes)


def poisson_uncovered(count: float, lam: float) -> float:
    return count * exp(-lam)


def full_coverage_H_threshold(n: int, blocks_lambda_per_H2: float, pairs: float) -> float:
    """Smallest H with pairs * exp(-k H^2) < 1 where lambda = k H^2."""
    return (log(pairs) / blocks_lambda_per_H2) ** 0.5
