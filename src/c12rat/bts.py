"""Lemma BT-S (smooth numbers in an arithmetic progression in a short interval; notes/c12rat_results.md §5).
DIAGNOSTIC: exact check of the explicit inequality

  #{n in (X0, X0+Z] : n = a mod D, P+(n) <= y}  <=  sum_{m1 in (W/y, W], P+(m1) <= y, (m1, D) = 1} (floor(Z/(m1 D)) + 1)

valid for gcd(a, D) = 1 and 1 <= W <= X0 (every y-smooth n > W has a divisor in (W/y, W]).
"""
from __future__ import annotations

from math import gcd

import numpy as np


def lhs(smooth: np.ndarray, X0: int, Z: int, a: int, D: int) -> int:
    start = X0 + 1 + ((a - (X0 + 1)) % D)
    idx = np.arange(start, X0 + Z + 1, D)
    return int(smooth[idx].sum())


def rhs(smooth: np.ndarray, X0: int, Z: int, D: int, y: int, W: int) -> int:
    lo = W // y + 1
    m = np.arange(lo, W + 1)
    m = m[smooth[m]]
    m = m[np.gcd(m, D) == 1]
    return int((Z // (m * D) + 1).sum())


def rhs_analytic(smooth: np.ndarray, Z: int, D: int, y: int, W: int) -> float:
    """(Z/D) sum 1/m1 + Psi(W) >= rhs (dropping the floors and the coprimality)."""
    lo = W // y + 1
    m = np.arange(lo, W + 1)
    m = m[smooth[m]]
    return float(Z / D * (1.0 / m).sum() + len(m))
