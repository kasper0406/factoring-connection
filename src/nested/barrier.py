"""Exact helper computations for notes/nested_results.md (barrier propositions).  DIAGNOSTIC / analysis only.

Functions here are elementary and exact (Fractions / integers); floats only for display.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from typing import Dict, List

from src.thmA import arith


def is_smooth(x: int, B: int) -> bool:
    """x != 0 is B-smooth (trial division; for tests only)."""
    x = abs(x)
    assert x > 0
    for p in arith.primes_upto(B):
        while x % p == 0:
            x //= p
        if x == 1:
            return True
    return x == 1


# ---------------------------------------------------------------- Theorem AF: Gamma law of the cofactor norm
def gamma_upper_tail(r: int, V: Fraction) -> Fraction:
    """Pr[Gamma(r,1) >= V] for integer shape r >= 1, as exp(-V) * sum_{j<r} V^j/j!  (returned WITHOUT the
    exp(-V) factor, i.e. the Poisson partial sum; multiply by exp(-V) in floating point for display)."""
    s, term = Fraction(0), Fraction(1)
    for j in range(r):
        if j:
            term = term * V / j
        s += term
    return s


def af_weight(r: int, X: Fraction, T: Fraction) -> float:
    """Upper bound factor in Cor AF: Pr_v[X e^{-v} <= T] for v ~ Gamma(r,1), i.e. Pr[v >= log(X/T)]."""
    import math
    V = math.log(float(X / T)) if X > T else 0.0
    if V <= 0:
        return 1.0
    s, term = 0.0, 1.0
    for j in range(r):
        if j:
            term *= V / j
        s += term
    return math.exp(-V) * s


# ---------------------------------------------------------------- Prop OD: same-field children vs parameters
def min_children(d: int, lam: Fraction) -> Fraction:
    """Size-model lower bound for the number of algebraic children of a node: (d/2)/lam (log|F| >= (d/2) log q,
    children <= q^lam)."""
    return Fraction(d, 2) / lam


def overdetermined(d: int, lam: Fraction, depth: int, params_extra: int = 3) -> bool:
    """True iff a parent+children configuration of the given depth (1: parent with children; 2: plus
    grandchildren) has more algebraic values than the field has parameters (d + params_extra)."""
    m = min_children(d, lam)
    nodes = sum(m ** j for j in range(depth + 1))
    return nodes > d + params_extra


# ---------------------------------------------------------------- Prop BT: backtracking cost exponent (model)
def backtrack_log_cost(m: Fraction, levels: int, log_inv_p: List[Fraction]) -> Fraction:
    """Model: log of the multiplicative search overhead prod_j (m * p_j^{-m}) = sum_j (log m + m log(1/p_j));
    log m is passed implicitly as 0 here (it is O(levels * log m), negligible); returns m * sum log(1/p_j)."""
    assert len(log_inv_p) == levels
    return m * sum(log_inv_p)


# ---------------------------------------------------------------- Theorem NR: the constant
def nr_constant(b: float, b1: float, s: float, delta: float) -> float:
    """c_NR = max(2b, (b + 1/b1)/(1 - 2^{-1/3}), (b + 2s)/(1 - tau_F^{1/3})), tau_F = delta*b1 < 1 (reduced-basis box,
    ||P|| <= 3 e^{sE} q; after proof critic issue 1)."""
    tau = delta * b1
    assert 0 < tau < 1
    return max(2 * b, (b + 1 / b1) / (1 - 2 ** (-1 / 3)), (b + 2 * s) / (1 - tau ** (1 / 3)))


def ad_heuristic_ok(b: float, s: float, delta: float, lam: float, kappa_F: float, m0: float) -> bool:
    """Standard heuristic compatibility of Hypothesis AD (Poisson model, reduced-basis box, worst-case sizes):
    2s > (1/3)[(kappa_F + m0)/(lam b) + delta/lam + delta s/(lam b)]."""
    rhs = ((kappa_F + m0) / (lam * b) + delta / lam + delta * s / (lam * b)) / 3
    return 2 * s > rhs
