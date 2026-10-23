"""multidescent size model and bookkeeping (exact rationals).  MODEL arithmetic, not proofs of supply.

Units: exponents of D = t/ln t (t = ln N), with the common factor ln t suppressed, as in recurse §3.4 and
escapeR2 §7.2 (chain sub-model).  A level-i target has size L[zeta_i] (log-size D^{zeta_i} ln t); the field has
degree d = D^delta and base-m heights log M = log H = t/(d+1) (exponent 1 - delta).

j-level chain (j = number of descent levels; the last level's values are B-smooth, B = L[alpha]):
  zeta_0 = 1 - alpha (split), zeta_j = alpha (factor base),
  rational value exponent 1 - delta, algebraic value exponent max(1 - delta, delta + zeta_i),
  smoothness over L[zeta_{i+1}] at cost L[alpha]:  value - zeta_{i+1} <= alpha.
Closed form: alpha_j = (j+1)/(3j+2), delta_j = j/(3j+2), zeta_i = 1 - alpha - i(3 alpha - 1)
(= escapeR2's (r+2)/(3r+5) with r = j - 1; codex multidescent_1 independently).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Optional


def alpha_j(j: int) -> Fr:
    assert j >= 1
    return Fr(j + 1, 3 * j + 2)


def chain_feasible(alpha: Fr, j: int) -> Optional[Dict]:
    """Exact feasibility of the j-level chain at exponent alpha; returns a witness (delta, zetas) or None.
    Constraints: delta >= 1 - 2 alpha (last level: 1 - delta <= alpha + zeta_j = 2 alpha), and
    delta + zeta_i - zeta_{i+1} <= alpha for i < j (algebraic values), zeta_i >= zeta_{i+1}, rational values
    1 - delta <= alpha + zeta_{i+1}.  Summing the algebraic constraints: j delta + (1 - 2 alpha) <= j alpha."""
    lo = 1 - 2 * alpha
    hi = alpha - (1 - 2 * alpha) / j
    if lo > hi or lo <= 0:
        return None
    delta = lo
    step = (1 - 2 * alpha) / j               # uniform descent from 1 - alpha to alpha; needs delta + step <= alpha
    zetas = [1 - alpha - i * step for i in range(j + 1)]
    assert zetas[-1] == alpha
    # verify every constraint exactly
    for i in range(j):
        if delta + zetas[i] - zetas[i + 1] > alpha:
            return None
        if 1 - delta > alpha + zetas[i + 1]:
            return None
        if max(1 - delta, delta + zetas[i]) - zetas[i + 1] > alpha:
            return None
    return {"alpha": alpha, "delta": delta, "zetas": zetas}


def min_alpha_grid(j: int, denom: int = 600) -> Fr:
    """Least alpha = k/denom at which the chain is feasible (brute-force check of the closed form)."""
    for k in range(denom // 3, denom // 2 + 1):
        a = Fr(k, denom)
        if chain_feasible(a, j) is not None:
            return a
    raise ValueError


def children_exponent(j: int) -> Fr:
    """D-exponent of the number of large algebraic children of a level-i node at alpha_j:
    (delta + zeta_i) - zeta_{i+1} = delta + 3 alpha - 1 = alpha (the same at every level)."""
    a = alpha_j(j)
    w = chain_feasible(a, j)
    return w["delta"] + w["zetas"][0] - w["zetas"][1]


def two_level_parameters() -> Dict[str, Fr]:
    a = alpha_j(2)
    w = chain_feasible(a, 2)
    return {"alpha": a, "delta": w["delta"], "zeta1": w["zetas"][1], "level1_alg_value": w["delta"] + w["zetas"][0],
            "level2_value": max(1 - w["delta"], w["delta"] + w["zetas"][1]), "children": children_exponent(2)}


# ------------------------------------------------------------------ union-bound (Prop UB) bookkeeping
def union_bound_precision(j: int) -> List[Fr]:
    """D-exponent of the relative-variance precision eta_i needed at depth i = 1..j-1 by the union-bound route
    (Prop UB): eta_1 = O(1/m), eta_i = O(1/m^i) with m = D^{alpha_j} (typical count); returned as the exponent -i*alpha_j."""
    a = alpha_j(j)
    return [-i * a for i in range(1, j)]


def local_factor_covariance_exponent(j: int) -> Fr:
    """Heuristic D-exponent of the field-to-field (local-factor) relative covariance of smooth counts:
    (xi(u)/ln B)^2 with xi ~ alpha ln D, ln B ~ D^alpha ln t  ->  D^{-2 alpha} (ln-factors dropped).  HEURISTIC."""
    return -2 * alpha_j(j)


def union_bound_margin(j: int) -> List[Fr]:
    """Exponent margin (needed precision exponent) - (covariance exponent) at each depth; a NEGATIVE margin means
    the unconditional local-factor covariance alone exceeds the needed precision, so the hypothesis must be
    stated differently; local-data conditioning does not repair it (notes §3.5).  Zero = borderline (constants decide)."""
    c = local_factor_covariance_exponent(j)
    return [p - c for p in union_bound_precision(j)]


def split_loss_ok(alpha: Fr, kappa: Fr, C: Fr) -> bool:
    """A per-child loss (ln N)^{-C} with m = D^alpha children and a split of 2k = 2 kappa D^alpha primes costs
    exp(O((2 kappa + 1) C D^alpha ln t)) = L[alpha, O(1)]: always within budget (navigate §5.4 bookkeeping)."""
    return alpha > 0 and kappa > 0 and C >= 0
