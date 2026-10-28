"""Exact helpers for notes/fruit_results.md (no floating point).

1. Descent exponent ladders (heuristic smoothness model, leading order in the exponent):
   - pure descent (free degree D^gamma, factor base L[alpha]):     alpha_k = (k+1)/(3k+2)
   - hybrid (descend INTO an LV field: degree ~ D^{1/3}, factor base L[1/3]):
                                                                 alpha_k = (k+2)/(3k+3)
   `chain_feasible` checks the constraint chains directly.
2. Prop IV-LV counting: for fixed coprime box pairs (a,b), (a',b') (non-proportional, b b' != 0),
   the number of m in [L1, L2) with (a-bm)(a'-b'm) a nonzero perfect square is at most
   2 * sum_{d | Delta, d squarefree, d > 0} floor(sqrt((|a| + |b| L2)/d)),  Delta = ab' - a'b.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import isqrt

import sympy


def exponent_pure(k: int) -> Fr:
    return Fr(k + 1, 3 * k + 2)


def exponent_hybrid(k: int) -> Fr:
    return Fr(k + 2, 3 * k + 3)


def chain_feasible(alpha: Fr, k: int, mode: str) -> bool:
    """Leading-order feasibility of a k-level descent at cost L[alpha].
    Exponents: split primes L[beta_1] need 1 - beta_1 <= alpha; a special-q step from L[beta_i]
    in a degree-D^gamma field has algebraic log-size D^{gamma+beta_i}, rational log-size
    D^{1-gamma}; children L[beta_{i+1}]-smooth cost L[max(gamma+beta_i, 1-gamma) - beta_{i+1}].
    The last children must be smooth over the factor base L[phi].
    mode 'pure': gamma and phi free (phi = alpha is optimal: LA costs L[2 phi]-scale... at the
    exponent level the factor base may be L[alpha]); mode 'hybrid': gamma = phi = 1/3."""
    if mode == "hybrid":
        gammas = [Fr(1, 3)]
        phi = Fr(1, 3)
    else:
        phi = alpha
        # gamma enters through 1-gamma-beta <= alpha (rational child) and gamma+beta-beta' <= alpha;
        # smallest admissible gamma is optimal: gamma = max(0, 1 - 2 alpha)  (last-level rational side)
        gammas = [max(Fr(0), 1 - 2 * alpha)]
    for gamma in gammas:
        # backward: the largest beta_k allowed so that the last level lands in L[phi]
        bmax = alpha + phi - gamma          # gamma + beta_k - phi <= alpha
        if 1 - gamma - phi > alpha:         # last rational child must be L[phi]-smooth
            continue
        for _ in range(k - 1):
            # gamma + beta_i - beta_{i+1} <= alpha  with beta_{i+1} <= bmax
            bmax = alpha + bmax - gamma
        if 1 - alpha <= bmax:              # split: beta_1 >= 1 - alpha
            return True
    return False


def squarefree_divisors(n: int):
    n = abs(n)
    ps = list(sympy.factorint(n).keys()) if n > 1 else []
    out = [1]
    for p in ps:
        out += [d * p for d in out]
    return out


def rational_square_bound(a, b, a2, b2, L1, L2) -> int:
    Delta = a * b2 - a2 * b
    assert Delta != 0 and b != 0 and b2 != 0
    top = abs(a) + abs(b) * L2
    return 2 * sum(isqrt(top // d) for d in squarefree_divisors(Delta))


def count_rational_square_m(a, b, a2, b2, L1, L2) -> int:
    cnt = 0
    for m in range(L1, L2):
        P = (a - b * m) * (a2 - b2 * m)
        if P > 0 and isqrt(P) ** 2 == P:
            cnt += 1
    return cnt
