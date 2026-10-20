"""Ledgers and constants for the rank-r volume theorem (agent rank3). Units as in resbarrier/heurnfs:
x = log N; D = delta (x/log x)^{1/3}; log B = beta t, entropy lambda = gamma t with t = x^{1/3} (log x)^{2/3};
log X = eta x^{2/3} (log x)^{1/3}. Constants are L_N[1/3, .] exponents (single N, RAM, omega = 2).

Corollary R3 (proved (model)): log X >= (r-1) x/(D-1) + (D/r)(lambda - log C_r)  =>  eta >= (r-1)/delta + delta gamma/r.
H (smoothness): gamma >= beta + eta/(3 beta). Cost >= max(gamma, 2 beta).
Scope H5_r (box error): r(D+1)(2E+1)^{r-1} <= B/4, i.e. (r-1) e <= beta, i.e. gamma <= r beta/(r-1).
Floating point is used only in the numeric cross-check `numeric_floor` (analysis only)."""
from __future__ import annotations

from fractions import Fraction as Fr

import sympy as sp


def unscoped_constant(r: int) -> sp.Expr:
    """Floor if the lattice-point error could be ignored: c^3 = 64 (r-1)/(9 r). r = 2: (32/9)^{1/3} = 1.52629 (SNFS);
    r = 3: (128/27)^{1/3} = 1.67989 (attained by SNFS-type 3-dim sieving in the size model); r -> oo: (64/9)^{1/3}."""
    return sp.cbrt(sp.Rational(64 * (r - 1), 9 * r))


def scoped_constant(r: int) -> sp.Expr:
    """Floor inside H5_r. r = 2: H5 (gamma <= 2 beta) is not binding: (32/9)^{1/3}. r >= 3: gamma = r beta/(r-1) binds:
    beta^{3/2} >= 2(r-1)/3, cost 2 beta = 2 (2(r-1)/3)^{2/3}; r = 3: 2.42285."""
    if r == 2:
        return unscoped_constant(2)
    return 2 * (sp.Rational(2 * (r - 1), 3)) ** sp.Rational(2, 3)


def numeric_floor(r: int, scoped: bool, grid: int = 4000) -> float:
    """Brute-force check of the two closed forms: minimise max(gamma, 2 beta) over beta, gamma, delta subject to
    gamma >= beta + ((r-1)/delta + delta gamma/r)/(3 beta) (optimal delta closed form) and, if scoped,
    gamma <= r beta/(r-1)."""
    best = float("inf")
    for i in range(1, grid):
        beta = 3.0 * i / grid
        gmax = r * beta / (r - 1) if scoped else 6.0
        # smallest feasible gamma: gamma - 2 sqrt((r-1) gamma / r)/(3 beta) >= beta; LHS increasing for gamma large
        lo, hi = beta, max(gmax, beta)
        def ok(g):
            return g - 2 * ((r - 1) * g / r) ** 0.5 / (3 * beta) >= beta
        if not ok(hi):
            continue
        for _ in range(60):
            mid = (lo + hi) / 2
            if ok(mid):
                hi = mid
            else:
                lo = mid
        best = min(best, max(hi, 2 * beta))
    return best


def liouville_height_floor(r: int, D: int) -> Fr:
    """Proposition R4 (proved): sum over the D embeddings of log max_j |L_i(u_j)| >= ((r-1)/(D-1)) x - O(D log r).
    In height units (log N): h_total >= (r-1)/(D-1)."""
    return Fr(r - 1, D - 1)


def gr_height_floor_lower_count(r: int, d: int) -> Fr:
    """resbarrier G_r (heuristic, Schmidt lower count, both sides degree d): (r-1)/(d+1)."""
    return Fr(r - 1, d + 1)


def aggregate_vs_liouville(x: float, r: int, D: int) -> dict:
    """Compare the aggregate ledger's norm size x/d + D log r (d = D/2) with the rigorous R4 floor (r-1)x/(D-1)."""
    import math
    return {"aggregate_logX": 2 * x / D + D * math.log(r), "R4_floor": (r - 1) * x / (D - 1) - D * math.log(12 * (r + 1))}
