"""Exact cost ledgers for the intermediate-k regimes (agent midk).

Units: x = log N. Heights h and the exponents a, b are in units of x (resp. of log X ~ x^a C^b with
C = log sieve region at the working L[alpha] scale, C = x^alpha). Constants are L_N[alpha, .].

Three families (see notes/midk_prereg.md):
  A  rank-r linked modules (G_r regime)     — parameter r,   fallback alpha 1/4, under G_r 1/3
  B  descent ladder                          — parameter j,   alpha_j = (j+1)/(3j+2), inf 1/3
  C  floor-sum ladder                        — deterministic N^{1/(j+2)} metric (not L[1/3])

The reduction log X >~ x^a C^b  =>  alpha = a/(2-b) is resbarrier.ledger.alpha_floor (reused).
"""
from __future__ import annotations

from fractions import Fraction as Fr

import sympy as sp

from src.resbarrier.ledger import alpha_floor, gr_min, kl_min

x, delta, c, kappa, r = sp.symbols("x delta c kappa r", positive=True)


# ---------------------------------------------------------------------------
# Family A: the rank-r aggregate linked-module ledger.
# A candidate Sum a_j (w_1j, w_2j), a_j in {0,+-1}, has per-side log-norm
#     log X = (height floor) + (degree term) = kappa*x/d + d*c,
# with degree-per-side d, c = log(E*sqrt(r)) >= 1 the per-degree box factor, and
#     kappa = 1       if the r-1 linking conditions are DEPENDENT (G_r fails favourably),
#             r - 1   if they are INDEPENDENT (G_r holds): each costs one height unit.
# The rigorous 2-plane resultant gives kappa >= 1 unconditionally (Bean/Thunder area on a
# rank-2 subplane), so the height floor x/d is rigorous.
# ---------------------------------------------------------------------------

def min_logX_over_d(kappa_val=1, c_val=c):
    """Exact min over d>0 of kappa*x/d + d*c.  = 2*sqrt(kappa*c*x) at d = sqrt(kappa*x/c)."""
    expr = kappa_val * x / delta + delta * c_val  # delta plays the role of d here
    dstar = sp.solve(sp.diff(expr, delta), delta)
    dstar = [s for s in dstar if s.is_positive][0]
    return sp.simplify(expr.subs(delta, dstar)), sp.simplify(dstar)


def amgm_floor_exponents(kappa_val=1):
    """log X_min = 2 sqrt(kappa c x) = x^{1/2} * (2 sqrt(kappa c)); so (a, b) = (1/2, 0) in x^a C^b
    because c = log(E sqrt r) is only poly-logarithmic in the number of candidates C.  Hence the
    fallback floor alpha = alpha_floor(1/2, 0) = 1/4, and under G_r (kappa = r-1 ~ C, so the box
    entropy enters log X with b = 1/2) alpha = alpha_floor(1/2, 1/2) = 1/3."""
    logX, dstar = min_logX_over_d(kappa_val)
    return {"logX_min": logX, "d_star": dstar}


def family_A_alpha(regime: str) -> Fr:
    """alpha for Family A. 'fallback' (G_r fails, kappa=1, b=0) -> 1/4; 'Gr' (kappa=r-1 ~ C, b=1/2) -> 1/3."""
    if regime == "fallback":
        return alpha_floor(Fr(1, 2), Fr(0))
    if regime == "Gr":
        return alpha_floor(Fr(1, 2), Fr(1, 2))
    raise ValueError(regime)


def aggregate_optimum_scaling():
    """The attaining L[1/4] scaling (resbarrier sec 4.2, codex resbarrier_2 #7): exact exponents of
    (log N) and (log log N). d = sqrt(x/log x), r ~ x^{1/4}(log x)^{3/4}, log X ~ sqrt(x log x),
    B = L[1/4], u = log X / log B ~ x^{1/4}. Returned as a dict of (exp_x, exp_loglog) pairs."""
    # exponents written as (power of log N, power of log log N)
    return {
        "d": (Fr(1, 2), Fr(-1, 2)),          # sqrt(x/log x)
        "r": (Fr(1, 4), Fr(3, 4)),           # x^{1/4}(log x)^{3/4}
        "logX": (Fr(1, 2), Fr(1, 2)),        # sqrt(x log x)
        "logB": (Fr(1, 4), Fr(3, 4)),        # L[1/4]
        "u": (Fr(1, 4), Fr(-1, 4)),          # x^{1/4}(log x)^{-1/4}  ~ x^{1/4}
        "alpha": Fr(1, 4),
    }


def gr_floor_monotone(d: int, rmax: int = 8):
    """Under G_r the height floor gr_min(r,d,d) is non-decreasing in r (so r=2 is optimal, alpha>=1/3).
    Returns the list of (r, gr_min) for 2<=r<=min(rmax,d)."""
    out = []
    for rr in range(2, min(rmax, d) + 1):
        out.append((rr, gr_min(rr, d, d)))
    return out


def kl_floor_monotone(d: int, rmax: int = 8):
    """With a rational side, the rigorous KL floor min(h1+h2) is non-decreasing in r."""
    out = []
    for rr in range(2, min(rmax, d) + 1):
        out.append((rr, kl_min(rr, d)))
    return out


# ---------------------------------------------------------------------------
# Family B: the descent ladder alpha_j = (j+1)/(3j+2).
# ---------------------------------------------------------------------------

def ladder_alpha(j: int) -> Fr:
    """multidescent/pc3: cost exponent of a j-level descent, conditional on C1/C2 (pair correlation)."""
    return Fr(j + 1, 3 * j + 2)


def ladder_table(jmax: int = 8):
    return [(j, ladder_alpha(j)) for j in range(1, jmax + 1)]


# ---------------------------------------------------------------------------
# Family C: the floor-sum ladder, deterministic N^{1/(j+2)} metric.
# ---------------------------------------------------------------------------

def floorsum_exponent(j: int) -> Fr:
    """Deterministic cost exponent N^{1/(j+2)} of the degree-j floor-sum ladder (floorhard/topblock)."""
    return Fr(1, j + 2)
