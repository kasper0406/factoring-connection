"""Exact size-model bookkeeping for special-form polynomials (notes/snfsrig_results.md section 1).

Conventions (as in notes/recurse_results.md section 3.3 and notes/lv5unif_results.md):
  t = log N, D = t / log t, Lambda_alpha = t^alpha (log t)^(1-alpha), B = exp(a Lambda_alpha),
  d ~ D^x (degree), k ~ D^y (split: k primes of size N^(1/(2k)), nu = 1/(4k)).
A one-level design has norm-product exponent
  theta = theta_coeff + theta_box,   theta_box = (d - 1) nu  (box radius V ~ sqrt(l) S),
and the coefficient term theta_coeff = (log||f|| + log||g||)/log N.  For every pair with a LINEAR side g and
deg f = d, Res(f, g) is a nonzero multiple of N, so (d + 1)||f|| ||g||^d >= N and theta_coeff >= 1/d - o(1)
(Lemma RES in the notes; this module only does the exponent arithmetic).

Everything is exact (fractions.Fraction); no floating point is used for any decision.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------------------------------------
# 1. One-level floor as a function of the coefficient exponent p (theta_coeff ~ d^{-p})
# ---------------------------------------------------------------------------------------------------------


def one_level_floor(p) -> Fr:
    """Least alpha for which some (x, y) satisfies
         (R1) y <= alpha,  (R2/R3) p*x >= 1 - 2 alpha  and  y >= x + 1 - 2 alpha,  x >= 0.
    Closed form: (1 + p)/(3p + 2) for p > 0 (p = 1 gives 2/5; p -> infinity gives 1/3); p = 0 gives 1/2."""
    p = Fr(p)
    if p < 0:
        raise ValueError("p >= 0")
    if p == 0:
        return Fr(1, 2)
    return (1 + p) / (3 * p + 2)


def one_level_feasible(alpha, p, x, y) -> bool:
    """Check the three exponent constraints exactly at a given (alpha, x, y)."""
    alpha, p, x, y = Fr(alpha), Fr(p), Fr(x), Fr(y)
    return (x >= 0 and y <= alpha and p * x >= 1 - 2 * alpha and y >= x + 1 - 2 * alpha)


def one_level_optimum(p) -> Tuple[Fr, Fr, Fr]:
    """(alpha, x, y) at the floor: x = (1 - 2 alpha)/p, y = alpha."""
    al = one_level_floor(p)
    if Fr(p) == 0:
        return al, Fr(0), al
    return al, (1 - 2 * al) / Fr(p), al


# ---------------------------------------------------------------------------------------------------------
# 2. Constants at alpha = 2/5: theta * D^{1/5} -> theta0(d0, kappa), d = d0 D^{1/5}, k = kappa D^{2/5}
# ---------------------------------------------------------------------------------------------------------

DESIGNS = ("U", "GNFS_heur", "SNFS")


def theta0(design: str, d0, kappa) -> Fr:
    """Leading constant of theta D^{1/5}.
    U          : Theorem U's rigorous design (random field; floors h >= max(2nu, mu-(d-2)nu, 1-d mu+nu)):
                 3/d0 if kappa >= d0^2/4, else 2/d0 + d0/(4 kappa)   [lv5unif params.theta0]
    GNFS_heur  : base-m pair, ||f|| ~ ||g|| ~ N^{1/(d+1)}, no randomisation floors: 2/d0 + d0/(4 kappa)
    SNFS       : ||f|| = N^{o(1)}, ||g|| ~ N^{1/d}: 1/d0 + d0/(4 kappa)  (the RES minimum 1/d is attained)."""
    d0, kappa = Fr(d0), Fr(kappa)
    if design == "U":
        if kappa >= d0 * d0 / 4:
            return Fr(3) / d0
        return Fr(2) / d0 + d0 / (4 * kappa)
    if design == "GNFS_heur":
        return Fr(2) / d0 + d0 / (4 * kappa)
    if design == "SNFS":
        return Fr(1) / d0 + d0 / (4 * kappa)
    raise ValueError(design)


def theta0_min_squared(design: str, kappa) -> Fr:
    """(min over d0 of theta0)^2, exact.  U and GNFS_heur: 2/kappa (at d0 = sqrt(8 kappa), Kern regime, which
    is admissible for U since kappa < d0^2/4 = 2 kappa); SNFS: 1/kappa (at d0 = 2 sqrt(kappa))."""
    kappa = Fr(kappa)
    if design in ("U", "GNFS_heur"):
        return Fr(2) / kappa
    if design == "SNFS":
        return Fr(1) / kappa
    raise ValueError(design)


def two_pow_2_5_bracket(den: int = 10 ** 6) -> Tuple[Fr, Fr]:
    lo_n = 1319507 * den // 10 ** 6
    while Fr(lo_n + 1, den) ** 5 < 4:
        lo_n += 1
    while Fr(lo_n, den) ** 5 > 4:
        lo_n -= 1
    return Fr(lo_n, den), Fr(lo_n + 1, den)


def cost_constant(design: str, a, d0, kappa, C_gamma) -> Fr:
    """Upper bracket of Theorem U's cost formula with theta0 of the given design:
         c = max(2a, a + s + max(a, sigma), (a + s)/(1 - 2^{-2/5})),  s = 2 kappa (1 + C_gamma), sigma = 2 theta0/(5a).
    For SNFS and GNFS_heur this is HYPOTHETICAL: it assumes a relation supply gamma >= (log N)^{-C_gamma} that is
    proved only for design U (and, for SNFS, is exactly what is missing; notes section 2)."""
    a, kappa, C_gamma = Fr(a), Fr(kappa), Fr(C_gamma)
    sg = 2 * theta0(design, d0, kappa) / (5 * a)
    s = 2 * kappa * (1 + C_gamma)
    lo_t, _hi_t = two_pow_2_5_bracket()
    fhi = lo_t / (lo_t - 1)              # upper bracket of 1/(1 - 2^{-2/5})
    return max(2 * a, a + s + max(a, sg), (a + s) * fhi)


def optimise(design: str, C_gamma, grid_a=None, grid_kappa=None, grid_d0=None) -> Dict:
    grid_a = grid_a or [Fr(i, 20) for i in range(1, 81)]
    grid_kappa = grid_kappa or [Fr(i, 400) for i in range(1, 81)]
    best = None
    for kappa in grid_kappa:
        d0s = grid_d0 or [Fr(i, 20) for i in range(1, 300) if kappa <= Fr(i, 20) ** 2 <= 32 * kappa]
        for d0 in d0s:
            for a in grid_a:
                c = cost_constant(design, a, d0, kappa, C_gamma)
                if best is None or c < best["c"]:
                    best = {"design": design, "a": a, "kappa": kappa, "d0": d0, "c": c,
                            "theta0": theta0(design, d0, kappa)}
    return best


# ---------------------------------------------------------------------------------------------------------
# 3. Multi-level descent in one field (recurse section 3.4): exponent count of levels
# ---------------------------------------------------------------------------------------------------------


def multilevel_levels(alpha) -> int:
    """Levels needed from beta_1 = 1 - alpha down to alpha with decrement 3 alpha - 1 (alpha > 1/3).
    Identical for GNFS and SNFS heights: the coefficient term enters as (log N)^{1 - delta} either way."""
    alpha = Fr(alpha)
    if alpha <= Fr(1, 3):
        raise ValueError("bounded number of levels needs alpha > 1/3 (alpha = 1/3: O(log log N) levels)")
    num, den = 1 - 2 * alpha, 3 * alpha - 1
    q = num / den
    return int(q) if q.denominator == 1 else int(q) + 1


# ---------------------------------------------------------------------------------------------------------
# 4. Prop R2 (no algebraic children) recomputed: GNFS premise (M1) vs SNFS premise (M1-S)
# ---------------------------------------------------------------------------------------------------------


def r2_window_nonempty_exponent(alpha) -> bool:
    """GNFS (M1): nested Prop R2's window (A, W) is non-empty iff Lambda^2 U^3 = o(t^2), Lambda ~ U ~ t^alpha up to
    log factors, i.e. iff 5 alpha < 2.  Non-empty window => no finite cost-L[alpha] tree without algebraic children."""
    return 5 * Fr(alpha) < 2


def r2s_contraction_exponent(alpha) -> Fr:
    """SNFS (M1-S): at a node with rational target e^z choose d = c' Xi/(z/2 + log S), Xi = Lambda U (algebraic value
    B-smooth within budget).  Rational cofactor ~ t/d ~ t z/(2 c' Xi); children y = cofactor/U, so
        y/z  ~  t/(2 c' Lambda U^2)  =  t^{1 - 3 alpha} (log t)^{3 alpha - 1} / const.
    Returns the t-exponent 1 - 3 alpha of the contraction ratio (< 0: polynomial contraction; = 0: constant)."""
    return 1 - 3 * Fr(alpha)


def r2s_floor() -> Fr:
    """Least alpha with a non-expanding rational-only descent in the SNFS model: 1/3 (attained with constant
    contraction, see r2s_alpha_third)."""
    return Fr(1, 3)


def r2s_alpha_third(b, s, phiF, phiR, b1) -> Dict:
    """alpha = 1/3, E = t^{1/3}(log t)^{2/3}, Q = t^{2/3}(log t)^{1/3}, D' = Q/E, t = D'Q.
    B = e^{bE}; box radius e^{sE} sqrt(q); algebraic smoothness u_F = phiF D' (cost e^{phiF E/3});
    rational-cofactor smoothness u_R = phiR D' (cost e^{phiR E/3}); split targets z_1 = b1 Q.
    Node with target z >= bE: d(z) = phiF b Q/(z/2 + sE) (so the algebraic value e^{d(z/2+sE)} = e^{phiF b Q} is
    B-smooth with u = phiF D'); rational cofactor t/d = D'(z/2 + sE)/(phiF b); children y = (z/2 + sE)/(phiF phiR b).
    Returns the exact contraction bound lam (y <= lam z for all z >= bE), the top degree d1 = 2 phiF b/b1, the
    heuristic per-target good-point exponent 2s - (phiF + phiR)/3 (must be > 0 for an all-target supply under the
    random-integer model), and the top-level rational-value exponent tau = 1/d1 (values N^{tau})."""
    b, s, phiF, phiR, b1 = map(Fr, (b, s, phiF, phiR, b1))
    lam = (Fr(1, 2) + s / b) / (phiF * phiR * b)
    d1 = 2 * phiF * b / b1
    margin = 2 * s - (phiF + phiR) / 3
    return {"lam": lam, "contracts": lam < 1, "d1": d1, "d1_ge_1": d1 >= 1,
            "good_point_margin": margin, "heuristic_ok": margin > 0,
            "tau_top": (1 / d1) if d1 >= 1 else None,
            "max_degree_over_Dprime": 2 * phiF}


# ---------------------------------------------------------------------------------------------------------
# 5. The size-model table (data; derivations in the notes)
# ---------------------------------------------------------------------------------------------------------

TABLE: List[Tuple[str, str, str, str]] = [
    # (design, one-level alpha, multi-level one-field alpha (algebraic children), rational-only alpha)
    ("GNFS-type (M1 holds; random fields available)", "2/5", "1/3 (endogenous algebraic children)", "2/5 (nested Prop R2)"),
    ("SNFS (||f|| = N^{o(1)}, ||g|| = N^{1/d}; fixed fields)", "2/5 (theta0 smaller by sqrt 2 vs heuristic GNFS)",
     "1/3", "1/3 (Prop R2-S, per-node degree)"),
]
