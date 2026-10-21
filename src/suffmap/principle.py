"""Smoothness-size principle and size-entropy laws (agent suffmap; notes/suffmap_results.md, section 1).

Everything here is exact (sympy / fractions). No floating point is used in any check.

Units. x = log N, F_a = x^a (log x)^{1-a}, L[a, c] = exp((c + o(1)) F_a).

HEURISTIC H: an examined candidate whose smoothness integers have total log-size lam is B-smooth with probability
rho(lam / log B) * L^{o(1)}, with -log rho(u) = (1 + o(1)) u log u; tests and generation cost L^{o(1)} each;
sparse linear algebra on B columns costs B^omega (omega = 2); dependencies are fruitful.

Theorem SS (sampling form; notes section 1.2).
  (a) Fixed size lam: cost exp((1 + o(1)) sqrt(2 lam log lam)) = L(X)^{sqrt 2} with X = e^lam (Pomerance).
  (b) In L-units: lam = s F_{2 alpha}  ==>  cost L[alpha, 2 sqrt(alpha s)], with L[alpha, 2 sqrt(alpha s)]
      candidates needed and B = L[alpha, sqrt(alpha s)].
  (c) Size-entropy law lam(C) (smallest achievable total log-size when e^C candidates are examined):
      cost = e^{C*}, C* = least C with C^2 >= 2 lam(C) log lam(C), for every omega <= 2.
      For lam = A x^a (log x)^{1-a-b} C^b (0 <= b < 2): alpha = a / (2 - b) and c^{2-b} = 4 alpha A.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, Tuple

import sympy as sp

X, LX = sp.symbols("x ell", positive=True)  # x = log N, ell = log x (treated as independent large symbol)


# ----------------------------------------------------------------------------- Theorem SS (b): fixed size
def sigma_needed(alpha) -> Fr:
    """Candidate log-size exponent sigma (lam = s F_sigma) needed for an L[alpha] sampling method: sigma = 2 alpha."""
    return 2 * Fr(alpha)


def c_sampling(alpha, s) -> sp.Expr:
    """Constant of L[alpha, c] for candidates of total log-size s F_{2 alpha}: c = 2 sqrt(alpha s)."""
    return sp.nsimplify(2) * sp.sqrt(sp.Rational(Fr(alpha).numerator, Fr(alpha).denominator)
                                     * sp.Rational(Fr(s).numerator, Fr(s).denominator))


def s_needed(alpha, c) -> sp.Expr:
    """Size coefficient s allowed for a target L[alpha, c]: s = c^2 / (4 alpha)."""
    a = sp.Rational(Fr(alpha).numerator, Fr(alpha).denominator)
    return sp.nsimplify(c) ** 2 / (4 * a)


def leading_coefficients(alpha, s, b) -> Dict[str, sp.Expr]:
    """Exact leading-order L[alpha]-coefficients for factor base B = L[alpha, b] and candidates s F_{2 alpha}.

    log B = b F_a; u = s F_{2a} / (b F_a) = (s/b) x^a (log x)^{-a}; log u = (a + o(1)) log x, so
    u log u = (a s / b) F_a. Trials = B / rho: coefficient b + a s / b. Linear algebra: 2 b.
    """
    a = sp.nsimplify(alpha)
    s = sp.nsimplify(s)
    b = sp.nsimplify(b)
    return {"trials": b + a * s / b, "linear_algebra": 2 * b, "rho_penalty": a * s / b}


def optimum_fixed_size(alpha, s) -> Dict[str, sp.Expr]:
    """Exact optimum of max(b + alpha s / b, 2 b): b* = sqrt(alpha s), both terms equal 2 sqrt(alpha s)."""
    a = sp.nsimplify(alpha)
    s = sp.nsimplify(s)
    bb = sp.Symbol("b", positive=True)
    trials = bb + a * s / bb
    bstar = sp.solve(sp.Eq(sp.diff(trials, bb), 0), bb)
    bstar = [r for r in bstar if r.is_positive][0]
    t = sp.simplify(trials.subs(bb, bstar))
    la = sp.simplify(2 * bstar)
    return {"b": sp.simplify(bstar), "trials": t, "linear_algebra": la, "c": sp.simplify(sp.Max(t, la))}


# ----------------------------------------------------------------------------- Theorem SS (c): size-entropy law
def exponent_from_law(a, b) -> Fr:
    """alpha for lam(C) = A x^a (log x)^{1-a-b} C^b: alpha = a / (2 - b)."""
    a, b = Fr(a), Fr(b)
    if not (0 <= b < 2):
        raise ValueError("need 0 <= b < 2")
    return a / (2 - b)


def log_power_of_law(a, b) -> Fr:
    """The (log x)-power g that makes the law L-homogeneous: g = 1 - a - b."""
    return 1 - Fr(a) - Fr(b)


def constant_from_law(A, a, b) -> sp.Expr:
    """c with c^{2-b} = 4 alpha A, alpha = a/(2-b): the L[alpha, c] constant of the law A x^a (log x)^g C^b."""
    al = exponent_from_law(a, b)
    alr = sp.Rational(al.numerator, al.denominator)
    bb = sp.Rational(Fr(b).numerator, Fr(b).denominator)
    return sp.simplify((4 * alr * sp.nsimplify(A)) ** (1 / (2 - bb)))


def verify_law_balance(A, a, b) -> sp.Expr:
    """Leading-order check that C = c F_alpha solves C^2 = 2 lam(C) log lam(C) for the law (exact, symbolic).

    Substitute x = X, log x = ell, keep log lam = (a + alpha b) log x (its leading term). Returns the ratio
    C^2 / (2 lam * (a + alpha b) * ell), which must simplify to 1.
    """
    al = exponent_from_law(a, b)
    alr = sp.Rational(al.numerator, al.denominator)
    ar = sp.Rational(Fr(a).numerator, Fr(a).denominator)
    br = sp.Rational(Fr(b).numerator, Fr(b).denominator)
    g = 1 - ar - br
    c = constant_from_law(A, a, b)
    Fa = X ** alr * LX ** (1 - alr)
    Cv = c * Fa
    lam = sp.nsimplify(A) * X ** ar * LX ** g * Cv ** br
    loglam_leading = (ar + alr * br) * LX
    return sp.simplify(sp.powsimp(Cv ** 2 / (2 * lam * loglam_leading), force=True))


# Catalogue of size-entropy laws. (A, a, b) with lam = A x^a (log x)^{1-a-b} C^b. Ledgers from
# notes/heurnfs_results.md (B1, B2), notes/resbarrier_results.md (4.2) and the pieces cited in the note.
LAWS: Dict[str, Tuple[sp.Expr, Fr, Fr]] = {
    # Dixon: random x^2 mod N, size N, independent of C.
    "Dixon": (sp.Integer(1), Fr(1), Fr(0)),
    # QS: (t + m)^2 - N ~ 2 t sqrt(N): size sqrt(N) (t = L[1/2] adds o(x)).
    "QS": (sp.Rational(1, 2), Fr(1), Fr(0)),
    # GNFS / base-m, any two-polynomial selection under G: log X >= 2x/D + D C/2 >= 2 sqrt(x C).
    "GNFS": (sp.Integer(2), Fr(1, 2), Fr(1, 2)),
    # SNFS (kappa = 1): log X >= x/D + D C/2 >= sqrt(2 x C).
    "SNFS": (sp.sqrt(2), Fr(1, 2), Fr(1, 2)),
    # sieving dimension t = 3 (heurnfs B2): 2(t-1)x/D + D C/t >= 2 sqrt(2(t-1) x C / t).
    "GNFS t=3": (2 * sp.sqrt(sp.Rational(4, 3)), Fr(1, 2), Fr(1, 2)),
    # hypothetical aggregate rank-r family (resbarrier 4.2, stylised): lam ~ A sqrt(x log x), A = 1 here.
    "aggregate rank-r (stylised, A=1)": (sp.Integer(1), Fr(1, 2), Fr(0)),
}


def law_table() -> Dict[str, Dict[str, object]]:
    out = {}
    for name, (A, a, b) in LAWS.items():
        out[name] = {"alpha": exponent_from_law(a, b), "c": constant_from_law(A, a, b),
                     "balance_ratio": verify_law_balance(A, a, b)}
    return out


# ----------------------------------------------------------------------------- targets (sigma-alpha table)
TARGETS = (Fr(1, 3), Fr(1, 4), Fr(1, 5), Fr(1, 10))


def sigma_alpha_table() -> Dict[Fr, Dict[str, object]]:
    """For each target alpha: sigma = 2 alpha, s allowed for c = 1 and for c = 1.9019-level (s = c^2/(4 alpha)),
    and the zero-entropy law exponent a (b = 0) and the b = 1/2 law exponent a = 3 alpha / 2."""
    out = {}
    for al in TARGETS:
        out[al] = {
            "sigma": sigma_needed(al),
            "s_for_c=1": s_needed(al, 1),
            "law_a_if_b=0": 2 * al,
            "law_a_if_b=1/2": Fr(3, 2) * al,
            "family_size_exponent": al,
        }
    return out


def polytime_size_bound() -> str:
    """Sampling reaches polynomial time iff C* = O(log x) iff lam = O((log x)^2 / log log x) (proved in section 1.2)."""
    return "e^{C*} <= x^{K(1+o(1))}  <=>  lam(C*) <= (1+o(1)) K^2 (log x)^2 / (4 log log x)  (poly(x) generation/tests)"


def quasi_poly_check(k) -> Fr:
    """For lam = (log x)^k: C* ~ sqrt(2 k) (log x)^{k/2} (log log x)^{1/2}; returns the (log x)-exponent k/2."""
    return Fr(k, 2)


# ----------------------------------------------------------------------------- listing regime (pilist Table 1.1(A))
def one_sided_cost(alpha, kappa, nu, tau) -> sp.Expr:
    """Congruence listing, one-sided: factor base L[alpha, kappa], window log-size nu F_{2 alpha}, lister cost
    L[alpha, tau] per output. Exponent coefficient max(kappa + tau + alpha nu / kappa, 2 kappa)."""
    a, k, n, t = (sp.nsimplify(v) for v in (alpha, kappa, nu, tau))
    return sp.Max(k + t + a * n / k, 2 * k)


def one_sided_optimum_free_lister(alpha, nu) -> sp.Expr:
    """tau = 0: optimum 2 sqrt(alpha nu), identical to Theorem SS(b) with s = nu (only the residue side counts)."""
    return optimum_fixed_size(alpha, nu)["c"]


def doubly_cost(kappa, tau) -> sp.Expr:
    """Doubly-split listing: max(kappa + tau, 2 kappa); infimum tau as kappa -> 0 (any alpha)."""
    k, t = sp.nsimplify(kappa), sp.nsimplify(tau)
    return sp.Max(k + t, 2 * k)


def single_interval_program() -> Dict[str, sp.Expr]:
    """Prop SI: smooth u constructed in single intervals [kN, kN + h], h = L[2/3, nu], y = L[1/3, kappa].

    Existence in one interval needs nu >= 2/(3 kappa) (u ~ N is y-smooth with probability L[2/3, -2/(3 kappa)]).
    The residue v = u - kN in [0, h] is then y-smooth with probability L[1/3, -nu/(3 kappa)]. Cost coefficient
    max(kappa + nu/(3 kappa), 2 kappa), minimised at nu = 2/(3 kappa), kappa^3 = 2/9: c = (16/9)^{1/3}.
    """
    k = sp.Symbol("kappa", positive=True)
    trials = k + sp.Rational(2, 9) / k ** 2
    kstar = sp.solve(sp.Eq(trials, 2 * k), k)
    kstar = [r for r in kstar if r.is_positive][0]
    # check the intersection lies left of the minimiser of `trials` (so it is the optimum of the max)
    kmin = [r for r in sp.solve(sp.Eq(sp.diff(trials, k), 0), k) if r.is_positive][0]
    return {"kappa": sp.simplify(kstar), "nu": sp.simplify(sp.Rational(2, 3) / kstar),
            "c": sp.simplify(2 * kstar), "c_cubed": sp.simplify((2 * kstar) ** 3),
            "intersection_left_of_min": bool(sp.simplify(kmin - kstar) > 0)}


def single_interval_exponent(alpha) -> Fr:
    """Single intervals at smoothness L[alpha]: the interval must have length >= L[1 - alpha] (one smooth u), so the
    residue has log-size ~ F_{1-alpha}, smooth at L[alpha] with probability L[1 - 2 alpha]^{-1}: exponent
    max(alpha, 1 - 2 alpha), minimised at alpha = 1/3."""
    a = Fr(alpha)
    return max(a, 1 - 2 * a)


def congruence_box_exponent(alpha) -> Fr:
    """Existence: an L[alpha]-doubly-smooth congruence family needs box (U/N) of size L[1 - alpha] (pilist 1.2)."""
    return 1 - Fr(alpha)


def nfs_doubly_box_exponent(alpha) -> Fr:
    """Existence for NFS-shaped doubly listing at factor base L[alpha]: log E ~ x^gamma with gamma >= 1 - 2 alpha.

    Norm log-size >= 2 x^{(1+gamma)/2}; u log u ~ x^{(1+gamma)/2 - alpha} log x <= 2 x^gamma needs
    (1 + gamma)/2 - alpha <= gamma.
    """
    return 1 - 2 * Fr(alpha)


if __name__ == "__main__":  # pragma: no cover
    import pprint
    print("fixed-size optimum, alpha=1/4, s=1:", optimum_fixed_size(Fr(1, 4), 1))
    pprint.pprint(law_table())
    pprint.pprint(sigma_alpha_table())
