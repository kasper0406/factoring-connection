"""sicheck model: exact (sympy) re-derivation of the Prop SI cost program and its consistency checks.

Units: x = log N, F_a = x^a (log x)^{1-a}, L[a, c] = exp((c + o(1)) F_a). Smoothness model (CEP): an integer of
log-size s F_b is L[a, k]-smooth with probability L[b - a, -(b - a) s / k] (u = s F_b / (k F_a), log u ~ (b - a) log x).
Everything returned is an exact sympy expression or Fraction; no floating point enters a check.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict

import sympy as sp

kap, nu, tau, al = sp.symbols("kappa nu tau alpha", positive=True)


def smooth_penalty(size_level, size_coef, base_level, base_coef):
    """CEP in L-units: an integer of log-size size_coef * F_{size_level} is L[base_level, base_coef]-smooth with
    probability L[size_level - base_level, -(size_level - base_level) * size_coef / base_coef].
    Returns (level, coefficient) of the probability (coefficient negative)."""
    lev = sp.nsimplify(size_level) - sp.nsimplify(base_level)
    return lev, -lev * sp.nsimplify(size_coef) / sp.nsimplify(base_coef)


def si_program() -> Dict[str, sp.Expr]:
    """Independent re-derivation of Prop SI.

    - u ~ kN has log-size x = 1 * F_1: smooth with probability L[2/3, -2/(3 kappa)]  -> existence threshold nu_0.
    - v <= h = L[2/3, nu]: log-size nu F_{2/3}: smooth with probability L[1/3, -nu/(3 kappa)].
    - relations needed: pi(y) = L[1/3, kappa]; oracle calls: L[1/3, kappa + nu/(3 kappa) + tau]; LA: L[1/3, 2 kappa].
    """
    lev_u, c_u = smooth_penalty(1, 1, sp.Rational(1, 3), kap)
    lev_v, c_v = smooth_penalty(sp.Rational(2, 3), nu, sp.Rational(1, 3), kap)
    assert lev_u == sp.Rational(2, 3) and lev_v == sp.Rational(1, 3)
    nu0 = sp.simplify(-c_u)                               # 2/(3 kappa)
    trials = kap - c_v + tau                              # kappa + nu/(3 kappa) + tau
    la = 2 * kap
    t0 = sp.simplify(trials.subs({nu: nu0, tau: 0}))      # kappa + 2/(9 kappa^2)
    ks = [r for r in sp.solve(sp.Eq(t0, la), kap) if r.is_positive][0]
    kmin = [r for r in sp.solve(sp.Eq(sp.diff(t0, kap), 0), kap) if r.is_positive][0]
    c = sp.simplify(2 * ks)
    return {"nu_threshold": nu0, "trials": sp.simplify(trials), "trials_at_threshold": t0,
            "kappa_star": sp.simplify(ks), "kappa_star_cubed": sp.simplify(ks ** 3),
            "kappa_min_trials": sp.simplify(kmin), "kappa_min_cubed": sp.simplify(kmin ** 3),
            "c": c, "c_cubed": sp.simplify(c ** 3), "nu_star": sp.simplify(nu0.subs(kap, ks)),
            "d_trials_d_nu": sp.diff(trials, nu)}


def si_constant_with_k_exponent(theta) -> sp.Expr:
    """Multipliers k ~ N^theta: u has log-size (1 + theta) x, threshold 2(1 + theta)/(3 kappa), crossing
    kappa^3 = 2(1 + theta)/9, so c^3 = 16(1 + theta)/9 (codex sicheck_1 section 2; re-derived here).
    Returns c^3."""
    th = sp.nsimplify(theta)
    lev_u, c_u = smooth_penalty(1, 1 + th, sp.Rational(1, 3), kap)
    t0 = kap + (-c_u) / (3 * kap)
    ks = [r for r in sp.solve(sp.Eq(t0, 2 * kap), kap) if r.is_positive][0]
    return sp.simplify((2 * ks) ** 3)


def si_optimum_is_global() -> bool:
    """max(kappa + nu/(3 kappa), 2 kappa) over nu >= 2/(3 kappa): trials increase in nu, so nu = threshold; then
    max(t0, 2 kappa) is minimised at the crossing because t0 is decreasing there (kappa* < kappa_min) and 2 kappa is
    increasing. Checked exactly: d trials / d nu = 1/(3 kappa) > 0, kappa*^3 = 2/9 < kappa_min^3 = 4/9."""
    P = si_program()
    return bool(sp.simplify(P["d_trials_d_nu"] - 1 / (3 * kap)) == 0
                and P["kappa_star_cubed"] == sp.Rational(2, 9) and P["kappa_min_cubed"] == sp.Rational(4, 9)
                and P["c_cubed"] == sp.Rational(16, 9))


def size_principle_consistency() -> Dict[str, sp.Expr]:
    """Theorem SS(b) (Pomerance's principle) with alpha = 1/3: random part of log-size s F_{2/3} gives
    c = 2 sqrt(s/3) at factor base L[1/3, sqrt(s/3)]. For Prop SI the random part is v with s = nu* = 2/(3 kappa*);
    SS(b)'s optimal base sqrt(s/3) must equal kappa* and its c must equal (16/9)^{1/3}."""
    P = si_program()
    s = P["nu_star"]
    b_opt = sp.sqrt(s / 3)
    c_ss = 2 * b_opt
    return {"s": s, "b_opt_minus_kappa_star": sp.simplify(b_opt - P["kappa_star"]),
            "c_ss_cubed": sp.simplify(c_ss ** 3)}


def s_of_constant(c) -> sp.Expr:
    """Random log-size coefficient s (in F_{2/3}) that a single-family alpha = 1/3 method with constant c can afford:
    s = 3 c^2 / 4."""
    return sp.Rational(3, 4) * sp.nsimplify(c) ** 2


def comparison_table() -> Dict[str, Dict[str, sp.Expr]]:
    """Constants and their 'random log-size budgets' s = 3c^2/4 (exact)."""
    consts = {
        "Prop SI": sp.Rational(16, 9) ** sp.Rational(1, 3),
        "Kirchner special / resbarrier k=1": sp.Rational(16, 9) ** sp.Rational(1, 3),
        "SNFS": sp.Rational(32, 9) ** sp.Rational(1, 3),
        "factory per N (Coppersmith 1993)": None,   # 1.6386 is not a closed form here; numeric only
        "GNFS": sp.Rational(64, 9) ** sp.Rational(1, 3),
    }
    out = {}
    for k, c in consts.items():
        if c is None:
            continue
        out[k] = {"c": c, "c_cubed": sp.simplify(c ** 3), "s": sp.simplify(s_of_constant(c))}
    return out


def single_interval_exponent_program():
    """Claim (iii) re-derived at smoothness level alpha (y = L[alpha, kappa]).

    Smooth density near kN: L[1 - alpha, -(1 - alpha)/kappa], so an interval must have log-length >= (1-alpha)/kappa
    * F_{1-alpha} (threshold). The residue then has log-size (1-alpha)/kappa F_{1-alpha} and is smooth with probability
    L[1 - 2 alpha, -(1 - 2 alpha)(1 - alpha)/kappa^2]. Calls: L[alpha] * L[1 - 2 alpha]; LA L[alpha].
    Exponent max(alpha, 1 - 2 alpha) for alpha <= 1/2."""
    lev_u, c_u = smooth_penalty(1, 1, al, kap)
    nu_thr = -c_u
    lev_v, c_v = smooth_penalty(lev_u, nu_thr, al, kap)
    return {"interval_level": sp.simplify(lev_u), "nu_threshold": sp.simplify(nu_thr),
            "residue_penalty_level": sp.simplify(lev_v), "residue_penalty_coef": sp.simplify(c_v)}


def single_interval_exponent(alpha: Fr) -> Fr:
    a = Fr(alpha)
    return max(a, 1 - 2 * a)


def below_threshold_never_helps() -> sp.Expr:
    """Intervals of log-length H below the threshold log(1/delta): calls per relation ~ (1/(h delta)) / rho_v. In
    L-units at alpha = 1/3 with h = L[2/3, nu], nu < nu0: extra factor L[2/3, nu0 - nu] for empty intervals, while the
    residue gain is only L[1/3, (nu0 - nu)/(3 kappa)]. Returns the level-2/3 coefficient of the net change, which is
    nu0 - nu > 0 (so lowering nu below nu0 costs an L[2/3] factor)."""
    nu0 = 2 / (3 * kap)
    return sp.simplify(nu0 - nu)


def oracle_cost_threshold() -> Dict[str, sp.Expr]:
    """With oracle cost L[1/3, tau] per call (and nu at threshold), the program is max(kappa + 2/(9kappa^2) + tau,
    2 kappa). Return the tau at which the optimised constant equals MNFS's 1.9019 is suffmap's 0.757 (numeric check
    only, outside the exact tests)."""
    t0 = kap + sp.Rational(2, 9) / kap ** 2
    kmin = sp.Rational(4, 9) ** sp.Rational(1, 3)
    return {"min_trials_no_tau": sp.simplify(t0.subs(kap, kmin)), "kappa_min": kmin}


if __name__ == "__main__":  # pragma: no cover
    import pprint
    pprint.pprint(si_program())
    print("optimum global:", si_optimum_is_global())
    pprint.pprint(size_principle_consistency())
    pprint.pprint({k: {kk: (vv, sp.N(vv, 6)) for kk, vv in v.items()} for k, v in comparison_table().items()})
    pprint.pprint(single_interval_exponent_program())
    pprint.pprint(oracle_cost_threshold())
