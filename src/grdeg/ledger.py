"""Exact exponent-level ledger for the resbarrier open item #1 (degeneracy-driven G_r failure).

DIAGNOSTIC code. No factoring, no secret inputs. Pure symbolic / rational arithmetic (sympy), floating point only in
printed summaries.

Deliverable #1 of agent grdeg: pin down exactly what G_r asserts, which claim must fail and in which direction to open
alpha < 1/3, and the precise (r, d, D, box-side) regime that realises L[1/4].

Setup (resbarrier sec 2.6 / sec 4.2; rank3 R4). Two nonlinear sides of degree d (D = 2d), rank-r linked module,
candidates a_j in {0, +-1} so C := log(#candidates) = r*log 3. x = log N. Write everything at the *leading exponent of
x* (logarithms absorbed), with

    d = x^delta,   r = x^rho,   log B = x^alpha,   C = x^rho  (since C ~ r*log 3),
    cost ~ max(#candidates, B^2) = L[alpha]  =>  rho = alpha.

The algebraic/N part of the smoothness target F = W W' has log-size
    log X  =  kappa_gen * x        (the "N-saving" / generator-height term)
            + (D / r) * C          (the archimedean entropy term = 2*d*log3 ~ x^delta).

The two competing values of kappa_gen (as an *exponent of x*, i.e. log X contribution x^{e}):
  * AGGREGATE (what a G_r failure would need): generator-height floor 1/d, independent of r  =>  x/d = x^{1-delta}.
  * R4 / G_r (proved for general position): floor (r-1)/(D-1) ~ (r/2)/d                       =>  x^{1+rho-delta}.

CEP smoothness: with u = log X / log B, a yield needs u*log u <= C, i.e. at the exponent level
    exponent(log X) - alpha  <=  alpha        (since u = x^{E - alpha}, u log u ~ x^{E-alpha}, C = x^alpha).

This module optimises delta to minimise the achievable alpha under each kappa_gen and returns the exact rational alpha.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict

import sympy as sp


def aggregate_alpha() -> Dict[str, object]:
    """The L[1/4] regime: kappa_gen = 1/d (floor independent of r).

    log X exponent E(delta) = max(1 - delta, delta)   [x/d term vs 2 d log3 entropy term].
    Smoothness feasibility: E - alpha <= alpha  with rho = alpha and r << d (rho < delta).
    Minimise alpha.
    """
    delta = sp.symbols("delta", positive=True)
    # entropy term exponent is delta; N term exponent is (1 - delta). log X exponent = max.
    # Best delta balances them: 1 - delta = delta => delta = 1/2, E = 1/2.
    E = sp.Max(1 - delta, delta)
    delta_star = sp.Rational(1, 2)
    E_star = E.subs(delta, delta_star)  # = 1/2
    # feasibility E - alpha <= alpha  => alpha >= E/2
    alpha = E_star / 2
    return {
        "name": "aggregate (G_r FAILS: generator floor 1/d)",
        "kappa_gen_exponent": "1/delta  (x/d, independent of r)",
        "delta_star": delta_star,
        "logX_exponent": E_star,
        "alpha": Fr(int(sp.numer(alpha)), int(sp.denom(alpha))),
        "rho(=alpha)": Fr(int(sp.numer(alpha)), int(sp.denom(alpha))),
        "r_vs_d": "rho=1/4 < delta=1/2, so r ~ x^{1/4} << d ~ x^{1/2}: valid",
    }


def r4_alpha() -> Dict[str, object]:
    """The L[1/3] floor: kappa_gen = (r-1)/(D-1) ~ (r/2)/d  (rank3 R4, proved in general position).

    log X exponent E(delta, rho) = max(1 + rho - delta, delta).
    rho = alpha. feasibility: E - alpha <= alpha => alpha >= E/2.
    Optimise delta (balance 1 + rho - delta = delta => delta = (1+rho)/2), then solve alpha >= E/2 with rho=alpha.
    """
    alpha = sp.symbols("alpha", positive=True)
    rho = alpha
    delta_star = (1 + rho) / 2
    E_star = delta_star  # balanced value = (1+alpha)/2
    # feasibility at optimum: E_star - alpha <= alpha  => (1+alpha)/2 - alpha <= alpha
    sol = sp.solve(sp.Eq((1 + alpha) / 2 - alpha, alpha), alpha)
    a = sol[0]
    return {
        "name": "R4 / G_r (PROVED general position)",
        "kappa_gen_exponent": "(r-1)/(D-1) ~ (r/2)/d  => x^{1+rho-delta}",
        "delta_star": sp.nsimplify(delta_star.subs(alpha, a)),
        "logX_exponent": sp.nsimplify(E_star.subs(alpha, a)),
        "alpha": Fr(int(sp.numer(a)), int(sp.denom(a))),
        "note": "the r-dependence of the N-saving exactly cancels the 1/r entropy gain",
    }


def crossover_floor(kappa_gen_exp_in_rho):
    """General: given the N-saving exponent as 1 + c*rho - delta for a constant c in [0,1]
    (c=0 is the aggregate floor 1/d; c=1 is R4), return the achievable alpha.

    c interpolates how strongly the generator-height floor grows with r. c=0 -> 1/4, c=1 -> 1/3.
    """
    c = sp.nsimplify(kappa_gen_exp_in_rho)
    alpha = sp.symbols("alpha", positive=True)
    rho = alpha
    # N exponent = 1 + c*rho - delta ; entropy exponent = delta ; balance -> delta=(1+c*rho)/2
    E = (1 + c * rho) / 2
    sol = sp.solve(sp.Eq(E - alpha, alpha), alpha)
    a = sol[0]
    return Fr(int(sp.numer(a)), int(sp.denom(a)))


def summary() -> str:
    agg = aggregate_alpha()
    r4 = r4_alpha()
    lines = [
        "=== grdeg ledger: the G_r failure boundary (deliverable #1) ===",
        "",
        "G_r asserts: the r-1 linking conditions are INDEPENDENT, which (rank3 R1) forces",
        "  v_p(Delta_F) >= r(r-1)(D-2)!/(D-r)! at every p | N, hence the generator-height / N-saving",
        "  exponent kappa_gen = (r-1)/(D-1) ~ (r/2)/d.  This is the claim that must FAIL.",
        "",
        "Direction of failure to open alpha<1/3: the linking conditions must be *correlated*",
        "  (too MANY linked families -> point count exceeds Schmidt's generic bound), which by R1",
        "  is possible ONLY if Delta_F = 0, i.e. the module is frame-less degenerate.",
        "",
        f"AGGREGATE (G_r fails, floor 1/d):   alpha = {agg['alpha']}   at delta={agg['delta_star']} "
        f"(d~x^1/2, r~x^1/4, D=2d, a_j in {{0,+-1}})",
        f"R4 / G_r holds (floor (r-1)/(D-1)): alpha = {r4['alpha']}   ({r4['note']})",
        "",
        "Interpolation c in [0,1] (N-saving exponent 1 + c*rho - delta):",
    ]
    for c in [Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1)]:
        lines.append(f"  c={c}:  alpha = {crossover_floor(c)}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
