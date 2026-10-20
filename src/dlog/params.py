"""Exact cost constants for Theorem DL (notes/dlog_results.md §4).  Exact arithmetic (Fractions / sympy).

Cost model (exponents of L = L_p[1/2, 1]; B = L^a):
  rows:            a + max(s, theta/(2a))      s = shape-test exponent (1/sqrt2 with Lenstra-Pomerance on
                                               |u0|, |v0| <= sqrt(2p));  theta from Theorem A (src/thmA/params.py)
  exact integer LA: omega * a                  omega = 3 (Dixon 1982 p-adic lifting, dense), any omega >= 2 allowed
  c_DL(theta, s, omega) = min_a max(a + s, a + theta/(2a), omega*a).
For theta <= 1/2, s = 1/sqrt2, omega = 3:  c_DL = (1 + theta)/sqrt2 (Theorem A's constant), attained at a = theta/sqrt2.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List

import sympy as sp

from src.thmA import params as tp

SQRT2 = sp.sqrt(2)


def _as_sp(x):
    if isinstance(x, Fr):
        return sp.Rational(x.numerator, x.denominator)
    return sp.nsimplify(x) if not isinstance(x, sp.Basic) else x


def c_dl(theta, s=None, omega=3):
    """Exact min_a max(a + s, a + theta/(2a), omega*a) over a > 0 (sympy expression) and the minimiser."""
    th = _as_sp(theta)
    s = SQRT2 / 2 if s is None else _as_sp(s)
    om = _as_sp(omega)
    a = sp.symbols("a", positive=True)
    fs = [a + s, a + th / (2 * a), om * a]
    cands = {sp.sqrt(th / 2), th / (2 * s), sp.sqrt(th / (2 * (om - 1))), s / (om - 1)}
    best = None
    for c in cands:
        c = sp.nsimplify(c) if not c.free_symbols else c
        if not (c > 0):
            continue
        val = sp.Max(*[f.subs(a, c) for f in fs])
        val = sp.simplify(val)
        if best is None or sp.simplify(val - best[0]) < 0:
            best = (val, c)
    return best


def c_thmA(theta) -> sp.Expr:
    return (1 + _as_sp(theta)) / SQRT2


def la_never_binds(theta, omega=3) -> bool:
    """omega*theta/sqrt2 <= (1+theta)/sqrt2  <=>  (omega-1) theta <= 1."""
    th = _as_sp(theta)
    return bool((_as_sp(omega) - 1) * th <= 1)


def table(ds=(3, 5, 7, 9, 11, 19), omegas=(3,)) -> List[Dict]:
    out = []
    for d in ds:
        kmin = tp.min_k(d)
        for k in sorted({kmin, 50, 1000} - {None}) if kmin else (50, 1000):
            mu, th = tp.best_mu(d, k)
            row = {"d": d, "k": k, "mu": mu, "theta": th}
            for om in omegas:
                val, a = c_dl(th, omega=om)
                row[f"cDL_w{om}"] = val
                row[f"cDL_w{om}_float"] = float(val)
            row["cA"] = float(c_thmA(th))
            row["la_binds_w3"] = not la_never_binds(th, 3)
            out.append(row)
        inf_th = sp.Rational(3, d + 1)
        out.append({"d": d, "k": "inf", "theta": inf_th, "cDL_w3": c_dl(inf_th)[0],
                    "cDL_w3_float": float(c_dl(inf_th)[0])})
    return out


if __name__ == "__main__":
    for r in table(omegas=(3, sp.Rational(237, 100))):
        print(r)
