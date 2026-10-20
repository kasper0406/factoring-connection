"""Exact constants for the recursive shape test (agent recurse; notes/recurse_results.md §1.4).

Theory (notes/recurse_results.md Theorem RA, Lemma F):
  with smoothness bound B = L^a and a recursive shape test, the expected time is L^{c(a)+o(1)},
      c(a) = max( a + theta/(2a),  (2 + sqrt2) a )          [ >= 2a, so linear algebra never binds ]
  minimised at a* = sqrt((sqrt2 - 1) theta / 2), with
      c*^2 = (1 + sqrt2) theta.
  Fixed-depth recursion (LP at the bottom):  c_0 = 1,  c_{j+1} = (theta/c_j + c_j)/sqrt2  (valid while
  c_j^2 >= theta), so c_1 = (1 + theta)/sqrt2 is Theorem A.
All decisions below are EXACT: sympy algebraic numbers, or equivalent rational inequalities.
theta(d, k, mu, eta0) is imported from src/thmA/params.py (agent thmA), not re-derived.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List

import sympy as sp

from src.thmA import params as tpar

SQ2 = sp.sqrt(2)


def cstar_sq(th) -> sp.Expr:
    """c*^2 = (1 + sqrt2) theta (exact)."""
    return sp.nsimplify(1 + SQ2) * sp.Rational(th.numerator, th.denominator) if isinstance(th, Fr) \
        else (1 + SQ2) * th


def cstar(th) -> sp.Expr:
    return sp.sqrt(cstar_sq(th))


def astar(th) -> sp.Expr:
    t = sp.Rational(th.numerator, th.denominator) if isinstance(th, Fr) else th
    return sp.sqrt((SQ2 - 1) * t / 2)


def c_of_a(a, th) -> sp.Expr:
    """Fixpoint exponent of the self-recursive algorithm with B = L^a (exact)."""
    t = sp.Rational(th.numerator, th.denominator) if isinstance(th, Fr) else th
    a = sp.nsimplify(a)
    return sp.Max(a + t / (2 * a), (2 + SQ2) * a)


# ---------------------------------------------------------------- rational certificates (no radicals)
def rec_lt_1(th: Fr) -> bool:
    """c* < 1  <=>  (1+sqrt2) theta < 1  <=>  theta < sqrt2 - 1  <=>  (1 + theta)^2 < 2."""
    return th >= 0 and (1 + th) ** 2 < 2


def rec_beats_A(th: Fr) -> bool:
    """c* < (1+theta)/sqrt2  <=>  2(1+sqrt2)theta < (1+theta)^2  <=>  2 sqrt2 theta < 1 + theta^2
    <=> 8 theta^2 < (1 + theta^2)^2 (both sides positive).  Equivalent to theta < sqrt2 - 1 for 0<theta<1."""
    return 8 * th * th < (1 + th * th) ** 2


def la_slack(th: Fr) -> bool:
    """2 a* < c*  <=> 2 < 2 + sqrt2 : always true (ratio 2a*/c* = 2 - sqrt2)."""
    return True


def bs_ok(th: Fr) -> bool:
    """Prop BS cost condition B <= S^2, i.e. a*^2 <= theta/2  <=> (sqrt2-1)/2 <= 1/2 : always."""
    return th > 0


def fixed_depth(th, depth: int) -> List[sp.Expr]:
    """c_0 = 1 (Lenstra-Pomerance), c_{j+1} = g(c_j), g(c) = min_a max(2a, a + c/sqrt2, a + theta/(2a))
    = (theta/c + c)/sqrt2 when c^2 >= theta (else sqrt(2 theta))."""
    t = sp.Rational(th.numerator, th.denominator) if isinstance(th, Fr) else th
    cs = [sp.Integer(1)]
    for _ in range(depth):
        c = cs[-1]
        if sp.simplify(c ** 2 - t) >= 0:
            cs.append(sp.nsimplify((t / c + c) / SQ2))
        else:
            cs.append(sp.sqrt(2 * t))
    return cs


def table(ds=(7, 9, 11, 13, 19, 29, 49, 99), ks=(50, 1000)) -> List[Dict]:
    rows = []
    for d in ds:
        kmin = None
        for k in range(1, 20000):
            if Fr(k) <= Fr(d + 1, 4):
                continue
            mu, th = tpar.best_mu(d, k)
            if rec_lt_1(th) and tpar.admissible(d, k, mu)["theta_le_1"]:
                kmin = k
                break
        mu0, th0 = tpar.best_mu(d, kmin)
        row = {"d": d, "kmin": kmin, "theta_kmin": th0,
               "c_rec_kmin": cstar(th0), "c_A_kmin": (1 + sp.Rational(th0.numerator, th0.denominator)) / SQ2}
        for k in ks:
            _, th = tpar.best_mu(d, k)
            row[f"c_rec_k{k}"] = cstar(th)
            row[f"c_A_k{k}"] = (1 + sp.Rational(th.numerator, th.denominator)) / SQ2
        th_inf = Fr(3, d + 1)
        row["theta_inf"] = th_inf
        row["c_rec_inf"] = cstar(th_inf)
        row["c_A_inf"] = (1 + sp.Rational(3, d + 1)) / SQ2
        row["c_noshape_inf"] = sp.sqrt(2 * sp.Rational(3, d + 1))
        rows.append(row)
    return rows


def main():
    print("c*^2 = (1+sqrt2) theta;  a* = sqrt((sqrt2-1) theta/2);  2a*/c* = 2 - sqrt2 =", float(2 - SQ2))
    hdr = ["d", "kmin(c<1)", "theta@kmin", "c_rec@kmin", "c_A@kmin", "c_rec k=50", "c_A k=50",
           "c_rec k=1000", "c_rec inf", "c_A inf", "sqrt(2theta) inf"]
    print(" | ".join(hdr))
    for r in table():
        print(" | ".join(str(x) for x in [r["d"], r["kmin"], r["theta_kmin"],
              f"{float(r['c_rec_kmin']):.5f}", f"{float(r['c_A_kmin']):.5f}",
              f"{float(r['c_rec_k50']):.5f}", f"{float(r['c_A_k50']):.5f}",
              f"{float(r['c_rec_k1000']):.5f}", f"{float(r['c_rec_inf']):.5f}",
              f"{float(r['c_A_inf']):.5f}", f"{float(r['c_noshape_inf']):.5f}"]))
    print("(c_rec > c_A only when theta >= sqrt2 - 1; then use LP / Theorem A: the best constant is the min.)")
    print("fixed-depth iteration at d = 9, theta = 3/10 (infimum):")
    cs = fixed_depth(Fr(3, 10), 8)
    print("  ", [f"{float(c):.5f}" for c in cs], " limit", f"{float(cstar(Fr(3, 10))):.5f}")


if __name__ == "__main__":
    main()
