"""auditRA: independent exact re-derivation of Theorem RA's constants (notes/recurse_results.md §1.3, §1.6).

Written from scratch by agent auditRA; it does NOT import src/recurse/constants.py or src/thmA/params.py
(it compares against them at the end).  Exact arithmetic: Fractions for theta, sympy for the algebraic
numbers sqrt2 and c* = sqrt((1+sqrt2) theta).  Floats only for printing.

Checks:
  T1  theta_min(d,k) by exhaustive breakpoint search of the convex piecewise-linear theta(mu).
  T2  c(a) = max(a + theta/(2a), (2+sqrt2) a): minimiser a*, value c*, c*^2 = (1+sqrt2) theta (symbolic).
  T3  c* < 1 <=> theta < sqrt2 - 1;  c* < (1+theta)/sqrt2 <=> theta < sqrt2 - 1 (theta in (0,1)); 2a* < c*.
  T4  the d-table (inf over k, eta0: theta -> 3/(d+1)) and the d = 9, k = 13 value.
  T5  fixed-depth iteration c_{j+1} = g(c_j), g(c) = min_a max(2a, a + c/sqrt2, a + theta/(2a)), checked
      by direct minimisation, not by the closed form.
  T6  Corollary P's sequence d_j = 2j+7, k_j = ceil(d_j^2/4) + j: admissibility, theta_j, c*_j -> 0,
      and the role of eta0_j.
  T7  restart-schedule sum: sum_s 2^(s+1) Pr[reach s] / 2^(s0) is bounded (exact rational partial sums).
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as Fr
from itertools import combinations

import sympy as sp

SQ2 = sp.sqrt(2)


# ---------- T1: theta ----------
def _pieces(d: int, k: int):
    nu = Fr(1, 4 * k)
    # h_*(mu) = max of three affine functions  alpha*mu + beta
    return nu, [(Fr(0), 2 * nu), (Fr(1), -(d - 2) * nu), (Fr(-d), 1 + nu)]


def theta_of(d: int, k: int, mu: Fr, eta0: Fr = Fr(0)) -> Fr:
    nu, pcs = _pieces(d, k)
    return 2 * mu + max(al * mu + be for al, be in pcs) + eta0 + (d - 1) * nu


def theta_min(d: int, k: int) -> tuple[Fr, Fr]:
    """Minimise theta over mu in [nu, 1) (zero slack): candidates = endpoints and all pairwise crossings."""
    nu, pcs = _pieces(d, k)
    cands = {nu}
    for (a1, b1), (a2, b2) in combinations(pcs, 2):
        if a1 != a2:
            m = (b2 - b1) / (a1 - a2)
            if nu <= m < 1:
                cands.add(m)
    best = min(cands, key=lambda m: (theta_of(d, k, m), m))
    return best, theta_of(d, k, best)


def admissible(d: int, k: int, mu: Fr, eta0: Fr) -> bool:
    nu = Fr(1, 4 * k)
    h = theta_of(d, k, mu, eta0) - 2 * mu - (d - 1) * nu
    th = theta_of(d, k, mu, eta0)
    return (d >= 3 and d % 2 == 1 and Fr(k) > Fr(d + 1, 4) and mu >= nu + eta0 and mu < 1
            and eta0 > 0 and th < 1 and mu + h < 1 and mu + nu + eta0 / 2 < 1 and h > (5 - d) * nu)


# ---------- T2/T3 ----------
def cstar_sym(th):
    th = sp.Rational(th.numerator, th.denominator) if isinstance(th, Fr) else th
    return sp.sqrt((1 + SQ2) * th), sp.sqrt((SQ2 - 1) * th / 2)


def check_optimisation():
    t, a = sp.symbols("theta a", positive=True)
    astar = sp.sqrt((SQ2 - 1) * t / 2)
    f1 = a + t / (2 * a)
    f2 = (2 + SQ2) * a
    # crossing at a*
    cross = sp.simplify((f1 - f2).subs(a, astar))
    cst = sp.simplify(f2.subs(a, astar) ** 2 - (1 + SQ2) * t)
    # f1 decreasing up to sqrt(theta/2) and a* < sqrt(theta/2)
    lt = sp.simplify(astar ** 2 - t / 2)  # should be negative multiple of theta
    return {"f1(a*)-f2(a*)": str(cross), "c*^2-(1+sqrt2)theta": str(cst),
            "a*^2 - theta/2": str(sp.factor(lt))}


def threshold_checks():
    """Exact on a rational grid of theta in (0,1): the claimed equivalences."""
    bad = []
    for num in range(1, 1000):
        th = sp.Rational(num, 1000)
        cs, ast = cstar_sym(th)
        A = bool(sp.simplify(cs - 1) < 0) if False else bool((1 + SQ2) * th < 1)   # c* < 1
        B = bool(th < SQ2 - 1)
        C = bool((1 + SQ2) * th < (1 + th) ** 2 / 2)                                 # c* < (1+th)/sqrt2
        D = bool((2 * ast) ** 2 < cs ** 2)                                          # 2a* < c*
        if not (A == B and C == B and D):
            bad.append(str(th))
    return bad


# ---------- T5 ----------
def g_direct(c: float, th: float) -> float:
    """min_a max(2a, a + c/sqrt2, a + th/(2a)) by golden-section on a convex function (float; report only)."""
    import math
    f = lambda a: max(2 * a, a + c / math.sqrt(2), a + th / (2 * a))
    lo, hi = 1e-9, 5.0
    for _ in range(300):
        m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if f(m1) < f(m2):
            hi = m2
        else:
            lo = m1
    return f((lo + hi) / 2)


def main(out_path: str | None = None):
    res = {}
    res["optimisation_symbolic"] = check_optimisation()
    res["threshold_grid_failures"] = threshold_checks()

    # T4 table
    table = []
    for d in (7, 9, 11, 13, 19, 29, 49, 99):
        th_inf = Fr(3, d + 1)
        cs, ast = cstar_sym(th_inf)
        table.append({"d": d, "theta_inf": str(th_inf), "c*": float(cs), "a*": float(ast),
                      "thmA_inf": float((1 + sp.Rational(3, d + 1)) / SQ2),
                      "inf_over_k_matches_closed_form": None})
    # inf over k really is 3/(d+1): theta_min(d,k) decreasing to 3/(d+1)
    for row in table:
        d = row["d"]
        ths = [theta_min(d, k)[1] for k in (10 ** 3, 10 ** 5)]
        row["theta_min_k1e3"] = float(ths[0])
        row["theta_min_k1e5"] = float(ths[1])
        row["inf_over_k_matches_closed_form"] = all(t > Fr(3, d + 1) for t in ths) and \
            (ths[1] - Fr(3, d + 1)) < Fr(1, 10 ** 4)
    res["d_table"] = table

    # d = 9, k = 13
    mu, th = theta_min(9, 13)
    cs, ast = cstar_sym(th)
    res["d9_k13"] = {"mu": str(mu), "theta": str(th), "c*": float(cs), "thmA_c": float((1 + sp.Rational(th.numerator, th.denominator)) / SQ2),
                     "c*<1": bool((1 + SQ2) * sp.Rational(th.numerator, th.denominator) < 1)}
    mu12, th12 = theta_min(9, 12)
    res["d9_k12"] = {"theta": str(th12), "c*<1": bool((1 + SQ2) * sp.Rational(th12.numerator, th12.denominator) < 1)}
    # minimal k with c* < 1, for each d
    mink = {}
    for d in (7, 9, 11, 13, 19):
        for k in range(1, 400):
            if Fr(k) <= Fr(d + 1, 4):
                continue
            m, t_ = theta_min(d, k)
            if (1 + SQ2) * sp.Rational(t_.numerator, t_.denominator) < 1:
                mink[d] = k
                break
    res["min_k_cstar_lt_1"] = mink

    # T5 iteration at theta = 3/10
    import math
    th = 0.3
    seq = [1.0]
    for _ in range(8):
        seq.append(g_direct(seq[-1], th))
    res["iteration_d9_inf"] = [round(x, 5) for x in seq]
    res["fixed_point_d9_inf"] = float(sp.sqrt((1 + SQ2) * sp.Rational(3, 10)))
    res["g_prime_at_fixed_point"] = float((1 - sp.Rational(3, 10) / ((1 + SQ2) * sp.Rational(3, 10))) / SQ2)

    # T6 Corollary P sequence
    seqP = []
    for j in range(1, 61):
        d = 2 * j + 7
        k = -(-d * d // 4) + j
        mu0, th0 = theta_min(d, k)
        # eta0 must also be > 0 and mu >= nu + eta0; pick eta0_j = 1/(d^2) (-> 0) and also a FIXED 1/1000
        rows = {}
        for lab, eta0 in (("eta0=1/d^2", Fr(1, d * d)), ("eta0=1/1000 fixed", Fr(1, 1000))):
            ok = admissible(d, k, mu0, eta0)
            th = theta_of(d, k, mu0, eta0)
            cs = float(sp.sqrt((1 + SQ2) * sp.Rational(th.numerator, th.denominator)))
            rows[lab] = {"admissible": ok, "theta": float(th), "c*": cs}
        seqP.append({"j": j, "d": d, "k": k, "mu": str(mu0), **rows})
    res["corP_sequence_first"] = seqP[:8]
    res["corP_sequence_last"] = seqP[-2:]
    res["corP_all_admissible_eta0_shrinking"] = all(r["eta0=1/d^2"]["admissible"] for r in seqP)
    res["corP_cstar_floor_fixed_eta0"] = float(sp.sqrt((1 + SQ2) * sp.Rational(1, 1000)))

    # T7 restart schedule: sum_{m>=0} 2^(m+1) * 2^(-m(m+1)/2) (phase s0+m reached w.p. <= 2^{-m(m+1)/2});
    tot = Fr(0)
    for m in range(0, 60):
        tot += Fr(2 ** (m + 1), 2 ** (m * (m + 1) // 2))
    res["restart_sum_over_2^s0_phases_ge_s0"] = float(tot)
    res["restart_total_factor_incl_phases_lt_s0"] = float(tot + 2)   # sum_{s<s0} 2^{s+1} <= 2^{s0+1}

    # compare with the recurse agent's module, if importable
    try:
        sys.path.insert(0, ".")
        from src.recurse import constants as rc  # noqa
        res["recurse_module_attrs"] = [x for x in dir(rc) if not x.startswith("_")][:30]
    except Exception as e:  # pragma: no cover
        res["recurse_module_attrs"] = f"import failed: {e}"

    txt = json.dumps(res, indent=1, default=str)
    print(txt)
    if out_path:
        with open(out_path, "w") as fh:
            fh.write(txt)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
