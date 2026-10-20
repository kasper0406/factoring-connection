"""Exact cost bookkeeping for Theorem DL25 (rigorous L_p[2/5] discrete logs; notes/dlog25_results.md §6).

All arithmetic is exact (Fractions); 2^{2/5} and 2^{-2/5} are bracketed by rationals (src/lv5unif/params.py).

Cost model, exponents of L := L_p[2/5, 1] = exp(Lambda), Lambda = (log p)^{2/5}(log log p)^{3/5}:
  parameters (a, kappa, d0) as in Theorem U (notes/lv5unif_results.md §1) with N := p;
  s     = 2 kappa (1 + C_gamma)       1/Pr[row per attempt] = L^{s + o(1)}       (Theorem K-U with Lemma S_p'')
  sigma = 2 theta0(d0, kappa)/(5a)   special-q batch per attempt T_SQ = L^{max(a, sigma) + o(1)}
  M_rows = L^{a + o(1)} rows; post-processing (Prop XK-U) = L^{3a + o(1)}
  shape test = FACTOR_U on |u0|, |v0| <= sqrt(2p):  L_{sqrt(2p)}[2/5, c_sub] = L_p[2/5, 2^{-2/5} c_sub + o(1)]
  ==>  c_DL = max(3a,  a + s + max(a, sigma),  a + s + 2^{-2/5} c_sub)                         (6.1)
  where c_sub = c_U(a', kappa', d0') is Theorem U's constant for ANY admissible subroutine tuple:
        c_U = max(2a, a + s + max(a, sigma), (a + s)/(1 - 2^{-2/5}))                            (Theorem U, (10.3))
Facts verified here exactly (and proved in the note, Lemma C25):
  (i)  with the same tuple for both levels, c_DL <= max(3a, c_U);
  (ii) for any c* and any tuple, max(3a, T1, a+s+2^{-2/5}c*) < c*  implies  c_U(tuple) < c*;
       hence inf c_DL >= inf c_U when the subroutine uses c* = inf c_U, so the two infima coincide whenever 3a does not
       bind at a (near-)optimal tuple.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, Tuple

from src.lv5unif import params as up

theta0 = up.theta0
sigma = up.sigma


def two_pow_minus_2_5_bracket() -> Tuple[Fr, Fr]:
    """Rationals lo < 2^{-2/5} < hi (from the exact bracket of 2^{2/5})."""
    lo_t, hi_t = up.two_pow_2_5_bracket()
    return 1 / hi_t, 1 / lo_t


def s_of(kappa, C_gamma) -> Fr:
    return 2 * Fr(kappa) * (1 + Fr(C_gamma))


def terms(a, d0, kappa, C_gamma, c_sub: Tuple[Fr, Fr]) -> Dict[str, Tuple[Fr, Fr]]:
    """Exact brackets of the three terms of (6.1); c_sub is a bracket (lo, hi) of the subroutine constant."""
    a, d0, kappa = Fr(a), Fr(d0), Fr(kappa)
    s = s_of(kappa, C_gamma)
    sg = sigma(a, d0, kappa)
    m_lo, m_hi = two_pow_minus_2_5_bracket()
    la = (3 * a, 3 * a)
    sq = (a + s + max(a, sg),) * 2
    sh = (a + s + m_lo * c_sub[0], a + s + m_hi * c_sub[1])
    return {"LA_3a": la, "specialq": sq, "shape": sh}


def c_dl(a, d0, kappa, C_gamma, c_sub: Tuple[Fr, Fr]) -> Tuple[Fr, Fr]:
    t = terms(a, d0, kappa, C_gamma, c_sub)
    return max(v[0] for v in t.values()), max(v[1] for v in t.values())


def binding(a, d0, kappa, C_gamma, c_sub) -> str:
    t = terms(a, d0, kappa, C_gamma, c_sub)
    return max(t, key=lambda k: t[k][1])


def c_u(a, d0, kappa, C_gamma) -> Tuple[Fr, Fr]:
    return up.cost_constant(a, d0, kappa, C_gamma)


def default_grids():
    grid_a = [Fr(i, 20) for i in range(1, 41)]
    grid_kappa = [Fr(1, 2 ** j) for j in range(0, 12)] + [Fr(i, 200) for i in range(1, 41)]
    return grid_a, grid_kappa


def _d0_cands(kappa):
    out = set()
    for num in range(1, 400):
        d0 = Fr(num, 20)
        if Fr(1, 4) * 8 * kappa <= d0 * d0 <= 4 * 8 * kappa:
            out.add(d0)
    return out


def optimise_u(C_gamma, grid_a=None, grid_kappa=None) -> Dict:
    ga, gk = default_grids()
    return up.optimise(Fr(C_gamma), grid_a or ga, grid_kappa or gk)


def optimise_dl(C_gamma, c_sub: Tuple[Fr, Fr], grid_a=None, grid_kappa=None) -> Dict:
    ga, gk = default_grids()
    grid_a, grid_kappa = grid_a or ga, grid_kappa or gk
    best = None
    for kappa in grid_kappa:
        for d0 in _d0_cands(kappa):
            for a in grid_a:
                lo, hi = c_dl(a, d0, kappa, C_gamma, c_sub)
                if best is None or hi < best["hi"]:
                    best = {"a": a, "kappa": kappa, "d0": d0, "lo": lo, "hi": hi,
                            "binding": binding(a, d0, kappa, C_gamma, c_sub)}
    return best


def same_tuple_dominated(a, d0, kappa, C_gamma) -> bool:
    """Lemma C25(i): with the subroutine at the same tuple, c_DL <= max(3a, c_U) (checked on the upper brackets;
    the inequality a + s + 2^{-2/5} c_U <= c_U holds because (a+s) <= (1 - 2^{-2/5}) c_U)."""
    cu = c_u(a, d0, kappa, C_gamma)
    lo, hi = c_dl(a, d0, kappa, C_gamma, cu)
    # use the exact identity on the 'recursion' term rather than the (slightly loose) brackets:
    s = s_of(kappa, C_gamma)
    m_lo, m_hi = two_pow_minus_2_5_bracket()
    f_lo, f_hi = up.recursion_factor_bracket()      # brackets of 1/(1 - 2^{-2/5})
    # a+s+2^{-2/5} c_U <= c_U  <=>  (a+s) <= (1-2^{-2/5}) c_U ; c_U >= (a+s)/(1-2^{-2/5}) by definition (lower bracket)
    return (a + s) * f_lo <= cu[1] and lo <= max(3 * Fr(a), cu[1]) + (m_hi - m_lo) * cu[1]


def xk_cost_ledger() -> list:
    """Prop XK-U cost ledger (exponent of L, extra factors).  Every extra factor is poly(d, log p) = L^{o(1)}
    because d = O(D^{1/5}) <= log p; see notes/dlog25_results.md §4."""
    return [
        ("exponent matrix E (Lemma XK0-U, incl. r-maximal orders at <= O(log p) bad r <= B)", "a", "poly(d, log p) per entry; B^{1+o(1)} bad-prime scan"),
        ("rank profile mod random prime q (dense elimination)", "3a", "M * C_col * rho"),
        ("Dixon p-adic lifting, d right-hand sides", "3a", "d * rho^3 * polylog"),
        ("verification x^T E = 0", "2a", "nnz(E) * beta_x, beta_x = O(rho log log p)"),
        ("unit elimination (Lemma XK4-U): roots of f, logs of <= 2kM atoms, LLL in dim <= 2d", "2a", "poly(d) * (M * P + kappa^2), P, kappa = O(d beta_x + d^2 + d log p)"),
        ("W(x*), A_*, B_* (modular exponentiations with L^{a}-bit exponents)", "2a", "M * beta_x * polylog"),
    ]


if __name__ == "__main__":
    import json
    out = {}
    for Cg in (0, 1, 5, 20):
        bu = optimise_u(Cg)
        csub = (bu["lo"], bu["hi"])
        bd = optimise_dl(Cg, csub)
        out[str(Cg)] = {"U": {k: str(v) for k, v in bu.items()}, "DL": {k: str(v) for k, v in bd.items()},
                        "U_hi_float": float(bu["hi"]), "DL_hi_float": float(bd["hi"]), "DL_lo_float": float(bd["lo"])}
        print(Cg, "U <=", float(bu["hi"]), bu["a"], bu["kappa"], bu["d0"],
              "| DL in [", float(bd["lo"]), ",", float(bd["hi"]), "]", bd["a"], bd["kappa"], bd["d0"], bd["binding"])
    print(json.dumps(out, indent=1))
