"""mf2L: leading-order exponent bookkeeping for the input substitution of notes/mf2L_results.md §4
(Drappeau's Lemme 5  ->  Harper's published large-values bound (4.1) + Hyp SC).  Exact rationals only.

A quantity is represented by the leading term of its logarithm, c * D^e * lam^f, as the pair (e, f); lam = log log N,
D = log N / lam.  "A << B" holds at leading order iff (e_A, f_A) < (e_B, f_B) lexicographically.  This checks exponent
arithmetic only; it is not an independent proof (same status as src/thmAB/params.py, which it imports).

Design: c12alg's 2L_DW(alpha), alpha in [3/8, 2/5).  Four character families (notes/mf2L_results.md §1):
  child_alg  (thmAB Lemma EN-A, FAM-A):  y = B,  values log x ~ D^{2a} lam, moduli Delta'/delta <= 2 S2^2
  child_rat  (c12rat Lemma EN, FAM):     y = B,  values D^{2a} lam, moduli <= 4 Q2 S2^2  (Q2 <= B2)
  par_alg    (thmAB §9 parents):         y = B2, values D^{2-3a} lam, moduli <= 2 V1^2 / q
  par_rat    (c12rat RAT-2):             y = B,  values D^{2a} lam, moduli <= 2 V1^2
Each family has a loss exponent (the log of the family normaliser 𝔇 is ~ D^loss lam^2) and a u-exponent.
"""
from __future__ import annotations

from fractions import Fraction as F
from typing import Dict, List, Tuple

from src.thmAB.params import theta, theta_par

Term = Tuple[F, F]


def families(al: F) -> Dict[str, Dict[str, Term]]:
    one, two, zero = F(1), F(2), F(0)
    th, thp = theta(al), theta_par(al)
    return {
        "child_alg": dict(y=(al, one), x=(2 * al, one), Q=(al, one), loss=(th, two), u=(al, zero)),
        "child_rat": dict(y=(al, one), x=(2 * al, one), Q=(max(al, 4 * al - 1), one), loss=(th, two), u=(al, zero)),
        "par_alg": dict(y=(4 * al - 1, one), x=(2 - 3 * al, one), Q=(1 - al, one), loss=(thp, two),
                        u=(3 - 7 * al, zero)),
        # c12rat's parent family uses the child thinness threshold kappa0 = e^{3 C_T D^theta lam^2} in h_*,
        # and Xi_p = e^{O(D^{1-2a} lam)}; theta >= 1 - 2a on [3/8, 2/5), so the loss is (theta, 2).
        "par_rat": dict(y=(al, one), x=(max(2 * al, 1 - al), one), Q=(1 - al, one),
                        loss=(max(th, 1 - 2 * al), two), u=(al, zero)),
    }


def kappa(al: F) -> F:
    """Smallest admissible Hyp SC exponent: conductors up to exp(u^kappa), saving exp(-u^kappa), any kappa > this."""
    return max(f["loss"][0] / f["u"][0] for f in families(al).values())


def comparisons(al: F) -> List[Tuple[str, Term, Term]]:
    out = []
    for name, f in families(al).items():
        out += [
            (f"{name}: Harper (4.1) needs y <= x^(1/100)", f["y"], f["x"]),
            (f"{name}: Farey points 4Q^2 <= x^(1/1000) (moduli are small conductors)", f["Q"], f["x"]),
            (f"{name}: R1 = e^(C D^loss lam^2) is y^o(1)", f["loss"], f["y"]),
            (f"{name}: Hyp SC exponent kappa < 1 (loss << u)", (f["loss"][0], F(0)), (f["u"][0], F(0))),
            # Lemma SC-D regime: (log x)^{1/2} <= log y (strictly, at leading order), so u <= log y
            (f"{name}: SC-D regime (log x)^(1/2) << log y", (f["x"][0] / 2, f["x"][1] / 2), f["y"]),
        ]
    return out


def informational(al: F) -> Dict[str, str]:
    """Comparisons that are EXPECTED to fail; they document why the residual input is what it is."""
    out = {}
    for name, f in families(al).items():
        # Siegel-Walfisz range exp(c sqrt(log y)) versus the needed conductor range exp(D^loss lam^2)
        sw = (f["y"][0] / 2, f["y"][1] / 2)
        out[f"{name}: needed conductors inside Siegel-Walfisz range"] = "yes" if f["loss"] < sw else "no"
        # lambda-gap: pointwise reach e^{u} vs classical large sieve threshold rho^{-1/2} = e^{(u log u)/2}
        out[f"{name}: pointwise e^-u + classical LS cover all conductors"] = (
            "yes" if (f["u"][0], F(0)) >= (f["u"][0], F(1)) else "no (lambda-gap; closed by Harper (4.1))")
    return out


def margins(al: F) -> Dict[str, str]:
    return {n: ("ok" if a < b else ("borderline" if a == b else "FAIL")) for n, a, b in comparisons(al)}


def all_ok(al: F) -> bool:
    return all(v == "ok" for v in margins(al).values())


if __name__ == "__main__":  # pragma: no cover
    for al in [F(3, 8), F(19, 50), F(7, 18), F(39, 100), F(399, 1000)]:
        print(al, "kappa_min =", kappa(al), "all ok:", all_ok(al))
        for k, v in informational(al).items():
            print("   ", k, "->", v)
