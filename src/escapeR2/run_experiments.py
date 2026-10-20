"""Run E1-E3 of notes/escapeR2_prereg.md.  Usage:
  .venv/bin/python -m src.escapeR2.run_experiments experiments/raw
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations

from src.escapeR2 import resplit, sizemodel as SM, transfer as TR


def e1() -> dict:
    out = {"a": [], "b": None, "c": {}}
    # (a) product groups with coordinate point-stabiliser subgroups
    for sizes in ([3, 3], [3, 4], [3, 3, 3]):
        G, blocks = TR.direct_product(sizes)
        Hs = [TR.stabiliser(G, b[0]) for b in blocks]
        # H_j = Stab(point of block j) contains every other factor, so K_j = Q(root of f_j)
        for D in TR.cyclic_subgroups(G):
            rep = TR.transfer_report(G, Hs, D)
            rep2 = TR.transfer_report(G, Hs[:2], D)
            out["a"].append({"sizes": sizes, "D_order": len(D), **rep, "tau2_pair12": rep2["tau2"]})
    # (b) V4: the three quadratic subfields of a biquadratic field
    G, Hs = TR.v4()
    D = TR.closure([TR.identity(4)], 4)
    out["b"] = {**TR.transfer_report(G, Hs, D), "pairs_tau2": [TR.transfer_report(G, [Hs[i], Hs[j]], D)["tau2"]
                                                              for i, j in combinations(range(3), 2)]}
    # (c) census
    for name, G in TR.census_groups().items():
        subs = TR.subgroups_two_generated(G)
        Ds = TR.cyclic_subgroups(G)
        n_pairs = n_triples = viol_bound = viol_tau2 = pos_tau = crit = crit_pos = 0
        for D in Ds:
            for H1, H2 in combinations(subs, 2):
                for Hs in ([H1, H2], [H2, H1]):
                    rep = TR.transfer_report(G, Hs, D)
                    n_pairs += 1
                    viol_tau2 += rep["tau2"] != 0
                    viol_bound += rep["tau"] > rep["bound"]
            for H1, H2, H3 in combinations(subs, 3):
                for Hs in ([H1, H2, H3], [H2, H1, H3], [H3, H1, H2]):
                    rep = TR.transfer_report(G, Hs, D)
                    n_triples += 1
                    viol_bound += rep["tau"] > rep["bound"]
                    pos_tau += rep["tau"] > 0
                    if rep["trivial_meet"]:
                        crit += 1
                        crit_pos += rep["tau"] > 0
        out["c"][name] = {"order": len(G), "subgroups": len(subs), "D_classes": len(Ds), "pair_cases": n_pairs,
                          "triple_cases": n_triples, "P1_violations": viol_bound, "P3_violations": viol_tau2,
                          "triples_tau_pos": pos_tau, "triples_trivial_meet": crit,
                          "triples_trivial_meet_tau_pos": crit_pos}
    return out


def e2() -> dict:
    out = {"Q1": [], "Q2": None, "Q3": [], "Q4": None}
    tiny = Fr(1, 10 ** 6)
    for r in range(0, 9):
        a = SM.closed_no_budget(r)
        out["Q1"].append({"r": r, "closed": str(a), "feasible_at_closed": SM.chain_feasible(a, r, None),
                          "infeasible_below": not SM.chain_feasible(a - tiny, r, None)})
    a1, r1 = SM.best_floor(Fr(1), rmax=8)
    out["Q2"] = {"alpha_star": float(a1), "argmin_r": r1, "exact_check_2/5": SM.chain_feasible(Fr(2, 5), 0, Fr(1)),
                 "below_infeasible_all_r": all(not SM.chain_feasible(Fr(2, 5) - tiny, r, Fr(1)) for r in range(9))}
    for C in [Fr(1), Fr(5, 4), Fr(3, 2), Fr(2), Fr(5, 2), Fr(8, 3), Fr(3), Fr(4), Fr(6), Fr(10)]:
        a, r = SM.best_floor(C, rmax=12)
        closed = min(SM.closed_budget(rr, C) for rr in range(13))
        out["Q3"].append({"C": str(C), "alpha_star": float(a), "argmin_r": r, "closed_min": str(closed),
                          "closed_float": float(closed), "agree_1e-9": abs(a - closed) < Fr(1, 10 ** 9),
                          "above_1/3": a > Fr(1, 3)})
    # Q4: alpha = 3/8, r = 1, delta = 1/4, zeta_1 = 1/2:  v_0 = max(3/4, 1/4 + 5/8) = 7/8, children exponent 7/8 - 1/2
    v0 = max(Fr(3, 4), Fr(1, 4) + Fr(5, 8))
    out["Q4"] = {"v0": str(v0), "children_exponent": str(v0 - Fr(1, 2)), "d_exponent": "1/4",
                 "budget_C_needed": str((v0 - Fr(1, 2)) / Fr(1, 4))}
    return out


def e3() -> dict:
    return {"seeds": [resplit.run_seed(s) for s in range(1, 6)]}


def main(outdir: str) -> None:
    for name, fn in (("e1", e1), ("e2", e2), ("e3", e3)):
        t = time.time()
        res = fn()
        res["seconds"] = round(time.time() - t, 1)
        with open(f"{outdir}/escapeR2_{name}.json", "w") as fh:
            json.dump(res, fh, indent=1, default=str)
        print(name, "done", res["seconds"], "s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "experiments/raw")
