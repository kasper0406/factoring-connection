"""Verdicts for the pre-registered e2 / e2s / e3 (notes/coverage38_prereg.md).
Run: PYTHONPATH=. .venv/bin/python -m src.coverage38.analyze"""
from __future__ import annotations

import json
import os
from math import exp

import numpy as np

RAW = "experiments/raw"
ANA = "experiments/analysis"


def _load(name):
    path = os.path.join(RAW, f"coverage38_{name}.jsonl")
    if not os.path.exists(path):
        return []
    return [json.loads(l) for l in open(path)]


def e2():
    rows = _load("e2")
    complete = len(rows) == 14
    out = {"complete": complete, "cells": [], "D1": {}, "D2": complete, "D3": complete,
           "D4": complete}
    for r in rows:
        d = r["hstar_all"] - r["hstar_hard"]
        d3 = abs(r["hstar_hard"] - r["pred_hard"])
        ok4 = 0.5 <= r["tail_enrich_easy"] <= 2.5 if r["G"] == 4 else None
        out["D2"] &= d <= 3
        out["D3"] &= d3 <= 4
        if ok4 is not None:
            out["D4"] &= ok4
        out["cells"].append({k: r[k] for k in ("n", "J", "G", "frac_easy", "hstar_all",
                                               "hstar_hard", "pred_all", "pred_hard",
                                               "tail_enrich_easy")})
    for n in sorted({r["n"] for r in rows}):
        f = {r["G"]: r["frac_easy"] for r in rows if r["n"] == n}
        out["D1"][n] = {"ratio": f[16] / f[4], "pass": 0.10 <= f[16] / f[4] <= 0.35}
    return out


def e2s():
    rows = _load("e2s")
    out = {"complete": len(rows) == 4, "cells": [], "S1": len(rows) == 4, "S2": {}}
    for r in rows:
        lo, hi = (0.8, 1.25) if r["G"] == 4 else (0.7, 1.4)
        out["S1"] &= lo <= r["R"] <= hi
        out["cells"].append({k: r[k] for k in ("n", "G", "frac_easy", "Href", "unc_frac_all",
                                               "unc_frac_easy", "unc_frac_hard", "R", "z",
                                               "hstar_all", "hstar_hard")})
    for n in sorted({r["n"] for r in rows}):
        f = {r["G"]: r["frac_easy"] for r in rows if r["n"] == n}
        out["S2"][n] = {"ratio": f[16] / f[4], "pass": 0.10 <= f[16] / f[4] <= 0.35}
    return out


def e3():
    rows = _load("e3")
    # completeness (post-run fix after GPT-6.1-sol review): 3 q x 2 H, 1000 real, 100 surrogates
    complete = (len(rows) == 6 and len({(r["q"], r["H"]) for r in rows}) == 6
                and all(len(r["real"]) >= 990 and len(r["surrogate"]) == 100 for r in rows))
    out = {"complete": complete, "cells": [], "O1": complete, "O2": complete,
           "O3": complete, "O4": complete, "O5": {}}
    pool = {}
    for r in rows:
        real, sur = r["real"], r["surrogate"]
        unc = [x for x in real if x["hits"] == 0]
        cov = [x for x in real if x["hits"] > 0]
        med = lambda L, k: float(np.median([x[k] for x in L]))
        o1 = min(x["E_abs"] for x in real) / r["lam_q"]
        o2 = med(real, "E_abs") / r["E_pred"]
        o3 = {k: med(unc, k) / med(cov, k) for k in ("E_abs", "Mlow")}
        o4 = {k: med(sur, k) / med(real, k) for k in ("E_abs", "Mlow", "Mall")}
        frac_hoeff = float(np.mean([x["Mall"] <= r["hoeff"] for x in sur]))
        out["O1"] &= o1 >= 5
        out["O2"] &= 0.8 <= o2 <= 1.25
        out["O3"] &= all(0.9 <= v <= 1.1 for v in o3.values())
        out["O4"] &= all(0.8 <= v <= 1.25 for v in o4.values()) and frac_hoeff >= 0.5
        pool.setdefault(r["H"], [0, 0, 0.0])
        pool[r["H"]][0] += len(unc)
        pool[r["H"]][1] += len(real)
        pool[r["H"]][2] += len(real) * exp(-r["lam_q"])
        out["cells"].append({"q": r["q"], "H": r["H"], "k0": r["k0"], "lam_q": r["lam_q"],
                             "J": r["J"], "sqrtJ": r["sqrtJ"], "n_unc": len(unc),
                             "min_Eabs_over_lam": o1, "Eabs_over_Epred": o2,
                             "unc_over_cov": o3, "sur_over_real": o4,
                             "frac_sur_below_hoeff": frac_hoeff,
                             "median_Mall_real_over_sqrtJ": med(real, "Mall") / r["sqrtJ"],
                             "median_Mlow_real_over_sqrtJ": med(real, "Mlow") / r["sqrtJ"]})
    for H, (u, t, e) in pool.items():
        out["O5"][H] = {"obs": u, "exp": e, "ratio": u / e, "pass": 0.5 <= u / e <= 2}
    return out


if __name__ == "__main__":
    res = {"e2": e2(), "e2s": e2s(), "e3": e3()}
    os.makedirs(ANA, exist_ok=True)
    with open(os.path.join(ANA, "coverage38_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
