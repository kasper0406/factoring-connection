"""Verdicts for the pre-registered predictions (notes/cov38b_prereg.md + addendum A).
Run: PYTHONPATH=. .venv/bin/python -m src.cov38b.analyze"""
from __future__ import annotations

import json
import os
from math import log, isqrt

import numpy as np

RAW = "experiments/raw"
OUT = "experiments/analysis/cov38b_summary.json"


def _load(name):
    with open(os.path.join(RAW, f"cov38b_{name}.jsonl")) as f:
        return [json.loads(x) for x in f]


def e1():
    R = _load("e1")
    pr = [r for r in R if not r["all_residues"]]
    al = {r["q"]: r for r in R if r["all_residues"]}
    assert len(pr) == 15 and len(al) == 9, "grid incomplete"
    p1a = all(r["hstar"] <= r["thue"] for r in pr)
    inb = [r["H05"] - 1 <= r["hstar"] <= r["H95"] + 1 for r in pr]
    p1b = sum(inb) >= 12
    high = sum(r["hstar"] > r["H95"] + 1 for r in pr)
    p1c = all(r["hstar"] < r["thue"] for r in pr if r["Q0"] >= 10 ** 7)
    d = [al[q]["hstar"] - next(r["hstar"] for r in pr if r["q"] == q) for q in al]
    p1d = all(0 <= x <= 3 for x in d)
    rows = [{k: r[k] for k in ("Q0", "q", "Y", "J", "thue", "n", "hstar", "H05", "H50", "H95", "c_emp")}
            for r in pr]
    return {"P1a": p1a, "P1b": p1b, "P1b_in_band": int(sum(inb)), "P1b_high_side": high, "P1c": p1c,
            "P1d": p1d, "all_minus_prime": d, "rows": rows}


def e2():
    R = _load("e2")
    out = {}
    kap = {r["q"]: r for r in R if r["part"] == "kappa"}
    assert len(kap) == 3
    out["P2a"] = all(0.5 <= r["kappa_over_sqrtJlnq"] <= 2 for r in kap.values())
    out["kappa_ratio"] = {q: r["kappa_over_sqrtJlnq"] for q, r in kap.items()}
    W = [r for r in R if r["part"] == "wiener"]
    # exact expected Wiener grid (addendum A: k0 = q // (8 Hp^3) >= 16), both shapes
    exp_grid = {(q, Hp, sh) for q in kap for Hp in (1, 2, 4, 8, 16, 32) if q // (8 * Hp ** 3) >= 16
                for sh in ("sharp", "tri")}
    assert {(r["q"], r["Hp"], r["shape"]) for r in W} == exp_grid and len(W) == len(exp_grid), "E2 grid incomplete"
    ok_ratio, slopes, ok_lb, lb_ratio = True, {}, True, {}
    for q in kap:
        for sh in ("sharp", "tri"):
            rs = sorted([r for r in W if r["q"] == q and r["shape"] == sh], key=lambda r: r["Hp"])
            for r in rs:
                if r["Hp"] >= 2 and not (0.2 <= r["W"] / r["Hp"] <= 5):
                    ok_ratio = False
                if r["LB"] > r["W"] + 1 + 1e-9:
                    ok_lb = False
                if sh == "tri" and r["Hp"] >= 2:
                    lb_ratio[f"{q}:{r['Hp']}"] = r["LB"] / r["Hp"]
                    if not (0.1 <= r["LB"] / r["Hp"] <= 2):
                        ok_lb = False
            if sh == "tri":
                pts = [(log(r["Hp"]), log(r["W"])) for r in rs if r["Hp"] >= 4]
                if len(pts) >= 2:
                    x, y = np.array(pts).T
                    slopes[q] = float(np.polyfit(x, y, 1)[0])
    out["P2b"] = ok_ratio and all(0.7 <= s <= 1.15 for s in slopes.values())
    out["P2b_slopes"] = slopes
    out["P2c"] = ok_lb
    out["LB_over_Hp"] = lb_ratio
    het = {r["q"]: r for r in R if r["part"] == "HET"}
    e1rows = {r["q"]: r for r in _load("e1") if not r["all_residues"]}
    p2d = True
    for q, r in het.items():
        for k in ("HET_farey", "HET_arc"):
            v = r[k] if r[k] is not None else 10 ** 9
            if v < r["sqrtY"]:
                p2d = False
    q7 = max(het)
    p2d = p2d and e1rows[q7]["hstar"] < het[q7]["sqrtY"]
    out["P2d"] = p2d
    out["HET"] = {q: {"sqrtY": r["sqrtY"], "farey": r["HET_farey"], "arc": r["HET_arc"],
                      "hstar_E1": e1rows[q]["hstar"]} for q, r in het.items()}
    return out


def e3():
    R = _load("e3")
    assert sorted(r["n"] for r in R) == [48, 64, 80]
    out = {"P3a": True, "P3b": True, "slopes": {}, "frac64": {}}
    for r in R:
        G = np.array([8, 16, 32, 64])
        f = np.array([r[f"frac_{g}"] for g in G])
        s = float(np.polyfit(np.log(G), np.log(f), 1)[0])
        out["slopes"][r["n"]] = s
        out["frac64"][r["n"]] = r["frac_64"]
        out["P3a"] &= -1.5 <= s <= -1.0
        out["P3b"] &= r["frac_64"] <= 0.02
    return out


if __name__ == "__main__":
    res = {"e1": e1(), "e2": e2(), "e3": e3()}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1, default=str)
    for k, v in res.items():
        print(k, {kk: vv for kk, vv in v.items() if kk != "rows"})
