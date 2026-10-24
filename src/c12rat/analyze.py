"""Analysis of c12rat E1-E3 against the pre-registered predictions (notes/c12rat_prereg.md).  DIAGNOSTIC.

Usage: python -m src.c12rat.analyze {E1|E2|E3} raw.jsonl summary.json
"""
from __future__ import annotations

import json
import sys

import numpy as np


def load(path):
    with open(path) as f:
        return [json.loads(l) for l in f if l.strip()]


def E1(rows):
    out = {}
    shp = [r for r in rows if r["kind"] == "shape"]
    p1a = True
    p1b = {}
    for r in shp:
        for s in ("0.00390625", "0.015625", "0.0625", "0.25"):
            if r["F_" + s] > r["bound_" + s] + 1e-12:
                p1a = False
        if r["q"] == 100003:
            for s in ("0.015625", "0.0625", "0.25"):
                p1b[s] = r["F_" + s] / float(s)
    out["shape"] = shp
    out["P1a_held"] = p1a
    out["P1b_ratios"] = p1b
    out["P1b_held"] = all(0.85 <= v <= 1.05 for v in p1b.values())
    th = [r for r in rows if r["kind"] == "thin"]
    mult_ok = all(r["mult_max"] <= r["mult_bound"] for r in th)
    tf_ok = all(r["f_%g" % K] <= r["fb_%g" % K] + 1e-12 for r in th for K in (2, 8, 32))
    out["P1c_mult_held"] = mult_ok
    out["P1c_TF_held"] = tf_ok
    out["P1c_held"] = mult_ok and tf_ok
    Ef = {K: float(np.mean([r["f_%g" % K] for r in th])) for K in (2, 8, 32)}
    out["E_f"] = Ef
    out["P1d_ratios"] = [Ef[8] / Ef[2], Ef[32] / Ef[8]]
    out["P1d_held"] = all(0.15 <= x <= 0.45 for x in out["P1d_ratios"])
    sk = [r for r in th if r["sigma"] < 1 / 16]
    share = sum(r["f_32"] for r in sk) / max(sum(r["f_32"] for r in th), 1e-300)
    out["n_labels"] = len(th)
    out["n_skew"] = len(sk)
    out["P1e_share"] = share
    out["P1e_held"] = 0.3 <= share <= 0.9
    # descriptive: thin fraction times K for regular labels
    reg = [r for r in th if r["sigma"] >= 1 / 4]
    out["regular_fK"] = {K: float(np.mean([r["f_%g" % K] * K for r in reg])) for K in (2, 8, 32)}
    return out


def E2(rows):
    L = np.array([r["lhs"] for r in rows], float)
    R = np.array([r["rhs"] for r in rows], float)
    RA = np.array([r["rhs_an"] for r in rows], float)
    ex = np.array([r["expect"] for r in rows], float)
    out = dict(n=len(rows))
    out["P2a_violations"] = int((L > R).sum())
    out["P2a_violations_analytic"] = int((L > RA + 1e-9).sum())
    out["P2a_held"] = out["P2a_violations"] == 0
    ok = ex > 0
    ratio = L[ok] / ex[ok]
    out["P2b_median_count_over_expect"] = float(np.median(ratio))
    out["P2b_held"] = 0.8 <= out["P2b_median_count_over_expect"] <= 1.25
    pos = L > 0
    out["median_loss_rhs_over_count"] = float(np.median(R[pos] / L[pos]))
    out["q10_q90_loss"] = [float(np.quantile(R[pos] / L[pos], 0.1)), float(np.quantile(R[pos] / L[pos], 0.9))]
    return out


def E3(rows):
    out = {}
    allr = [r for r in rows if r["law"] == "all"]
    prim = [r for r in rows if r["law"] == "prim"]
    out["P3a_held"] = all(abs(r["C_nu"]) <= r["B1"] + 1e-12 for r in rows)
    out["P3b_max_absR_all"] = float(max(abs(r["R_lab"]) for r in allr))
    out["P3b_held"] = out["P3b_max_absR_all"] <= 0.015
    dif = [abs(r["R_lab"] - r["R_prim"]) for r in prim]
    out["P3c_max_absdiff_prim"] = float(max(dif))
    out["P3c_held"] = out["P3c_max_absdiff_prim"] <= 0.02
    rat = [r["sd_pair_R"] / max(abs(r["R_lab"]), 1e-4) for r in allr]
    out["P3d_median_ratio"] = float(np.median(rat))
    out["P3d_held"] = out["P3d_median_ratio"] >= 3
    sk = [r["phi_ex"] for r in allr if r["kind"] == "skew"]
    un = [r["phi_ex"] for r in allr if r["kind"] == "uniform"]
    out["P3e_diff"] = float(np.mean(sk) - np.mean(un))
    out["P3e_held"] = out["P3e_diff"] >= -0.02
    keys = ("q", "r", "kind", "sigma", "npts", "m", "phi_ex", "C_nu", "B1", "R_lab", "R_prim", "absR", "thin_ex",
            "sd_pair_R", "sd_pair_cov")
    out["labels"] = [{k: r.get(k) for k in keys} | {"law": r["law"]} for r in rows]
    for law, rr in (("all", allr), ("prim", prim)):
        out["mean_" + law] = {k: float(np.mean([r[k] for r in rr])) for k in ("phi_ex", "C_nu", "B1", "R_lab", "thin_ex")}
    return out


if __name__ == "__main__":
    which, path, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    res = {"E1": E1, "E2": E2, "E3": E3}[which](load(path))
    with open(outp, "w") as f:
        json.dump(res, f, indent=1, default=float)
    show = {k: v for k, v in res.items() if k not in ("labels", "shape")}
    print(json.dumps(show, indent=1, default=float))
