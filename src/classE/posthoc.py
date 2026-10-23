"""HARNESS, post hoc (NOT pre-registered): corrected planted-control expectation, effective sample size and top
statistics for the classE battery (notes/classE_results.md section 5, deviation D1).

The registered V1 expectation rho*sqrt(n) assumes within-stratum target variance 1/4.  For a planted
x = y w.p. (1+rho)/2, E[U] = rho * sum_h S_yy,h and V ~ (1/4) sum_h S_yy,h, so the correct expectation is
z_exp = 2 rho sqrt(sum_h S_yy,h).  For tb under the top-12-bit floor, S_yy,h < n_h/4 because the size of N predicts tb.
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from src.classE.harness import MENU, load, Z_THRESH
from src.eps.analysis import strat_z


def fast_strat_z(x: np.ndarray, y: np.ndarray, h: np.ndarray) -> dict:
    """Same statistic as src/eps/analysis.py::strat_z (vectorised by grouping; strata with n_h < 2 dropped)."""
    _, inv, cnt = np.unique(h, return_inverse=True, return_counts=True)
    keep = cnt[inv] >= 2
    x, y, inv = x[keep], y[keep], inv[keep]
    G = cnt.size
    nh = np.bincount(inv, minlength=G).astype(float)
    mx = np.bincount(inv, x, G) / np.maximum(nh, 1)
    my = np.bincount(inv, y, G) / np.maximum(nh, 1)
    xd, yd = x - mx[inv], y - my[inv]
    U = float(xd @ yd)
    sxx = np.bincount(inv, xd * xd, G)
    syy_ = np.bincount(inv, yd * yd, G)
    ok = nh >= 2
    V = float((sxx[ok] * syy_[ok] / (nh[ok] - 1)).sum())
    if V <= 1e-12 * max(1.0, float(sxx.sum())):
        return {"z": 0.0, "degenerate": True}
    return {"z": U / math.sqrt(V), "degenerate": False}


def syy(y: np.ndarray, h: np.ndarray) -> float:
    _, inv, cnt = np.unique(h, return_inverse=True, return_counts=True)
    G = cnt.size
    my = np.bincount(inv, y, G) / cnt
    yd = y - my[inv]
    return float((yd * yd)[cnt[inv] >= 2].sum())


def main(ns):
    out = {}
    for n in ns:
        path = f"experiments/raw/classE_main_n{n}.jsonl.gz"
        d = load(path)
        import gzip
        Ns = []
        with gzip.open(path, "rt") as fh:
            fh.readline()
            for line in fh:
                Ns.append(int(json.loads(line)["N"]))
        d["Ntop16"] = np.array([x >> (n - 16) for x in Ns], dtype=np.int64)
        d["Njoint"] = np.array([((x >> (n - 8)) << 4) | (x % 16) for x in Ns], dtype=np.int64)
        res = {}
        for tgt, key in (("eps", "Nmod"), ("tb", "Ntop"), ("tb_top16", None), ("tb_joint", None)):
            if tgt == "eps":
                y = (d["eps"] == 1).astype(float)
            else:
                y = d["tb"].astype(float)
            if key is not None:
                h = d[key]
            elif tgt == "tb_top16":
                h = d["Ntop16"]
            else:
                h = d["Njoint"]
            S = syy(y, h)
            neff = 4 * S
            zs = []
            for j, (name, grp, role) in enumerate(MENU):
                if role.startswith("menu"):
                    z = fast_strat_z(d["X"][:, j], y, h)
                    if not z["degenerate"]:
                        zs.append((abs(z["z"]), name, round(z["z"], 2)))
            zs.sort(reverse=True)
            res[tgt] = {"sum_Syy": S, "planted_z_expected_corrected": 2 * 0.05 * math.sqrt(S),
                        "effective_n": neff, "mde_rho_effective": (Z_THRESH + 0.8416) / math.sqrt(neff),
                        "top5": [(nm, z) for _, nm, z in zs[:5]],
                        "n_menu_abs_z_ge_6": sum(1 for a, _, _ in zs if a >= 6)}
        res["path_hit_instances"] = sum(1 for x in d["path_hits"] if x)
        res["path_hits_examples"] = [x for x in d["path_hits"] if x][:5]
        res["nonunit_instances"] = sum(1 for x in d["nonunits"] if x)
        out[n] = res
        print(n, json.dumps(res, indent=None), flush=True)
    with open("experiments/analysis/classE_posthoc.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or [48, 64, 128])
