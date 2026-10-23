"""Prediction table for E1/E2 (notes/navigate_prereg.md).  .venv/bin/python -m src.navigate.analyze"""
from __future__ import annotations

import json
import os
from collections import defaultdict

E1 = "experiments/raw/navigate_E1.jsonl"
E2 = "experiments/raw/navigate_E2.jsonl"
OUT = "experiments/analysis/navigate_summary.json"


def e1():
    agg = defaultdict(lambda: defaultdict(float))
    for line in open(E1):
        r = json.loads(line)
        b = r["bits"]
        for alg in ("S0", "G1", "J2"):
            d = r[alg]
            A = agg[(b, alg)]
            for k in ("outputs", "distinct", "positive", "tests", "splits", "null_pairs", "null_succ"):
                A[k] += d.get(k, 0) or 0
            if d.get("mean_tested_logsize") is not None:
                A["ls_sum"] += d["mean_tested_logsize"]
                A["ls_n"] += 1
    out = {}
    for b in sorted({k[0] for k in agg}):
        S0, G1, J2 = agg[(b, "S0")], agg[(b, "G1")], agg[(b, "J2")]
        Y = {k: v["outputs"] / v["tests"] for k, v in (("S0", S0), ("G1", G1), ("J2", J2))}
        Yd = {k: v["distinct"] / v["tests"] for k, v in (("S0", S0), ("G1", G1), ("J2", J2))}

        def R(A):
            # successes per null-tested split / null success rate (same splits)
            return (A["outputs"] / A["null_pairs"]) / (A["null_succ"] / A["null_pairs"]) if A["null_succ"] else None
        row = {
            "Y": Y, "Y_distinct": Yd,
            "J2/S0": Y["J2"] / Y["S0"], "G1/S0": Y["G1"] / Y["S0"],
            "J2/S0_distinct": Yd["J2"] / Yd["S0"], "G1/S0_distinct": Yd["G1"] / Yd["S0"],
            "R_S0": R(S0), "R_G1": R(G1),
            "mean_logsize_S0": S0["ls_sum"] / S0["ls_n"], "mean_logsize_G1": G1["ls_sum"] / G1["ls_n"],
            "splits_per_output": {k: v["splits"] / max(1, v["outputs"]) for k, v in (("S0", S0), ("G1", G1), ("J2", J2))},
            "positive_frac": {k: v["positive"] / max(1, v["outputs"]) for k, v in (("S0", S0), ("G1", G1), ("J2", J2))},
            "outputs": {k: v["outputs"] for k, v in (("S0", S0), ("G1", G1), ("J2", J2))},
        }
        row["P1a"] = 0.85 <= row["R_S0"] <= 1.15
        row["P1b"] = 0.85 <= row["J2/S0"] <= 1.15
        row["P1c"] = row["G1/S0"] >= 1.15
        row["P1d"] = 0.85 <= row["R_G1"] <= 1.15
        out[b] = row
    return out


def e2():
    cells = []
    meta = []
    for line in open(E2):
        r = json.loads(line)
        meta.append({k: r[k] for k in ("seed", "y", "N", "phi", "H", "H_has_nonsquares", "vectors")})
        for c in r["cells"]:
            cells.append({"seed": r["seed"], "y": r["y"], **c})
    fails = defaultdict(list)
    for c in cells:
        mu = c["mu"]
        if 0.1 <= mu <= 5 and abs(c["h"] - c["poisson_h"]) > 0.05:
            fails["P2a"].append(c)
        if mu >= 0.5 and not (0.8 <= c["D"] <= 1.6):
            fails["P2b"].append(c)
        if 0.2 <= mu <= 10 and not (0.2 <= c["f75"] <= 0.5):
            fails["P2c"].append(c)
        if "h_sq" in c and 0.1 <= mu <= 5 and abs(c["h_sq"] - c["h_nonsq"]) > 0.05:
            fails["P2d"].append(c)
    profiles = {}
    for line in open(E2):
        r = json.loads(line)
        profiles[f"y{r['y']}_s{r['seed']}"] = r["walk_profile"]
    return {"meta": meta, "cells": cells,
            "P2a": not fails["P2a"], "P2b": not fails["P2b"], "P2c": not fails["P2c"], "P2d": not fails["P2d"],
            "fails": {k: [(c["seed"], c["y"], c["beta"]) for c in v] for k, v in fails.items()},
            "profiles": profiles}


def main():
    s = {"E1": e1(), "E2": e2()}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(s, open(OUT, "w"), indent=1, default=str)
    for b, row in s["E1"].items():
        print(b, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()
                  if k in ("J2/S0", "G1/S0", "R_S0", "R_G1", "mean_logsize_S0", "mean_logsize_G1", "P1a", "P1b", "P1c",
                           "P1d", "J2/S0_distinct", "G1/S0_distinct")})
        print("   Y", {k: round(v, 4) for k, v in row["Y"].items()}, "splits/out",
              {k: round(v, 1) for k, v in row["splits_per_output"].items()}, "pos", row["positive_frac"])
    E = s["E2"]
    for c in E["cells"]:
        print("E2", c["seed"], c["y"], c["beta"], "mu=%.3f h=%.3f pois=%.3f D=%.3f f75=%.3f max=%d" %
              (c["mu"], c["h"], c["poisson_h"], c["D"], c["f75"], c["max"]),
              ("sq %.3f nonsq %.3f" % (c["h_sq"], c["h_nonsq"])) if "h_sq" in c else "")
    print({k: E[k] for k in ("P2a", "P2b", "P2c", "P2d")}, E["fails"])


if __name__ == "__main__":
    main()
