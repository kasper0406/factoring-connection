"""Scoring of the registered predictions P1-P5 (notes/crtdecode_prereg.md). Floats only in summaries."""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")
B = {"C1": 3000, "C2": 6000}
DECS = ["DT", "DT8", "PH-1.25", "PH-1.5", "RT"]


def e1():
    rows = [json.loads(l) for l in open(os.path.join(RAW, "crtdecode_E1.jsonl"))]
    out = {}
    for cfg in ("C1", "C2"):
        R = [r for r in rows if r["cfg"] == cfg]
        S = [r for r in R if r["n_smooth"] > 0]
        res = {"windows": len(R), "with_smooth": len(S),
               "mean_smooth_per_window": sum(r["n_smooth"] for r in R) / len(R)}
        for d in DECS:
            costs = [r[d]["cost"] if r[d]["found"] else B[cfg] for r in S]
            res[d] = {"capped_mean": sum(costs) / len(costs), "success": sum(r[d]["found"] for r in S) / len(S),
                      "tests": sum(r[d]["tests"] for r in S), "pointers": sum(r[d]["pointers"] for r in S)}
        res["SV_mean"] = sum(r["SV"] for r in S) / len(S)
        # P2: hits at test time, PH (both) and DT8 pooled
        g = defaultdict(lambda: [0, 0])
        for r in S + [r for r in R if r["n_smooth"] == 0]:
            for d in ("PH-1.25", "PH-1.5", "DT8"):
                for k, (n, s) in r[d]["by_hits"].items():
                    key = "1" if k == "1" else ">=2"
                    g[key][0] += n
                    g[key][1] += s
        res["P2_groups"] = {k: {"tested": v[0], "smooth": v[1], "rate": v[1] / v[0] if v[0] else None}
                            for k, v in g.items()}
        dt = res["DT"]["capped_mean"]
        res["P1"] = {d: res[d]["capped_mean"] / dt for d in ("PH-1.25", "PH-1.5", "DT8")}
        res["P1_holds"] = all(res[d]["capped_mean"] >= 0.8 * dt for d in ("PH-1.25", "PH-1.5"))
        a, b = g["1"], g[">=2"]
        res["P2_holds"] = (a[0] >= 30 and b[0] >= 30 and b[1] / b[0] >= 2 * a[1] / a[0]) if (a[0] and b[0]) else None
        res["P3_holds"] = res["DT"]["capped_mean"] < res["RT"]["capped_mean"] < res["SV_mean"]
        out[cfg] = res
    return out


def e2():
    rows = [json.loads(l) for l in open(os.path.join(RAW, "crtdecode_E2.jsonl"))]
    out = {}
    for cfg in ("C1", "C2"):
        R = [r for r in rows if r["cfg"] == cfg]
        pi = R[0]["pi_y"]
        out[cfg] = {"windows": len(R), "pi_y": pi, "exact": sum(r["exact"] for r in R),
                    "depth_max": max(r["depth"] for r in R), "depth_mean": sum(r["depth"] for r in R) / len(R),
                    "groups": R[0]["groups"], "bound": 3 * pi + 1,
                    "P4_holds": all(r["exact"] for r in R), "P5_holds": all(r["depth"] <= 3 * pi + 1 for r in R),
                    "leaf_over_h_mean": sum(r["leaf_incidences"] / r["h"] for r in R) / len(R),
                    "leaf_eq_sieve": all(r["leaf_incidences"] == r["sieve_incidences"] for r in R),
                    "output_ok": all((r["output"] is not None) == (r["n_smooth"] > 0) for r in R)}
    return out


if __name__ == "__main__":
    print(json.dumps({"E1": e1(), "E2": e2()}, indent=1))
