"""Score the pre-registered predictions P1-P9 (notes/abcsunit_prereg.md) from the raw JSONL."""
from __future__ import annotations

import json
import math
import os
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
DESIGN = {(24, 1.25, "9/8", 0), (24, 1.25, "5/4", 0), (24, 1.5, "9/8", 0), (24, 1.5, "5/4", 0)}


def load(name):
    with open(os.path.join(RAW, name)) as f:
        return [json.loads(l) for l in f if l.strip()]


def key(r):
    return (r["n"], r["kappa"], r["K"], r["seed"])


def frac_X_over_N(r):
    return r["X"] // r["N"]


def score(rows, exclude_design=False):
    rows = [r for r in rows if "skipped" not in r]
    if exclude_design:
        rows = [r for r in rows if key(r) not in DESIGN]
    out = {}
    s1 = [r for r in rows if r["M1"] >= 100]
    bad1 = [(key(r), round(r["R1"], 3)) for r in s1 if not (0.6 <= r["R1"] <= 1.5)]
    med1 = statistics.median(r["R1"] for r in s1)
    out["P1"] = (not bad1 and 0.8 <= med1 <= 1.2, {"n_runs": len(s1), "median_R1": round(med1, 3),
                                                     "out_of_band": bad1,
                                                     "min_R1": round(min(r["R1"] for r in s1), 3),
                                                     "max_R1": round(max(r["R1"] for r in s1), 3)})
    s1p = [r for r in rows if r["K"] == "9/8"]
    bad1p = [(key(r), round(r["R0"], 3)) for r in s1p if not r["R0"] < 0.85]
    out["P1'"] = (not bad1p, {"n_runs": len(s1p), "violations": bad1p,
                              "max_R0": round(max(r["R0"] for r in s1p), 3)})
    bad1b = [(key(r), r["R1b"] - r["R1"]) for r in rows if abs(r["R1b"] - r["R1"]) > 0.05]
    out["P1b"] = (not bad1b, {"max_abs_diff": max(abs(r["R1b"] - r["R1"]) for r in rows), "violations": bad1b})
    out["P2"] = (all(r["H_eq_phi"] for r in rows), {"n_runs": len(rows)})
    out["P3"] = (all(r["propAB_holds"] for r in rows), {})
    over = [key(r) for r in rows if r["rank"] > r["rank_pred"]]
    full = [r for r in rows if r["forest_size"] >= r["R_target"]]
    miss = [(key(r), r["rank"], r["rank_pred"]) for r in full if r["rank"] != r["rank_pred"]]
    out["P4"] = (not over and not miss, {"runs_with_R_relations": len(full), "rank_misses": miss,
                                         "rank_exceeds_bound": over,
                                         "runs_short_of_R": [(key(r), r["forest_size"], r["R_target"], r["rank"],
                                                              r["rank_pred"]) for r in rows
                                                             if r["forest_size"] < r["R_target"]]})
    reach = [r for r in rows if r["rank"] == r["rank_pred"]]
    tf = sum(r["fruitful"] for r in reach)
    tt = sum(r["deps_tested"] for r in reach)
    zbad = [(key(r), r["z_image_dim"], r["fruitful"]) for r in reach if r["z_image_dim"] != 2 or r["fruitful"] == 0]
    out["P5"] = (0.40 <= tf / tt <= 0.60 and not zbad, {"runs": len(reach), "pooled_fruitful": round(tf / tt, 3),
                                                       "fruitful": tf, "tested": tt, "violations": zbad})
    s6 = [r for r in rows if r["M1_smoothk"] >= 100 and frac_X_over_N(r) > r["y"]]
    bad6 = [(key(r), round(r["R3"], 3)) for r in s6 if not (0.6 <= r["R3"] <= 1.5)]
    out["P6"] = (not bad6, {"n_runs": len(s6), "violations": bad6,
                            "R3_values": [(key(r), round(r["R3"], 3)) for r in s6]})
    bad7 = [key(r) for r in rows if r["coprime_pairs"] and r["q_gt1eps"] / r["coprime_pairs"] > 1e-3]
    v7 = [key(r) for r in rows if r["abcv_violations"]]
    out["P7"] = (not bad7 and not v7, {"frac_violations": bad7, "abcv_violations": v7,
                                       "total_q_gt1": sum(r["q_gt1"] for r in rows),
                                       "total_q_gt1eps": sum(r["q_gt1eps"] for r in rows),
                                       "max_q": round(max(r["qmax"] for r in rows), 4),
                                       "coprime_pairs_total": sum(r["coprime_pairs"] for r in rows)})
    s8 = [r for r in rows if r["M1_smoothk"] >= 100]
    bad8 = [(key(r), r["triply"], r["pi_y"]) for r in s8 if r["triply"] < r["pi_y"]]
    # upper half (added after review): registered loose bound 3*7^(2(s+2)+3) and the proved 3*7^(2s+3)
    up_loose = [key(r) for r in rows if math.log10(max(r["triply"], 1)) > r["evertse_bound_log10"]]
    up_tight = [key(r) for r in rows
                if math.log10(max(r["triply"], 1)) > math.log10(3) + (2 * r["pi_y"] + 3) * math.log10(7)]
    bad8 = bad8 + [("upper", k) for k in up_loose + up_tight]
    out["P8"] = (not bad8, {"n_runs": len(s8), "violations": bad8,
                            "min_triply_over_pi": round(min(r["triply"] / r["pi_y"] for r in s8), 2)})
    return out


def score_E4(rows):
    ok = True
    by = {}
    for r in rows:
        by.setdefault(r["nprimes"], []).append(r)
    for n, rs in by.items():
        rs.sort(key=lambda r: r["Hbits"])
        c = [r["count"] for r in rs]
        ok &= all(a <= b for a, b in zip(c, c[1:])) and all(x <= rs[0]["vkm"] for x in c)
        ok &= rs[-1]["Hbits"] == 40 and rs[-1]["count"] == rs[-1]["vkm"]
    return ok


if __name__ == "__main__":
    rows = load("abcsunit_E123.jsonl")
    for excl in (False, True):
        print("=== exclude design runs" if excl else "=== all runs")
        for k, (ok, info) in score(rows, excl).items():
            print(k, "HOLDS" if ok else "FAILS", json.dumps(info, default=str)[:900])
    print("P9", "HOLDS" if score_E4(load("abcsunit_E4.jsonl")) else "FAILS")
