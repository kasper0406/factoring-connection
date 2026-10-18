"""Summaries of the step-4 experiments (E1-E3) with null expectations.

  .venv/bin/python -m src.cycloharvey.analyze

Null for ACCIDENTAL hits: a (baby, giant) pair collides modulo a prime l
by chance with probability ~1/l (E2: projections lambda; E1: ~1/ord_l(alpha)
<= ... we use 1/l as the order-of-magnitude null), so the expected number of
accidental collisions is ~ baby * giant * (1/p + 1/q).  N bits only (the
null uses sqrt(N) ~ p ~ q, so no factor is needed).
"""
from __future__ import annotations

import collections
import json
import math
import os

import numpy as np

from src.common.metrics import read_jsonl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
OUT = os.path.join(ROOT, "experiments", "analysis", "cycloharvey_experiments.json")


def main():
    out = {}
    e1 = read_jsonl(os.path.join(RAW, "cycloharvey_e1_harvey.jsonl"))
    by = collections.defaultdict(list)
    for r in e1:
        by[r["bits"]].append(r)
    t1, xs, ys = {}, [], []
    for b, rs in sorted(by.items()):
        cost = np.mean([r["m"] + r["model_giant"] for r in rs])
        mul = np.mean([r["counts_full"].get("rn_mul", 0) + r["counts_full"].get("rn_mul_pow", 0) for r in rs])
        null = np.mean([2 * r["m"] * r["model_giant"] / 2 ** (b / 2) for r in rs])
        t1[b] = {"n": len(rs), "success": sum(r["success"] for r in rs),
                 "class": dict(collections.Counter(r["class"] for r in rs)),
                 "status": dict(collections.Counter(r["status"] for r in rs)),
                 "impl_giant_equals_model": sum(r["impl_giant"] == r["model_giant"] for r in rs),
                 "mean_list_cost": float(cost), "mean_mults": float(mul),
                 "mean_gcd": float(np.mean([r["counts_full"].get("gcd", 0) for r in rs])),
                 "null_accidental_per_instance": float(null)}
        xs.append(b)
        ys.append(math.log2(cost))
    out["E1"] = t1
    out["E1_slope_listcost_per_bit"] = float(np.polyfit(xs, ys, 1)[0])
    out["E1_slope_listcost_per_bit_40plus"] = float(np.polyfit(xs[2:], ys[2:], 1)[0])

    e2 = read_jsonl(os.path.join(RAW, "cycloharvey_e2_cyclo.jsonl"))
    by = collections.defaultdict(list)
    for r in e2:
        by[(r["bits"], r["m"])].append(r)
    t2 = {}
    for (b, m), rs in sorted(by.items()):
        t2[f"{b}|m={m}"] = {
            "n": len(rs), "success": sum(r["success"] for r in rs),
            "class": dict(collections.Counter(r["class"] for r in rs)),
            "D": rs[0]["D"],
            "giant_ratio_vs_Fp": float(sum(r["giant"] or 0 for r in rs) / sum(r["harvey_giant"] or 1 for r in rs)),
            "R_mul_ratio_vs_Fp_mul": float(sum(r["R_mul"] for r in rs) / sum(r["harvey_rn_mul"] for r in rs)),
            "base_ops_lower_ratio": float(sum(r["R_base_ops_lower"] for r in rs) / sum(r["harvey_rn_mul"] for r in rs)),
            "null_accidental_per_instance": float(np.mean([2 * r["mb"] * (r["giant"] or 0) / 2 ** (b / 2) for r in rs])),
        }
    out["E2"] = t2

    e3 = read_jsonl(os.path.join(RAW, "cycloharvey_e3_onesided.jsonl"))
    t3 = collections.defaultdict(lambda: [0, 0, 0])
    for r in e3:
        k = f"{r['F']}|p==q mod m:{r['p_eq_q_mod_m']}|m={r['m']}"
        t3[k][0] += 1
        t3[k][1] += r["identity_mod_p"]
        t3[k][2] += r["also_mod_q"]
    out["E3"] = {k: {"n": v[0], "identity_mod_p": v[1], "two_sided": v[2]} for k, v in sorted(t3.items())}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
