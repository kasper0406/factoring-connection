"""uprime: analysis of E1 (hybrid BV toy) and E2 (mixed acceptance, conditional child mean).  Statistics only
(floats allowed for summaries and bootstrap CIs, AGENT_RULES rule 1).

Run: .venv/bin/python -m src.uprime.analyze E1 <e1.jsonl>     |    ... E2 <e2.jsonl>
"""
from __future__ import annotations

import json
import math
import random
import sys
from typing import Dict, List

from src.uprime import hybrid_bv


def analyze_e1(path: str) -> Dict:
    s = hybrid_bv.summarize(path)
    cells = s["cells"]
    s["P1a"] = all(c["bad_F1"] == 0 for c in cells)
    s["P1b"] = all(0.8 <= c["H_ratio"] <= 1.25 for c in cells)
    s["P1c"] = all(0.8 <= c["medV_ratio"] <= 1.25 for c in cells)
    s["P1d"] = all(0.3 <= c["medz_F1"] <= 1.5 for c in cells)
    return s


def _boot(groups: List[List[float]], stat, B: int = 2000, seed: int = 7):
    rng = random.Random(seed)
    vals = []
    for _ in range(B):
        sample = [g for _ in groups for g in [groups[rng.randrange(len(groups))]]]
        vals.append(stat(sample))
    vals.sort()
    return vals[int(0.025 * B)], vals[int(0.975 * B) - 1]


def analyze_e2(path: str) -> Dict:
    rows = [json.loads(l) for l in open(path)]
    hdr = rows[0]
    S2 = hdr["cfg"]["S2"]
    boxes = [r for r in rows[1:] if r.get("kind") == "box" and r.get("n", 0) > 0]
    kids = [r for r in rows[1:] if r.get("kind") == "child"]
    both = sum(b["both"] for b in boxes)
    pred = sum(b["rat"] * b["alg2"] / b["n"] for b in boxes)
    removed = sum(b["removed_sq"] + b["removed_disc"] for b in boxes)
    out = {"boxes": len(boxes), "points": sum(b["n"] for b in boxes), "rat": sum(b["rat"] for b in boxes),
           "alg2": sum(b["alg2"] for b in boxes), "both": both, "acc": sum(b["acc"] for b in boxes),
           "removed_sq": sum(b["removed_sq"] for b in boxes), "removed_disc": sum(b["removed_disc"] for b in boxes),
           "removed_N": sum(b["removed_N"] for b in boxes), "sqdiv_any_among_alg2": sum(b["sqdiv_any"] for b in boxes),
           "product_ratio": both / pred if pred else float("nan"),
           "removed_frac": removed / both if both else float("nan")}
    # product-rule ratio CI by field clusters
    fields: Dict = {}
    for b in boxes:
        fields.setdefault((b["seed"], b["field"]), []).append(b)
    fl = list(fields.values())

    def pr(sample):
        num = sum(b["both"] for g in sample for b in g)
        den = sum(b["rat"] * b["alg2"] / b["n"] for g in sample for b in g)
        return num / den if den else float("nan")

    out["product_ratio_CI_fields"] = _boot(fl, pr)
    # conditional child mean vs matched null
    kk = [k for k in kids if k["null"] is not None]
    out["children"] = len(kids)
    out["children_with_null"] = len(kk)
    mc = sum(k["cond"]["good"] for k in kk) / len(kk)
    mn = sum(k["null"]["good"] for k in kk) / len(kk)
    out["mean_cond"], out["mean_null"], out["R_CM"] = mc, mn, mc / mn
    parents: Dict = {}
    for k in kk:
        parents.setdefault((k["seed"], k["field"], k["ell"]), []).append(k)
    pl = list(parents.values())

    def rcm(sample):
        num = sum(k["cond"]["good"] for g in sample for k in g)
        den = sum(k["null"]["good"] for g in sample for k in g)
        return num / den if den else float("nan")

    out["R_CM_CI_parents"] = _boot(pl, rcm)
    fk: Dict = {}
    for k in kk:
        fk.setdefault((k["seed"], k["field"]), []).append(k)
    out["R_CM_CI_fields"] = _boot(list(fk.values()), rcm)
    # zero probability vs Poisson
    z_c = sum(1 for k in kk if k["cond"]["good"] == 0) / len(kk)
    out["P0_cond"], out["P0_poisson_cond"] = z_c, math.exp(-mc)
    out["P0_ratio"] = z_c / math.exp(-mc)
    out["irregular_frac"] = sum(1 for k in kids if k["cond"]["irregular"]) / len(kids)
    out["irregular_bound"] = 3 / S2 ** 2
    # predictions
    out["P2a"] = 0.8 <= out["product_ratio"] <= 1.3
    out["P2b"] = out["removed_frac"] <= 0.02
    lo, hi = out["R_CM_CI_parents"]
    out["P2c"] = (0.8 <= out["R_CM"] <= 1.5) and lo >= 0.6 and hi <= 2.0
    out["P2d"] = out["irregular_frac"] <= out["irregular_bound"]
    out["P2e"] = 0.8 <= out["P0_ratio"] <= 1.3
    return out


if __name__ == "__main__":
    which, path = sys.argv[1], sys.argv[2]
    res = analyze_e1(path) if which == "E1" else analyze_e2(path)
    print(json.dumps(res, indent=1, default=str))
