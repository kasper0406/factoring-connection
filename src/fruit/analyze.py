"""ANALYSIS for E1-E3 exactly as pre-registered in notes/fruit_prereg.md (commit 9f11355).

Usage: python -m src.fruit.analyze <e1.jsonl> <e2.jsonl> [summary.json]
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from math import sqrt


def z_binom(k, n, p):
    if n == 0:
        return 0.0
    return (k - n * p) / sqrt(n * p * (1 - p))


def load(path):
    return [json.loads(l) for l in open(path)]


def e1(recs):
    out = {}
    base = [r for r in recs if r.get("h") == 0]
    nfr = sum(sum(1 for b in r["base_uniform"] if b) for r in base)
    nsm = sum(sum(1 for b in r["base_uniform"] if b is not None) for r in base)
    out["P1_base_fruitful"] = [nfr, nsm]
    hyb = [r for r in recs if r.get("h", 0) > 0 and "r" in r]
    out["accidental"] = sum(1 for r in recs if "accidental_factor" in r)
    cells = defaultdict(lambda: [0, 0])
    pooled = defaultdict(lambda: [0, 0])
    for r in hyb:
        if r.get("forced") is None:
            continue
        c = cells[(r["lev"], r["h"])]
        c[0] += int(r["forced"])
        c[1] += 1
        pooled[r["h"]][0] += int(r["forced"])
        pooled[r["h"]][1] += 1
    out["P2_cells"] = {f"{k[0]}|h={k[1]}": [v[0], v[1], round(z_binom(v[0], v[1], .5), 2)] for k, v in sorted(cells.items())}
    out["P2_pooled"] = {f"h={h}": [v[0], v[1], round(z_binom(v[0], v[1], .5), 2)] for h, v in sorted(pooled.items())}
    obs = exp = var = 0
    for r in hyb:
        if "uniform" not in r:
            continue
        dK, rr = r["dimK"], r["r"]
        pu = (2 ** dK - 2 ** (dK - rr)) / (2 ** dK - 1)
        p = pu / 2
        obs += int(r["uniform"])
        exp += p
        var += p * (1 - p)
    out["P3_uniform"] = [obs, round(exp, 1), round((obs - exp) / sqrt(var), 2) if var else None]
    viol = sum(1 for r in hyb if not (max(0, r["h"] - r["deficit"]) <= r["r"] <= r["h"]))
    out["P4_violations"] = viol
    absorb = defaultdict(lambda: [0, 0])
    for r in hyb:
        if r["h"] == 1:
            a = absorb[(r["bits"], r["lev"])]
            a[0] += int(r["r"] >= 1)
            a[1] += 1
    out["P5_absorb_h1"] = {f"{k[0]}b|{k[1]}": [v[0], v[1], round(v[0] / v[1], 3)] for k, v in sorted(absorb.items())}
    lev_pool = defaultdict(lambda: [0, 0])
    for (b, lev), v in absorb.items():
        lev_pool[lev][0] += v[0]
        lev_pool[lev][1] += v[1]
    out["P5_absorb_h1_pooled"] = {k: [v[0], v[1], round(v[0] / v[1], 3)] for k, v in lev_pool.items()}
    # r >= 1 by h and level, and deficits
    rpos = defaultdict(lambda: [0, 0])
    dA = defaultdict(list)
    newc = defaultdict(list)
    for r in hyb:
        rpos[(r["lev"], r["h"])][0] += int(r["r"] >= 1)
        rpos[(r["lev"], r["h"])][1] += 1
        if r["h"] == 1:
            dA[(r["bits"], r["lev"])].append(r["deficitA"])
            newc[(r["bits"], r["lev"])].append(r["newcols"])
    out["r_pos"] = {f"{k[0]}|h={k[1]}": round(v[0] / v[1], 3) for k, v in sorted(rpos.items())}
    out["deficitA_median"] = {f"{k[0]}b|{k[1]}": sorted(v)[len(v) // 2] for k, v in sorted(dA.items())}
    out["newcols_h1_mean"] = {f"{k[0]}b|{k[1]}": round(sum(v) / len(v), 2) for k, v in sorted(newc.items())}
    # absorption conditional on the descent row adding no new column
    cond = defaultdict(lambda: [0, 0])
    for r in hyb:
        if r["h"] == 1 and r["newcols"] == 0:
            cond[r["lev"]][0] += int(r["r"] >= 1)
            cond[r["lev"]][1] += 1
    out["absorb_h1_given_no_newcol"] = {k: [v[0], v[1]] for k, v in cond.items()}
    return out


def e2(recs):
    out = {}
    out["accidental"] = sum(1 for r in recs if "accidental_factor" in r)
    hyb = [r for r in recs if r.get("kind") == "E2" and r.get("h", 0) > 0]
    cells = defaultdict(lambda: [0, 0])
    pooled = defaultdict(lambda: [0, 0])
    nonsq = att = 0
    for r in hyb:
        for mode in ("forced", "uniform"):
            v = r.get(mode)
            if v is None:
                continue
            att += 1
            nonsq += v == "nonsquare"
        v = r.get("forced")
        if v is None or v == "nonsquare":
            continue
        grp = "tri" if r["omega"] == 3 else "semi"
        c = cells[(r["bits"], grp, r["h"])]
        c[0] += int(v)
        c[1] += 1
        if r["h"] == 1:
            pooled[grp][0] += int(v)
            pooled[grp][1] += 1
    pexp = {"semi": .5, "tri": .75}
    out["Q1_cells"] = {f"{k[0]}b|{k[1]}|h={k[2]}": [v[0], v[1], round(z_binom(v[0], v[1], pexp[k[1]]), 2)]
                       for k, v in sorted(cells.items())}
    out["Q1_pooled_h1"] = {g: [v[0], v[1], round(z_binom(v[0], v[1], pexp[g]), 2)] for g, v in pooled.items()}
    out["Q2_nonsquare"] = [nonsq, att]
    obs = lo = var = 0
    for r in hyb:
        v = r.get("uniform")
        if v is None or v == "nonsquare":
            continue
        p = (1 - 2 ** (-r["r"])) * (1 - 2 ** (1 - r["omega"]))
        obs += int(v)
        lo += p
        var += p * (1 - p)
    out["Q3_uniform"] = [obs, round(lo, 1), round((obs - lo) / sqrt(var), 2) if var else None]
    out["Q4_violations"] = sum(1 for r in hyb if not (max(0, r["h"] - r["deficit"]) <= r["r"] <= r["h"]))
    sb = [r["sqrt_branch"] for r in hyb if r.get("sqrt_branch") is not None]
    out["Q5_sqrt_branch"] = [sum(1 for x in sb if x[0]), len(sb), min((x[1] for x in sb), default=None)]
    ab = defaultdict(lambda: [0, 0])
    for r in hyb:
        if r["h"] == 1:
            ab[r["bits"]][0] += int(r["r"] >= 1)
            ab[r["bits"]][1] += 1
    out["Q6_absorb_h1"] = {f"{k}b": [v[0], v[1], round(v[0] / v[1], 3)] for k, v in ab.items()}
    rpos = defaultdict(lambda: [0, 0])
    for r in hyb:
        rpos[r["h"]][0] += int(r["r"] >= 1)
        rpos[r["h"]][1] += 1
    out["Q6_rpos_by_h"] = {h: round(v[0] / v[1], 3) for h, v in rpos.items()}
    out["Q6_newcols_h1_mean"] = round(sum(r["newcols"] for r in hyb if r["h"] == 1) / max(1, sum(1 for r in hyb if r["h"] == 1)), 2)
    out["Q6_deficitA_median"] = sorted(r["deficitA"] for r in hyb if r["h"] == 1)[len([r for r in hyb if r["h"] == 1]) // 2]
    base0 = [r for r in recs if r.get("kind") == "E2" and r.get("h") == 0 and r.get("base_uniform") is not None]
    out["base_only_uniform_fruitful"] = [sum(int(r["base_uniform"]) for r in base0), len(base0)]
    for kind in ("E3", "E3rand"):
        e3 = [r for r in recs if r.get("kind") == kind and "full_dup_pairs" in r]
        out[f"R1_{kind}_full_dup_pairs"] = sum(r["full_dup_pairs"] for r in e3)
        by = defaultdict(lambda: [0, 0, 0])
        for r in e3:
            by[r["bits"]][0] += r["rat_sq_members"]
            by[r["bits"]][1] += r["M"]
            by[r["bits"]][2] += 1
        out[f"R2_{kind}_ratsq_members"] = {f"{k}b": [v[0], v[1], round(v[0] / v[1], 4), v[2]] for k, v in by.items()}
    secs = {}
    for r in recs:
        if "secs_instance" in r:
            secs[(r["N"], r["bits"])] = r["secs_instance"]
    out["secs_total_instances"] = round(sum(secs.values()), 1)
    return out


def main():
    res = {"E1": e1(load(sys.argv[1])), "E2": e2(load(sys.argv[2]))}
    s = json.dumps(res, indent=1)
    print(s)
    if len(sys.argv) > 3:
        open(sys.argv[3], "w").write(s)


if __name__ == "__main__":
    main()
