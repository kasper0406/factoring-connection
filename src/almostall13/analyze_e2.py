"""E2 analysis (notes/almostall13_prereg.md): Q1-Q4.

Usage: .venv/bin/python -m src.almostall13.analyze_e2 experiments/raw/almostall13_e2.jsonl \
           experiments/raw/almostall13_e2_summary.json
Floating point is used only for the summary statistics (rates, variances, CIs).
"""
from __future__ import annotations

import json
import random
import sys
from math import sqrt
from typing import Dict, List


def mean(xs):
    return sum(xs) / len(xs)


def var(xs):
    m = mean(xs)
    return sum((x - m) ** 2 for x in xs) / (len(xs) - 1)


def sd_between(ps: List[float], n: int) -> float:
    v = var(ps) - mean([p * (1 - p) for p in ps]) / (n - 1)
    return sqrt(max(0.0, v))


def anova_F(groups: List[List[float]]) -> float:
    allx = [x for g in groups for x in g]
    gm = mean(allx)
    k, n = len(groups), len(allx)
    ssb = sum(len(g) * (mean(g) - gm) ** 2 for g in groups)
    ssw = sum((x - mean(g)) ** 2 for g in groups for x in g)
    return (ssb / (k - 1)) / (ssw / (n - k))


def jensen(rows: List[Dict]) -> Dict:
    pooled = {}
    for k in ("a", "r"):
        num = sum(c["base_d1"] for r in rows for p in r["parents"] for c in p["children"] if c["kind"] == k)
        den = sum(1 for r in rows for p in r["parents"] for c in p["children"] if c["kind"] == k)
        pooled[k] = num / den if den else 0.0
    obs = exp_pool = exp_field = 0.0
    nodes = 0
    for r in rows:
        fr = {}
        for k in ("a", "r"):
            xs = [c["base_d1"] for p in r["parents"] for c in p["children"] if c["kind"] == k]
            fr[k] = mean(xs) if xs else 0.0
        for p in r["parents"]:
            if len(p["children"]) < 2:
                continue
            nodes += 1
            obs += all(c["d1"] for c in p["children"])
            e1 = e2 = 1.0
            for c in p["children"]:
                e1 *= pooled[c["kind"]]
                e2 *= fr[c["kind"]]
            exp_pool += e1
            exp_field += e2
    return {"nodes": nodes, "obs": obs, "exp_pool": exp_pool, "exp_field": exp_field, "pooled_rates": pooled,
            "J_pool": obs / exp_pool if exp_pool else None, "J_field": obs / exp_field if exp_field else None}


def corr(xs, ys):
    mx, my = mean(xs), mean(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = sqrt(sum((x - mx) ** 2 for x in xs))
    sy = sqrt(sum((y - my) ** 2 for y in ys))
    return sxy / (sx * sy) if sx and sy else None


def main():
    rows = [json.loads(l) for l in open(sys.argv[1])]
    fix = [r for r in rows if r["regime"] == "FIX"]
    rnd = [r for r in rows if r["regime"] == "RND"]
    n = rows[0]["n_uniform"]
    out: Dict = {"fields_FIX": len(fix), "fields_RND": len(rnd), "n_uniform": n}
    for key in ("u_alg", "u_rat"):
        pf = [r[key] / n for r in fix]
        pr = [r[key] / n for r in rnd]
        groups = {}
        for r in fix:
            groups.setdefault(r["N_seed"], []).append(r[key] / n)
        gl = [groups[s] for s in sorted(groups)]
        rng = random.Random(f"almostall13|E2|boot|{key}")
        diffs = []
        for _ in range(2000):
            bs_r = [pr[rng.randrange(len(pr))] for _ in range(len(pr))]
            bs_f = []
            for _ in range(len(gl)):
                g = gl[rng.randrange(len(gl))]
                bs_f.extend(g[rng.randrange(len(g))] for _ in range(len(g)))
            diffs.append(mean(bs_r) - mean(bs_f))
        diffs.sort()
        Fobs = anova_F(gl)
        rngp = random.Random(f"almostall13|E2|perm|{key}")
        flat = [x for g in gl for x in g]
        sizes = [len(g) for g in gl]
        ge = 0
        for _ in range(2000):
            rngp.shuffle(flat)
            gs, i = [], 0
            for sz in sizes:
                gs.append(flat[i:i + sz]); i += sz
            ge += anova_F(gs) >= Fobs
        sdr, sdf = sd_between(pr, n), sd_between(pf, n)
        out[key] = {"mean_FIX": mean(pf), "mean_RND": mean(pr), "diff": mean(pr) - mean(pf),
                    "ci95": [diffs[49], diffs[1949]], "sd_btw_FIX": sdf, "sd_btw_RND": sdr,
                    "sd_ratio": sdr / sdf if sdf else None, "range_FIX": [min(pf), max(pf)],
                    "range_RND": [min(pr), max(pr)], "group_means_FIX": [mean(g) for g in gl],
                    "anova_F": Fobs, "anova_perm_p": (1 + ge) / 2001,
                    "corr_roots50_FIX": corr(pf, [r["roots50"] for r in fix]),
                    "corr_roots50_RND": corr(pr, [r["roots50"] for r in rnd])}
    out["jensen_FIX"] = jensen(fix)
    out["jensen_RND"] = jensen(rnd)
    out["time_total_s"] = sum(r["time"] for r in rows)
    a = out["u_alg"]
    out["verdict"] = {
        "Q1": abs(a["diff"]) <= 0.05 and -0.08 <= a["ci95"][0] and a["ci95"][1] <= 0.08,
        "Q2": a["sd_ratio"] is not None and 0.67 <= a["sd_ratio"] <= 1.5,
        "Q3": a["anova_perm_p"] > 0.01,
        "Q4": all(out[j]["J_pool"] is not None and out[j]["J_pool"] > 1.2 and 0.75 <= out[j]["J_field"] <= 1.33
                  for j in ("jensen_FIX", "jensen_RND")),
    }
    json.dump(out, open(sys.argv[2], "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
