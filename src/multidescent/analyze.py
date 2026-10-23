"""Analysis of multidescent E1-E3 (notes/multidescent_prereg.md).  Floating point only for summary statistics.

Run:  .venv/bin/python -m src.multidescent.analyze <E1|E2|E3> <in.jsonl> [<in2.jsonl> ...]
"""
from __future__ import annotations

import json
import math
import random
import statistics as st
import sys
from collections import defaultdict
from typing import Dict, List


def _load(paths: List[str]) -> List[Dict]:
    rows = []
    for p in paths:
        with open(p) as fh:
            rows.extend(json.loads(line) for line in fh if line.strip())
    return rows


def _boot_ci(groups: List[List[float]], stat, B: int = 2000, seed: int = 1):
    """Cluster bootstrap (resample groups) of a statistic of the pooled values."""
    rng = random.Random(seed)
    vals = []
    for _ in range(B):
        pick = [groups[rng.randrange(len(groups))] for _ in groups]
        vals.append(stat([x for g in pick for x in g]))
    vals.sort()
    return vals[int(0.025 * B)], vals[int(0.975 * B) - 1]


# ------------------------------------------------------------------ E1
def analyze_e1(rows: List[Dict], B: int = 2 ** 11, B2: int = 2 ** 16) -> Dict:
    nodes = [r for r in rows if not r.get("summary")]
    summ = [r for r in rows if r.get("summary")]
    lo, hi = int(math.log2(B)), int(math.log2(B2))
    out = {"nodes": len(nodes), "fields": len(summ)}
    from src.multidescent.toy import primes_upto
    primes = [p for p in primes_upto(B2) if p > B]
    for kind in ("alg", "rat"):
        res = {}
        for j in range(lo, hi):
            wlo, whi = 2 ** j, 2 ** (j + 1)
            # 4 sub-buckets of equal log-width
            edges = [wlo * 2 ** (k / 4) for k in range(5)]
            emp = [0] * 4
            ref_unif = [0.0] * 4      # uniform on available labels (per field), pooled over fields
            ref_w = [0.0] * 4         # labels weighted by 1/q
            for r in nodes:
                for q in r[kind]:
                    if wlo < q <= whi:
                        k = min(3, int(4 * math.log2(q / wlo)))
                        emp[k] += 1
            for s in summ:
                nroots = {int(q): k for q, k in s["nroots"].items()}
                for q in primes:
                    if wlo < q <= whi:
                        k = min(3, int(4 * math.log2(q / wlo)))
                        lab = nroots.get(q, 0) if kind == "alg" else 1
                        ref_unif[k] += lab * s["nodes"]
                        ref_w[k] += lab * s["nodes"] / q
            tot = sum(emp)
            if tot == 0:
                continue
            su, sw = sum(ref_unif), sum(ref_w)
            res[f"2^{j}"] = {
                "children": tot,
                "ratio_vs_uniform_labels": [round((emp[k] / tot) / (ref_unif[k] / su), 3) if ref_unif[k] else None
                                            for k in range(4)],
                "ratio_vs_1/q_labels": [round((emp[k] / tot) / (ref_w[k] / sw), 3) if ref_w[k] else None
                                        for k in range(4)],
            }
        out[kind] = res
    return out


# ------------------------------------------------------------------ E2
def analyze_e2(rows: List[Dict]) -> Dict:
    groups = defaultdict(list)
    for r in rows:
        groups[(r["seed"], r["field"], r["ell"])].append(r)
    G = list(groups.values())
    diffs = [[x["p_child"] - x["p_unif"] for x in g] for g in G]
    zc = [[1.0 if x["p_child"] == 0 else 0.0 for x in g] for g in G]
    zu = [[1.0 if x["p_unif"] == 0 else 0.0 for x in g] for g in G]
    zd = [[(1.0 if x["p_child"] == 0 else 0.0) - (1.0 if x["p_unif"] == 0 else 0.0) for x in g] for g in G]
    allr = [x for g in G for x in g]
    mean = lambda v: sum(v) / len(v)
    return {
        "children": len(allr), "parents": len(G),
        "mean_p_child": round(mean([x["p_child"] for x in allr]), 4),
        "mean_p_unif": round(mean([x["p_unif"] for x in allr]), 4),
        "diff": round(mean([d for g in diffs for d in g]), 4),
        "diff_CI": [round(v, 4) for v in _boot_ci(diffs, mean)],
        "zero_frac_child": round(mean([z for g in zc for z in g]), 4),
        "zero_frac_unif": round(mean([z for g in zu for z in g]), 4),
        "zero_diff_CI": [round(v, 4) for v in _boot_ci(zd, mean)],
        "mean_log2_q_child": round(mean([math.log2(x["q"]) for x in allr]), 4),
        "mean_log2_q_unif": round(mean([math.log2(x["qu"]) for x in allr]), 4),
    }


# ------------------------------------------------------------------ E3
def analyze_e3(rows: List[Dict]) -> Dict:
    out = {"parents": len(rows), "viol": sum(r["viol"] for r in rows)}
    strata = {"m<=d-1": [0.0, 0.0, 0], "m>d-1": [0.0, 0.0, 0]}
    nstrata = {"m<=d-1": [0.0, 0.0, 0], "m>d-1": [0.0, 0.0, 0]}
    ex_f, ex_l, pois = [], [], []
    per_parent = []
    for r in rows:
        d = r["d"]
        cs = list(r["counts"].values())
        if not cs or not cs[0]:
            continue
        R = len(cs[0])
        good = [[1 if c > 0 else 0 for c in col] for col in cs]
        J = sum(all(good[i][h] for i in range(len(cs))) for h in range(R)) / R
        Pi = 1.0
        for g in good:
            Pi *= sum(g) / R
        key = "m<=d-1" if len(cs) <= d - 1 else "m>d-1"
        strata[key][0] += J
        strata[key][1] += Pi
        strata[key][2] += 1
        ncs = [c for c in r.get("null_counts", {}).values()]
        Jn = Pin = None
        if ncs and len(ncs) == len(cs) and all(len(c) == R for c in ncs):
            ng = [[1 if c > 0 else 0 for c in col] for col in ncs]
            Jn = sum(all(ng[i][h] for i in range(len(ncs))) for h in range(R)) / R
            Pin = 1.0
            for g in ng:
                Pin *= sum(g) / R
            nstrata[key][0] += Jn
            nstrata[key][1] += Pin
            nstrata[key][2] += 1
        per_parent.append((key, J, Pi, Jn, Pin))
        for col in cs:
            mu = sum(col) / len(col)
            if mu >= 0.5:
                var = sum((c - mu) ** 2 for c in col) / (len(col) - 1)
                ex_f.append(var / mu ** 2 - 1 / mu)
                p0 = sum(1 for c in col if c == 0) / len(col)
                pois.append((p0, math.exp(-mu)))
        for col in r["line"].values():
            if len(col) >= 10:
                mu = sum(col) / len(col)
                if mu >= 0.5:
                    var = sum((c - mu) ** 2 for c in col) / (len(col) - 1)
                    ex_l.append(var / mu ** 2 - 1 / mu)
    for k, (J, Pi, n) in strata.items():
        out[k] = {"parents": n, "sumJ": round(J, 4), "sumPi": round(Pi, 4),
                  "R_joint": round(J / Pi, 4) if Pi > 0 else None}
    # cluster bootstrap over parents for R_joint per stratum
    for k, (J, Pi, n) in nstrata.items():
        out[k]["null_R_joint"] = round(J / Pi, 4) if Pi > 0 else None
    for k in strata:
        grp = [[(J, Pi)] for (kk, J, Pi, _, _) in per_parent if kk == k]
        if len(grp) >= 5:
            def stat(v):
                s1 = sum(a for a, _ in v)
                s2 = sum(b for _, b in v)
                return s1 / s2 if s2 > 0 else float("nan")
            out[k]["R_joint_CI"] = [round(x, 3) for x in _boot_ci(grp, stat)]
    for k in strata:
        grp = [[(J, Pi, Jn, Pin)] for (kk, J, Pi, Jn, Pin) in per_parent if kk == k and Jn is not None]
        if len(grp) >= 5:
            def stat2(v):
                a = sum(x[0] for x in v) / max(1e-12, sum(x[1] for x in v))
                b = sum(x[2] for x in v) / max(1e-12, sum(x[3] for x in v))
                return a / b if b > 0 else float("nan")
            out[k]["alg_over_null"] = round(stat2([x for g in grp for x in g]), 4)
            out[k]["alg_over_null_CI"] = [round(x, 3) for x in _boot_ci(grp, stat2)]
    if out["m<=d-1"].get("R_joint") and out["m>d-1"].get("R_joint"):
        out["ratio_over/under"] = round(out["m>d-1"]["R_joint"] / out["m<=d-1"]["R_joint"], 4)
    out["fiber_excess_relvar_median"] = round(st.median(ex_f), 4) if ex_f else None
    out["fiber_excess_n"] = len(ex_f)
    out["line_excess_relvar_median"] = round(st.median(ex_l), 4) if ex_l else None
    out["line_excess_n"] = len(ex_l)
    if pois:
        rat = [a / b for a, b in pois if b > 0.02]
        out["P0_over_poisson_median"] = round(st.median(rat), 3) if rat else None
        out["P0_over_poisson_n"] = len(rat)
    return out


def main() -> None:
    exp, paths = sys.argv[1], sys.argv[2:]
    rows = _load(paths)
    res = {"E1": analyze_e1, "E2": analyze_e2, "E3": analyze_e3}[exp](rows)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
