"""HARNESS analysis of the signsplit runs (uses factors only for diagnostics).

  .venv/bin/python -m src.signsplit.analyze
Writes experiments/analysis/signsplit_analysis.json.
"""
from __future__ import annotations

import json
from collections import defaultdict
from math import log2
from statistics import median

import numpy as np

from src.common.instances import balanced_semiprime

RAW2 = "experiments/raw/signsplit_runs_rerun.jsonl"
OUT = "experiments/analysis/signsplit_analysis.json"
FROZEN = {32: (97, 97, 272), 40: (267, 269, 1080), 48: (779, 789, 4308),
          56: (2297, 2349, 17220), 64: (6797, 6871, 68868)}


def main():
    rows = [json.loads(l) for l in open(RAW2)]
    P = {}
    for r in rows:
        key = (r["bits"], r["seed"])
        if key not in P:
            P[key] = balanced_semiprime(*key).labels["p"]
    out = {}
    # post-hoc: value-based V2 hit classification + deduplicated V2 counts
    from .run import identity_class_v2
    from .hh import v2_dedup_counts
    from multiprocessing import Pool
    v2rows = [r for r in rows if r["variant"] == "V2"]
    for r in v2rows:
        if r.get("operational"):
            r["class_id"] = identity_class_v2({"factor": 1 if r["success"] else None,
                                               "status": r["status"], "hit": r.get("hit")},
                                              r["N"], P[(r["bits"], r["seed"])], r["M"],
                                              r["m0"], r["tR"], r["g"])
    with Pool(24) as pool:
        dd = pool.starmap(v2_dedup_counts, [(r["N"], r["M"], r["m0"], r["tR"], r["g"]) for r in v2rows])
    for r, d in zip(v2rows, dd):
        rows.append({**r, "variant": "V2d", "work": d["work_dedup"], "class_id": None,
                     "right_distinct": d["right_distinct"]})

    # best over menu per instance
    best = defaultdict(lambda: {"V0": None, "V3": None, "V2": None, "V2d": None})
    for r in rows:
        d = best[(r["bits"], r["seed"])]
        if d[r["variant"]] is None if r["variant"] in d else False:
            d[r["variant"]] = r["work"]
        elif r["variant"] in d and r["work"] < d[r["variant"]]:
            d[r["variant"]] = r["work"]
    bits = sorted({b for b, _ in best})
    med = {v: [median(best[(b, s)][v] for s in range(20)) for b in bits] for v in ("V0", "V3", "V2", "V2d")}
    out["median_best_over_V0"] = {v: dict(zip(bits, [median(best[(b, s)][v] / best[(b, s)]["V0"] for s in range(20)) for b in bits])) for v in ("V3", "V2", "V2d")}
    out["median_best_work"] = {v: dict(zip(bits, med[v])) for v in med}
    out["slope_log2work_per_bit"] = {v: float(np.polyfit(bits, [log2(x) for x in med[v]], 1)[0])
                                    for v in med}
    out["slope_upper_half"] = {v: float(np.polyfit(bits[2:], [log2(x) for x in med[v][2:]], 1)[0])
                               for v in med}
    out["frozen_vs_measured"] = {b: {"V0": (FROZEN[b][0], med["V0"][i]),
                                     "V3": (FROZEN[b][1], med["V3"][i]),
                                     "V2": (FROZEN[b][2], med["V2"][i]),
                                     "max_rel_dev": max(abs(med[v][i] / FROZEN[b][k] - 1)
                                                        for k, v in enumerate(("V0", "V3", "V2")))}
                                 for i, b in enumerate(bits)}

    # V1 vs V0 construction savings (matched)
    sav = defaultdict(list)
    for r in rows:
        if r["variant"] != "V1":
            continue
        r0 = next(x for x in rows if x["variant"] == "V0" and x["bits"] == r["bits"]
                  and x["seed"] == r["seed"] and x["M"] == r["M"])
        if r0["status"] == r["status"] == "collision" or True:
            sav[(r["bits"], r["M"])].append((
                r["counts"].get("rn_mul", 0) / max(1, r0["counts"].get("rn_mul", 1)),
                r["counts"].get("lattice_red", 0) / max(1, r0["counts"].get("lattice_red", 1)),
                r["work"] == r0["work"]))
    out["V1_over_V0"] = {f"{b}/{M}": {"median_rn_mul_ratio": median(x[0] for x in v),
                                      "median_lattice_red_ratio": median(x[1] for x in v),
                                      "work_identical": all(x[2] for x in v)}
                         for (b, M), v in sorted(sav.items())}

    # coefficient inflation: |c| bits of the true-branch vector, V3 vs V0 (E1, same m0)
    infl = defaultdict(list)
    for r in rows:
        if r["variant"] == "V3" and r["exp"] == "E1":
            r0 = next(x for x in rows if x["variant"] == "V0" and x["bits"] == r["bits"]
                      and x["seed"] == r["seed"] and x["M"] == r["M"])
            MA = 1
            for m in r["anti"]:
                MA *= m
            infl[(r["M"], tuple(r["anti"]))].append(
                (max(r["abs_c_bits"], r["abs_d_bits"]) - max(r0["abs_c_bits"], r0["abs_d_bits"]),
                 log2(MA)))
    out["V3_coeff_inflation_bits"] = {f"{M}/{A}": {"median_extra_bits": median(x[0] for x in v),
                                                   "predicted_log2_MA": v[0][1]}
                                      for (M, A), v in sorted(infl.items())}

    # hit classification and the order-wrap null expectation
    cls = defaultdict(lambda: defaultdict(int))
    wrap = defaultdict(list)
    for r in rows:
        if r.get("class_id") is None:
            continue
        p = P[(r["bits"], r["seed"])]
        if r["variant"] == "V2":
            span = r["work"] and (r["lambda"] * 2 + 1) * r["M"] ** 2
        else:
            span = (2 * r["lambda"] + 1) * r["M"] ** 2
        cls[(r["bits"], r["variant"])][r["class_id"]] += 1
        wrap[(r["variant"], "span>p" if span > p else "span<=p")].append(r["class_id"] == "ACCIDENTAL_ID")
    # null model for accidental (mod-p only) hits: #pairs tested / p
    nb = defaultdict(list)
    for r in rows:
        if r.get("class_id") in (None, "NONE") or r["variant"] == "V2d":
            continue
        p = P[(r["bits"], r["seed"])]
        if r["variant"] == "V2":
            pairs = r["left"] * r["right"]
        else:
            pairs = r["giant"] * r["baby"]
        x = pairs / p
        b = "<0.01" if x < 0.01 else "0.01-0.1" if x < 0.1 else "0.1-1" if x < 1 else ">=1"
        nb[(r["variant"] if r["variant"] != "V1" else "V0", b)].append(
            (r["class_id"] == "ACCIDENTAL_ID", x))
    out["accidental_rate_by_pairs_over_p"] = {
        f"{v}/{b}": {"n": len(x), "accidental_rate": sum(a for a, _ in x) / len(x),
                     "median_pairs_over_p": median(y for _, y in x)}
        for (v, b), x in sorted(nb.items())}
    # sharper null: Pr[some tested pair collides mod p or mod q by order], using
    # E[1/ord] = c_n/(n) with c_n = sum_{d|n} phi(d)/d for a uniform base (harness:
    # factors p-1, q-1 with sympy).
    from math import exp
    from sympy import divisors, totient
    cache = {}

    def cfun(n):
        if n not in cache:
            cache[n] = sum(totient(d) / d for d in divisors(n))
        return cache[n]
    agg = defaultdict(lambda: [0, 0.0, 0])
    for r in rows:
        if r.get("class_id") in (None, "NONE") or r["variant"] == "V2d":
            continue
        p = P[(r["bits"], r["seed"])]
        q = r["N"] // p
        pairs = r["left"] * r["right"] if r["variant"] == "V2" else r["giant"] * r["baby"]
        lam_ = pairs * (float(cfun(p - 1)) / (p - 1) + float(cfun(q - 1)) / (q - 1))
        a = agg[r["variant"]]
        a[0] += 1
        a[1] += 1 - exp(-lam_)
        a[2] += r["class_id"] == "ACCIDENTAL_ID"
    out["accidental_null_vs_observed"] = {v: {"n": a[0], "expected_P_any_accidental_pair": a[1] / a[0],
                                              "observed_accidental_first_hit": a[2] / a[0]}
                                          for v, a in sorted(agg.items())}
    out["hit_classes"] = {f"{b}/{v}": dict(c) for (b, v), c in sorted(cls.items())}
    out["accidental_rate_by_window_span"] = {f"{v}/{k}": {"n": len(x), "accidental_rate": sum(x) / len(x)}
                                             for (v, k), x in sorted(wrap.items())}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=1, default=str))
