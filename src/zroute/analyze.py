"""zroute analysis: aggregates experiments/raw/zroute_e{1,2,3,4}.jsonl against the pre-registered predictions.

Usage: PYTHONPATH=. .venv/bin/python -m src.zroute.analyze > experiments/analysis/zroute_summary.json
"""
from __future__ import annotations

import json
import os
import statistics as st
from collections import defaultdict

RAW = "experiments/raw"


def load(name):
    p = os.path.join(RAW, f"zroute_{name}.jsonl")
    if not os.path.exists(p):
        return []
    return [json.loads(l) for l in open(p) if l.strip()]


def e1(recs):
    out = {"instances": len(recs)}
    stages = [s for r in recs for s in r.get("stages", []) if s.get("M")]
    out["stages"] = len(stages)
    out["P1a_verify_fail"] = sum(s["verify_fail"] for s in stages)
    out["P1a_nonunit"] = sum(s["nonunit"] for s in stages)
    out["P1a_undecided"] = sum(s["undecided"] for s in stages)
    eq = sum(1 for s in stages if s["rank_out"] == s["rank_raw_rational"])
    out["P1b_rank_out_eq_raw"] = f"{eq}/{len(stages)}"
    out["P1c_PIN"] = f"{sum(1 for r in recs if r.get('PIN'))}/{len(recs)}"
    pin_st = [s for s in stages if s["PIN"]]
    out["P1d_miller_given_PIN"] = f"{sum(1 for s in pin_st if s['miller_ok'])}/{len(pin_st)}"
    out["P1d_max_base"] = max((s["miller_base"] or 0) for s in pin_st) if pin_st else None
    out["unit_rank_seen"] = sorted({s["unit_rank_seen"] for s in stages})
    out["factored_any_stage"] = f"{sum(1 for r in recs if r.get('factored_any_stage'))}/{len(recs)}"
    nonpin_ok = [s for s in stages if not s["PIN"]]
    out["miller_ok_without_PIN"] = f"{sum(1 for s in nonpin_ok if s['miller_ok'])}/{len(nonpin_ok)}"
    by = defaultdict(list)
    for r in recs:
        last = [s for s in r.get("stages", []) if s.get("M")]
        if last:
            s = last[-1]
            by[(r["fam"], r["bits"])].append({"n": r["n"], "M": s["M"], "C_alg": s["C_alg"], "rank": s["rank_out"],
                                              "PIN": s["PIN"], "mult": s["mult"], "secs": r["secs"],
                                              "rss": r["rss_mb"], "unpinned": len(s["unpinned"])})
    out["cells"] = {f"{k[0]}{k[1]}": v for k, v in by.items()}
    return out


def e2(recs):
    out = {}
    for d in ("orbit", "random"):
        rs = [r for r in recs if r["design"] == d]
        ftf = [r["rows_to_full"] for r in rs]
        out[d] = {"instances": len(rs), "rows_to_full": ftf, "n": [r["n"] for r in rs],
                  "PIN": sum(1 for r in rs if r["PIN"]), "miller_ok": sum(1 for r in rs if r["miller_ok"]),
                  "miller_ok_given_PIN": sum(1 for r in rs if r["PIN"] and r["miller_ok"])}
    rnd = [r for r in recs if r["design"] == "random"]
    orb = [r for r in recs if r["design"] == "orbit"]
    out["P2a_random_full_at_n"] = f"{sum(1 for r in rnd if r['rows_to_full'] == r['n'])}/{len(rnd)}"
    out["P2b_orbit_gt_3n"] = f"{sum(1 for r in orb if r['rows_to_full'] is None or r['rows_to_full'] > 3 * r['n'])}/{len(orb)}"
    out["orbit_ratio_rows_to_full_over_n"] = [round(r["rows_to_full"] / r["n"], 2) if r["rows_to_full"] else None
                                              for r in orb]
    return out


def e3(recs):
    o = [r["orbit_hits"] for r in recs]
    ra = [r["rand_hits"] for r in recs]
    J = recs[0]["J"] if recs else 1
    def stats(x):
        med = st.median(x)
        return {"mean": st.mean(x), "median": med, "min": min(x), "max": max(x), "var": st.pvariance(x),
                "min_over_median": min(x) / med, "var_over_mean": st.pvariance(x) / st.mean(x)}
    return {"instances": len(recs), "orbit": stats(o), "random": stats(ra),
            "P3a_mean_ratio": st.mean(o) / st.mean(ra), "delta_random": st.mean(ra) / J,
            "P3d_orbit_cycles": sum(1 for r in recs if r["cycle_at"] is not None),
            "corr_orbit_random": _corr(o, ra)}


def _corr(x, y):
    mx, my = st.mean(x), st.mean(y)
    sx = sum((a - mx) ** 2 for a in x) ** 0.5
    sy = sum((b - my) ** 2 for b in y) ** 0.5
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def e4(recs, delta):
    out = {"instances": len(recs)}
    comp = [r for r in recs if r["complete"]]
    out["P4a_dominant_and_full"] = f"{sum(1 for r in comp if r['diag_dominant'] and r['rank'] == r['n'])}/{len(comp)}"
    out["P4b_miller_given_complete"] = f"{sum(1 for r in comp if r['miller_ok'])}/{len(comp)}"
    by = defaultdict(list)
    for r in recs:
        by[r["bits"]].append(r)
    out["P4c_complete_by_bits"] = {b: f"{sum(1 for r in rs if r['complete'])}/{len(rs)}" for b, rs in by.items()}
    out["tries_mean_by_bits"] = {b: round(st.mean(r["tries_mean"] for r in rs if r["tries_mean"]), 2)
                                 for b, rs in by.items()}
    out["tries_max_by_bits"] = {b: max(r["tries_max"] for r in rs if r["tries_max"]) for b, rs in by.items()}
    if delta and 40 in by:
        m = st.mean(r["tries_mean"] for r in by[40] if r["tries_mean"])
        out["P4d_window"] = [0.5 / delta, 2 / delta]
        out["P4d_tries_mean_40"] = m
        out["P4d_pass"] = 0.5 / delta <= m <= 2 / delta
    return out


def main():
    res = {}
    r3 = load("e3")
    res["e3"] = e3(r3) if r3 else None
    delta = res["e3"]["delta_random"] if r3 else None
    for k, f in (("e1", e1), ("e2", e2)):
        r = load(k)
        res[k] = f(r) if r else None
    r4 = load("e4")
    res["e4"] = e4(r4, delta) if r4 else None
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
