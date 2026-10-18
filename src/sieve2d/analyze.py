"""sieve2d: evaluation of the frozen predictions P1-P6 (notes/sieve2d_prereg.md).

  PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments summarize

Reads experiments/raw/sieve2d_{exh,samp}.jsonl and writes
experiments/analysis/sieve2d_summary.json.  The decision rules are applied exactly
as frozen; floats are used only for ratios and the (pre-registered) least-squares fit.
"""
from __future__ import annotations

import json
import os
from math import log, sqrt

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
OUT = os.path.join(ROOT, "experiments", "analysis", "sieve2d_summary.json")


def load(name):
    with open(os.path.join(RAW, name)) as fh:
        return [json.loads(l) for l in fh if l.strip()]


def ratio_cells(cells, lo, hi):
    sel = [c for c in cells if c["E_model"] >= 30]
    bad = [(c["n"], c["kind"], c["family"], c["block"], c["H"], c["uncovered"],
            round(c["E_model"], 1), round(c["uncovered"] / c["E_model"], 3))
           for c in sel if not (lo <= c["uncovered"] / c["E_model"] <= hi)]
    so = sum(c["uncovered"] for c in sel)
    se = sum(c["E_model"] for c in sel)
    rs = [c["uncovered"] / c["E_model"] for c in sel]
    return sel, bad, so, se, rs


def wls_slope(x, y, w):
    x, y, w = map(np.asarray, (x, y, w))
    W = w.sum()
    mx, my = (w * x).sum() / W, (w * y).sum() / W
    return float((w * (x - mx) * (y - my)).sum() / (w * (x - mx) ** 2).sum())


def p1(samp, exh):
    out = {}
    sel, bad, so, se, rs = ratio_cells(samp, 0.85, 1.15)
    pooled = so / se
    out["E2"] = {"cells": len(sel), "violations": bad, "pooled": pooled,
                 "ratio_min": min(rs), "ratio_max": max(rs),
                 "pass": (not bad) and 0.95 <= pooled <= 1.05}
    for kind in ("prime", "int"):
        s2, b2, so2, se2, rs2 = ratio_cells([c for c in samp if c["kind"] == kind], .85, 1.15)
        out["E2_" + kind] = {"cells": len(s2), "violations": len(b2), "pooled": so2 / se2}
    sel, bad, so, se, rs = ratio_cells(exh, 0.67, 1.5)
    out["E1_descriptive"] = {"cells": len(sel), "violations": bad, "pooled": so / se,
                             "ratio_min": min(rs), "ratio_max": max(rs),
                             "pass": not bad}
    out["pass"] = out["E2"]["pass"]
    return out


def p2(samp):
    b0 = [c for c in samp if c["block"] == 0 and c["E_model"] < 5]
    so = sum(c["uncovered"] for c in b0)
    se = sum(c["E_model"] for c in b0)
    thr = se + 3 * sqrt(se) + 3
    partA = so <= thr
    big = [c for c in samp if c["lam"] >= 6]
    viol = [(c["n"], c["kind"], c["block"], c["H"], round(c["lam"], 2), c["uncovered"],
             c["pairs"]) for c in big
            if not (c["uncovered"] / c["pairs"] <= 0.1 * 3 / c["lam"] ** 2)]
    return {"A_cells": len(b0), "A_sum_obs": so, "A_sum_E": se, "A_threshold": thr,
            "A_pass": partA, "B_cells": len(big), "B_violations": viol,
            "B_pass": not viol, "pass": partA and not viol}


def p3(samp):
    cells = [c for c in samp if c["block"] == 0 and c["kind"] == "prime"
             and c["uncovered"] >= 10]
    y = [log(c["uncovered"] / c["pairs"]) for c in cells]
    w = [c["uncovered"] for c in cells]
    s_eff = wls_slope([c["lam_eff"] for c in cells], y, w)
    s_nom = wls_slope([c["lam"] for c in cells], y, w)
    return {"cells": len(cells), "slope_lam_eff": s_eff, "slope_lam": s_nom,
            "pass_eff": -1.15 <= s_eff <= -0.85, "pass_nom": -1.3 <= s_nom <= -0.85,
            "pass": (-1.15 <= s_eff <= -0.85) and (-1.3 <= s_nom <= -0.85)}


def p4(samp):
    key = lambda c: (c["n"], c["block"], c["H"])
    pr = {key(c): c for c in samp if c["kind"] == "prime" and c["E_model"] >= 30}
    it = {key(c): c for c in samp if c["kind"] == "int" and c["E_model"] >= 30}
    ks = sorted(set(pr) & set(it))
    ri = sum(it[k]["uncovered"] for k in ks) / sum(it[k]["E_model"] for k in ks)
    rp = sum(pr[k]["uncovered"] for k in ks) / sum(pr[k]["E_model"] for k in ks)
    return {"cells": len(ks), "int_pooled": ri, "prime_pooled": rp, "ratio": ri / rp,
            "pass": 0.9 <= ri / rp <= 1.1}


def p5(samp):
    cells = [c for c in samp if c["lam"] >= 1]
    rows, bad = [], []
    for c in cells:
        rv = c["var_M"] / c["model_var"]
        r4 = c["m4_M"] / c["model_m4"]
        rq = c["Qstat"] / c["Qmodel"] if c["kind"] == "int" else None
        ok = 0.9 <= rv <= 1.1 and 0.8 <= r4 <= 1.2 and (rq is None or 0.8 <= rq <= 1.2)
        rows.append((rv, r4, rq))
        if not ok:
            bad.append((c["n"], c["kind"], c["block"], c["H"], round(c["lam"], 2),
                        round(rv, 3), round(r4, 3), None if rq is None else round(rq, 3)))
    rvs = [r[0] for r in rows]
    r4s = [r[1] for r in rows]
    rqs = [r[2] for r in rows if r[2] is not None]
    return {"cells": len(cells), "var_ratio_range": [min(rvs), max(rvs)],
            "m4_ratio_range": [min(r4s), max(r4s)], "Q_ratio_range": [min(rqs), max(rqs)],
            "violations": bad, "pass": not bad}


def p6(exh):
    cells = [c for c in exh if c["n"] == 30 and c["family"] == "thm" and c["block"] == 0
             and c["kind"] == "prime" and 0.01 <= c["uncovered"] / c["pairs"] <= 0.5]
    rows = [(c["H"], round(c["uncovered"] / c["pairs"], 4), c["qdisp"].get("nq"),
             c["qdisp"].get("D")) for c in cells]
    ok = all(r[3] is not None and 0.5 <= r[3] <= 2 for r in rows)
    return {"cells": len(rows), "rows": rows, "pass": ok and len(rows) > 0}


def tables(samp, exh):
    t2 = [(c["n"], c["kind"], c["block"], c["H"], c["Y"], round(c["lam"], 3),
           round(c["lam_eff"], 3), c["pairs"], c["uncovered"], round(c["E_model"], 2),
           round(c["var_M"] / c["model_var"], 3), round(c["m4_M"] / c["model_m4"], 3))
          for c in samp]
    t1 = [(c["n"], c["kind"], c["family"], c["block"], c["H"], c["Y"], round(c["lam"], 3),
           c["pairs"], c["uncovered"], round(c["E_model"], 2)) for c in exh]
    return t1, t2


def summarize():
    samp = load("sieve2d_samp.jsonl")
    exh = load("sieve2d_exh.jsonl")
    res = {"P1": p1(samp, exh), "P2": p2(samp), "P3": p3(samp), "P4": p4(samp),
           "P5": p5(samp), "P6": p6(exh)}
    res["decision"] = {
        "model_approx_exp": res["P1"]["pass"] and res["P2"]["pass"] and res["P3"]["pass"],
        "primes_vs_integers_lt10pct": res["P4"]["pass"],
        "fourth_moment_main_term_accurate": res["P5"]["pass"]}
    t1, t2 = tables(samp, exh)
    res["table_E1"] = t1
    res["table_E2"] = t2
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1)
    for k in ("P1", "P2", "P3", "P4", "P5", "P6"):
        r = dict(res[k])
        r.pop("rows", None)
        print(k, json.dumps(r)[:1500])
    print("decision", res["decision"])
    return res


if __name__ == "__main__":
    summarize()
