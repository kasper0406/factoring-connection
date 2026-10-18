"""Summaries for `unbalanced` experiments (floating point only for fits / display)."""
from __future__ import annotations

import json
import os
from collections import defaultdict
from math import exp
from statistics import mean, median

import numpy as np

from src.common.metrics import read_jsonl

RAW = "experiments/raw"
ANA = "experiments/analysis/unbalanced_summary.json"


def _slope(xs, ys):
    if len(xs) < 2:
        return None
    return float(np.polyfit(np.array(xs, float), np.log2(np.array(ys, float)), 1)[0])


def _bad(R):
    bad = defaultdict(int)
    ncov = 0
    for r in R:
        for b in r["blocks"]:
            lm = b.get("lemma")
            if lm:
                ncov += lm["n_cover"]
                for k, v in lm["bad"].items():
                    bad[k] += v
    return dict(bad), ncov


def summarize():
    out = {}
    p = f"{RAW}/unbalanced_lemma.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        bad, ncov = _bad(R)
        out["lemma"] = {"runs": len(R), "incorrect": sum(not r["correct"] for r in R),
                        "success": sum(r["success"] for r in R), "lemma_bad": bad,
                        "covering_targets": ncov,
                        "certified": sum(r["certified"] for r in R),
                        "certified_not_success": sum(r["certified"] and not r["success"] for r in R)}
    pw, pc = f"{RAW}/unbalanced_oper1_W1.jsonl", f"{RAW}/unbalanced_oper1_C1.jsonl"
    if os.path.exists(pw):
        RW = read_jsonl(pw)
        RC = read_jsonl(pc) if os.path.exists(pc) else []
        cmap = {(r["n"], r["j"], r["seed"]): r for r in RC}
        by = defaultdict(list)
        for r in RW:
            by[(r["n"], r["j"])].append(r)
        tab = {}
        for (n, j), rs in sorted(by.items()):
            cs = [cmap[(n, j, r["seed"])] for r in rs if (n, j, r["seed"]) in cmap]
            bad, ncov = _bad(rs)
            cls = defaultdict(int)
            for r in rs:
                if r["first_block"] is not None:
                    cls[r["blocks"][r["first_block"]]["class"]] += 1
            ccls = defaultdict(int)
            for r in cs:
                if r["first_block"] is not None:
                    ccls[r["blocks"][r["first_block"]]["class"]] += 1
            ent = {"runs": len(rs), "W1_success": sum(r["success"] for r in rs),
                   "W1_incorrect": sum(not r["correct"] for r in rs),
                   "W1_certified": sum(r["certified"] for r in rs),
                   "W1_cert_not_success": sum(r["certified"] and not r["success"] for r in rs),
                   "W1_lemma_bad": bad, "W1_cover_targets": ncov,
                   "W1_med_E_first": median([r["E_first"] for r in rs if r["E_first"]]) if any(r["E_first"] for r in rs) else None,
                   "W1_med_E_full": median([r["E_full"] for r in rs]),
                   "W1_med_E_block0": median([r["E_block0"] for r in rs]),
                   "W1_mean_coverP_block0": mean([r["blocks"][0].get("cover_p", 0) for r in rs]),
                   "W1_first_class": dict(cls)}
            if cs:
                ent.update({"C1_success": sum(r["success"] for r in cs),
                            "C1_incorrect": sum(not r["correct"] for r in cs),
                            "C1_med_E_full": median([r["E_full"] for r in cs]),
                            "C1_med_E_block0": median([r["E_block0"] for r in cs]),
                            "C1_mean_coverP_block0": mean([r["blocks"][0].get("cover_p", 0) for r in cs]),
                            "C1_first_class": dict(ccls),
                            "ratio_E_full_W1_over_C1": median([r["E_full"] for r in rs]) / median([r["E_full"] for r in cs])})
            tab[f"{n},{j}"] = ent
        out["oper1"] = tab
    p = f"{RAW}/unbalanced_oper2.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        by = defaultdict(list)
        for r in R:
            by[(r["n"], r["j"])].append(r)
        tab = {}
        for (n, j), rs in sorted(by.items()):
            bad, ncov = _bad(rs)
            cls = defaultdict(int)
            for r in rs:
                if r["first_block"] is not None:
                    cls[r["blocks"][r["first_block"]]["class"]] += 1
            tab[f"{n},{j}"] = {"runs": len(rs), "success": sum(r["success"] for r in rs),
                               "incorrect": sum(not r["correct"] for r in rs),
                               "certified": sum(r["certified"] for r in rs),
                               "cert_not_success": sum(r["certified"] and not r["success"] for r in rs),
                               "lemma_bad": bad, "cover_targets": ncov,
                               "Y": rs[0]["design_summary"]["Y"], "Delta": rs[0]["design_summary"]["Delta"],
                               "med_E_first": median([r["E_first"] for r in rs if r["E_first"]]) if any(r["E_first"] for r in rs) else None,
                               "med_E_block0": median([r["E_block0"] for r in rs]),
                               "med_E_full": median([r["E_full"] for r in rs]),
                               "med_list_block0": median([r["list_block0"] for r in rs]),
                               "first_class": dict(cls)}
        slopes = {}
        for n in sorted({k[0] for k in by}):
            js = sorted(j for (nn, j) in by if nn == n)
            slopes[n] = {q: _slope(js, [tab[f"{n},{j}"][q] for j in js])
                         for q in ("med_E_block0", "med_E_full", "med_list_block0")}
        out["oper2"] = tab
        out["oper2_slopes_vs_j"] = slopes
    p = f"{RAW}/unbalanced_exc.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        by = defaultdict(list)
        for r in R:
            by[(r["n"], r["j"])].append(r)
        tab = {}
        for (n, j), rs in sorted(by.items()):
            lam0 = mean(r["lam"][0] for r in rs)
            tab[f"{n},{j}"] = {"runs": len(rs), "certified": sum(r["certified"] for r in rs) / len(rs),
                               "block0_uncovered": sum(r["blocks"][0]["n_cover"] == 0 for r in rs) / len(rs),
                               "poisson_block0": exp(-lam0), "mean_lam0": lam0,
                               "block0_clean": sum(r["blocks"][0]["clean"] for r in rs) / len(rs),
                               "all_uncovered": sum(all(b["n_cover"] == 0 for b in r["blocks"]) for r in rs),
                               "full_over_cover": [sum(b["n_full"] for r in rs for b in r["blocks"]),
                                                   sum(b["n_cover"] for r in rs for b in r["blocks"])]}
        out["exc"] = tab
    os.makedirs(os.path.dirname(ANA), exist_ok=True)
    with open(ANA, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1)[:12000])


# ---------------------------------------------------------------------------
# Decision rules R-a .. R-f of notes/unbalanced_prereg.md §4 (added after the freeze;
# evaluation code only -- it reads the frozen raw records and the frozen config).
# ---------------------------------------------------------------------------
def decisions():
    cfg = json.load(open("experiments/configs/unbalanced_experiments.json"))
    pred = cfg["predictions"]
    out = {}
    # R-a
    bad_tot, inc, cns = defaultdict(int), 0, 0
    runs = 0
    for f in ["unbalanced_lemma", "unbalanced_oper1_W1", "unbalanced_oper2"]:
        p = f"{RAW}/{f}.jsonl"
        if not os.path.exists(p):
            continue
        R = read_jsonl(p)
        runs += len(R)
        b, _ = _bad(R)
        for k, v in b.items():
            bad_tot[k] += v
        inc += sum(not r["correct"] for r in R)
        cns += sum(r["certified"] and not r["success"] for r in R)
    pc = f"{RAW}/unbalanced_oper1_C1.jsonl"
    RC = read_jsonl(pc) if os.path.exists(pc) else []
    inc_c1 = sum(not r["correct"] for r in RC)
    out["R-a"] = {"runs": runs, "incorrect": inc, "incorrect_C1": inc_c1, "lemma_bad": dict(bad_tot),
                  "certified_not_success": cns,
                  "holds": inc == 0 and inc_c1 == 0 and sum(bad_tot.values()) == 0 and cns == 0}
    # R-b
    RW = read_jsonl(f"{RAW}/unbalanced_oper1_W1.jsonl") if os.path.exists(f"{RAW}/unbalanced_oper1_W1.jsonl") else []
    byW = defaultdict(list)
    for r in RW:
        byW[(r["n"], r["j"])].append(r)
    succ = {f"{n},{j}": sum(r["success"] for r in rs) for (n, j), rs in sorted(byW.items())}
    RE = read_jsonl(f"{RAW}/unbalanced_exc.jsonl")
    byE = defaultdict(list)
    for r in RE:
        byE[(r["n"], r["j"])].append(r)
    exc = {}
    for (n, j), rs in sorted(byE.items()):
        lam0 = mean(r["lam"][0] for r in rs)
        unc = sum(r["blocks"][0]["n_cover"] == 0 for r in rs) / len(rs)
        cert = sum(r["certified"] for r in rs) / len(rs)
        exc[f"{n},{j}"] = {"certified": cert, "uncovered0": unc, "bound": exp(-lam0) + 0.04,
                           "ok": cert >= 0.99 and unc <= exp(-lam0) + 0.04}
    out["R-b"] = {"W1_success": succ, "exc": exc,
                  "holds": all(v >= 23 for v in succ.values()) and all(e["ok"] for e in exc.values())}
    # R-c, R-d
    byC = defaultdict(list)
    for r in RC:
        byC[(r["n"], r["j"])].append(r)
    ratio, rd = {}, {}
    for n in sorted({k[0] for k in byW}):
        js = sorted(j for (nn, j) in byW if nn == n)
        if (n, 0) in byW and (n, 0) in byC:
            ratio[n] = median([r["E_full"] for r in byW[(n, 0)]]) / median([r["E_full"] for r in byC[(n, 0)]])
        J1 = max(js)
        if (n, 0) in byC and (n, J1) in byC:
            cov = lambda rs: mean(r["blocks"][0].get("cover_p", 0) + r["blocks"][0].get("cover_q", 0) for r in rs)
            c0, cJ = cov(byC[(n, 0)]), cov(byC[(n, J1)])
            sC = sum(r["success"] for r in byC[(n, J1)])
            sW = sum(r["success"] for r in byW[(n, J1)])
            rd[n] = {"C1_cov_j0": c0, "C1_cov_J1": cJ, "ratio": cJ / c0 if c0 else None,
                     "primary_ok": cJ <= 0.3 * c0, "C1_succ_J1": sC, "W1_succ_J1": sW,
                     "secondary_ok": sC < sW}
    J1s = {int(k): v["J1"] for k, v in pred.items() if k != "regime2"}
    rc_each = {n: (2 <= v <= J1s.get(n, 6) + 2) for n, v in ratio.items()}
    drift = ratio[80] / ratio[48] if 80 in ratio and 48 in ratio else None
    out["R-c"] = {"ratio_W1_over_C1_j0": ratio, "each_ok": rc_each, "drift_80_over_48": drift,
                  "holds": all(rc_each.values()) and drift is not None and 0.5 <= drift <= 2}
    out["R-d"] = {"per_n": rd, "holds_primary": all(v["primary_ok"] for v in rd.values()),
                  "holds_secondary": all(v["secondary_ok"] for v in rd.values())}
    # R-e, R-f
    p2 = f"{RAW}/unbalanced_oper2.jsonl"
    if os.path.exists(p2):
        R2 = read_jsonl(p2)
        by2 = defaultdict(list)
        for r in R2:
            by2[(r["n"], r["j"])].append(r)
        slopes, calib, s2 = {}, {}, {}
        for n in sorted({k[0] for k in by2}):
            js = sorted(j for (nn, j) in by2 if nn == n)
            if len(js) >= 2:
                slopes[n] = _slope(js, [median([r["E_block0"] for r in by2[(n, j)]]) for j in js])
        for (n, j), rs in sorted(by2.items()):
            s2[f"{n},{j}"] = sum(r["success"] for r in rs)
            exp_l = pred["regime2"][f"{n},{j}"]["expected_block0_list"]
            med = median([r["list_block0"] for r in rs])
            calib[f"{n},{j}"] = {"median_list": med, "expected": exp_l, "ratio": med / exp_l,
                                 "ok": exp_l / 1.5 <= med <= 1.5 * exp_l}
        e_ok = all(0.08 <= slopes.get(n, -1) <= 0.24 for n in (80, 96))
        refuted = all(slopes.get(n, 1) <= 0.04 for n in (80, 96))
        out["R-e"] = {"slopes_E_block0": slopes, "success": s2,
                      "slope_ok": e_ok, "success_ok": all(v >= 14 for v in s2.values()),
                      "refuted_R16": refuted}
        out["R-f"] = {"calibration": calib, "holds": all(c["ok"] for c in calib.values())}
    path = "experiments/analysis/unbalanced_decisions.json"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    return out
