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
