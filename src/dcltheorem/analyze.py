"""Summaries for dcltheorem experiments (floating point only for fits)."""
from __future__ import annotations

import json
import os
from collections import defaultdict
from statistics import median

import numpy as np

from src.common.metrics import read_jsonl

RAW = "experiments/raw"
ANA = "experiments/analysis/dcltheorem_summary.json"


def _slope(xs, ys):
    xs, ys = np.array(xs, float), np.log2(np.array(ys, float))
    return float(np.polyfit(xs, ys, 1)[0]) if len(xs) >= 2 else None


def summarize():
    out = {}
    fr = (48, 80)
    p = f"{RAW}/dcltheorem_lemma.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        bad = defaultdict(int)
        for r in R:
            for b in r["blocks"]:
                for k, v in (b.get("lemma") or {}).get("bad", {}).items():
                    bad[k] += v
        out["lemma"] = {"runs": len(R), "incorrect": sum(not r["correct"] for r in R),
                        "lemma_bad": dict(bad),
                        "certified": sum(r["certified"] for r in R),
                        "certified_not_success": sum(r["certified"] and not r["success"] for r in R),
                        "success": sum(r["success"] for r in R),
                        "covering_targets_clean": sum(b["lemma"]["n_cover"] for r in R for b in r["blocks"] if b.get("lemma") and b["lemma"]["clean"])}
    p = f"{RAW}/dcltheorem_oper.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        by = defaultdict(list)
        for r in R:
            by[(r["design"], r["n"])].append(r)
        tab = {}
        for (d, n), rs in sorted(by.items()):
            bad = defaultdict(int)
            cls = defaultdict(int)
            for r in rs:
                for b in r["blocks"]:
                    for k, v in (b.get("lemma") or {}).get("bad", {}).items():
                        bad[k] += v
                fb = r["first_block"]
                if fb is not None:
                    cls[r["blocks"][fb]["class"]] += 1
            ok = [r for r in rs if r["success"]]
            tab[f"{d}@{n}"] = {
                "runs": len(rs), "success": len(ok), "incorrect": sum(not r["correct"] for r in rs),
                "certified": sum(r["certified"] for r in rs),
                "certified_not_success": sum(r["certified"] and not r["success"] for r in rs),
                "block0_success": sum(bool(r["blocks"][0].get("factor")) for r in rs),
                "first_block_hist": {str(k): sum(r["first_block"] == k for r in rs) for k in range(4)},
                "lemma_bad": dict(bad), "first_class": dict(cls),
                "acc_null_block0_median": median(r["blocks"][0].get("acc_null", 0) for r in rs),
                "med_E_first": median(r["E_first"] for r in ok) if ok else None,
                "med_E_block0": median(r["E_block0"] for r in rs),
                "med_E_full": median(r["E_full"] for r in rs),
                "med_list_block0": median(r["list_block0"] for r in rs),
                "med_list_full": median(r["list_full"] for r in rs),
                "med_wall": median(r["wall_s"] for r in rs),
                "H": rs[0]["design_summary"]["H"], "Y": rs[0]["design_summary"]["Y"],
                "lam": rs[0]["design_summary"]["lam"]}
        out["oper"] = tab
    p = f"{RAW}/dcltheorem_hh.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        by = defaultdict(list)
        for r in R:
            by[r["n"]].append(r)
        out["hh"] = {str(n): {"runs": len(rs), "success": sum(r["success"] for r in rs),
                              "status": {s: sum(r["status"] == s for r in rs) for s in {r["status"] for r in rs}},
                              "med_E_collision": median([r["E"] for r in rs if r["status"] == "collision"] or [0]),
                              "med_E_all": median(r["E"] for r in rs),
                              "list": median(r["list"] for r in rs),
                              "med_wall": median(r["wall_s"] for r in rs)}
                     for n, rs in sorted(by.items())}
    # slopes
    sl = {}
    if "oper" in out:
        for d in ("C1", "C2", "C3"):
            for key in ("med_E_first", "med_E_block0", "med_E_full", "med_list_block0", "med_list_full"):
                pts = [(int(k.split("@")[1]), v[key]) for k, v in out["oper"].items()
                       if k.startswith(d + "@") and v[key] and fr[0] <= int(k.split("@")[1]) <= fr[1]]
                if len(pts) >= 2:
                    sl[f"{d}:{key}"] = _slope(*zip(*pts))
    if "hh" in out:
        for key in ("med_E_collision", "med_E_all", "list"):
            pts = [(int(n), v[key]) for n, v in out["hh"].items() if v[key] and fr[0] <= int(n) <= fr[1]]
            sl[f"HH:{key}"] = _slope(*zip(*pts))
    out["slopes_48_80"] = sl
    if "oper" in out and "hh" in out:
        out["ratio_C1_first_over_HH"] = {
            n: out["oper"][f"C1@{n}"]["med_E_first"] / out["hh"][n]["med_E_collision"]
            for n in out["hh"] if f"C1@{n}" in out["oper"] and out["hh"][n]["med_E_collision"]}
    p = f"{RAW}/dcltheorem_exc.jsonl"
    if os.path.exists(p):
        R = read_jsonl(p)
        by = defaultdict(list)
        for r in R:
            by[r["n"]].append(r)
        ex = {}
        for n, rs in sorted(by.items()):
            ncov = sum(b["n_cover"] for r in rs for b in r["blocks"])
            nfull = sum(b["n_full"] for r in rs for b in r["blocks"])
            ex[str(n)] = {"pairs": len(rs), "certified": sum(r["certified"] for r in rs) / len(rs),
                          "block0_clean": sum(r["blocks"][0]["clean"] for r in rs) / len(rs),
                          "block0_cert": sum(r["blocks"][0]["cert"] for r in rs) / len(rs),
                          "block0_uncovered": sum(r["blocks"][0]["n_cover"] == 0 for r in rs) / len(rs),
                          "all_blocks_uncovered": sum(all(b["n_cover"] == 0 for b in r["blocks"]) for r in rs) / len(rs),
                          "full_over_cover": nfull / max(1, ncov),
                          "ord_q2_bits_median": median(r["ord_q2_bits"] for r in rs),
                          "q_bits_median": median(r["q_bits"] for r in rs)}
        out["exc"] = ex
    os.makedirs(os.path.dirname(ANA), exist_ok=True)
    with open(ANA, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
