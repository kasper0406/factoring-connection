"""Evaluate the pre-registered predictions (notes/dsr_prereg.md) on experiments/raw/dsr_e{1,2}.jsonl.
  PYTHONPATH=. .venv/bin/python -m src.dsr.analyze > experiments/analysis/dsr_summary.json
"""
from __future__ import annotations

import json
import math
from statistics import median

import numpy as np


def load(path):
    return [json.loads(line) for line in open(path)]


def e1(rows):
    out = {}
    rA = {s: [r["A1"] / np.mean(r["ctrl_A1"]) for r in rows if r["set"] == s] for s in ("S_can", "S_bal")}
    rI = {s: [r["beta"] / np.mean(r["ctrl_beta"]) for r in rows if r["set"] == s] for s in ("S_can", "S_bal")}
    out["P1a_rA_range"] = {s: [round(min(v), 3), round(max(v), 3)] for s, v in rA.items()}
    out["P1a_pass"] = all(0.7 <= x <= 1.3 for v in rA.values() for x in v)
    pr = [r for r in rows if r["fam"] == "prime" and r["set"] == "S_can"]
    slope = float(np.polyfit(np.log([r["size"] for r in pr]), np.log([r["A1"] for r in pr]), 1)[0])
    out["P1b_slope_S_can"] = round(slope, 4)
    out["P1b_pass"] = 0.4 <= slope <= 0.6
    frac = sum(x > 1 for x in rI["S_can"]) / len(rI["S_can"])
    out["P1c_frac_rinf_gt1"] = frac
    out["P1c_pass"] = bool(frac >= 0.75)
    out["rinf_range"] = {s: [round(min(v), 2), round(max(v), 2)] for s, v in rI.items()}
    out["parseval_max_err"] = max(max(r["parseval_rel_err"], r["ctrl_parseval"]) for r in rows)
    return out


def e2(rows):
    out = {}
    key = lambda r: (r["bits"], r["seed"], r["set"])  # noqa: E731
    by = {}
    for r in rows:
        by.setdefault(key(r), {})[r["list"]] = r
    nus = [r["nu"] for r in rows if r["list"] == "R_G" and r["nu"] is not None]
    nR = sum(r["list"] == "R_G" for r in rows)
    out["P2a_frac_in_range"] = sum(0.8 <= x <= 1.5 for x in nus) / nR
    out["P2a_pass"] = out["P2a_frac_in_range"] >= 0.9
    out["P2a_nu_range"] = [round(min(nus), 3), round(max(nus), 3)]

    def ratio(lst, setname):
        v = []
        for k, d in by.items():
            if k[2] != setname or lst not in d:
                continue
            a, b = d[lst]["m_star"], d["R_G"]["m_star"]
            v.append(math.inf if a is None else a / b)
        return v

    rw = ratio("R_walk", "S_can")
    out["P2b_median_ratio"] = round(median(rw), 3)
    out["P2b_pass"] = 0.7 <= median(rw) <= 1.5
    cens = {f"{lst}/{s}": sum(r["m_star"] is None for r in rows if r["list"] == lst and r["set"] == s)
            for lst, s in (("D_box", "S_can"), ("D_slide", "S_can"), ("D_box", "S_bal"))}
    out["P2c_censored_counts"] = cens
    out["P2c_pass"] = sum(cens.values()) == 0
    rb = ratio("D_box", "S_can")
    out["P2d_median_ratio"] = median(rb)
    out["P2d_pass"] = median(rb) >= 1.5
    rs = ratio("D_slide", "S_can")
    out["P2e_median_ratio"] = median(rs)
    out["P2e_pass"] = 0.7 <= median(rs) <= 1.5
    nb = {b: [r["nu"] if r["nu"] is not None else math.inf for r in rows
              if r["list"] == "D_box" and r["set"] == "S_can" and r["bits"] == b] for b in (20, 22, 24)}
    out["P2f_median_nu_D_box"] = {b: median(v) for b, v in nb.items()}
    out["P2f_pass"] = median(nb[24]) <= 2 * median(nb[20])
    # descriptive: DSR targets only
    desc = {}
    for lst in ("R_G", "R_walk", "D_box", "D_slide", "D_orbit"):
        for s in ("S_can", "S_bal"):
            v = [r for r in rows if r["list"] == lst and r["set"] == s]
            if not v:
                continue
            desc[f"{lst}/{s}"] = {
                b: {"median_T*": median([r["T_star"] if r["T_star"] is not None else math.inf
                                         for r in v if r["bits"] == b]),
                    "unhit_total": sum(r["targets_unhit"] for r in v if r["bits"] == b),
                    "mean_first_hit": round(float(np.mean([r["target_hits_mean"] for r in v if r["bits"] == b])), 2),
                    "median_delta": round(median([r["delta"] for r in v if r["bits"] == b]), 4)}
                for b in (20, 22, 24)}
    out["descriptive_DSR_targets"] = desc
    return out


def main():
    res = {"E1": e1(load("experiments/raw/dsr_e1.jsonl")), "E2": e2(load("experiments/raw/dsr_e2.jsonl"))}
    print(json.dumps(res, indent=1, default=lambda x: "inf" if x == math.inf else str(x)))


if __name__ == "__main__":
    main()
