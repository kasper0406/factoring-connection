"""Experiment: scaling of the three interval-product baselines (owner factC).

HARNESS script (instance generation knows factors; the algorithms do not).
Usage:
  .venv/bin/python -m src.trackC.exp_baselines --config experiments/configs/trackC_baselines.json
Writes JSONL to experiments/raw/trackC_baselines.jsonl and a fit summary to
experiments/analysis/trackC_baselines_fit.json.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

from src.common.instances import balanced_semiprime
from src.common.metrics import Metrics, append_jsonl
from src.trackC.baselines import bgs_interval, naive_interval, ps_interval

FN = {"naive": naive_interval, "ps": ps_interval, "bgs": bgs_interval}


def fit(xs, ys):
    lx = [math.log2(x) for x in xs]
    ly = [math.log2(max(y, 1e-12)) for y in ys]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    sxx = sum((a - mx) ** 2 for a in lx)
    sxy = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    return sxy / sxx if sxx else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = json.load(open(args.config))
    out = cfg["out"]
    if os.path.exists(out):
        os.remove(out)
    rows = []
    for bits in cfg["bits"]:
        inst = balanced_semiprime(bits, cfg["seed"])
        N = inst.N
        p = inst.labels["factors"][0]          # harness only: to check K < p
        for method, backend, logK_max in cfg["runs"]:
            for e in range(cfg["logK_min"], logK_max + 1):
                K = 1 << e
                if K >= p:
                    continue
                m = Metrics()
                t0 = time.perf_counter()
                if method == "naive":
                    r = FN[method](N, 1, K, m)
                else:
                    r = FN[method](N, 1, K, m, backend=backend)
                dt = time.perf_counter() - t0
                rec = {"exp": "trackC_baselines", "bits": bits, "seed": cfg["seed"], "method": method,
                       "backend": backend, "K": K, "log2K": e, "wall": dt,
                       "nonunit": r["nonunit"], "aux_hit": r["aux_hit"],
                       "counts": r["counts"], "maxima": r["maxima"]}
                append_jsonl(out, rec)
                rows.append(rec)
                print(bits, method, backend, e, f"{dt:.3f}s", r["counts"], flush=True)
    # fits
    summary = {}
    for bits in cfg["bits"]:
        for method, backend, _ in cfg["runs"]:
            sel = [r for r in rows if r["bits"] == bits and r["method"] == method
                   and r["backend"] == backend and r["K"] >= cfg["fit_min_K"]]
            if len(sel) < 3:
                continue
            Ks = [r["K"] for r in sel]
            key = f"{bits}:{method}:{backend}"
            summary[key] = {
                "K_range": [min(Ks), max(Ks)],
                "slope_wall": fit(Ks, [r["wall"] for r in sel]),
                "slope_rn_mul": fit(Ks, [r["counts"].get("rn_mul", 0) + 1 for r in sel]),
                "slope_pmul_model": fit(Ks, [r["counts"].get("pmul_model", 0) + 1 for r in sel])
                if method != "naive" else None,
                "slope_pmul": fit(Ks, [r["counts"].get("pmul", 0) + 1 for r in sel])
                if method != "naive" else None,
                "max_stored_coeffs": max(r["maxima"].get("stored_coeffs", 0) for r in sel),
            }
    os.makedirs(os.path.dirname(cfg["fit_out"]), exist_ok=True)
    json.dump(summary, open(cfg["fit_out"], "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
