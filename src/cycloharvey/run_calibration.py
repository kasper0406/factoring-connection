"""Run the calibration gate: exact model costs for the five calibrated
algorithms over a range of N sizes, plus exact rational LP exponents.

  .venv/bin/python -m src.cycloharvey.run_calibration

Writes experiments/raw/cycloharvey_calibration.jsonl and
experiments/analysis/cycloharvey_calibration.json.  N is a size parameter
only (N = 2^n + 2^(n/2) + 1-ish odd number); no factorisation is involved.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from src.common.metrics import append_jsonl
from .calibration import (best_harvey, best_hh, best_hittmeir, best_lehman,
                          cost_bsgs_sum, exponent_lp)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw", "cycloharvey_calibration.jsonl")
ANA = os.path.join(ROOT, "experiments", "analysis", "cycloharvey_calibration.json")
CFG = os.path.join(ROOT, "experiments", "configs", "cycloharvey_calibration.json")

CONFIG = {
    "sizes": {"bsgs_sum": list(range(40, 161, 8)),
              "lehman": list(range(30, 67, 4)),
              "hittmeir": list(range(40, 97, 8)),
              "harvey": list(range(40, 121, 8)),
              "hh": list(range(40, 129, 8))},
    "caps": "r, eta <= 2^22 (sizes chosen so the optimum is below the cap)",
    "N_of_n": "N = 2^n + 2^(n//2) + 1 (size parameter only)",
    "fit": "least squares of log2(cost) on log2(N) over the upper half of sizes",
}


def N_of(n: int) -> int:
    return (1 << n) + (1 << (n // 2)) + 1


def summarize():
    """Recompute the summary from the raw JSONL (used when the run is cut)."""
    from src.common.metrics import read_jsonl
    recs = read_jsonl(RAW)
    summary = {}
    for alg in CONFIG["sizes"]:
        rs = sorted((r for r in recs if r["alg"] == alg and r["n"] in CONFIG["sizes"][alg]),
                    key=lambda r: r["n"])
        xs = np.array([r["n"] for r in rs], float)
        ys = np.array([np.log2(r["cost"]) for r in rs], float)
        h = len(xs) // 2
        lp, arg = exponent_lp(alg)
        summary[alg] = {"lp_exponent": str(lp), "lp_args": {k: str(v) for k, v in arg.items()},
                        "sizes": [int(x) for x in xs],
                        "fit_slope_all": float(np.polyfit(xs, ys, 1)[0]),
                        "fit_slope_upper_half": float(np.polyfit(xs[h:], ys[h:], 1)[0]),
                        "last_local_slope": float((ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))}
    with open(ANA, "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


def main():
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    with open(CFG, "w") as f:
        json.dump(CONFIG, f, indent=1)
    if os.path.exists(RAW):
        os.remove(RAW)
    fns = {"bsgs_sum": cost_bsgs_sum, "lehman": best_lehman,
           "hittmeir": best_hittmeir, "harvey": best_harvey, "hh": best_hh}
    summary = {}
    for alg, sizes in CONFIG["sizes"].items():
        xs, ys = [], []
        for n in sizes:
            t = time.time()
            rec = fns[alg](N_of(n))
            rec = {k: (int(v) if isinstance(v, (np.integer,)) else v) for k, v in rec.items()}
            rec.update({"n": n, "wall_s": round(time.time() - t, 3)})
            append_jsonl(RAW, rec)
            xs.append(n)
            ys.append(np.log2(rec["cost"]))
            print(alg, n, rec["cost"], rec["wall_s"], flush=True)
        xs, ys = np.array(xs, float), np.array(ys, float)
        h = len(xs) // 2
        slope_all = float(np.polyfit(xs, ys, 1)[0])
        slope_hi = float(np.polyfit(xs[h:], ys[h:], 1)[0])
        local = [float((ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])) for i in range(len(xs) - 1)]
        lp, arg = exponent_lp(alg)
        summary[alg] = {"lp_exponent": str(lp), "lp_args": {k: str(v) for k, v in arg.items()},
                        "fit_slope_all": slope_all, "fit_slope_upper_half": slope_hi,
                        "last_local_slope": local[-1]}
    with open(ANA, "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "summarize":
        summarize()
    else:
        main()
