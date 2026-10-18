"""E1 driver: exact minimal-height rank-r GAP covers (see notes/divcover_prereg.md §2).

Run from the project root:
    .venv/bin/python -m src.divcover.run_exact [--jobs 30]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

from src.divcover.core import c_min_height, random_model_log_height, shapes_of_rank, _binary

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CFG = os.path.join(ROOT, "experiments/configs/divcover_exact.json")


def run_cell(n, L, r, mode, Hmax, timeout):
    shapes = shapes_of_rank(r, L)
    best, arg, censored, nparams = None, None, [], 0
    t0 = time.time()
    cur_max = Hmax
    for sh in shapes:
        try:
            res = c_min_height(n, sh, cur_max, mode, timeout=timeout)
        except subprocess.TimeoutExpired:
            censored.append(list(sh))
            continue
        nparams += res["nparams"]
        if res["best"] is not None and (best is None or res["best"] < best):
            best = res["best"]
            arg = {"shape": list(sh), "b": res["b"], "c": list(res["c"])}
            cur_max = best - 1  # still exact: later shapes only need to beat best
    return {
        "n": n, "L": L, "rank": r, "mode": mode, "Hmax": Hmax, "H_star": best,
        "argmin": arg, "timed_out_shapes": censored, "nparams": nparams,
        "wall_s": round(time.time() - t0, 3),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=30)
    args = ap.parse_args()
    cfg = json.load(open(CFG))
    _binary()  # compile once
    cache = {}
    cells = []
    for n in cfg["n_grid"]:
        for L in cfg["L_grid"]:
            if L > cfg["L_max_frac_of_n"] * n:
                continue
            for r in cfg["ranks"]:
                shapes = shapes_of_rank(r, L)
                if not shapes:
                    continue
                lh = random_model_log_height(n, shapes, cache)
                if math.exp(lh) > cfg["cap_pred_height"][str(r)]:
                    continue
                for mode in cfg["modes"]:
                    Hmax = math.ceil(cfg["Hmax_factor"][mode] * math.exp(lh))
                    cells.append((n, L, r, mode, Hmax, lh))
    print(f"{len(cells)} cells", flush=True)
    out = os.path.join(ROOT, cfg["raw"])
    with ProcessPoolExecutor(max_workers=args.jobs) as ex, open(out, "w") as fh:
        futs = {ex.submit(run_cell, n, L, r, mode, Hmax, cfg["timeout_s"]): (lh,)
                for (n, L, r, mode, Hmax, lh) in cells}
        for f in as_completed(futs):
            rec = f.result()
            rec["logH_pred_prime"] = futs[f][0]
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(rec["n"], rec["L"], rec["rank"], rec["mode"], rec["H_star"],
                  round(math.exp(rec["logH_pred_prime"]), 1), rec["wall_s"], flush=True)


if __name__ == "__main__":
    main()
