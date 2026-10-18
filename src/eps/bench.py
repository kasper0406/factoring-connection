"""Per-instance cost benchmark on the 'bench' seed range (never reused).  Computes statistics only; it
does not compute or look at eps / QR_2 labels."""
from __future__ import annotations

import argparse
import json
import time
from multiprocessing import Pool

import numpy as np

from src.common.instances import balanced_semiprime
from src.eps.harness import SEED_BASE
from src.eps.operational import MENU_NAMES, compute_stats


def _one(args):
    n, seed = args
    inst = balanced_semiprime(n, seed)
    if inst.N % 4 != 1:
        return None
    t0 = time.perf_counter()
    st, aux = compute_stats(inst.N)
    return {"t": time.perf_counter() - t0, "tg": aux["t"], "x": [st[k] for k in MENU_NAMES],
            "first": {k: v for k, v in aux.items() if k.startswith("first_")}}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[32, 48, 64])
    ap.add_argument("--seeds", type=int, default=400)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="experiments/raw/eps_bench.jsonl")
    a = ap.parse_args(argv)
    res = {}
    for n in a.n:
        jobs = [(n, SEED_BASE["bench"] + j) for j in range(a.seeds)]
        with Pool(a.workers) as pool:
            rr = [r for r in pool.map(_one, jobs) if r is not None]
        t = np.array([r["t"] for r in rr])
        X = np.array([r["x"] for r in rr])
        groups = {g: float(np.mean([r["tg"][g] for r in rr])) for g in rr[0]["tg"]}
        res[n] = {"instances": len(rr), "mean_s": float(t.mean()), "p95_s": float(np.percentile(t, 95)),
                  "group_mean_s": groups,
                  "n_distinct": {k: int(len(np.unique(X[:, j]))) for j, k in enumerate(MENU_NAMES)},
                  "pc_rate": {k: float(np.mean([r["first"][k] > 0 for r in rr])) for k in rr[0]["first"]},
                  "pc_rate_T64": {k: float(np.mean([0 < r["first"][k] <= 64 for r in rr])) for k in rr[0]["first"]}}
        print(json.dumps({"n": n, **{k: v for k, v in res[n].items() if k != "n_distinct"}}), flush=True)
    with open(a.out, "w") as fh:
        for n, v in res.items():
            fh.write(json.dumps({"n": n, **v}) + "\n")


if __name__ == "__main__":
    main()
