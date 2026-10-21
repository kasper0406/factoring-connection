"""E1 (pre-registered, notes/schnorrcvp_prereg.md): exact fac-relation counts in boxes vs heuristic H_S / H_S'.

Usage: .venv/bin/python -m src.schnorrcvp.run_e1 [--quick]
Writes experiments/raw/schnorrcvp_E1.jsonl.
"""
from __future__ import annotations

import json
import sys
import time
from math import gcd

from src.common.instances import balanced_semiprime

from .core import first_primes
from .enum_rel import relations_in_box
from .psi import expected_in_box, expected_primitive_in_box

GRID = [  # (bits, n, V, W)
    (24, 10, 10**3, 10**5), (28, 10, 10**3, 10**5), (32, 10, 10**3, 10**5),
    (24, 15, 10**3, 10**5), (32, 15, 10**3, 10**5), (40, 15, 10**3, 10**5),
    (24, 25, 10**3, 2 * 10**4), (32, 25, 10**3, 2 * 10**4), (40, 25, 10**3, 2 * 10**4),
]
SEEDS = range(10)


def main(quick: bool = False):
    out = open("experiments/raw/schnorrcvp_E1.jsonl", "w")
    grid = GRID[:2] if quick else GRID
    seeds = range(2) if quick else SEEDS
    for bits, n, V, W in grid:
        P = first_primes(n)
        for seed in seeds:
            N = balanced_semiprime(bits, seed).N
            t0 = time.time()
            r = relations_in_box(N, P, V, W)
            prim = [(u, v, w) for (u, v, w) in r["relations"] if gcd(v, w) == 1]
            t1 = time.time()
            naive = expected_in_box(N, P, V, W, K=2000, h=0.04)
            refined = expected_primitive_in_box(N, P, V, W, K=2000, h=0.04)
            rec = {"bits": bits, "n": n, "V": V, "W": W, "seed": seed, "N": str(N),
                   "obs_all": r["count"], "obs_primitive": len(prim),
                   "pred_naive": float(naive["E"]), "pred_naive_lo": float(naive["E_lo"]),
                   "pred_naive_hi": float(naive["E_hi"]), "pred_primitive": float(refined["E"]),
                   "n_v": r["n_v"], "n_w": r["n_w"], "primitive_list": [[str(u), str(v), str(w)] for u, v, w in prim],
                   "t_enum": round(t1 - t0, 2), "t_model": round(time.time() - t1, 2)}
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(bits, n, seed, rec["obs_all"], rec["obs_primitive"], round(rec["pred_naive"], 2),
                  round(rec["pred_primitive"], 2), rec["t_enum"], rec["t_model"], flush=True)


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
