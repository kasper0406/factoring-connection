"""HARNESS: Track A (a) circuit-section growth experiment.

Sections the hash-consed circuit of F_N along several paths and records unique
reachable nodes per level.  Paths:
  diag_q      bits of the support index q        (diagnostic: uses factors)
  diag_pp     bits of the support index p*p      (diagnostic: uses factors)
  op_zero     all-zero bits                      (operational)
  op_random   seeded random bits                 (operational)
  op_exact    bits chosen by ExactScalarSelector (operational)

Usage:
  .venv/bin/python -m src.trackA.run_circuits --bits 16 24 32 --seeds 2 --levels 9 \
      --out experiments/raw/trackA_circuits.jsonl
"""
from __future__ import annotations

import argparse
import random
import time

from ..circuit_sections import run_path
from ..common import instances as I
from ..common.metrics import append_jsonl
from .selectors import ExactScalarSelector, descend


def bits_of(i, k):
    return [(i >> t) & 1 for t in range(k)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 24, 32])
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--levels", type=int, default=9)
    ap.add_argument("--node-cap", type=int, default=2_500_000)
    ap.add_argument("--out", default="experiments/raw/trackA_circuits.jsonl")
    a = ap.parse_args()
    for bits in a.bits:
        for seed in range(a.seeds):
            inst = I.balanced_semiprime(bits, seed)
            N = inst.N
            p, q = inst.factors
            k = a.levels
            rng = random.Random(f"circ|{N}|{seed}")
            r_exact, _, _, _ = descend(N, ExactScalarSelector(), seed=seed, max_level=k)
            paths = {"diag_q": bits_of(q, k), "diag_pp": bits_of(p * p, k),
                     "op_zero": [0] * k, "op_random": [rng.getrandbits(1) for _ in range(k)],
                     "op_exact": bits_of(r_exact, k)}
            for kar in (False, True):
                for name, path in paths.items():
                    t0 = time.perf_counter()
                    rows = run_path(N, path, node_cap=a.node_cap, karatsuba=kar)
                    rec = {"experiment": "trackA_circuits", "N": N, "bits": bits, "seed": seed,
                           "path": name, "karatsuba": kar, "path_bits": path, "rows": rows,
                           "wall_s": time.perf_counter() - t0}
                    append_jsonl(a.out, rec)
                    print(bits, seed, name, "kar" if kar else "naive",
                          [r_["nodes"] for r_ in rows], flush=True)


if __name__ == "__main__":
    main()
