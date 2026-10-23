"""POST HOC (not pre-registered): E2's counts in the sparse regime (K = pi(y) >> number of letters).

Written after E2's P2a-P2c failed, to test the explanation that E2's toy cells were in a *dense* regime (each of
K = 6-8 primes used with many exponents), where C_beta is a low-frequency function of the coset of the relation
lattice, while NAV at L[1/3] is sparse (K >> letters).  Same code as e2_counts.run_instance; only y and beta differ.

    .venv/bin/python -m src.navigate.posthoc_e2
"""
from __future__ import annotations

import json
import time

from src.navigate import e2_counts as E

CONFIGS = {512: [0.8, 0.9, 1.0, 1.1], 2048: [0.8, 0.9, 1.0]}
OUT = "experiments/raw/navigate_E2_posthoc.jsonl"


def main():
    E.LO = 0.5
    with open(OUT, "w") as fh:
        for y, betas in CONFIGS.items():
            for seed in E.SEEDS:
                t0 = time.time()
                r = E.run_instance(seed, y, betas)
                r["posthoc"] = True
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                print(y, seed, r["vectors"], r["H"], r["phi"], round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
