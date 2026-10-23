"""POST HOC (not pre-registered): fixed-window flatness of smooth x smooth products mod N.

Replaces E2-M4, which the codex adversarial round showed to be uninformative (its lists were consecutive smooth
numbers whose products never wrap mod N, so only the random target dilation was measured). Here A is a set of
random y-smooth integers in [N^0.55, 2 N^0.55] (products wrap ~N^0.1 times), and we count pairs (a, b), a < b, with
a*b mod N in the FIXED window [0, W] (the Prop SI window near kN). Expectation under flatness: C(|A|,2)*(W+1)/N.
Run: .venv/bin/python -m src.smoothbuild.posthoc_m4 -> experiments/raw/smoothbuild_E2_posthoc_m4.jsonl
"""
from __future__ import annotations

import json
import random

from src.common.instances import balanced_semiprime

from .e2_class import is_smooth

OUT = "experiments/raw/smoothbuild_E2_posthoc_m4.jsonl"
LIST = 1500


def main() -> None:
    with open(OUT, "w") as fh:
        for bits in (40, 48, 56, 64):
            y = 1 << (bits // 4)
            hits = 0
            expct = 0.0
            for seed in range(6):
                N = balanced_semiprime(bits, seed).N
                rng = random.Random(f"smoothbuild-M4post|{bits}|{seed}")
                lo = int(N ** 0.55)
                A = set()
                while len(A) < LIST:
                    z = rng.randrange(lo, 2 * lo)
                    if is_smooth(z, y):
                        A.add(z)
                A = sorted(A)
                pairs = LIST * (LIST - 1) // 2
                W = 400 * N // pairs
                h = 0
                for i, a in enumerate(A):
                    for b in A[i + 1:]:
                        if a * b % N <= W:
                            h += 1
                hits += h
                expct += pairs * (W + 1) / N
            row = {"bits": bits, "y": y, "hits": hits, "expected": expct, "O/E": hits / expct}
            fh.write(json.dumps(row) + "\n")
            print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
