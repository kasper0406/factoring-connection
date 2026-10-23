"""E1 (pre-registered, notes/navigate_prereg.md): toy navigation algorithms vs sampling. HARNESS.

The harness draws r, passes only a = r^2 mod N to the operational algorithms of core.py, verifies every output exactly,
and runs the size-matched coprime null for S0 and G1.

    .venv/bin/python -m src.navigate.e1_walks
"""
from __future__ import annotations

import json
import math
import os
import random
import time

from src.common.instances import balanced_semiprime
from src.navigate.core import (SmoothTester, primes_upto, reduced_key, run_g1, run_j2, run_s0, verify_relation)

BITS = (40, 48, 56)
SEEDS = range(6)
TARGETS = 100
BUDGET = 400
OUT = "experiments/raw/navigate_E1.jsonl"


class NullPairs:
    """Size-matched coprime null: for each tested split, an independent coprime pair of the same bit lengths."""

    def __init__(self, fb, seed):
        self.rng = random.Random(seed)
        self.tester = SmoothTester(fb)
        self.pairs = 0
        self.succ = 0

    def __call__(self, au: int, av: int):
        bu, bv = au.bit_length(), av.bit_length()
        while True:
            U = self.rng.randrange(1 << (bu - 1), 1 << bu) if bu > 1 else 1
            V = self.rng.randrange(1 << (bv - 1), 1 << bv) if bv > 1 else 1
            if math.gcd(U, V) == 1:
                break
        self.pairs += 1
        if self.tester.is_smooth(U) and self.tester.is_smooth(V):
            self.succ += 1


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        for bits in BITS:
            y = int(2 ** (bits / 4))
            fb = primes_upto(y)
            for seed in SEEDS:
                inst = balanced_semiprime(bits, seed)
                N = inst.N
                assert all(N % p for p in fb)
                hrng = random.Random(f"navigate-E1-targets|{bits}|{seed}")
                t0 = time.time()
                for ti in range(TARGETS):
                    while True:
                        r = hrng.randrange(1, N)
                        if math.gcd(r, N) == 1:
                            break
                    a = r * r % N
                    rec = {"bits": bits, "seed": seed, "N": N, "y": y, "target": ti}
                    for name, fn, use_null in (("S0", run_s0, True), ("G1", run_g1, True), ("J2", run_j2, False)):
                        nul = NullPairs(fb, f"navigate-E1-null|{name}|{bits}|{seed}|{ti}") if use_null else None
                        R = fn(N, a, fb, f"{seed}|{ti}", budget=BUDGET, null_cb=nul) if use_null else \
                            fn(N, a, fb, f"{seed}|{ti}", budget=BUDGET)
                        outs = R["outputs"]
                        ok = [verify_relation(u, v, a, N, fb) for u, v in outs]
                        assert all(ok), (name, bits, seed, ti)
                        dist = len({reduced_key(u, v) for u, v in outs})
                        pos = sum(1 for u, v in outs if u > 0 and v > 0)
                        d = {"outputs": len(outs), "distinct": dist, "positive": pos, "tests": R["tests"],
                             "splits": R["splits"]}
                        if use_null:
                            d["null_pairs"] = nul.pairs
                            d["null_succ"] = nul.succ
                        if "tested_logsizes" in R:
                            ls = R["tested_logsizes"]
                            d["mean_tested_logsize"] = sum(ls) / len(ls) if ls else None
                        rec[name] = d
                    fh.write(json.dumps(rec) + "\n")
                print(bits, seed, round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
