"""E2 (HARNESS): residue-class construction vs sampling at matched size.

Pre-registered in notes/smoothbuild_prereg.md (E2). Run:
    .venv/bin/python -m src.smoothbuild.e2_class
Writes experiments/raw/smoothbuild_E2.jsonl (one row per (bits, N-seed, method)).

Constructors (M0, M1, M2, M4) see only (N, t). The smoothness test of the small pieces uses flint factorisation of
the PIECE (never of N); this file is harness code and is not imported by operational code.
"""
from __future__ import annotations

import json
import os
import random
import time
from typing import List

import flint

from src.common.instances import balanced_semiprime

from .core import gauss_split, primes_upto

OUT = "experiments/raw/smoothbuild_E2.jsonl"
SIZES = [40, 48, 56, 64]
SEEDS = list(range(6))
N_TARGETS = 20_000
N_TARGETS_M2 = 2_000
K_MULT = 16
LIST_SIZE = 1024
M4_TARGETS = 3


def is_smooth(z: int, y: int) -> bool:
    z = abs(z)
    if z == 0:
        return False
    if z == 1:
        return True
    return max(int(p) for p, _e in flint.fmpz(z).factor()) <= y


def same_bitlength(rng: random.Random, z: int) -> int:
    b = abs(z).bit_length()
    if b <= 1:
        return 1
    return rng.randrange(1 << (b - 1), 1 << b)


def smooth_list(lo: int, y: int, count: int) -> List[int]:
    out = []
    z = lo
    while len(out) < count:
        if is_smooth(z, y):
            out.append(z)
        z += 1
    return out


def run_one(bits: int, seed: int, fh) -> None:
    inst = balanced_semiprime(bits, seed)
    N = inst.N
    y = 1 << (bits // 4)
    rng = random.Random(f"smoothbuild-E2|{bits}|{seed}")
    nrng = random.Random(f"smoothbuild-E2-null|{bits}|{seed}")
    mults = primes_upto(60)[:K_MULT]

    def target() -> int:
        while True:
            r = rng.randrange(2, N - 1)
            t = r * r % N
            if t:
                return t

    # M0 and M1 on the same targets, with nulls
    m0 = m0_null = m1 = m1_null = 0
    for _ in range(N_TARGETS):
        t = target()
        m0 += is_smooth(t, y)
        m0_null += is_smooth(nrng.randrange(1, N), y)
        u, v = gauss_split(t, N)
        m1 += is_smooth(u, y) and is_smooth(v, y)
        m1_null += is_smooth(same_bitlength(nrng, u), y) and is_smooth(same_bitlength(nrng, v), y)
    rows = [
        {"method": "M0", "trials": N_TARGETS, "hits": m0, "null_hits": m0_null, "tests_per_trial": 1},
        {"method": "M1", "trials": N_TARGETS, "hits": m1, "null_hits": m1_null, "tests_per_trial": 2},
    ]
    # M2: K multipliers per target
    split_hits = target_hits = 0
    for _ in range(N_TARGETS_M2):
        t = target()
        any_hit = False
        for s in mults:
            u, v = gauss_split(t * s % N, N)
            h = is_smooth(u, y) and is_smooth(v, y)
            split_hits += h
            any_hit = any_hit or h
        target_hits += any_hit
    rows.append({"method": "M2", "trials": N_TARGETS_M2 * K_MULT, "hits": split_hits,
                 "targets": N_TARGETS_M2, "target_hits": target_hits, "K": K_MULT, "tests_per_trial": 2})
    # M4: flatness of products of two smooth lists in a window mod N
    lo = 1 << (bits // 2 - 2)
    A = smooth_list(lo, y, LIST_SIZE)
    W = -(-400 * N // (LIST_SIZE * LIST_SIZE))
    hits = 0
    expected_num = 0
    for _ in range(M4_TARGETS):
        t = target()
        ti = pow(t, -1, N)
        for a in A:
            c = a * ti % N
            for b in A:
                w = c * b % N
                if 1 <= w <= W:
                    hits += 1
        expected_num += LIST_SIZE * LIST_SIZE * W
    rows.append({"method": "M4", "trials": M4_TARGETS * LIST_SIZE * LIST_SIZE, "hits": hits,
                 "expected": expected_num / N, "W": W, "list_lo": lo})
    for r in rows:
        r.update({"exp": "E2", "bits": bits, "seed": seed, "N": N, "y": y})
        fh.write(json.dumps(r) + "\n")
    fh.flush()


def main() -> None:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    t0 = time.time()
    with open(OUT, "w") as fh:
        for bits in SIZES:
            for seed in SEEDS:
                run_one(bits, seed, fh)
                print(bits, seed, round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
