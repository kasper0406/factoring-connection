"""HARNESS: empirical threshold (lattice dimension vs known bits) for formats (a)/(b).

For each (n, seed, fmt, d): scan k downward from a surely-sufficient value until two
consecutive failures; every attempt is one JSONL record.  fmt:
  low   : p mod 2^k                       (decode_low_bits)
  oddM  : p mod M, M a public random odd k-bit integer (decode_residue)
  high  : top k bits of p as an interval   (decode_interval)
Oracle (knows p) and decoder (sees N + hint) are separate calls.

Run:  .venv/bin/python -m src.trackE.exp_threshold experiments/configs/trackE_threshold.json
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from multiprocessing import Pool

from src.common.instances import balanced_semiprime
from src.common.metrics import Metrics, append_jsonl
from src.trackE import oracles
from src.trackE.coppersmith import (balanced_promise, best_m, decode_interval, decode_low_bits,
                                    decode_residue)


def heuristic_k(n: int, d: int, N: int, fmt: str, bits_p: int) -> int:
    L = balanced_promise(N)[0]
    m, ly = best_m(N.bit_length(), math.log2(L), d, hermite=1.02)
    if fmt == "high":   # W = 2^{bits_p - k} - 1, Y ~ W/2
        return max(1, math.ceil(bits_p - 1 - ly))
    # low/oddM: Y ~ sqrt(N)/(4 M)
    return max(1, math.ceil(N.bit_length() / 2 - 2 - ly))


def attempt(N, inst, fmt, k, d, m):
    mt = Metrics()
    if fmt == "low":
        h = oracles.oracle_low_bits(inst, k)                      # oracle stage
        t0 = time.perf_counter()
        g = decode_low_bits(N, h["k"], h["r"], d=d, m=m, metrics=mt)  # decoder stage
    elif fmt == "oddM":
        rng = random.Random(f"trackE-oddM|{N.bit_length()}|{k}")   # public, factor-independent
        M = rng.getrandbits(k) | (1 << (k - 1)) | 1
        h = oracles.oracle_residue(inst, M)
        t0 = time.perf_counter()
        g = decode_residue(N, h["M"], h["r"], d=d, m=m, metrics=mt)
    elif fmt == "high":
        h = oracles.oracle_high_bits(inst, k)
        t0 = time.perf_counter()
        g = decode_interval(N, h["A"], h["W"], d=d, m=m, metrics=mt)
    else:
        raise ValueError(fmt)
    dt = time.perf_counter() - t0
    ok = g is not None and N % g == 0 and 1 < g < N
    return mt.to_record(success=ok, decode_s=round(dt, 4))


def task(args):
    n, seed, fmt, d, out = args
    inst = balanced_semiprime(n, seed)
    N = inst.N
    L = balanced_promise(N)[0]
    m = best_m(N.bit_length(), math.log2(L), d, hermite=1.02)[0]
    bits_p = inst.labels["p"].bit_length()   # only used to size the scan start (harness)
    k = min(n // 2 - 1, heuristic_k(n, d, N, fmt, bits_p) + 4)
    recs, fails, best_ok = [], 0, None
    # go up until a success (bounded), then down until 2 consecutive failures
    while k < n // 2:
        r = attempt(N, inst, fmt, k, d, m)
        r.update(n=n, seed=seed, fmt=fmt, d=d, m_used=m, k=k)
        recs.append(r)
        if r["success"]:
            break
        k += 2
    if recs and recs[-1]["success"]:
        k -= 1
        while k >= 1 and fails < 2:
            r = attempt(N, inst, fmt, k, d, m)
            r.update(n=n, seed=seed, fmt=fmt, d=d, m_used=m, k=k)
            recs.append(r)
            fails = 0 if r["success"] else fails + 1
            k -= 1
    for r in recs:
        append_jsonl(out, r)
    return (n, seed, fmt, d, len(recs))


def main(cfg_path: str):
    cfg = json.load(open(cfg_path))
    out = cfg["out"]
    jobs = []
    for n_str, dims in cfg["dims"].items():
        n = int(n_str)
        for seed in range(cfg["seeds"]):
            for fmt in cfg["formats"]:
                for d in dims:
                    jobs.append((n, seed, fmt, d, out))
    jobs.sort(key=lambda j: -(j[0] * j[3] ** 2))   # biggest first
    with Pool(cfg.get("workers", 16)) as pool:
        for res in pool.imap_unordered(task, jobs):
            print(res, flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
