"""Registered experiments E1–E3 of notes/siconstruct_prereg.md.  Writes experiments/raw/siconstruct_E{1,2,3}.jsonl.

Usage: .venv/bin/python -m src.siconstruct.run_exp E1|E2|E3
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time

from src.common.instances import balanced_semiprime

from .constructors import aim2d_construct, brute_interval, mitm_construct
from .core import ht_local_density, lpath_params, primes_upto, random_smooth_list
from .theory import mitm_expected_hits, mitm_list_size

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f"siconstruct_{name}.jsonl"), "w")


def e1():
    bits = 48
    f = _out("E1")
    for seed in (0, 1, 2):
        N = balanced_semiprime(bits, seed).N
        pp = lpath_params(bits)
        y, h = pp["y"], pp["h"]
        P = primes_upto(y)
        for K in (1, 4, 16, 64, 256):
            rng = random.Random(f"E1|{seed}|{K}")
            A0 = 1 << 20
            B0 = -(-K * N // A0)
            A = random_smooth_list(A0, 2 * A0, P, 4096, rng)
            B = random_smooth_list(B0, 2 * B0, P, 4096, rng)
            t = time.time()
            r = mitm_construct(N, h, A, B)
            pred = float(mitm_expected_hits(len(A), len(B), h, N))
            rec = {"exp": "E1", "bits": bits, "seed": seed, "N": N, "y": y, "h": h, "K": K,
                   "nA": len(A), "nB": len(B), "hits": r["hits"], "outputs": len(r["outputs"]),
                   "pred_hits": pred, "ratio": r["hits"] / pred, "queries": r["queries"], "points": r["points"],
                   "gen": r["gen"], "secs": time.time() - t,
                   "k_values": sorted({k for _, (k, _) in r["outputs"]})[:50]}
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print(bits, seed, K, r["hits"], round(pred, 1), round(r["hits"] / pred, 3), flush=True)
    f.close()


def e2():
    f = _out("E2")
    R0 = 64
    for bits in (32, 40, 48, 56, 64, 72, 80):
        pp = lpath_params(bits)
        y, h = pp["y"], pp["h"]
        P = primes_upto(y)
        for seed in (0, 1):
            N = balanced_semiprime(bits, seed).N
            rng = random.Random(f"E2|{bits}|{seed}")
            M = mitm_list_size(N, h, R0)
            A0 = 1 << ((bits + 1) // 2)
            B0 = -(-N // A0)
            t = time.time()
            A = random_smooth_list(A0, 2 * A0, P, M, rng)
            B = random_smooth_list(B0, 2 * B0, P, M, rng)
            tg = time.time() - t
            r = mitm_construct(N, h, A, B)
            pred = float(mitm_expected_hits(M, M, h, N))
            cost = r["gen"] + r["queries"] + r["hits"]
            rec = {"exp": "E2", "method": "mitm", "bits": bits, "seed": seed, "N": N, "y": y, "h": h,
                   "M": M, "hits": r["hits"], "outputs": len(r["outputs"]), "pred_hits": pred,
                   "queries": r["queries"], "points": r["points"], "gen": r["gen"], "cost": cost,
                   "cost_per_output": cost / max(1, len(r["outputs"])),
                   "model_per_output": 2 * math.sqrt(N / (R0 * (h + 1))),
                   "secs_gen": tg, "secs": time.time() - t}
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print("mitm", bits, seed, M, r["hits"], round(pred, 1), round(cost / max(1, len(r["outputs"])), 1),
                  flush=True)
            # brute-force baseline (measured for bits <= 56; H-T model for all)
            dens_model = ht_local_density(float(N), y, P)
            rec2 = {"exp": "E2", "method": "brute", "bits": bits, "seed": seed, "N": N, "y": y, "h": h,
                    "ht_density": dens_model, "ht_per_output": 1 / dens_model}
            if bits <= 56:
                W = min(h + 1, 1 << 26)
                tot = 0
                found = 0
                t = time.time()
                for k in range(1, 9):
                    rb = brute_interval(N, k, h, P, W)
                    tot += rb["cost"]
                    found += len(rb["outputs"])
                rec2.update({"sieved": tot, "found": found, "measured_density": found / tot,
                             "measured_per_output": tot / max(1, found), "secs": time.time() - t})
            f.write(json.dumps(rec2) + "\n")
            f.flush()
            print("brute", bits, seed, rec2.get("found"), rec2.get("sieved"), round(1 / dens_model), flush=True)
    f.close()


def e3():
    f = _out("E3")
    K = 16
    for bits in (40, 48, 56):
        pp = lpath_params(bits)
        y, h = pp["y"], pp["h"]
        P = primes_upto(y)
        lh = h.bit_length() - 1
        for seed in (0, 1):
            N = balanced_semiprime(bits, seed).N
            for S0exp in (lh // 2, lh - 4):
                rng = random.Random(f"E3|{bits}|{seed}|{S0exp}")
                S0 = 1 << S0exp
                S = random_smooth_list(S0, 2 * S0, P, 64, rng)
                t = time.time()
                r = aim2d_construct(N, h, K, S, P, max_per_sk=2048)
                # prediction: sum of H-T local densities at the tested cofactor sizes (bucketed by bit length)
                cache = {}
                pred = 0.0
                for bl, cnt in r["cofactor_bits"].items():
                    if bl not in cache:
                        cache[bl] = ht_local_density(float(2 ** bl) * 0.75, y, P)
                    pred += cnt * cache[bl]
                dN = ht_local_density(float(N) * K / 2, y, P)
                nout = len(r["outputs"])
                rec = {"exp": "E3", "bits": bits, "seed": seed, "N": N, "y": y, "h": h, "K": K, "S0exp": S0exp,
                       "nS": len(S), "tested": r["tested"], "outputs": nout, "pred_outputs": pred,
                       "ratio": nout / pred if pred else None, "cof_bits_min": r["cofactor_bits_min"],
                       "per_output_tests": r["tested"] / max(1, nout), "brute_per_output_model": 1 / dN,
                       "secs": time.time() - t}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(bits, seed, S0exp, r["tested"], nout, round(pred, 1), flush=True)
    f.close()


if __name__ == "__main__":
    {"E1": e1, "E2": e2, "E3": e3}[sys.argv[1]]()
