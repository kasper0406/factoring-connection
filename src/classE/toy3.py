"""HARNESS (knows the factors): exhaustive-population toy test of the SHIFTED d = 3 Fourier-Dedekind sums.

Population P_n: all N = pq, p < q primes of n/2 bits with the top two bits set, N = 1 (mod 4) (n in {20, 22}).
Statistics (OPERATIONAL part `toy3_stats(N)`, N only): for each registered pair (a, b) of structured arguments,
S3 = N^3 s_1(a, b, 1; N) (an integer) and Z3 = floor(N s_1(a, b, 1; N)) (the non-aggregate integer content):
  S3 mod 7, bitlen|S3|, Z3 mod 7, sign(Z3).
Targets eps and tb exactly as src/classE/harness.py.  Floors: eps -> N mod 2^8, tb -> top 6 bits of N.
Registered in notes/classE_prereg.md section 6.  The O(N) evaluator restricts this to toy sizes.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from math import gcd
from multiprocessing import Pool
from typing import Dict, List

import numpy as np

from src.common.instances import is_prime
from src.eps.analysis import strat_z
from src.classE.fd import fd_scaled
from src.classE.structure import arguments

PAIRS = (("isqrt", "cbrt2"), ("fermres", "euler3"), ("isqrt", "invsqrt"), ("hash2", "isqrt"))
STATS = ("S3mod7", "S3bits", "Z3mod7", "Z3sign")
NAMES = [f"{a}_{b}_{s}" for a, b in PAIRS for s in STATS]
Z_TOY = 4.5


def toy3_stats(N: int) -> Dict[str, int]:
    args = arguments(N)
    out = {}
    for a_name, b_name in PAIRS:
        a, b = args[a_name], args[b_name]
        if a <= 0 or b <= 0 or gcd(a * b, N) != 1:
            for s in STATS:
                out[f"{a_name}_{b_name}_{s}"] = 0
            continue
        S3 = fd_scaled(1, [a, b, 1], N)            # N^3 s
        Z3 = S3 // (N * N)                         # floor(N s)
        out[f"{a_name}_{b_name}_S3mod7"] = S3 % 7
        out[f"{a_name}_{b_name}_S3bits"] = abs(S3).bit_length()
        out[f"{a_name}_{b_name}_Z3mod7"] = Z3 % 7
        out[f"{a_name}_{b_name}_Z3sign"] = (Z3 > 0) - (Z3 < 0)
    return out


def population(n: int) -> List[tuple]:
    h = n // 2
    lo, hi = 3 << (h - 2), 1 << h
    ps = [x for x in range(lo | 1, hi, 2) if is_prime(x)]
    return [(p, q) for i, p in enumerate(ps) for q in ps[i + 1:] if (p * q) % 4 == 1]


def _work(pq):
    p, q = pq
    N = p * q
    t0 = time.perf_counter()
    st = toy3_stats(N)
    return {"N": N, "p": p, "q": q, "eps": 1 if p % 4 == 1 else -1, "tb": int(8 * p >= 7 * q),
            "x": [st[k] for k in NAMES], "t": time.perf_counter() - t0}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[20, 22])
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args(argv)
    for n in a.n:
        pop = population(n)
        with Pool(a.workers) as pool:
            recs = pool.map(_work, pop, chunksize=8)
        os.makedirs("experiments/raw", exist_ok=True)
        with open(f"experiments/raw/classE_toy3_n{n}.jsonl", "w") as fh:
            for r in recs:
                fh.write(json.dumps({k: (v if k != "t" else round(v, 4)) for k, v in r.items()}) + "\n")
        X = np.array([r["x"] for r in recs], dtype=np.float64)
        res = {"n": n, "instances": len(recs), "mean_s": float(np.mean([r["t"] for r in recs])), "stats": {}}
        for tgt, key in (("eps", lambda r: r["N"] % 256), ("tb", lambda r: r["N"] >> (n - 6))):
            y = np.array([(r["eps"] == 1) if tgt == "eps" else r["tb"] for r in recs], dtype=np.float64)
            h = np.array([key(r) for r in recs])
            zs = {}
            for j, nm in enumerate(NAMES):
                z = strat_z(X[:, j], y, h)
                zs[nm] = {"z": z["z"], "degenerate": z["degenerate"]}
            nondeg = [v["z"] for v in zs.values() if not v["degenerate"]]
            res["stats"][tgt] = {"z": zs, "max_abs_z": float(np.max(np.abs(nondeg))) if nondeg else None,
                                 "mean_z2": float(np.mean(np.square(nondeg))) if nondeg else None,
                                 "hits": [k for k, v in zs.items() if not v["degenerate"] and abs(v["z"]) >= Z_TOY]}
            print(n, tgt, "instances", len(recs), "max|z|", res["stats"][tgt]["max_abs_z"], "mean z2",
                  res["stats"][tgt]["mean_z2"], "hits", res["stats"][tgt]["hits"], flush=True)
        os.makedirs("experiments/analysis", exist_ok=True)
        with open(f"experiments/analysis/classE_toy3_n{n}.json", "w") as fh:
            json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
