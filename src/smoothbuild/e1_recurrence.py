"""E1 (HARNESS): recurrence of large primes in pre-factored polynomial families vs random integers.

Pre-registered in notes/smoothbuild_prereg.md (E1). Run:
    .venv/bin/python -m src.smoothbuild.e1_recurrence
Writes experiments/raw/smoothbuild_E1.jsonl.
"""
from __future__ import annotations

import json
import os
import random
import time
from math import comb
from typing import Dict, List

import flint

from .core import FAMILIES, GROUPS, burnside_rank, inv_square_sum_fixed, primes_upto, root_list

A = 1 << 15
P1, P2 = 1 << 16, 1 << 22
STARTS = [(1 << 16) + j * (1 << 15) for j in range(5)]
SCALE = 160
OUT = "experiments/raw/smoothbuild_E1.jsonl"


def family_counts(coeffs: List[int], primes: List[int]) -> Dict[str, object]:
    """For each replicate start a0: S = sum_p C(c_p, 2) with c_p from root sieving, plus sum_p rho(rho-1)/p^2."""
    S = [0] * len(STARTS)
    arr = 0  # fixed-point sum of rho(rho-1)/p^2
    for p in primes:
        roots = root_list(coeffs, p)
        k = len(roots)
        if k >= 2:
            arr += k * (k - 1) * ((1 << SCALE) // (p * p))
        if k == 0:
            continue
        for j, a0 in enumerate(STARTS):
            c = sum(1 for r in roots if (r - a0) % p < A)
            if c >= 2:
                S[j] += comb(c, 2)
    return {"S": S, "arr_fixed": arr}


def random_counts(primes_set_lo: int, primes_set_hi: int, seed: int) -> int:
    rng = random.Random(f"smoothbuild-E1-R|{seed}")
    cnt: Dict[int, int] = {}
    for _ in range(A):
        z = rng.randrange(1 << 63, 1 << 64)
        for q, _e in flint.fmpz(z).factor():
            q = int(q)
            if primes_set_lo < q <= primes_set_hi:
                cnt[q] = cnt.get(q, 0) + 1
    return sum(comb(c, 2) for c in cnt.values() if c >= 2)


def main() -> None:
    t0 = time.time()
    primes = [p for p in primes_upto(P2) if p > P1]
    inv2 = inv_square_sum_fixed(primes, SCALE)
    one = 1 << SCALE
    # control expectation: C(A,2) * sum_p (m_p / 2^63)^2
    er_fixed = 0
    for p in primes:
        m = ((1 << 64) - 1) // p - ((1 << 63) - 1) // p
        er_fixed += (m * m << SCALE) // (1 << 126) // 1
    pairs = comb(A, 2)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        for name, (coeffs, gname) in FAMILIES.items():
            order, orbits, rank = burnside_rank(GROUPS[gname], len(GROUPS[gname][0]))
            res = family_counts(coeffs, primes)
            E = (rank - 1) * pairs * inv2 / one
            E_arr = (A * A / 2) * res["arr_fixed"] / one
            for j, s in enumerate(res["S"]):
                row = {"exp": "E1", "family": name, "group": gname, "order": order, "rank": rank,
                       "replicate": j, "a0": STARTS[j], "A": A, "P1": P1, "P2": P2, "S": s,
                       "E_chebotarev": E, "E_arrangement": E_arr,
                       "ratio_cheb": s / E, "ratio_arr": s / E_arr}
                fh.write(json.dumps(row) + "\n")
            print(name, res["S"], round(E, 1), round(E_arr, 1), round(time.time() - t0, 1), flush=True)
        ER = pairs * er_fixed / one
        for seed in range(5):
            s = random_counts(P1, P2, seed)
            row = {"exp": "E1", "family": "R:uniform[2^63,2^64)", "group": None, "rank": None,
                   "replicate": seed, "A": A, "P1": P1, "P2": P2, "S": s, "E_chebotarev": ER,
                   "ratio_cheb": s / ER}
            fh.write(json.dumps(row) + "\n")
            print("R", seed, s, round(ER, 1), round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
