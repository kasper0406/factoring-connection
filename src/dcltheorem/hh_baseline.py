"""Baseline: Harvey--Hittmeir 2022 Alg. 4.3 (signsplit's operational V0, imported
unchanged), with its stage 3 (Harvey Alg. 4.1) replaced at call time by
stage3.harvey_find_collisions -- the same algorithm over fmpz_mod_poly (works
for N >= 2^64) with every scalar multiplication charged.  OPERATIONAL (sees N,
public M, m0, seed only).  Parameter choice: min over M in the primorials of
signsplit.accounting.best_shared (list count), as in signsplit.

Counter mapping to metric E: signsplit's Ring charges rn_mul (square-and-multiply
exact) and gcd/inv; stage 3 charges mulN/gcdN/invN/poly_cost.
"""
from __future__ import annotations

from typing import Dict

import src.signsplit.hh as sshh
from src.signsplit.accounting import best_shared

from .stage3 import harvey_find_collisions

PRIMORIALS = (1, 2, 6, 30, 210, 2310, 30030)


def hh_params(N: int) -> Dict:
    best = None
    for M in PRIMORIALS:
        b = best_shared(N, M, "V0")
        if best is None or b["work"] < best["work"]:
            best = b
    return best


def hh_run(N: int, params: Dict, seed: int = 0) -> Dict:
    orig = sshh.find_collisions
    sshh.find_collisions = harvey_find_collisions
    try:
        r = sshh.hh_search(N, params["M"], params["m0"], "V0", seed=seed)
    finally:
        sshh.find_collisions = orig
    return r


def hh_E(counts: Dict[str, int], n: int) -> int:
    lg = max(1, (n - 1).bit_length())
    return (counts.get("rn_mul", 0) + counts.get("mulN", 0) + counts.get("poly_cost", 0)
            + lg * (counts.get("gcd", 0) + counts.get("gcdN", 0) + counts.get("inv", 0)
                    + counts.get("invN", 0)))
