"""OPERATIONAL: TW* plus small-prime descent (count M of one fixed curve -> factor), for N in X_E.

Inputs: N, the short Weierstrass coefficients (a, b) of the fixed curve, M = #E(Z/N), a bound B, an rng.
Reuses ecount's Z/N group law (Lemma ZN: every non-unit met is a proper factor).
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Tuple

import gmpy2

from src.ecount.zn_curve import Split, ec_add, ec_mul, f_eval


def small_prime_parts(M: int, B: int) -> List[Tuple[int, int]]:
    """(l, e) with l <= B prime (trial division; l is tested prime by trial division of l itself) and l^e || M."""
    out = []
    m = M
    for l in range(2, B + 1):
        if any(l % r == 0 for r in range(2, int(l ** 0.5) + 1)):
            continue
        e = 0
        while m % l == 0:
            m //= l
            e += 1
        if e:
            out.append((l, e))
    return out


def descent_trial(N: int, a: int, b: int, M: int, parts: List[Tuple[int, int]],
                  rng: random.Random) -> Tuple[Optional[int], str]:
    x0 = rng.randrange(N)
    delta = f_eval(x0, a, b, N)
    g = gcd(delta, N)
    if 1 < g < N:
        return g, "delta_gcd"
    if g == N:
        return None, "delta_zero"
    ad = a * delta * delta % N
    P = (delta * x0 % N, delta * delta % N)
    try:
        if ec_mul(M, P, ad, N) is not None:
            return None, "affine"
        for (l, e) in parts:
            R = ec_mul(M // l ** e, P, ad, N)
            for _ in range(e):
                R = ec_mul(l, R, ad, N)
    except Split as s:
        return s.factor, "split"
    return None, "killed"


def xe_factor(N: int, a: int, b: int, M: int, rng: random.Random, B: int = 1000,
              max_trials: int = 20) -> Dict:
    parts = small_prime_parts(M, B)
    outcomes: Dict[str, int] = {}
    for t in range(1, max_trials + 1):
        fac, out = descent_trial(N, a, b, M, parts, rng)
        outcomes[out] = outcomes.get(out, 0) + 1
        if fac is not None:
            assert 1 < fac < N and N % fac == 0
            return {"factor": fac, "trials": t, "outcomes": outcomes}
    return {"factor": None, "trials": max_trials, "outcomes": outcomes}
