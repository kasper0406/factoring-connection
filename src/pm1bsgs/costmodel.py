"""OPERATIONAL-SAFE exact cost formulas (eq_ops) mirroring operational.py.

Depends only on public parameters (N, B1, B2, m, J, K, seed).  Tested against
the instrumented counts in tests/test_pm1bsgs.py.  Backtracking work (after a
gcd == N) is not included; it is reported separately by the runs.

Prime counts and maximal prime gaps up to SIEVE_LIMIT come from an exact numpy
sieve (cached on disk under experiments/analysis/pm1bsgs_sieve_cache.json).
"""
from __future__ import annotations

import json
import os
from functools import lru_cache
from typing import Dict, Optional, Tuple

import numpy as np

from .operational import (PRIME2_CHUNK, STAGE1_CHUNK, birthday_schedule,
                          mcost, pow_cost, rem_cost)

SIEVE_LIMIT = 1 << 30
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_CACHE = os.path.join(_ROOT, "experiments", "analysis", "pm1bsgs_sieve_cache.json")
_PRIMES: Optional[np.ndarray] = None
_PRIMES_LIMIT = 0


def _primes_np(limit: int) -> np.ndarray:
    """All primes <= limit (possibly more: callers filter)."""
    global _PRIMES, _PRIMES_LIMIT
    if _PRIMES is not None and _PRIMES_LIMIT >= limit:
        return _PRIMES
    limit = max(limit, 1 << 16)
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    s[4::2] = False
    for i in range(3, int(limit ** 0.5) + 1, 2):
        if s[i]:
            s[i * i::2 * i] = False
    _PRIMES = np.nonzero(s)[0].astype(np.int64)
    _PRIMES_LIMIT = limit
    return _PRIMES


def _popcount(a: np.ndarray) -> np.ndarray:
    a = a.astype(np.uint64)
    c = np.zeros(a.shape, dtype=np.int64)
    while np.any(a):
        c += (a & np.uint64(1)).astype(np.int64)
        a >>= np.uint64(1)
    return c


def _load_cache() -> Dict:
    try:
        with open(_CACHE) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _save_cache(d: Dict) -> None:
    os.makedirs(os.path.dirname(_CACHE), exist_ok=True)
    with open(_CACHE, "w") as f:
        json.dump(d, f, indent=0, sort_keys=True)


_MEM: Dict = {}


def _cached(key: str, fn):
    if key in _MEM:
        return _MEM[key]
    if not _MEM:
        _MEM.update(_load_cache())
        if key in _MEM:
            return _MEM[key]
    v = fn()
    _MEM[key] = v
    _save_cache(_MEM)
    return v


# ---------------------------------------------------------------- stage 1
def stage1_cost(B1: int) -> int:
    """rn_mul + gcd of operational.stage1 (no backtracking)."""
    if B1 <= 1:
        return 0
    def f():
        P = _primes_np(B1)
        P = P[P <= B1]
        mult = np.ones(len(P), dtype=np.int64)
        small = P[P * P <= B1]
        for idx, l in enumerate(small.tolist()):
            e, pk = 0, l
            while pk <= B1:
                e += 1
                pk *= l
            mult[idx] = e
        pc = _popcount(P)
        bits = _bitlen_minus1(P)
        cost_each = bits + pc - 1
        cost_each[P <= 1] = 0
        steps = int(mult.sum())
        return int((cost_each * mult).sum()) + (steps + STAGE1_CHUNK - 1) // STAGE1_CHUNK
    return _cached(f"s1|{B1}", f)


def _bitlen_minus1(a: np.ndarray) -> np.ndarray:
    a = a.astype(np.uint64)
    r = np.zeros(a.shape, dtype=np.int64)
    for sh in (32, 16, 8, 4, 2, 1):
        mask = a >= (np.uint64(1) << np.uint64(sh))
        r[mask] += sh
        a = np.where(mask, a >> np.uint64(sh), a)
    return r


# ---------------------------------------------------------------- trees
@lru_cache(maxsize=None)
def _tree_levels(n: int) -> Tuple[Tuple[int, ...], ...]:
    level = [1] * n
    levels = [tuple(level)]
    while len(level) > 1:
        nxt = [level[i] + level[i + 1] for i in range(0, len(level) - 1, 2)]
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
        levels.append(tuple(level))
    return tuple(levels)


@lru_cache(maxsize=None)
def tree_cost(n: int) -> int:
    levels = _tree_levels(n)
    c = 0
    for lv in levels[:-1]:
        for i in range(0, len(lv) - 1, 2):
            c += mcost(max(lv[i], lv[i + 1]))
    return c


@lru_cache(maxsize=None)
def remtree_cost(degF: int, n: int) -> int:
    levels = _tree_levels(n)
    c = 0
    top = levels[-1][0]
    dR = degF
    if dR >= top:
        c += rem_cost(dR, top)
        dR = top - 1
    cur = [dR]
    for k in range(len(levels) - 2, -1, -1):
        lv = levels[k]
        nxt = []
        for idx, dg in enumerate(lv):
            dPR = cur[idx // 2]
            single = (idx == len(lv) - 1) and (len(lv) % 2 == 1)
            if single or dPR < dg:
                nxt.append(dPR)
            else:
                c += rem_cost(dPR, dg)
                nxt.append(dg - 1)
        cur = nxt
    return c


# ---------------------------------------------------------------- stage 2
def bsgs_stage2_cost(mbaby: int, J: int, mode: str, N: int, block: Optional[int] = None) -> int:
    block = block or mbaby
    c = pow_cost(N) if mode == "N" else 0
    c += mbaby - 1
    c += tree_cost(mbaby)
    c += pow_cost(mbaby)
    c += J - 1
    for c0 in range(0, J, block):
        L = min(block, J - c0)
        c += tree_cost(L) + remtree_cost(mbaby, L) + L + 1
    return c


def prime_count_and_gap(B1: int, B2: int) -> Tuple[int, int, int]:
    """(#primes in (B1,B2], max gap between consecutive such primes, first prime)."""
    assert B2 <= SIEVE_LIMIT
    def f():
        P = _primes_np(B2)
        lo = int(np.searchsorted(P, B1, side="right"))
        hi = int(np.searchsorted(P, B2, side="right"))
        if hi <= lo:
            return [0, 0, 0]
        seg = P[lo:hi]
        g = int(np.diff(seg).max()) if hi - lo > 1 else 0
        return [hi - lo, g, int(seg[0])]
    return tuple(_cached(f"pg|{B1}|{B2}", f))


def prime_stage2_cost(B1: int, B2: int) -> int:
    npr, maxgap, r0 = prime_count_and_gap(B1, B2)
    if npr == 0:
        return 0
    c = pow_cost(r0) + 1                          # x = b^r0, b2 = b^2
    if npr > 1:
        c += max(0, (maxgap - 2) // 2)            # gap table b^4 .. b^maxgap
    c += npr + (npr - 1)                          # accumulate + step
    c += (npr + PRIME2_CHUNK - 1) // PRIME2_CHUNK
    return c


def birthday_stage2_cost(K: int, N: Optional[int] = None, seed: Optional[int] = None) -> int:
    if N is not None and seed is not None:
        c = (K - 1) + sum(birthday_schedule(N, seed, K))
    else:
        c = (K - 1) + (K - 1) // 2                # nominal: half the bits are 1
    if K < 2:
        return c
    c += tree_cost(K) + K + remtree_cost(K - 1, K) + K + 1
    return c


def base_cost() -> int:
    return 1  # gcd(a, N)
