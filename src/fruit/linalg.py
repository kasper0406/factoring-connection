"""OPERATIONAL: GF(2) linear algebra for hybrid relation matrices (rows = int bitmasks).

Lemma HY (notes/fruit_results.md §1): for a base matrix A (M rows) and h hidden-root rows D,
the projection pi of the left kernel K of [A; D] onto the D-coordinates has rank
    r = h - (rank[A; D] - rank A) = h - rank(D mod rowspace A),
so a uniform element of K uses a hidden-root row with probability exactly 1 - 2^{-r}.
"""
from __future__ import annotations

import random
from typing import List, Optional, Tuple

from ..nfsfruit.gf2 import left_kernel, rank_of


def active_columns(rows: List[int]) -> int:
    acc = 0
    for r in rows:
        acc |= r
    return bin(acc).count("1")


def hybrid_rank(base: List[int], desc: List[int]) -> dict:
    """Exact ranks for Lemma HY.  r = rank of pi: K -> F_2^h."""
    rA = rank_of(base)
    rAD = rank_of(list(base) + list(desc))
    h = len(desc)
    cact = active_columns(list(base) + list(desc))
    cA = active_columns(list(base))
    return dict(rankA=rA, rankAD=rAD, h=h, r=h - (rAD - rA), cact=cact, cA=cA,
                deficit=cact - rA, deficitA=cA - rA, newcols=cact - cA)


def kernel_basis(rows: List[int]) -> List[int]:
    return left_kernel(rows)[1]


def sample_uniform(kern: List[int], rng: random.Random) -> int:
    """Uniform element of span(kern) (zero allowed)."""
    S = 0
    for b in kern:
        if rng.getrandbits(1):
            S ^= b
    return S


def sample_nonzero(kern: List[int], rng: random.Random) -> Optional[int]:
    if not kern:
        return None
    while True:
        S = sample_uniform(kern, rng)
        if S:
            return S


def sample_forced(kern: List[int], M: int, rng: random.Random, tries: int = 400) -> Optional[int]:
    """Uniform element of K minus ker(pi), pi = coordinates >= M (the hidden-root rows).
    Selection uses only the kernel (public data), never a hidden root: transcript-safe."""
    if not kern or not any(b >> M for b in kern):
        return None
    for _ in range(tries):
        S = sample_uniform(kern, rng)
        if S >> M:
            return S
    return None
