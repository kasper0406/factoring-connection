"""OPERATIONAL (heuristic) Galois decoding: trace NUMBERS of many fixed curves -> residues of p (or s) mod l.

Inputs: N and trace integers a_p(E_i) (ordered, one prime) or a_N(E_i) (products).  The decoder only uses the
Haar/Chebotarev likelihood of galois.py; it never sees p, q or group orders and never factors.  Correctness of the
output is a HEURISTIC (independent full-image Frobenius classes for a fixed prime across many curves).
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
from math import isqrt
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

from src.tracedata.galois import class_counts, prod_likelihood


@lru_cache(maxsize=None)
def _logn(l: int) -> np.ndarray:
    n = class_counts(l).astype(float)
    with np.errstate(divide="ignore"):
        return np.log(n)          # [t, d]; column d = 0 is -inf (unused)


@lru_cache(maxsize=None)
def _logL(l: int) -> np.ndarray:
    L = prod_likelihood(l)
    with np.errstate(divide="ignore"):
        return np.log(L)          # [d, dq, u]


def ll_ordered(l: int, residues: Sequence[int]) -> np.ndarray:
    """Log-likelihood of d = p mod l (index d = 1..l-1; entry 0 = -inf) from trace residues a_p(E_i) mod l."""
    c = np.bincount(np.asarray(residues, dtype=np.int64) % l, minlength=l).astype(float)
    ll = c @ _logn(l)
    ll[0] = -np.inf
    return ll


def ll_product(l: int, Nmod: int, residues: Sequence[int]) -> np.ndarray:
    """Log-likelihood of d = p mod l (symmetric under d <-> N/d) from a_N(E_i) mod l."""
    c = np.bincount(np.asarray(residues, dtype=np.int64) % l, minlength=l).astype(float)
    LL = _logL(l)
    ll = np.full(l, -np.inf)
    for d in range(1, l):
        dq = Nmod * pow(d, -1, l) % l
        ll[d] = float(c @ LL[d, dq, :])
    return ll


def rank_p(l: int, residues: Sequence[int]) -> List[int]:
    ll = ll_ordered(l, residues)
    return [int(d) for d in np.argsort(-ll) if d != 0]


def rank_s(l: int, Nmod: int, residues: Sequence[int]) -> List[int]:
    """Ranked candidate values of s mod l (each s from its best d)."""
    ll = ll_product(l, Nmod, residues)
    best = {}
    for d in range(1, l):
        s = (d + Nmod * pow(d, -1, l)) % l
        if s not in best or ll[d] > best[s]:
            best[s] = ll[d]
    return [s for s, _ in sorted(best.items(), key=lambda kv: -kv[1])]


def _crt(res: Sequence[int], mods: Sequence[int]) -> Tuple[int, int]:
    x, m = 0, 1
    for r, l in zip(res, mods):
        t = ((r - x) * pow(m, -1, l)) % l
        x, m = x + m * t, m * l
    return x, m


def factor_from_ordered(N: int, traces: Sequence[int], ells: Sequence[int], top: int = 2,
                        max_lifts: int = 4) -> Optional[int]:
    """List-decode p mod l for each l, CRT with p odd, test divisibility of the first max_lifts lifts.
    Heuristic; needs prod(ells) * 2 * max_lifts > p.  Toy sizes."""
    lists = [rank_p(l, [t % l for t in traces])[:top] for l in ells]
    for combo in product(*lists):
        x, m = _crt(list(combo) + [1], list(ells) + [2])
        for k in range(max_lifts):
            p = x + k * m
            if 1 < p < N and N % p == 0:
                return p
    return None


def factor_from_products(N: int, aN_list: Sequence[int], ells: Sequence[int], top: int = 3,
                         max_lifts: int = 4) -> Optional[int]:
    """List-decode s mod l, CRT with s even, solve x^2 - s x + N = 0 over Z for the first max_lifts lifts."""
    lists = [rank_s(l, N % l, [a % l for a in aN_list])[:top] for l in ells]
    for combo in product(*lists):
        x, m = _crt(list(combo) + [0], list(ells) + [2])
        for k in range(max_lifts):
            s = x + k * m
            D = s * s - 4 * N
            if D < 0:
                continue
            r = isqrt(D)
            if r * r == D and (s - r) % 2 == 0:
                p = (s - r) // 2
                if 1 < p < N and N % p == 0:
                    return p
    return None
