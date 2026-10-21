"""Prop GLM (notes/factorbits_results.md §4.3): a single-bit function of the integer alone that is as hard
as factoring, via Goldreich-Levin with the GL randomness embedded in a cofactor.

    f(M) = parity( P+(M) AND M )          (P+ = largest prime factor; AND = bitwise)

For balanced N = pq (q the larger factor, L = bitlen(q)) and r in {0,1}^L, put k = r * N^{-1} mod 2^L and
M = N k.  If P+(k) < q (all but O(1/L) of r) then P+(M) = q and f(M) = <q, r> mod 2.  A predictor for f
with advantage eps on the distribution of M therefore gives a GL oracle for x = q.

`reduce_factor` is OPERATIONAL given a predictor (it sees N, queries the predictor, checks candidates by
division).  `planted_predictor` is HARNESS: it computes f from the factorisation and flips it on a
pseudorandom (1/2 - eps)-fraction of inputs, fixed per input.
"""
from __future__ import annotations

import hashlib
import random
from typing import Callable, Optional

import numpy as np


def f_bit(M: int, Pplus: int) -> int:
    return bin(Pplus & M).count("1") & 1


def fwht(a: np.ndarray) -> np.ndarray:
    a = a.astype(np.int64).copy()
    h = 1
    while h < len(a):
        a = a.reshape(-1, 2 * h)
        x, y = a[:, :h].copy(), a[:, h:].copy()
        a[:, :h], a[:, h:] = x + y, x - y
        a = a.reshape(-1)
        h *= 2
    return a


def reduce_factor(N: int, L: int, pred: Callable[[int], int], l: int, seed: int = 0) -> Optional[tuple]:
    """Goldreich-Levin (Rackoff's pairwise-independent version, all guesses at once by a Walsh-Hadamard
    transform).  Makes L * 2^l predictor queries; returns (q, queries) or None."""
    rng = random.Random(seed)
    mask = (1 << L) - 1
    Ninv = pow(N, -1, 1 << L)
    s = [rng.getrandbits(L) for _ in range(l)]
    # r_J for all subsets J (index J as an l-bit integer)
    rJ = np.zeros(1 << l, dtype=object)
    for J in range(1, 1 << l):
        low = J & -J
        rJ[J] = rJ[J ^ low] ^ s[low.bit_length() - 1]
    queries = 0
    bits_sign = np.zeros((L, 1 << l), dtype=np.int64)
    for i in range(L):
        vals = np.empty(1 << l, dtype=np.int64)
        for J in range(1 << l):
            r = int(rJ[J]) ^ (1 << i)
            k = (r * Ninv) & mask
            queries += 1
            vals[J] = pred(N * k) if k else 0
        bits_sign[i] = fwht(1 - 2 * vals)          # sum_J (-1)^{vals_J + <sigma, J>}
    # candidate for guess sigma: bit i = 1 iff the transform is negative
    for sigma in range(1 << l):
        q = 0
        for i in range(L):
            if bits_sign[i, sigma] < 0:
                q |= 1 << i
        if 1 < q < N and N % q == 0:
            return q, queries
    return None


def planted_predictor(p: int, q: int, eps: float, key: str) -> Callable[[int], int]:
    """HARNESS.  True f(M) for M = p q k, computed from the known factors and a trial-division
    factorisation of the cofactor k, flipped when a keyed hash of M falls below 1/2 - eps."""
    N = p * q
    thr = int((0.5 - eps) * 2 ** 64)

    def Pplus_small(k: int) -> int:
        best, d = 1, 2
        while d * d <= k:
            while k % d == 0:
                best, k = max(best, d), k // d
            d += 1
        return max(best, k) if k > 1 else best

    def pred(M: int) -> int:
        k = M // N
        P = max(q, Pplus_small(k))
        bit = f_bit(M, P)
        h = int.from_bytes(hashlib.sha256(f"{key}|{M}".encode()).digest()[:8], "big")
        return bit ^ (h < thr)
    return pred
