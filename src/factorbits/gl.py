"""Prop GLM (notes/factorbits_results.md §4.3): a single-bit function of the integer alone that is as hard
as factoring, via Goldreich-Levin with the GL randomness embedded in a cofactor.

    f(M) = parity( P+(M) AND M )          (P+ = largest prime factor; AND = bitwise)

For balanced N = pq (q the larger factor, L = bitlen(q), so q >= 2^(L-1)) and r in {0,1}^(L-1), put
k = r * N^{-1} mod 2^(L-1) and M = N k.  Then k < 2^(L-1) <= q, so P+(M) = q ALWAYS, M = r (mod 2^(L-1)), and
f(M) = <q mod 2^(L-1), r> XOR bit_{L-1}(M).  A predictor for f with advantage eps on the distribution of M
therefore gives the GL oracle B(r) = pred(M) XOR bit_{L-1}(M) for x = q mod 2^(L-1), with the same advantage
(no loss).  q = 2^(L-1) + x.  (L = ceil(n/2) under the equal-length convention; otherwise try both values.)

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
    transform; the J = 0 query r = e_i is included, which is harmless).  Decodes x = q mod 2^(L-1) with
    (L-1) * 2^l predictor queries; returns (q, queries) or None."""
    rng = random.Random(seed)
    L1 = L - 1
    mask = (1 << L1) - 1
    Ninv = pow(N, -1, 1 << L1)
    s = [rng.getrandbits(L1) for _ in range(l)]
    rJ = [0] * (1 << l)
    for J in range(1, 1 << l):
        low = J & -J
        rJ[J] = rJ[J ^ low] ^ s[low.bit_length() - 1]
    queries = 0
    bits_sign = np.zeros((L1, 1 << l), dtype=np.int64)
    for i in range(L1):
        vals = np.empty(1 << l, dtype=np.int64)
        for J in range(1 << l):
            r = rJ[J] ^ (1 << i)
            k = (r * Ninv) & mask
            queries += 1
            if k == 0:
                vals[J] = 0
                continue
            M = N * k
            vals[J] = pred(M) ^ ((M >> L1) & 1)
        bits_sign[i] = fwht(1 - 2 * vals)          # sum_J (-1)^{vals_J + <sigma, J>}
    for sigma in range(1 << l):
        x = 0
        for i in range(L1):
            if bits_sign[i, sigma] < 0:
                x |= 1 << i
        q = x | (1 << L1)
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
