"""sicheck HARNESS / DIAGNOSTIC helpers. May use p, q (instance labels). Never imported by src/sicheck/core.py."""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from src.common.instances import balanced_semiprime


def is_blum_diag(p: int, q: int) -> bool:
    return p % 4 == 3 and q % 4 == 3


def semiprimes_by_type_diag(nbits: int, count: int, blum: bool, start_seed: int = 0) -> List[Tuple[int, int, int, int]]:
    """First `count` balanced semiprimes (seed order) of the requested Blum type: (seed, N, p, q)."""
    out = []
    s = start_seed
    while len(out) < count:
        inst = balanced_semiprime(nbits, s)
        p, q = inst.labels["factors"]
        if is_blum_diag(p, q) == blum:
            out.append((s, inst.N, p, q))
        s += 1
    return out


def smooth_numbers_upto(X: int, primes) -> np.ndarray:
    """All positive y-smooth integers <= X, sorted (exact; int64 container). Copied pattern from abcsunit.core."""
    vals = np.array([1], dtype=np.int64)
    for p in sorted(primes):
        pieces = [vals]
        cur = vals
        while True:
            cur = cur[cur <= X // p] * p
            if cur.size == 0:
                break
            pieces.append(cur)
        vals = np.concatenate(pieces)
    vals.sort()
    return vals


class NullDensity:
    """Size-matched null for H_S: the exact density of y-smooth integers in a window around v.

    p(v) = #{smooth n in [v - D, v + D]} / (2D + 1), D = max(200, v // 20). Parity-split versions count only even
    (resp. odd) integers in the window. Exact counts (searchsorted on the exact sorted list)."""

    def __init__(self, hmax: int, primes):
        self.vals = smooth_numbers_upto(int(hmax * 1.06) + 400, primes)
        self.even = self.vals[self.vals % 2 == 0]
        self.odd = self.vals[self.vals % 2 == 1]

    @staticmethod
    def _win(v: int):
        D = max(200, v // 20)
        return max(1, v - D), v + D

    def p(self, v: int) -> float:
        a, b = self._win(v)
        c = np.searchsorted(self.vals, b, side="right") - np.searchsorted(self.vals, a, side="left")
        return c / (b - a + 1)

    def p_parity(self, v: int, parity: int) -> float:
        a, b = self._win(v)
        arr = self.even if parity == 0 else self.odd
        c = np.searchsorted(arr, b, side="right") - np.searchsorted(arr, a, side="left")
        # number of integers of the given parity in [a, b]
        n_par = (b - a + 1) // 2 + (1 if ((b - a + 1) % 2 == 1 and a % 2 == parity) else 0)
        return c / n_par
