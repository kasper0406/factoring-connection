"""Exact enumeration of all fac-relations in a box (HARNESS / measurement code for E1).

Counts every triple (u, v, w) with v P-smooth, 1 <= v <= V, 1 <= |w| <= W, u = vN + w > 0 P-smooth and
|w| P-smooth.  Exact: numpy int64 window sieve (all values < 2^62 asserted), every hit re-verified by
Python-int trial division.  Uses only N and the public factor base.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Sequence

import numpy as np

from .core import factor_over


def smooth_upto(X: int, P: Sequence[int]) -> List[int]:
    """All P-smooth integers in [1, X], sorted (exact generation by multiplication)."""
    out = [1]
    for p in P:
        new = []
        for m in out:
            x = m * p
            while x <= X:
                new.append(x)
                x *= p
        out.extend(new)
    return sorted(out)


def smooth_mask_window(lo: int, length: int, P: Sequence[int]) -> np.ndarray:
    """Boolean mask of P-smooth integers in [lo, lo + length) (lo >= 1), exact int64 sieve."""
    assert lo >= 1 and lo + length < 2 ** 62
    arr = np.arange(lo, lo + length, dtype=np.int64)
    for p in P:
        pk = p
        while pk < lo + length:
            start = (-lo) % pk
            if start < length:
                arr[start::pk] //= p
            pk *= p
            if pk > 2 ** 62:
                break
    return arr == 1


def relations_in_box(N: int, P: Sequence[int], V: int, W: int) -> Dict:
    vs = smooth_upto(V, P)
    smooth_w = np.zeros(W + 1, dtype=bool)
    smooth_w[1:] = smooth_mask_window(1, W, P)
    n_w = int(smooth_w.sum())
    rels = []
    for v in vs:
        lo = v * N - W
        assert lo >= 1
        mask = smooth_mask_window(lo, 2 * W + 1, P)
        idx = np.nonzero(mask)[0]
        for i in idx:
            w = int(i) - W
            if w == 0 or not smooth_w[abs(w)]:
                continue
            u = v * N + w
            assert factor_over(u, P) is not None and factor_over(w, P) is not None
            rels.append((u, v, w))
    coprime = [r for r in rels if gcd(r[0], r[1]) == 1]
    return {"n_v": len(vs), "n_w": n_w, "relations": rels, "count": len(rels),
            "count_coprime_uv": len(coprime)}
