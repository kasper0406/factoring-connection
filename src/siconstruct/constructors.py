"""Toy constructors for y-smooth u in I_k = [kN, kN + h] (agent `siconstruct`).  Operational: see only N, h, K
and N-oblivious smooth lists.  Every output is verified exactly (u = kN + w, 0 <= w <= h, u y-smooth).

* brute_interval: sieve a window of I_k (the baseline; cost = integers sieved).
* mitm_construct: multiplicative birthday u = a*b, a in A, b in B (sorted lists of N-oblivious smooth numbers),
  ONE pair of lists serving ALL multipliers k at once (the amortisation question).
* aim2d_construct: one constructed smooth factor s, the multiplier k and offset w aimed so that s | kN + w,
  and the cofactor m = (kN + w)/s TESTED for smoothness (the "2-D aim"; lattice sieving with a smooth special-q).
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter
from typing import Dict, List, Optional, Sequence

from .core import is_smooth, sieve_window


def brute_interval(N: int, k: int, h: int, P: Sequence[int], W: int) -> Dict:
    """Sieve [kN, kN + min(W, h+1)).  Cost = integers sieved."""
    W = min(W, h + 1)
    found, cost = sieve_window(k * N, W, P)
    return {"k": k, "outputs": found, "cost": cost}


def mitm_construct(N: int, h: int, A: Sequence[int], B: Sequence[int]) -> Dict:
    """All pairs (a, b) in A x B with a*b in [kN, kN + h] for some k >= 1.

    For each a and each k whose b-interval [ceil(kN/a), floor((kN + h)/a)] meets [min B, max B], bisect the
    sorted B.  Counters:
      gen      = |A| + |B|                 (list generation + sorting, L^{o(1)} each)
      queries  = number of (a, k) pairs bisected (the K-dependent part of this implementation)
      points   = number of integers b in those intervals (cost of an ideal Bohr-set enumerator + hash lookup)
      hits     = pairs found;  outputs = distinct u.
    """
    Bs = sorted(B)
    bmin, bmax = Bs[0], Bs[-1]
    queries = 0
    points = 0
    hits = []
    for a in A:
        k_lo = max(1, (a * bmin - h) // N)
        k_hi = (a * bmax) // N
        for k in range(k_lo, k_hi + 1):
            kN = k * N
            lo = -(-kN // a)
            hi = (kN + h) // a
            if hi < bmin or lo > bmax or hi < lo:
                continue
            queries += 1
            points += min(hi, bmax) - max(lo, bmin) + 1
            i = bisect_left(Bs, lo)
            j = bisect_right(Bs, hi)
            for b in Bs[i:j]:
                w = a * b - kN
                assert 0 <= w <= h
                hits.append((a * b, k, w))
    outs = {u: (k, w) for (u, k, w) in hits}
    return {"gen": len(A) + len(B), "queries": queries, "points": points, "hits": len(hits),
            "outputs": sorted(outs.items())}


def aim2d_construct(N: int, h: int, K: int, S: Sequence[int], P: Sequence[int],
                    max_per_sk: Optional[int] = None) -> Dict:
    """For each constructed smooth s and k <= K, every w in [0, h] with s | kN + w (w = (-kN mod s) + j s), at most
    `max_per_sk` of them per (s, k); test the cofactor m = (kN + w)/s for P-smoothness.
    Returns outputs and a histogram of the tested cofactor bit lengths."""
    tested = 0
    outs = {}
    cof_bits: Counter = Counter()
    for s in S:
        for k in range(1, K + 1):
            kN = k * N
            w = (-kN) % s
            j = 0
            while w <= h and (max_per_sk is None or j < max_per_sk):
                m = (kN + w) // s
                tested += 1
                j += 1
                cof_bits[m.bit_length()] += 1
                if is_smooth(m, P):
                    outs[kN + w] = (k, w, s)
                w += s
    return {"gen": len(S), "loop": len(S) * K, "tested": tested, "outputs": sorted(outs.items()),
            "cofactor_bits_min": min(cof_bits) if cof_bits else None,
            "cofactor_bits": dict(cof_bits)}
