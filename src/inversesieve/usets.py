"""Exact uncovered sets.

  U(N) = {x in [a0, b0] : x mod l not in R_l(N) for every band prime l}     (OPERATIONAL-SAFE: sees N only)
  partners(p) = {q in (p, 2p] : l does not divide Bp - Aq for all (A,B), l} (HARNESS: uses p)

Both are computed by successive filtering with numpy int64 (exact integer arithmetic; all values < 2^40).
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from .design import candidate_interval, forbidden_mask


def U_of_N(N: int, primes: Sequence[int], V, interval: Tuple[int, int] = None,
           order: str = "asc") -> np.ndarray:
    a0, b0 = candidate_interval(N) if interval is None else interval
    x = np.arange(a0, b0 + 1, dtype=np.int64)
    for l in (sorted(primes) if order == "asc" else list(primes)):
        if x.size == 0:
            break
        m = forbidden_mask(N, l, V)
        x = x[~m[x % l]]
    return x


def U_bruteforce(N: int, primes, V, interval=None) -> List[int]:
    """Reference: direct test of l | B x^2 - A N (slow; tests only)."""
    a0, b0 = candidate_interval(N) if interval is None else interval
    out = []
    for x in range(a0, b0 + 1):
        if not any((B * x * x - A * N) % l == 0 for l in primes for (A, B) in V):
            out.append(x)
    return out


def is_pair_uncovered(p: int, q: int, primes, V) -> bool:
    """HARNESS: (p, q) uncovered iff no l | Bp - Aq."""
    return not any((B * p - A * q) % l == 0 for l in primes for (A, B) in V)


def uncovered_partners(p: int, primes, V, qlo: int, qhi: int) -> np.ndarray:
    """HARNESS: integers q in [qlo, qhi] with l not dividing Bp - Aq for all weights and band primes.
    For l not dividing p,B,A: l | Bp - Aq <=> q = B p A^{-1} (mod l)."""
    q = np.arange(qlo, qhi + 1, dtype=np.int64)
    for l in sorted(primes):
        if q.size == 0:
            break
        m = np.zeros(l, dtype=bool)
        for (A, B) in V:
            if A % l == 0:
                if (B * p) % l == 0:
                    m[:] = True
                continue
            m[(B * p * pow(A % l, -1, l)) % l] = True
        q = q[~m[q % l]]
    return q
