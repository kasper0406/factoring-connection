"""Lehman / stationary-point side of the dichotomy (notes/coverage38_results.md §4).

* lehman_factor(N, h0, X0): OPERATIONAL (sees only N).  For weights 1 <= A <= B <= h0
  (gcd 1) it scans S = ceil(2 sqrt(ABN)) .. + X0 and tests whether S^2 - 4ABN = D^2.
  The true S = Bp + Aq satisfies S^2 - 4ABN = (Bp - Aq)^2, and
  S - 2 sqrt(ABN) <= (Bp - Aq)^2 / (4 sqrt(ABN)).  Cost: sum over weights of X0 + 1
  isqrt tests (the BSGS-shared version costs ~ h0^2 + sqrt(h0^2 X0)).
* lehman_score(p, q): HARNESS.  L(p,q) = min over 1 <= B <= Bmax, A in {floor, ceil}(Bp/q)
  of max(B^2, delta^2 / (4 sqrt(ABN))), delta = Bp - Aq.  Dirichlet gives L <= C P^{2/5};
  "Lehman-easy at gain G" means L <= P^{2/5} / G.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Optional

import numpy as np


def _proper(d: int, N: int) -> Optional[int]:
    return d if 1 < d < N and N % d == 0 else None


def lehman_factor(N: int, h0: int, X0: int) -> Optional[int]:
    """OPERATIONAL.  Returns a proper divisor or None."""
    if N < 4:
        return None
    if N % 2 == 0:
        return 2
    for B in range(1, h0 + 1):
        for A in range(1, B + 1):
            if gcd(A, B) != 1:
                continue
            k4 = 4 * A * B * N
            s0 = isqrt(k4)
            if s0 * s0 < k4:
                s0 += 1
            for S in range(s0, s0 + X0 + 1):
                d2 = S * S - k4
                D = isqrt(d2)
                if D * D == d2:
                    for c in ((S + D) // 2, (S - D) // 2):
                        r = _proper(gcd(c, N), N)
                        if r:
                            return r
    return None


def lehman_score_exact_le(p: int, q: int, X0: int, Bmax: int) -> bool:
    """HARNESS, exact: is there B <= min(Bmax, isqrt(X0)) and A >= 1 with
    delta^4 <= 16 X0^2 A B N (i.e. delta^2/(4 sqrt(ABN)) <= X0)?"""
    N = p * q
    for B in range(1, min(Bmax, isqrt(X0)) + 1):
        a0 = (B * p) // q
        for A in (a0, a0 + 1):
            if A < 1:
                continue
            d = B * p - A * q
            if d ** 4 <= 16 * X0 * X0 * A * B * N:
                return True
    return False


def lehman_score_vec(P: np.ndarray, Q: np.ndarray, Bmax: int) -> np.ndarray:
    """HARNESS.  L(p,q) as float64 (a classification statistic; thresholds are
    re-checked exactly near ties by lehman_classify)."""
    N = P.astype(np.float64) * Q.astype(np.float64)
    best = np.full(len(P), np.inf)
    for B in range(1, Bmax + 1):
        a0 = (B * P) // Q
        for A in (a0, a0 + 1):
            ok = A >= 1
            d = (B * P - A * Q).astype(np.float64)
            Af = np.maximum(A, 1).astype(np.float64)
            w = d * d / (4.0 * np.sqrt(Af * B * N))
            val = np.maximum(float(B * B), w)
            val = np.where(ok, val, np.inf)
            np.minimum(best, val, out=best)
    return best


def lehman_classify(P: np.ndarray, Q: np.ndarray, X0: int, Bmax: int,
                    score: np.ndarray = None) -> np.ndarray:
    """Exact boolean 'L(p,q) <= X0' (float score + exact re-check within 1e-6)."""
    if score is None:
        score = lehman_score_vec(P, Q, Bmax)
    easy = score <= X0
    near = np.nonzero(np.abs(score - X0) <= 1e-6 * max(1, X0))[0]
    for i in near:
        easy[i] = lehman_score_exact_le(int(P[i]), int(Q[i]), X0, Bmax)
    return easy
