"""genusleak continued fractions / infrastructure baby steps (OPERATIONAL): sees only kN.

The continued fraction of sqrt(kN) walks the principal cycle of reduced forms of discriminant 4kN:
A_{i-1}^2 - kN B_{i-1}^2 = (-1)^i Q_i. Recurrences follow Gower-Wagstaff, Math. Comp. 77 (2008) 551-588 (SQUFOF),
with P0 = isqrt(kN), Q0 = 1, Q1 = kN - P0^2.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

__all__ = ["cf_period", "reverse_cycle", "squfof_scan", "midpoint_split"]


def _is_square(x: int) -> Optional[int]:
    if x < 0:
        return None
    r = isqrt(x)
    return r if r * r == x else None


def cf_period(M: int, max_steps: int = 10 ** 8) -> Tuple[List[int], List[int]]:
    """Return (P, Q) lists over one full period of the CF of sqrt(M) (M not a square):
    Q[0] = 1, ..., Q[L] = 1 with L the period; P[i] pairs with Q[i+1] as in Gower-Wagstaff
    (P[i] = b_i Q[i] - P[i-1]). P has length L+1 with P[0] = P0."""
    P0 = isqrt(M)
    if P0 * P0 == M:
        raise ValueError("square")
    P = [P0]
    Q = [1, M - P0 * P0]
    i = 1
    while True:
        b = (P0 + P[i - 1]) // Q[i]
        Pi = b * Q[i] - P[i - 1]
        Qn = Q[i - 1] + b * (P[i - 1] - Pi)
        P.append(Pi)
        Q.append(Qn)
        if Q[i] == 1:  # period closes: Q_L = 1 with L = i
            return P[: i + 1], Q[: i + 1]
        i += 1
        if i > max_steps:
            raise RuntimeError("period too long")


def reverse_cycle(M: int, N: int, Pprev: int, r: int, max_steps: int = 10 ** 8) -> Tuple[Optional[int], int]:
    """SQUFOF reverse cycle from the square-root form (Q = r, P = Pprev).  Returns (nontrivial divisor of N or
    None, steps)."""
    P0 = isqrt(M)
    b0 = (P0 - Pprev) // r
    Pp = b0 * r + Pprev
    Qa, Qb = r, (M - Pp * Pp) // r
    steps = 0
    while True:
        b = (P0 + Pp) // Qb
        Pn = b * Qb - Pp
        Qn = Qa + b * (Pp - Pn)
        steps += 1
        if Pn == Pp:
            g = gcd(N, Pn)
            if 1 < g < N:
                return g, steps
            g = gcd(N, Qb)
            if 1 < g < N:
                return g, steps
            return None, steps
        Pp, Qa, Qb = Pn, Qb, Qn
        if steps > max_steps:
            return None, steps


def squfof_scan(N: int, k: int = 1) -> Dict:
    """Walk one full period of sqrt(kN). For every even i with Q_i a square, run the reverse cycle and record
    whether it yields a nontrivial factor ('proper square').  Returns period length, the list of proper square
    indices, the first proper index and the reverse-cycle steps used for it."""
    M = k * N
    P, Q = cf_period(M)
    L = len(Q) - 1
    proper, improper = [], []
    rev_steps_first = None
    for i in range(2, L + 1, 2):
        r = _is_square(Q[i])
        if r is None or r == 1:
            continue
        g, st = reverse_cycle(M, N, P[i - 1], r, max_steps=4 * L + 10)
        if g is not None:
            proper.append(i)
            if rev_steps_first is None:
                rev_steps_first = st
        else:
            improper.append(i)
    return {"L": L, "proper": proper, "improper": improper, "first": proper[0] if proper else None,
            "rev_steps_first": rev_steps_first}


def midpoint_split(N: int, k: int = 1) -> Dict:
    """The symmetry point of the principal period of sqrt(kN) (distance ~ R/2): for even L the ambiguous form
    with P_i = P_{i-1}; return gcd-derived nontrivial divisor of N if any."""
    M = k * N
    P, Q = cf_period(M)
    L = len(Q) - 1
    if L % 2:
        return {"L": L, "parity": "odd", "split": None}
    for i in range(1, L + 1):
        if P[i] == P[i - 1]:
            for x in (P[i], Q[i]):
                g = gcd(N, x)
                if 1 < g < N:
                    return {"L": L, "parity": "even", "index": i, "split": g, "Q": Q[i]}
            return {"L": L, "parity": "even", "index": i, "split": None, "Q": Q[i]}
    return {"L": L, "parity": "even", "index": None, "split": None}
