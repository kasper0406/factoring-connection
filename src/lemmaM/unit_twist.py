"""OPERATIONAL (exponential-time diagnostic algorithm) + HARNESS for E2: the unit-twist
involution of the continued-fraction relation set of sqrt(N).

Relations: A_k^2 - N B_k^2 = (-1)^{k+1} Q_{k+1}, i.e. A_k^2 = (-1)^{k+1} Q_{k+1} (mod N).
After one period P the fundamental unit eps = A_{P-1} + B_{P-1} sqrt(N) acts:
  Q_{k+P} = Q_k  and  A_{k+P} = a * A_k (mod N) with a = A_{P-1} mod N, a^2 = (-1)^P.
For a^2 = 1, a != +-1 (mod N) the map k -> k+P is an x-preserving, beta-flipping
involution on the period-doubled relation set (the Dixon epsilon), but it displaces the
relation by one full period (height factor eps ~ exp(regulator)).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional


def cf_period(N: int, cap: int) -> Optional[Dict]:
    """Continued fraction of sqrt(N) through two periods (N not a square).
    Returns dict(P, a, sign, checks) or None if the period exceeds cap."""
    a0 = isqrt(N)
    if a0 * a0 == N:
        return None
    m, d, a = 0, 1, a0
    A_prev, A = 1, a0 % N          # A_{-1}, A_0 (mod N)
    Qs: List[int] = []
    As: List[int] = [A]
    P = None
    k = 0
    while True:
        m = d * a - m
        d = (N - m * m) // d
        a = (a0 + m) // d
        k += 1
        Qs.append(d)              # Q_k (k >= 1)
        if P is None and d == 1:
            P = k
        A_prev, A = A, (a * A + A_prev) % N
        As.append(A)              # As[k] = A_k mod N
        if P is not None and k >= 2 * P:
            break
        if k > 2 * cap:
            return None
    unit_a = As[P - 1]
    sign = (-1) ** P              # norm of the fundamental unit
    ok_Q = all(Qs[i] == Qs[i + P] for i in range(P))
    ok_A = all(As[i + P] == unit_a * As[i] % N for i in range(P))
    ok_unit = (unit_a * unit_a - sign) % N == 0
    return dict(P=P, a=unit_a, sign=sign, ok_Q=ok_Q, ok_A=ok_A, ok_unit=ok_unit)


def classify_unit(N: int, a: int, sign: int, fs: List[int]) -> str:
    """HARNESS: 'nontrivial' if a^2 = 1 and a != +-1 mod N; 'trivial'; 'norm-1'."""
    if sign == -1:
        return "norm-1"
    g = gcd(a - 1, N)
    return "nontrivial" if 1 < g < N else "trivial"
