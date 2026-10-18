"""Counts of reduced primitive forms of discriminant -4N with small first coefficient.  OPERATIONAL: N only.

Form (a, b, c), b^2 - 4ac = -4N, reduced: |b| <= a <= c, and b >= 0 if |b| = a or a = c; primitive:
gcd(a, b, c) = 1.  b must be even, b = 2 beta, beta^2 = -N (mod a).  For a <= A_MAX << sqrt(N) we have
c = (beta^2 + N)/a > a, so only the |b| = a boundary needs the sign rule.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Tuple

A_MAX = 50
_PRIMES = {x for x in range(2, A_MAX + 1) if all(x % d for d in range(2, int(x ** 0.5) + 1))}


def small_forms(N: int, a_max: int = A_MAX) -> List[Tuple[int, int, int]]:
    out = []
    for a in range(1, a_max + 1):
        # b in (-a, a], b even: beta in (-a/2, a/2]
        lo = -((a - 1) // 2) if a % 2 else -(a // 2) + 1
        hi = a // 2
        for beta in range(lo, hi + 1):
            b = 2 * beta
            num = beta * beta + N
            if num % a:
                continue
            c = num // a
            if c < a:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            out.append((a, b, c))
    return out


def form_features(N: int) -> Dict[str, int]:
    F = small_forms(N)
    nonp = [f for f in F if f[0] >= 2]
    return {
        "cl_cnt": len(F),
        "cl_cnt_odd": sum(1 for f in F if f[0] % 2),
        "cl_cnt_prime": sum(1 for f in F if f[0] in _PRIMES),
        "cl_cnt_le10": sum(1 for f in F if f[0] <= 10),
        "cl_cnt_le25": sum(1 for f in F if f[0] <= 25),
        "cl_min_a": min((f[0] for f in nonp if f[0] > 2), default=A_MAX + 1),
    }
