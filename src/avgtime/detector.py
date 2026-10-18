"""OPERATIONAL: small-order detector (notes/avgtime_results.md, Lemma S).

Sees only N and the public bound T.  Finds a proper divisor of N whenever the
smallest order d_min = min_{r | N prime} ord_r(2) satisfies d_min <= T and the
primes r with ord_r(2) = d_min do not exhaust N (for N = pq: ord_p 2 != ord_q 2).

Baby-step giant-step on the exponent:
  m = ceil(sqrt T), babies 2^i (0 <= i < m), f(X) = prod (X - 2^i),
  giants 2^{jm} (1 <= j <= ceil(T/m));  G_j = gcd(f(2^{jm}), N).
  j* = first j with G_j > 1;  then gcd(2^k - 1, N) for k in ((j*-1)m, j*m].
Block j has G_j > 1 iff ((j-1)m, jm] contains a positive multiple of some
ord_r(2), so j* = ceil(d_min/m) and k = d_min is scanned; gcd(2^{d_min}-1, N)
is the product of the r with ord_r(2) | d_min, i.e. with ord_r(2) = d_min.
Cost: O(sqrt T) ring operations plus the product/remainder trees.
Not new: elementary BSGS order search (cf. Pollard 1974, Pollard--Strassen).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Optional

from src.common.metrics import Metrics
from src.dcltheorem.stage3 import eval_baby_poly_at_giants


def _proper(d: int, N: int) -> bool:
    return 1 < d < N and N % d == 0


def small_order_detector(N: int, T: int, met: Optional[Metrics] = None) -> Optional[int]:
    met = met or Metrics()
    if N < 4 or T < 1:
        return None
    if N % 2 == 0:
        return 2
    m = isqrt(T - 1) + 1 if T > 1 else 1          # ceil(sqrt T)
    J = -(-T // m)                                  # ceil(T/m); J*m >= T
    g = pow(2, m, N)
    giants, x = [], 1
    for _ in range(J):
        x = x * g % N
        met.inc("mulN")
        giants.append(x)
    vals = eval_baby_poly_at_giants(N, 2, m, giants, met)
    for j, v in enumerate(vals, start=1):
        met.inc("gcdN")
        G = gcd(v, N)
        if G == 1:
            continue
        # first block with a hit: scan k in ((j-1)m, jm]
        k0 = (j - 1) * m + 1
        y = pow(2, k0, N)
        for k in range(k0, j * m + 1):
            met.inc("gcdN")
            d = gcd(y - 1, N)
            if _proper(d, N):
                return d
            y = y * 2 % N
            met.inc("mulN")
        return None          # every prime of N has order d_min: detector fails
    return None
