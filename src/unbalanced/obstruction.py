"""DIAGNOSTIC: the curvature obstruction (Proposition O, notes/unbalanced_results.md §5)
and the exponent profile of §4.  Exact integer arithmetic for the coverage counts;
floats only for the display profile.

Model.  A *Lehman-linear target* tau = (alpha, beta, c0, M, r, g0, W) assigns to an
integer candidate p == r (mod M) the real value
    g(p) = (alpha p + beta N / p + c0) / M^2
and a giant with a baby list of W consecutive exponents *structurally covers* p iff
g0 <= g(p) < g0 + W.  Proposition O: for beta != 0 the covered candidates in
[x1, x2] number at most 2 sqrt(W x2^3 / (|beta| N)) + 2.

covered_count() is the brute-force count (exact rationals: g(p) >= g0 <=>
alpha p^2 + beta N + c0 p >= g0 M^2 p for p > 0).
"""
from __future__ import annotations

from math import isqrt, log
from typing import Iterable, Tuple


def covered_count(N: int, alpha: int, beta: int, c0: int, M: int, r: int, g0: int, W: int,
                  x1: int, x2: int) -> int:
    """#{p in [x1, x2], p == r (mod M), p > 0 : g0 <= g(p) < g0 + W}, exactly."""
    M2 = M * M
    cnt = 0
    p = x1 + ((r - x1) % M)
    while p <= x2:
        if p > 0:
            v = alpha * p * p + beta * N + c0 * p          # = p * M^2 * g(p)
            if g0 * M2 * p <= v < (g0 + W) * M2 * p:
                cnt += 1
        p += M
    return cnt


def prop_o_bound_sq(N: int, beta: int, W: int, x2: int) -> Tuple[int, int]:
    """Return (S, 2) with bound = 2*sqrt(S) + 2, S = W x2^3 / (|beta| N) rounded UP
    (so the integer test  cnt <= 2*isqrt_ceil(S) + 2  is implied by the real bound)."""
    num = W * x2 ** 3
    den = abs(beta) * N
    S = -(-num // den)
    return S, 2


def prop_o_bound_int(N: int, beta: int, W: int, x2: int) -> int:
    S, _ = prop_o_bound_sq(N, beta, W, x2)
    s = isqrt(S)
    if s * s < S:
        s += 1
    return 2 * s + 2


def best_window_coverage(N: int, alpha: int, beta: int, c0: int, M: int, r: int, W: int,
                         x1: int, x2: int) -> int:
    """max over g0 of covered_count (scan g0 over the integer range of g on the class)."""
    M2 = M * M
    vals = []
    p = x1 + ((r - x1) % M)
    while p <= x2:
        if p > 0:
            # floor(g(p))
            vals.append((alpha * p * p + beta * N + c0 * p) // (M2 * p))
        p += M
    if not vals:
        return 0
    vals.sort()
    best, j = 0, 0
    for i in range(len(vals)):            # sliding window on sorted floor values
        while vals[i] - vals[j] >= W:
            j += 1
        best = max(best, i - j + 1)
    return best


# ---------------------------------------------------------------------------
# exponent profile (display; floats allowed)
# ---------------------------------------------------------------------------
def e_dcl(theta: float) -> float:
    """Regime I/II exponent at ratio N^theta (Theorems U, U_theta): (1+theta)/6."""
    return (1.0 + theta) / 6.0


def e_ps(theta: float) -> float:
    """Pollard-Strassen exponent at ratio N^theta: p = N^{(1-theta)/2}, cost p^{1/2}."""
    return (1.0 - theta) / 4.0


def e_profile(theta: float) -> float:
    return min(e_dcl(theta), e_ps(theta))


def frac_theta_below(t: float) -> float:
    """Fraction of prime pairs pq <= X with N^{1/3} <= p <= sqrt N whose ratio exponent
    theta = log(q/p)/log N is <= t, asymptotically (Mertens: density of p ~ N^a is
    proportional to 1/(a(1-a)) da, a = (1-theta)/2):  log2((1+t)/(1-t)) on [0, 1/3]."""
    return log((1 + t) / (1 - t)) / log(2.0)


def frac_covered_at(eta: float) -> float:
    """Fraction of that population factored at exponent <= 1/6 + eta by the profile."""
    t1 = min(1 / 3, 6 * eta)
    t2 = max(0.0, 1 / 3 - 4 * eta)
    if t1 >= t2:
        return 1.0
    return frac_theta_below(t1) + 1.0 - frac_theta_below(t2)


__all__ = ["covered_count", "prop_o_bound_int", "best_window_coverage", "best_set_coverage", "prop_o_set_bound_int", "e_profile",
           "frac_theta_below", "frac_covered_at"]


# ---------------------------------------------------------------------------
# Proposition O' (arbitrary baby exponent sets)
# ---------------------------------------------------------------------------
def best_set_coverage(N: int, alpha: int, beta: int, c0: int, M: int, r: int, m: int,
                      x1: int, x2: int) -> int:
    """Max over ALL m-element sets S of integers of #{x : floor(g(x)) in S}: the optimum
    for a single giant is the m most popular floor values."""
    M2 = M * M
    cnt = {}
    p = x1 + ((r - x1) % M)
    while p <= x2:
        if p > 0:
            v = (alpha * p * p + beta * N + c0 * p) // (M2 * p)
            cnt[v] = cnt.get(v, 0) + 1
        p += M
    return sum(sorted(cnt.values(), reverse=True)[:m])


def prop_o_set_bound_int(N: int, beta: int, m: int, x2: int) -> int:
    """2m + 8 sqrt(m x2^3 / (|beta| N)), rounded up (Proposition O', per giant)."""
    S = -(-(64 * m * x2 ** 3) // (abs(beta) * N))
    s = isqrt(S)
    if s * s < S:
        s += 1
    return 2 * m + s
