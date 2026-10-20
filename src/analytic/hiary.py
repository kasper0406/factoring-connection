"""Hiary's deterministic N^{1/3} factoring algorithm (arXiv 1408.2608, Algorithm 1).

OPERATIONAL: `hiary_factor(N)` sees only N. Exact integer arithmetic throughout.

Each loop iteration tests the whole block [x - H, x + H] for divisors of N in
O(log N) operations:
  * take the continued-fraction convergent b/q of N/x^2 with the largest q <= 4H
    (Lemma 3.1: N/x^2 = b/q + e2/(q q'), q <= 4H <= q', |e2| < 1);
  * set a = [qN/x] (nearest integer, ties down);
  * every divisor x + h of N with |h| <= H is an integer root of
    g(h) = (qN - a x) + (b x - a) h + b h^2   (Lemma 3.2).
The block half-width is H = floor((17N)^{-1/3} x), so blocks grow linearly in x and
each dyadic range [X, 2X] needs about (17N)^{1/3} ln 2 / 2 blocks; the total is
about (17N)^{1/3} ln(N^{1/6}/17^{1/3}) / 2 = 0.214 N^{1/3} ln N blocks.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, Iterator, List, Optional, Tuple

import gmpy2


def icbrt_floor(m: int) -> int:
    """floor(m^{1/3}) for m >= 0."""
    r = int(gmpy2.iroot(gmpy2.mpz(m), 3)[0])
    return r


def icbrt_ceil(m: int) -> int:
    r = icbrt_floor(m)
    return r if r ** 3 == m else r + 1


def convergents(num: int, den: int) -> Iterator[Tuple[int, int]]:
    """Continued-fraction convergents (b_j, q_j) of num/den, den > 0."""
    b0, q0, b1, q1 = 0, 1, 1, 0  # h_{-2}, k_{-2}, h_{-1}, k_{-1}
    while den:
        t, r = divmod(num, den)
        b0, q0, b1, q1 = b1, q1, t * b1 + b0, t * q1 + q0
        yield b1, q1
        num, den = den, r


def best_convergent(num: int, den: int, qmax: int) -> Tuple[int, int, int]:
    """Convergent b/q of num/den with the largest q <= qmax; also returns #CF steps."""
    best = None
    steps = 0
    for b, q in convergents(num, den):
        steps += 1
        if q > qmax:
            break
        best = (b, q)
    assert best is not None  # q_0 = 1 <= qmax
    return best[0], best[1], steps


def nearest_int(num: int, den: int) -> int:
    """[num/den], nearest integer, ties rounded down (Hiary's convention). den > 0."""
    # floor((2 num + den - 1) / (2 den)) rounds half down
    return (2 * num + den - 1) // (2 * den)


def integer_roots_quadratic(c0: int, c1: int, c2: int) -> List[int]:
    """All integer h with c0 + c1 h + c2 h^2 = 0 (polynomial not identically zero)."""
    if c2 == 0:
        if c1 == 0:
            return []
        return [-c0 // c1] if c0 % c1 == 0 else []
    disc = c1 * c1 - 4 * c2 * c0
    if disc < 0:
        return []
    s = isqrt(disc)
    if s * s != disc:
        return []
    out = []
    for num in (-c1 + s, -c1 - s):
        if num % (2 * c2) == 0:
            out.append(num // (2 * c2))
    return sorted(set(out))


def hiary_factor(N: int, *, trial_division: bool = True, stats: Optional[Dict[str, int]] = None,
                 lo: Optional[int] = None, hi: Optional[int] = None) -> Optional[int]:
    """Return a nontrivial factor of N, or None if N is prime (Hiary Algorithm 1).

    Optional [lo, hi] restricts the search to divisors in that window (the
    "partial information" variant of Hiary's Sec. 2): the trial-division and loop
    phases are clipped to the window.  stats (if given) receives counters:
    blocks, cf_steps, trial_divisions, quad_solves, divisor_tests.
    """
    if stats is None:
        stats = {}
    for k in ("blocks", "cf_steps", "trial_divisions", "quad_solves", "divisor_tests"):
        stats.setdefault(k, 0)
    if N < 4:
        return None
    sN = isqrt(N)
    x0 = min(icbrt_ceil(17 * N), sN)
    wlo = 2 if lo is None else max(2, lo)
    whi = sN if hi is None else min(sN, hi)
    # trial division on [wlo, min(x0, whi)]
    if trial_division:
        for k in range(wlo, min(x0, whi) + 1):
            stats["trial_divisions"] += 1
            if N % k == 0:
                stats["found_in"] = "trial"
                return k
    x = x0 + 2 if wlo <= x0 else wlo + 1
    H = 1 if wlo <= x0 else max(1, icbrt_floor(x ** 3 // (17 * N)))
    while x - H <= whi:
        stats["blocks"] += 1
        b, q, steps = best_convergent(N, x * x, 4 * H)
        stats["cf_steps"] += steps
        a = nearest_int(q * N, x)
        c0 = q * N - a * x
        c1 = b * x - a
        c2 = b
        stats["quad_solves"] += 1
        for h in integer_roots_quadratic(c0, c1, c2):
            if -H <= h <= H:
                d = x + h
                stats["divisor_tests"] += 1
                if 1 < d < N and N % d == 0:
                    stats["found_in"] = "block"
                    return d
        x = x + 2 * H + 1
        H = icbrt_floor(x ** 3 // (17 * N))
        if H < 1:
            H = 1
    return None


def predicted_blocks(N: int) -> float:
    """Continuum prediction for the number of loop blocks when N is prime (fit/plot only)."""
    import math
    c = (17 * N) ** (1.0 / 3.0)
    return 0.5 * c * math.log(math.sqrt(N) / c)
