"""Degree-j floor-sum ladder (reduction L_j; notes/analytic_results.md Sec. 3).

Claim (proved in the notes, cost model = calls to an exact degree-j floor-sum oracle):
if floor sums  F(P; h1, h2) = sum_{h1 <= h < h2} floor(P(h))  of rational degree-j polynomials
are available, all divisors of N in [2B/c... , sqrt N] are found with O(N^{1/(j+2)} log N)
blocks and O(N^{1/(j+2)} log^2 N) oracle calls, B = N^{1/(j+2)}.

Block [x, x+H), H = floor(c x / B).  P_j(h) = sum_{i<=j} (-1)^i N h^i / x^{i+1} is the Taylor
polynomial of f(h) = N/(x+h); the remainder R(h) = N (-h)^{j+1} / (x^{j+1} (x+h)) has the fixed
sign (-1)^{j+1} for h >= 0 and |R| <= rho = N H^{j+1} / x^{j+2} <= c^{j+1} B / x.
If x + h | N then {f(h)} = 0, hence
   j odd  (R >= 0):  {P(h)} in [1 - rho, 1)  or h = 0;
   j even (R <= 0):  {P(h)} in [0, rho].
These "candidates" are counted by a difference of two degree-j floor sums and located by
bisection; each candidate is tested by one division.  This replaces the precision-1/x
requirement (which gives only N^{1/(j+1)} blocks) by precision ~1/H plus candidate location.

Oracles:
  * j = 1: OPERATIONAL exact Euclid-type floor sum (`floor_sum_linear`), polylog time.
  * j >= 2: SIMULATED oracle: exact brute force with per-block prefix sums.  Its cost is NOT
    the model cost; only the number of oracle calls is reported (no polylog algorithm for
    exact degree-2 floor sums is known).
Everything is exact integer/rational arithmetic.  `ladder_factor(N, j)` sees only N.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import gmpy2


def floor_sum_linear(n: int, m: int, a: int, b: int) -> int:
    """sum_{i=0}^{n-1} floor((a i + b) / m), m > 0, any integers a, b (Euclid-like, O(log))."""
    if n <= 0:
        return 0
    ans = 0
    # normalise a, b into [0, m)
    if a < 0 or a >= m:
        qa, a = divmod(a, m)
        ans += qa * n * (n - 1) // 2
    if b < 0 or b >= m:
        qb, b = divmod(b, m)
        ans += qb * n
    while True:
        if a >= m:
            ans += (n - 1) * n // 2 * (a // m)
            a %= m
        if b >= m:
            ans += n * (b // m)
            b %= m
        y_max = a * n + b
        if y_max < m:
            break
        n, b = divmod(y_max, m)
        m, a = a, m
    return ans


def iroot_floor(m: int, k: int) -> int:
    return int(gmpy2.iroot(gmpy2.mpz(m), k)[0])


class _Block:
    """Degree-j data for one block; P(h) = num(h) / den with integer polynomial num."""

    def __init__(self, N: int, x: int, H: int, j: int):
        self.N, self.x, self.H, self.j = N, x, H, j
        self.den = x ** (j + 1)
        # num(h) = sum_i (-1)^i N h^i x^{j-i}
        self.coef = [(-1) ** i * N * x ** (j - i) for i in range(j + 1)]
        # rho = N H^{j+1} / x^{j+2} = (N H^{j+1}) / (x * den)  -> as a shift on the numerator scale
        self.rho = Fraction(N * H ** (j + 1), x ** (j + 2))
        self._prefix = None

    def num(self, h: int) -> int:
        s = 0
        for c in reversed(self.coef):
            s = s * h + c
        return s


def _shifted_floor(num: int, den: int, shift: Fraction) -> int:
    """floor(num/den + shift) exactly (integer arithmetic)."""
    u, v = shift.numerator, shift.denominator
    return (num * v + u * den) // (den * v)


class Oracle:
    """Counts calls; j = 1 exact Euclid, j >= 2 simulated by brute force with caching."""

    def __init__(self):
        self.calls = 0
        self.sim_terms = 0

    def floor_range(self, blk: _Block, shift: Fraction, h1: int, h2: int) -> int:
        """sum_{h1 <= h < h2} floor(P(h) + shift)."""
        self.calls += 1
        if h2 <= h1:
            return 0
        if blk.j == 1:
            # P(h) + s = (c0 + c1 h)/den + s;  with s = u/v:  ((c0 v + u den) + c1 v h) / (den v)
            u, v = shift.numerator, shift.denominator
            c0, c1 = blk.coef
            a = c1 * v
            b = c0 * v + u * blk.den
            m = blk.den * v
            # shift index to start at h1
            b2 = b + a * h1
            return floor_sum_linear(h2 - h1, m, a, b2)
        # simulated: cache per (block, shift) prefix sums
        key = shift
        if blk._prefix is None:
            blk._prefix = {}
        if key not in blk._prefix:
            pre = [0]
            for h in range(blk.H):
                self.sim_terms += 1
                pre.append(pre[-1] + _shifted_floor(blk.num(h), blk.den, shift))
            blk._prefix[key] = pre
        pre = blk._prefix[key]
        return pre[h2] - pre[h1]


def _candidate_count(orc: Oracle, blk: _Block, h1: int, h2: int) -> int:
    if blk.j % 2 == 1:
        # {P} >= 1 - rho  <=>  floor(P + rho) > floor(P)
        return orc.floor_range(blk, blk.rho, h1, h2) - orc.floor_range(blk, Fraction(0), h1, h2)
    # {P} <= rho  <=>  floor(P) > floor(P - rho) or {P} = rho exactly (measure zero; handled by
    # using the closed shift rho' = rho + tiny would change semantics) -> use P - rho - eps:
    eps = Fraction(1, 4 * blk.den * blk.rho.denominator)
    return orc.floor_range(blk, Fraction(0), h1, h2) - orc.floor_range(blk, -blk.rho - eps, h1, h2)


def _locate(orc: Oracle, blk: _Block, h1: int, h2: int, cnt: int, out: List[int]):
    if cnt <= 0:
        return
    if h2 - h1 == 1:
        out.append(h1)
        return
    mid = (h1 + h2) // 2
    cl = _candidate_count(orc, blk, h1, mid)
    _locate(orc, blk, h1, mid, cl, out)
    _locate(orc, blk, mid, h2, cnt - cl, out)


def ladder_factor(N: int, j: int, c: Fraction = Fraction(1, 2),
                  stats: Optional[Dict] = None) -> Optional[int]:
    """Smallest-found nontrivial factor of N via the degree-j ladder, or None (N prime)."""
    if stats is None:
        stats = {}
    for k in ("blocks", "candidates", "trial_divisions", "divisor_tests"):
        stats.setdefault(k, 0)
    orc = Oracle()
    sN = isqrt(N)
    B = iroot_floor(N, j + 2)
    # first x with H = floor(c x / B) >= 1
    x0 = min(sN + 1, -((-B * c.denominator) // c.numerator))
    for d in range(2, min(x0, sN + 1)):
        stats["trial_divisions"] += 1
        if N % d == 0:
            stats["found_in"] = "trial"
            stats["oracle_calls"] = orc.calls
            return d
    x = x0
    while x <= sN:
        H = max(1, (c.numerator * x) // (c.denominator * B))
        H = min(H, sN - x + 1)
        blk = _Block(N, x, H, j)
        stats["blocks"] += 1
        cands: List[int] = []
        tot = _candidate_count(orc, blk, 0, H)
        _locate(orc, blk, 0, H, tot, cands)
        stats["located"] = stats.get("located", 0) + len(cands)
        if j % 2 == 1 and 0 not in cands:
            cands.append(0)  # R(0) = 0: P(0) = N/x exactly
        stats["candidates"] += len(cands)
        for h in sorted(cands):
            stats["divisor_tests"] += 1
            d = x + h
            if 1 < d < N and N % d == 0:
                stats["found_in"] = "block"
                stats["oracle_calls"] = orc.calls
                stats["sim_terms"] = orc.sim_terms
                return d
        x += H
    stats["oracle_calls"] = orc.calls
    stats["sim_terms"] = orc.sim_terms
    return None


def predicted_blocks(N: int, j: int, c: float = 0.5) -> float:
    """Continuum count (B/c) ln(sqrt N / x0), x0 = B/c  (fits/plots only)."""
    import math
    B = N ** (1.0 / (j + 2))
    x0 = B / c
    return (B / c) * math.log(math.sqrt(N) / x0) if math.sqrt(N) > x0 else 0.0
