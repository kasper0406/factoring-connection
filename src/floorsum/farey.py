"""Exact floor sums of convex/concave rational polynomials by Farey (hull-type) dissection.

OPERATIONAL, exact integer arithmetic.  P(h) = num(h) / den with an integer polynomial num of any
degree j and den > 0.  Requirement: P'' has constant sign on the integer range (the ladder's
P_j on a block is convex for every j >= 2, since h < H << x).

Key exactness lemma (proved in notes/floorsum_results.md, Lemma R):
    if q >= 1, p, r are integers and  r <= q P(u) - p u < r + 1,  then
        floor(P(u)) = floor((p u + r) / q).
Proof: write q P(u) = n + theta with n = p u + r and 0 <= theta < 1; n = q*floor(n/q) + s with
0 <= s <= q-1, so s + theta < q and floor((n + theta)/q) = floor(n/q).

Algorithm (one "window" per Farey fraction):
  * Q = ceil(beta^{-1/3}), beta = max|P''|/2 on the range (curvature scale).
  * At the current position h, take p/q = the last continued-fraction convergent with q <= Q of
    the local slope s = P(h+1) - P(h); window length L = max(1, floor(1/(q Q beta))).
  * On the window, g(u) = q P(u) - p u is convex/concave; split at its (integer) extremum; on each
    monotone side the integer function G(u) = floor(g(u)) is monotone, so its constant runs are
    found by binary search.  Each run is summed exactly by the Euclid floor sum.
Cost: O((#windows + #runs) log H) exact evaluations; heuristically #windows + #runs ~ H beta^{1/3}
(the integer-hull / affine-arc-length count) -- this is exactly the convexity trap.
"""
from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from math import floor
from typing import Dict, List, Optional, Sequence, Tuple

import gmpy2

from src.analytic.ladder import floor_sum_linear


def poly_eval(coef: Sequence[int], h: int) -> int:
    """coef[i] is the coefficient of h^i."""
    s = 0
    for c in reversed(coef):
        s = s * h + c
    return s


def _convergent_le(num: int, den: int, Q: int) -> Tuple[int, int]:
    """Last continued-fraction convergent p/q of num/den (den > 0) with q <= Q."""
    p0, q0, p1, q1 = 0, 1, 1, 0
    a, b = num, den
    while b != 0:
        t = a // b
        p2, q2 = t * p1 + p0, t * q1 + q0
        if q2 > Q:
            break
        p0, q0, p1, q1 = p1, q1, p2, q2
        a, b = b, a - t * b
    if q1 == 0:  # cannot happen for Q >= 1 (first convergent has q = 1)
        return num // den, 1
    return p1, q1


class Counters(dict):
    def bump(self, k, v=1):
        self[k] = self.get(k, 0) + v


class FloorSumConvex:
    """Run decomposition of h -> floor(num(h)/den) on [h_lo, h_hi); supports range queries.

    `sign` = +1 if P is convex on the range, -1 if concave (checked lazily by the caller).
    """

    def __init__(self, coef: Sequence[int], den: int, h_lo: int, h_hi: int,
                 counters: Optional[Counters] = None):
        assert den > 0 and h_hi >= h_lo
        self.coef = list(coef)
        while len(self.coef) > 1 and self.coef[-1] == 0:
            self.coef.pop()
        self.den = den
        self.lo, self.hi = h_lo, h_hi
        self.cnt = counters if counters is not None else Counters()
        # runs: (start, end, p, q, r) meaning floor(P(u)) = floor((p u + r)/q) on [start, end)
        self.runs: List[Tuple[int, int, int, int, int]] = []
        self._build()
        self.starts = [r[0] for r in self.runs]
        self.prefix = [0]
        for (a, b, p, q, r) in self.runs:
            self.prefix.append(self.prefix[-1] + self._run_sum(p, q, r, a, b))

    # -- exact helpers
    def _num(self, h: int) -> int:
        self.cnt.bump("evals")
        return poly_eval(self.coef, h)

    def _G(self, h: int, p: int, q: int) -> int:
        """floor(q P(h) - p h)."""
        return (q * self._num(h) - p * h * self.den) // self.den

    def _run_sum(self, p, q, r, a, b) -> int:
        """sum_{a <= u < b} floor((p u + r)/q)."""
        if b <= a:
            return 0
        self.cnt.bump("euclid_calls")
        return floor_sum_linear(b - a, q, p, p * a + r)

    def _beta(self) -> Fraction:
        """max |P''|/2 over the range, using second differences at the ends (P'' monotone on
        short ranges is not assumed: we take the max of |second differences| at both ends and
        the midpoint, which is an upper estimate used only for the window-length heuristic)."""
        lo, hi = self.lo, self.hi
        pts = {lo, max(lo, hi - 3), (lo + hi) // 2}
        m = Fraction(0)
        for t in pts:
            d2 = self._num(t + 2) - 2 * self._num(t + 1) + self._num(t)
            m = max(m, Fraction(abs(d2), 2 * self.den))
        return m

    def _check_convexity(self):
        """Necessary check of the precondition (constant sign of P''): the second differences at
        the two ends of the range must not have opposite signs.  Sufficient for degree <= 3 (P''
        linear); for degree >= 4 the caller must guarantee it (true for ladder blocks)."""
        lo, hi = self.lo, self.hi
        if hi - lo < 3:
            return
        d = [self._num(t + 2) - 2 * self._num(t + 1) + self._num(t) for t in (lo, hi - 3)]
        if (d[0] > 0 and d[1] < 0) or (d[0] < 0 and d[1] > 0):
            raise ValueError("FloorSumConvex: P'' changes sign on the range (split at the inflection)")

    # -- decomposition
    def _build(self):
        lo, hi = self.lo, self.hi
        if hi <= lo:
            return
        if len(self.coef) <= 2:  # linear: one run, exact Euclid
            c0 = self.coef[0]
            c1 = self.coef[1] if len(self.coef) > 1 else 0
            # floor((c1 h + c0)/den): p = c1, q = den, r = c0 with run form floor((p u + r)/q)
            self.runs.append((lo, hi, c1, self.den, c0))
            self.cnt.bump("windows")
            return
        self._check_convexity()
        beta = self._beta()
        if beta == 0:
            Q = hi - lo + 1
        else:
            # Q = ceil(beta^{-1/3})
            inv = 1 / beta
            Q = int(gmpy2.iroot(gmpy2.mpz(inv.numerator // inv.denominator + 1), 3)[0]) + 1
        Q = max(1, min(Q, hi - lo + 1))
        self.Q = Q
        h = lo
        while h < hi:
            self.cnt.bump("windows")
            s_num = self._num(h + 1) - self._num(h)  # slope numerator over den
            p, q = _convergent_le(s_num, self.den, Q)
            if beta == 0:
                L = hi - h
            else:
                Lf = 1 / (q * Q * beta)
                L = max(1, Lf.numerator // Lf.denominator)
            e = min(hi, h + L)
            self._window(h, e, p, q)
            h = e

    def _window(self, a: int, b: int, p: int, q: int):
        """Append runs for [a, b) using the line slope p/q (g(u) = q P(u) - p u convex/concave)."""
        G = lambda u: self._G(u, p, q)
        # integer extremum of g on [a, b): first u where the sign of g(u+1)-g(u) flips.
        # g convex: forward differences non-decreasing; concave: non-increasing.
        def dg(u):
            return (q * (self._num(u + 1) - self._num(u)) - p * self.den)
        if b - a <= 2:
            split = [a, b]
        else:
            d_first, d_last = dg(a), dg(b - 2)
            split = [a, b]
            if (d_first < 0 < d_last) or (d_first > 0 > d_last):
                # binary search for the first u in [a, b-2] with sign change
                neg_first = d_first < 0
                lo_, hi_ = a, b - 2  # invariant: dg(lo_) has the first sign; dg(hi_) not
                while hi_ - lo_ > 1:
                    mid = (lo_ + hi_) // 2
                    self.cnt.bump("bsearch")
                    if (dg(mid) < 0) == neg_first and dg(mid) != 0:
                        lo_ = mid
                    else:
                        hi_ = mid
                split = [a, hi_, b] if a < hi_ < b else [a, b]
        for i in range(len(split) - 1):
            self._monotone_runs(split[i], split[i + 1], p, q, G)

    def _monotone_runs(self, a: int, b: int, p: int, q: int, G):
        """[a, b): G is monotone there (g is monotone up to one integer step at the split)."""
        u = a
        while u < b:
            r = G(u)
            # largest v in [u, b-1] with G(v) == r  (G monotone => {G == r} is an interval)
            if G(b - 1) == r:
                v = b - 1
            else:
                lo_, hi_ = u, b - 1  # G(lo_) == r, G(hi_) != r
                while hi_ - lo_ > 1:
                    mid = (lo_ + hi_) // 2
                    self.cnt.bump("bsearch")
                    if G(mid) == r:
                        lo_ = mid
                    else:
                        hi_ = mid
                v = lo_
            self.runs.append((u, v + 1, p, q, r))
            self.cnt.bump("runs")
            u = v + 1

    # -- queries
    def prefix_sum(self, t: int) -> int:
        """sum_{lo <= h < t} floor(P(h))."""
        if t <= self.lo:
            return 0
        if t >= self.hi:
            return self.prefix[-1]
        i = bisect_right(self.starts, t) - 1
        a, b, p, q, r = self.runs[i]
        return self.prefix[i] + self._run_sum(p, q, r, a, t)

    def range_sum(self, h1: int, h2: int) -> int:
        self.cnt.bump("queries")
        if h2 <= h1:
            return 0
        return self.prefix_sum(h2) - self.prefix_sum(h1)

    def total(self) -> int:
        return self.prefix[-1]


def floor_sum_poly(coef: Sequence[int], den: int, h1: int, h2: int,
                   counters: Optional[Counters] = None) -> int:
    """sum_{h1 <= h < h2} floor(poly(h)/den) for P convex or concave on the range (exact)."""
    return FloorSumConvex(coef, den, h1, h2, counters).total()


def floor_sum_brute(coef: Sequence[int], den: int, h1: int, h2: int) -> int:
    return sum(poly_eval(coef, h) // den for h in range(h1, h2))
