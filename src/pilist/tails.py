"""Adaptive elimination with Theta(u) tail primes (agent pilist, Task 3; smoothenum section 3.7 open item).

Model H_tail (notes/pilist_results.md, Prop TAIL). Band columns are the primes of (B^{1-eps}, B]. Each
relation carries <= j *forced* band entries, chosen adaptively from the history (special-q on tail primes,
lattice index <= the box area U forces j = O(1)), and k *random* band entries; conditional on everything
generated before it, every random entry equals any fixed prime with probability <= p.

Prop TAIL (first moment, PROVED in the model): the probability that some nonempty subset of the first R
relations cancels on the band columns (mod 2) is at most
    sum_{s=1}^{R} C(R, s) sum_{w=0}^{j s} C(k s, w) * M(k s - w) * p^{(k s - w)/2},
M(n) = (n - 1)!! for even n (perfect matchings), 0 for odd n.  `first_moment_log` evaluates it exactly in
log space; `r_threshold` returns the largest R with bound <= 1/2.

Simulation (HARNESS + strategies): `simulate` runs uniform band columns (p = 1/n) with
- 'nonadaptive': j = 0;
- 'greedy': j forced entries chosen *before* the random part, among columns seen exactly once so far
  (the special-q cycle-closing strategy of smoothenum (f));
- 'oracle': ILLEGAL in the model (forced entries chosen *after* seeing the random part, to cancel its
  residual against the current basis) - shows what the model assumption excludes.
It returns R*, the number of relations at the first band dependency (exact F2 elimination, bitsets).
"""
from __future__ import annotations

import random
from math import comb, exp, lgamma, log
from typing import Dict, List


def _log_matchings(n: int) -> float:
    """log((n-1)!!) for even n >= 0; -inf for odd n."""
    if n % 2:
        return float("-inf")
    if n == 0:
        return 0.0
    # (n-1)!! = n! / (2^{n/2} (n/2)!)
    return lgamma(n + 1) - (n // 2) * log(2) - lgamma(n // 2 + 1)


def _logsumexp(xs: List[float]) -> float:
    xs = [x for x in xs if x != float("-inf")]
    if not xs:
        return float("-inf")
    m = max(xs)
    return m + log(sum(exp(x - m) for x in xs))


def _logcomb(n: int, r: int) -> float:
    return lgamma(n + 1) - lgamma(r + 1) - lgamma(n - r + 1)


_LF = [0.0]


def _lf(n: int):
    """numpy table of log(i!) for i <= n (grown on demand)."""
    import numpy as np
    global _LF
    if len(_LF) <= n:
        _LF = [0.0] * (n + 1)
        acc = 0.0
        for i in range(1, n + 1):
            acc += log(i)
            _LF[i] = acc
    return np.array(_LF[: n + 1])


def first_moment_log(R: int, k: int, j: int, p: float) -> float:
    """log of the Prop TAIL first-moment bound (exact up to float rounding of the log-sum)."""
    import numpy as np
    LF = _lf(k * R + 2)
    lp = log(p)
    out = []
    for s in range(1, R + 1):
        w = np.arange(0, j * s + 1)
        m = k * s - w
        ok = (m >= 0) & (m % 2 == 0)
        w, m = w[ok], m[ok]
        if len(w) == 0:
            continue
        lcomb = LF[k * s] - LF[w] - LF[k * s - w]
        lmatch = LF[m] - (m // 2) * log(2) - LF[m // 2]
        inner = lcomb + lmatch + (m / 2) * lp
        mx = inner.max()
        val = mx + log(np.exp(inner - mx).sum())
        out.append(LF[R] - LF[s] - LF[R - s] + val)
    arr = np.array(out)
    mx = arr.max()
    return float(mx + log(np.exp(arr - mx).sum()))


def r_threshold(k: int, j: int, p: float, target: float = 0.5, rmax: int = 10 ** 7) -> int:
    """Largest R with first-moment bound <= target (bisection; the bound is increasing in R)."""
    lo, hi = 1, 2
    lt = log(target)
    while hi < rmax and first_moment_log(hi, k, j, p) <= lt:
        lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if first_moment_log(mid, k, j, p) <= lt:
            lo = mid
        else:
            hi = mid
    return lo


class F2Basis:
    """Incremental F2 row space over bitset ints, pivot = lowest set bit."""

    def __init__(self):
        self.piv: Dict[int, int] = {}

    def reduce(self, v: int) -> int:
        while v:
            low = v & -v
            r = self.piv.get(low)
            if r is None:
                return v
            v ^= r
        return 0

    def add(self, v: int) -> bool:
        """Insert; returns True iff v was dependent (reduced to zero)."""
        v = self.reduce(v)
        if v == 0:
            return True
        self.piv[v & -v] = v
        return False


def simulate(n: int, k: int, j: int, strategy: str, seed: int, rmax: int = None) -> Dict:
    rng = random.Random(seed)
    rmax = rmax or 3 * n
    basis = F2Basis()
    counts = [0] * n
    once = set()          # columns seen exactly once (by total multiplicity parity of appearance count)
    for R in range(1, rmax + 1):
        forced: List[int] = []
        if strategy == "greedy" and j > 0 and once:
            pool = sorted(once)
            rng.shuffle(pool)
            forced = pool[:j]
        rand = [rng.randrange(n) for _ in range(k)]
        v = 0
        for c in rand:
            v ^= 1 << c
        if strategy == "oracle" and j > 0:
            res = basis.reduce(v)
            bits = []
            while res and len(bits) < j:
                low = res & -res
                bits.append(low.bit_length() - 1)
                res ^= low
            forced = bits
        for c in forced:
            v ^= 1 << c
        for c in rand + forced:
            counts[c] += 1
            if counts[c] == 1:
                once.add(c)
            elif c in once:
                once.discard(c)
        if v == 0 or basis.add(v):
            return {"n": n, "k": k, "j": j, "strategy": strategy, "seed": seed, "R_star": R}
    return {"n": n, "k": k, "j": j, "strategy": strategy, "seed": seed, "R_star": None}
