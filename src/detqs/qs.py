"""OPERATIONAL. Deterministic single-polynomial quadratic sieve (no randomness).

qs_run(N, y, T, dmax) sees only N and the public parameters.  Steps:
1. If N is a perfect square or a factor-base prime divides N, report it
   (classified ACCIDENTAL by the harness: trial division).
2. Exact sieve of Q_t = (m+t)^2 - N, t = 1..T (sieve.smooth_mask).
3. Relations in increasing t; each is added to an incremental F_2 eliminator.
   Every time the kernel dimension grows, the new kernel vector S gives
   X_S = prod(m+t_i) mod N, Y_S = prod l^{e_l(S)/2} mod N and g = gcd(X_S - Y_S, N).
   The run stops at the first proper g, or when dmax kernel vectors were tested.
Every nonunit met is output (rule 3): gcd(m+t, N) and gcd(Y_S, N) are checked.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import gmpy2

from .linalg import IncrementalF2
from .sieve import factor_base, smooth_relations


@dataclass
class QSResult:
    N: int
    y: int
    T: int
    fb_size: int
    Z: int                         # number of y-smooth Q_t, t <= T
    factor: Optional[int] = None
    how: str = "none"              # 'square' | 'trial' | 'unit' | 'dependency' | 'none'
    relations_used: int = 0
    deps_tested: int = 0
    dep_outcomes: List[int] = field(default_factory=list)   # 1 fruitful, 0 trivial
    dep_jacobi: List[int] = field(default_factory=list)     # Jacobi (X_S Y_S / N)
    dep_sizes: List[int] = field(default_factory=list)      # |S|
    first_fruitful: Optional[int] = None                    # 1-based kernel index


def qs_run(N: int, y: int, T: int, dmax: int = 30) -> QSResult:
    if N < 4:
        raise ValueError
    r = int(gmpy2.isqrt(N))
    if r * r == N:
        res = QSResult(N, y, T, 0, 0, factor=r, how="square")
        return res
    fb = factor_base(N, y)
    if fb.divides_N:
        l = fb.divides_N[0]
        return QSResult(N, y, T, len(fb.primes), 0, factor=l, how="trial")
    m, fb, ts, exps, Q = smooth_relations(N, y, T, fb)
    res = QSResult(N, y, T, len(fb.primes), len(ts))
    col = {l: j for j, l in enumerate(fb.primes)}
    elim = IncrementalF2()
    rel_x: List[int] = []
    rel_e: List[dict] = []
    for t, e in zip(ts, exps):
        x = m + t
        g = gmpy2.gcd(x, N)
        if 1 < g < N:
            res.factor, res.how = int(g), "unit"
            return res
        chk = 1
        for l, k in e.items():
            chk *= l ** k
        assert chk == int(Q[t - 1]), "exponent vector does not reproduce Q_t"
        vec = 0
        for l, k in e.items():
            if k & 1:
                vec |= 1 << col[l]
        rel_x.append(x)
        rel_e.append(e)
        combo = elim.add(vec)
        if combo is None:
            continue
        S = [i for i in range(len(rel_x)) if (combo >> i) & 1]
        X = 1
        tot: dict = {}
        for i in S:
            X = X * rel_x[i] % N
            for l, k in rel_e[i].items():
                tot[l] = tot.get(l, 0) + k
        Y = 1
        for l, k in tot.items():
            assert k % 2 == 0
            Y = Y * pow(l, k // 2, N) % N
        gy = gmpy2.gcd(Y, N)
        if 1 < gy < N:
            res.factor, res.how = int(gy), "unit"
            return res
        g = int(gmpy2.gcd(X - Y, N))
        res.deps_tested += 1
        res.dep_sizes.append(len(S))
        res.dep_jacobi.append(int(gmpy2.jacobi(X * Y % N, N)))
        fr = 1 if 1 < g < N else 0
        res.dep_outcomes.append(fr)
        if fr:
            res.factor, res.how = g, "dependency"
            res.first_fruitful = res.deps_tested
            res.relations_used = len(rel_x)
            return res
        if res.deps_tested >= dmax:
            break
    res.relations_used = len(rel_x)
    return res
