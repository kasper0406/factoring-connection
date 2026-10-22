"""HARNESS / DIAGNOSTIC: exact censuses E1 and E2 of notes/selfred_prereg.md (they use the factors).

E1: Bayes-optimal advantage of any function of N mod 2^j for predicting eps = (-1/X), ordered iid pairs from
    S_b = primes in [3*2^(b-2), 2^b) (diagonal included), against the proved bound and a random-label null.
E2: G_p = gcd_{q in (p,2p) prime} h(pq) for oblivious polynomials h of arbitrary height (Theorem SR1).
"""
from __future__ import annotations

import math
import random
from collections import Counter
from fractions import Fraction
from math import comb, gcd
from typing import Dict, List

import numpy as np
import sympy

from .rigidity import content, peval, pmul


# ============================================================================= E1
def prime_set(b: int) -> np.ndarray:
    return np.array(list(sympy.primerange(3 << (b - 2), 1 << b)), dtype=np.int64)


def _tv1_exact(cnt_hist: Counter, s: int, U: int, nnz: int) -> Fraction:
    """TV(law of X mod 2^j, uniform on the U units) with class counts given by histogram {count: #classes}."""
    tot = sum(h * abs(k * U - s) for k, h in cnt_hist.items()) + (U - nnz) * s
    return Fraction(tot, 2 * s * U)


def _tv2_exact(cnt_hist: Counter, s: int, U: int) -> Fraction:
    """TV(mu x mu, uniform on units^2) = sum over cells of max(0, mu(u)mu(v) - 1/U^2), exact."""
    ks = sorted(cnt_hist)
    tot = 0
    U2, s2 = U * U, s * s
    for k in ks:
        for l in ks:
            v = k * l * U2 - s2
            if v > 0:
                tot += cnt_hist[k] * cnt_hist[l] * v
    return Fraction(tot, s2 * U2)


_null_cache: Dict[int, float] = {}


def null_excess(k: int) -> float:
    """E[max(B, k-B)] - k/2 for B ~ Bin(k, 1/2) = m C(k, m) / 2^k with m = floor(k/2) + 1 (de Moivre's mean
    absolute deviation).  Exact (Fraction) for k <= 4000, log-gamma double precision above (deviation D1)."""
    if k in _null_cache:
        return _null_cache[k]
    if k == 0:
        v = 0.0
    else:
        m = k // 2 + 1
        if k <= 4000:
            v = float(Fraction(m * comb(k, m), 1 << k))
        else:
            lg = math.lgamma(k + 1) - math.lgamma(m + 1) - math.lgamma(k - m + 1)
            v = math.exp(math.log(m) + lg - k * math.log(2))
    _null_cache[k] = v
    return v


def e1_census(b: int, jmax: int = None) -> List[dict]:
    S = prime_set(b)
    s = len(S)
    n = 2 * b
    jmax = jmax or n
    eps = np.where(S % 4 == 1, 1, -1).astype(np.int8)
    N = np.multiply.outer(S, S).ravel()   # ordered pairs (X, Y), X-major
    E = (N % 4) == 1
    NE = N[E]
    epsE = np.repeat(eps, s)[E]
    nE = int(E.sum())
    rows = []
    for j in range(2, jmax + 1):
        U = 1 << (j - 1)
        cls = NE & ((1 << j) - 1)
        order = np.argsort(cls, kind="stable")
        c_sorted = cls[order]
        e_sorted = epsE[order]
        _, starts, counts = np.unique(c_sorted, return_index=True, return_counts=True)
        plus = np.add.reduceat((e_sorted > 0).astype(np.int64), starts)
        minus = counts - plus
        best = int(np.maximum(plus, minus).sum())
        Astar = Fraction(best, nE) - Fraction(1, 2)
        # null
        kc = Counter(counts.tolist())
        null = sum(h * null_excess(k) for k, h in kc.items()) / nE
        # TV of the marginal law of X mod 2^j
        xm = S & ((1 << j) - 1)
        xc = Counter(Counter(xm.tolist()).values())
        nnz = sum(xc.values())
        tv1 = _tv1_exact(xc, s, U, nnz)
        tv2 = _tv2_exact(xc, s, U)
        bound = (3 * tv2 / (1 - 2 * tv2)) if tv2 < Fraction(1, 2) else None
        rows.append(dict(b=b, n=n, j=j, s=s, pairs=s * s, nE=nE, classes=int(len(counts)),
                         Astar=float(Astar), Astar_num=Astar.numerator, Astar_den=Astar.denominator,
                         A_null=null, ratio=(float(Astar) / null if null > 0 else None),
                         TV1=float(tv1), TV2=float(tv2),
                         bound=(float(bound) if bound is not None else None),
                         P1_ok=(bound is None or Astar <= bound)))
    return rows


# ============================================================================= E2
def vl(x: int, l: int) -> int:
    if x == 0:
        return 10 ** 9
    v = 0
    while x % l == 0:
        x //= l
        v += 1
    return v


def unit_fixed_valuation(h: List[int], l: int, Emax: int = 12) -> int:
    """u_l(h) = min over z with l not dividing z of v_l(h(z)) (computed mod l^E, E increasing)."""
    for E in range(1, Emax + 1):
        mod = l ** E
        m = E
        for z in range(1, mod):
            if z % l:
                v = min(vl(peval(h, z) % mod, l), E)
                if v < m:
                    m = v
                    if m == 0:
                        return 0
        if m < E:
            return m
        if mod > 2 * 10 ** 6:
            break
    raise RuntimeError(f"u_l not resolved for l={l}")


def unit_fixed_divisor(h: List[int]) -> int:
    c = abs(content(h))
    D = len(h) - 1
    cand = set(int(l) for l in sympy.primerange(2, D + 2)) | set(sympy.factorint(c).keys())
    U = 1
    for l in cand:
        if c % l == 0 and l - 1 > D:
            u = vl(c, l)
        else:
            u = unit_fixed_valuation(h, l)
        U *= l ** u
    return U


def e2_polys(D: int, k: int) -> List[int]:
    rng = random.Random(f"selfred-e2|{D}|{k}")
    h = [rng.randint(-(1 << 64), 1 << 64) for _ in range(D + 1)]
    while h[-1] == 0:
        h[-1] = rng.randint(1, 1 << 64)
    c = 1
    for _ in range(rng.randint(1, 3)):
        c *= int(sympy.prevprime(rng.randint(3, 1 << 20)))
    h = [c * x for x in h]
    if k % 2 == 1 and D >= 2:
        planted = [0, 1]  # x
        for i in range(1, D + 1):
            planted = pmul(planted, [-i, 1])
        h = pmul(h, planted)
    return h


def e2_case(h: List[int], p: int) -> dict:
    Q = list(sympy.primerange(p + 1, 2 * p))
    G = 0
    for q in Q:
        v = peval(h, p * q)
        if v:
            G = gcd(G, v)
    G = abs(G)
    vp = vl(G, p)
    Gr = G // p ** vp
    c = abs(content(h))
    D = len(h) - 1
    L = 1 + Fraction(4 * D * p) / Fraction(len(Q)) / Fraction(math.log(p))  # display bound (float log)
    big = [l for l in sympy.factorint(Gr) if c % l]
    maxl = max(big) if big else 1
    hyp = len(Q) > D * (math.isqrt(p) + 2)
    return dict(p=p, D=D, Q=len(Q), vp=vp, Gr=Gr, c_bits=c.bit_length(), max_l_not_c=maxl,
                L=float(L), hyp=hyp, P1_ok=(not hyp) or maxl < float(L))


def e1_posthoc_dedup(b: int) -> List[dict]:
    """POST-HOC (not pre-registered; deviation D2): random-label null with one label per UNORDERED pair.

    The pre-registered null gave each ordered pair its own coin, but (X, Y) and (Y, X) give the same N and the
    same eps, so that null is too small by up to a factor 2 at large j.  Here A* is recomputed on unordered pairs
    (diagonal weight 1, off-diagonal weight 2; identical to the ordered A*), and the null gives one fair coin
    per unordered pair (weight approximated by the mean weight of the class population).
    """
    S = prime_set(b)
    s = len(S)
    eps = np.where(S % 4 == 1, 1, -1)
    iu = np.triu_indices(s)
    N = S[iu[0]] * S[iu[1]]
    w = np.where(iu[0] == iu[1], 1, 2)
    E = (N % 4) == 1
    N, w, e = N[E], w[E], eps[iu[0]][E]
    nE = int(w.sum())
    rows = []
    for j in range(2, 2 * b + 1):
        cls = N & ((1 << j) - 1)
        o = np.argsort(cls, kind="stable")
        _, st, cnt = np.unique(cls[o], return_index=True, return_counts=True)
        ww, ee = w[o], e[o]
        plus = np.add.reduceat(ww * (ee > 0), st)
        tot = np.add.reduceat(ww, st)
        Astar = float(np.maximum(plus, tot - plus).sum()) / nE - 0.5
        null = sum(h * null_excess(k) for k, h in Counter(cnt.tolist()).items()) * (tot.sum() / cnt.sum()) / nE
        rows.append(dict(b=b, j=j, Astar=Astar, A_null_dedup=float(null), ratio=Astar / float(null)))
    return rows
