"""HARNESS.  Coverage heights at the Theorem-W scaling Y = floor(q^{8/25}).

h(a) = min_{l in band} mu_l(a * qbar mod l), mu_l(t) = least height of a signed primitive weight
(A, B), B >= 1, A != 0, with A = t B (mod l); mu_l(0) = infinity (censored).  For a prime p in
(q/2, q), h(p) <= H iff the pair (p, q) is covered by (band, signed weights of height <= H).
"""
from __future__ import annotations

from math import isqrt, exp, log
from typing import Dict, List, Tuple

import numpy as np

from src.dchworst import coverage as cv
from src.dcltheorem.params import primes_upto

INF = 1 << 30


def y_of(q: int, num: int = 8, den: int = 25) -> int:
    """Y = floor(q^{num/den}) exactly (integer root)."""
    # largest Y with Y^den <= q^num
    lo, hi = 1, 1
    while hi ** den <= q ** num:
        hi *= 2
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** den <= q ** num:
            lo = mid
        else:
            hi = mid - 1
    return lo


def band(Y: int) -> List[int]:
    return cv.band(Y)


def thue_cap(l: int) -> int:
    return isqrt(l)


def mu_table(l: int) -> np.ndarray:
    """mu_l(t) for all t mod l; t = 0 -> INF.  Thue: mu_l(t) <= isqrt(l) for t != 0."""
    Hcap = isqrt(l) + 1
    tab = cv.hmin_table(l, Hcap, True)
    tab = tab.astype(np.int64)
    tab[1] = 1          # the diagonal weight (1, 1), excluded from dchworst's signed set
    tab[0] = INF
    assert tab[1:].max() <= isqrt(l), "Thue violated"
    return tab


def nu(tab: np.ndarray, Hmax: int) -> np.ndarray:
    """nu[H] = #{t != 0 : mu(t) <= H}, H = 0..Hmax."""
    v = tab[1:]
    c = np.bincount(np.minimum(v, Hmax + 1), minlength=Hmax + 2)[: Hmax + 1]
    return np.cumsum(c)


def survival(tabs: Dict[int, np.ndarray], Hmax: int) -> np.ndarray:
    """s[H] = prod_l (1 - nu_l(H)/(l-1)): independent-residue Pr[h > H] for a prime p."""
    s = np.ones(Hmax + 1)
    for l, tab in tabs.items():
        s *= 1.0 - nu(tab, Hmax) / (l - 1)
    return s


def h_quantile(s: np.ndarray, n: int, alpha: float) -> int:
    """Least H with exp(-n s(H)) >= alpha, i.e. Pr[max over n draws <= H] >= alpha."""
    for H in range(len(s)):
        if exp(-n * s[H]) >= alpha:
            return H
    return len(s)


def heights(a: np.ndarray, q: int, tabs: Dict[int, np.ndarray]) -> np.ndarray:
    """h(a) for an int64 array a (exact integer arithmetic)."""
    h = np.full(len(a), INF, dtype=np.int64)
    for l, tab in tabs.items():
        qi = pow(q % l, -1, l)
        t = (a % l) * qi % l
        np.minimum(h, tab[t], out=h)
    return h


def primes_in(lo: int, hi: int, chunk: int = 1 << 23):
    """Yield int64 arrays of the primes in (lo, hi) by a segmented sieve."""
    base = np.array(primes_upto(isqrt(hi) + 1), dtype=np.int64)
    start = lo + 1
    while start < hi:
        end = min(hi, start + chunk)          # [start, end)
        seg = np.ones(end - start, dtype=bool)
        for pr in base:
            pr = int(pr)
            if pr * pr >= end:
                break
            f = max(pr * pr, ((start + pr - 1) // pr) * pr)
            seg[f - start::pr] = False
        if start <= 1:
            seg[: 2 - start] = False
        yield np.nonzero(seg)[0].astype(np.int64) + start
        start = end


def worst_case(q: int, all_residues: bool = False) -> dict:
    """Exhaustive h* over primes p in (q/2, q) (or all a in [1, q-1])."""
    Y = y_of(q)
    L = band(Y)
    tabs = {l: mu_table(l) for l in L}
    Hmax = max(isqrt(l) for l in L) + 1
    hist = np.zeros(Hmax + 2, dtype=np.int64)
    n = 0
    argmax = None
    hstar = -1
    if all_residues:
        gen = (np.arange(s, min(q, s + (1 << 23)), dtype=np.int64) for s in range(1, q, 1 << 23))
    else:
        gen = primes_in(q // 2, q)
    for a in gen:
        h = heights(a, q, tabs)
        assert h.max() < INF
        hist += np.bincount(h, minlength=Hmax + 2)[: Hmax + 2]
        n += len(a)
        i = int(np.argmax(h))
        if h[i] > hstar:
            hstar, argmax = int(h[i]), int(a[i])
    s = survival(tabs, Hmax)
    lmin = min(L)
    return {"q": q, "Y": Y, "J": len(L), "lmin": lmin, "thue": isqrt(lmin), "n": n,
            "hstar": hstar, "argmax": argmax, "hist": hist.tolist(),
            "H05": h_quantile(s, n, 0.05), "H50": h_quantile(s, n, 0.5), "H95": h_quantile(s, n, 0.95),
            "c_emp": (0.5 * log(Y) - log(hstar)) / log(len(L)) if len(L) > 1 else None,
            "all_residues": all_residues}
