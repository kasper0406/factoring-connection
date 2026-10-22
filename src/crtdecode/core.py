"""Exact helpers for crtdecode.

Everything here is exact integer arithmetic (numpy uint64 for windows < 2^63, Python ints otherwise).
Nothing here sees a factor of any modulus N: the inputs are a window start X, its length h and the bound y.

- primes_upto, prime_powers_upto: factor-base data (primes_upto copied from sihard.core; owner rule 7).
- window_profile: exact valuation profile of [X, X + h] (smooth mask, incidence counts, cofactors).
- smooth_in_range: exact list of y-smooth integers in [lo, hi] (used for oblivious repertoires).
- is_smooth: exact trial division test of one integer.
- certified: the idealised unweighted Johnson certificate (Theorem J's class).
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

import numpy as np


@lru_cache(maxsize=None)
def primes_upto(n: int) -> Tuple[int, ...]:
    if n < 2:
        return ()
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return tuple(i for i in range(n + 1) if s[i])


def prime_powers_upto(y: int, bound: int) -> List[Tuple[int, int, int]]:
    """All (q, e, q^e) with q <= y prime, e >= 1, q^e <= bound."""
    out = []
    for q in primes_upto(y):
        m, e = q, 1
        while m <= bound:
            out.append((q, e, m))
            m *= q
            e += 1
    return out


def is_smooth(z: int, y: int) -> bool:
    """Exact: z >= 1 is y-smooth."""
    if z < 1:
        raise ValueError("z must be positive")
    for p in primes_upto(y):
        if p * p > z:
            break
        while z % p == 0:
            z //= p
    return z == 1 or z <= y


def window_profile(X: int, h: int, y: int) -> Dict[str, object]:
    """Exact valuation profile of the window [X, X + h].

    Returns: 'smooth' (bool array, length h + 1), 'cof' (cofactor after removing all y-smooth prime powers),
    'incidences' (total number of pairs (q^e, j) with q^e | X + j), 'first' (least smooth offset or -1).
    """
    if X < 1 or X + h >= (1 << 63):
        raise ValueError("window outside uint64-safe range")
    vals = np.arange(X, X + h + 1, dtype=np.uint64)
    cof = vals.copy()
    inc = 0
    top = X + h
    for (q, e, m) in prime_powers_upto(y, top):
        j0 = (-X) % m
        if j0 > h:
            continue
        idx = np.arange(j0, h + 1, m)
        inc += len(idx)
        cof[idx] //= np.uint64(q)
    smooth = cof == 1
    nz = np.flatnonzero(smooth)
    return {"smooth": smooth, "cof": cof, "incidences": int(inc), "first": int(nz[0]) if nz.size else -1}


def prefix_incidences(X: int, j_last: int, y: int) -> int:
    """Exact number of pairs (q^e, j), 0 <= j <= j_last, q^e | X + j (sieve cost of the prefix)."""
    tot = 0
    for (_q, _e, m) in prime_powers_upto(y, X + j_last):
        tot += (X + j_last) // m - (X - 1) // m
    return tot


def smooth_in_range(lo: int, hi: int, y: int) -> List[int]:
    """Exact list of y-smooth integers in [lo, hi] (hi < 2^63)."""
    n = hi - lo + 1
    cof = np.arange(lo, hi + 1, dtype=np.uint64)
    for (q, _e, m) in prime_powers_upto(y, hi):
        j0 = (-lo) % m
        if j0 < n:
            idx = np.arange(j0, n, m)
            cof[idx] //= np.uint64(q)
    return [lo + int(i) for i in np.flatnonzero(cof == 1)]


def certified(g: int, M: int, region_len: int) -> bool:
    """Idealised unweighted Johnson certificate (ln g)^2 >= ln|R| * ln M, decided EXACTLY.

    Supported exact cases (the only ones crtdecode uses):
    - region_len == M: certified iff g >= M (i.e. g = M, since g | M);
    - g == M: certified iff M >= region_len;
    - M = q^e and g = q^b (same prime q): (b ln q)^2 >= ln R * e ln q  <=>  q^(b*b) >= R^e.
    Anything else raises (no float comparison of logarithms).
    """
    if g <= 1:
        return False
    if region_len == M:
        return g >= M
    if g == M:
        return M >= region_len
    q = _prime_base(M)
    if q is not None and _prime_base(g) == q:
        e = _exp_of(M, q)
        b = _exp_of(g, q)
        return q ** (b * b) >= region_len ** e
    raise NotImplementedError("certified(): no exact rule for this (g, M, |R|)")


def _prime_base(m: int):
    if m < 2:
        return None
    for p in primes_upto(1 << 16):
        if m % p == 0:
            while m % p == 0:
                m //= p
            return p if m == 1 else None
    return None


def _exp_of(m: int, q: int) -> int:
    e = 0
    while m % q == 0:
        m //= q
        e += 1
    return e
