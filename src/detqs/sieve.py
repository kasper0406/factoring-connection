"""OPERATIONAL. Exact smoothness sieve for the deterministic single-polynomial QS.

Q_t(N) = (m + t)^2 - N, m = ceil(sqrt N), t = 1..T.  Factor base: all primes
l <= y such that x^2 = N (mod l) is solvable (l = 2, l | N, or (N/l) = 1).
Smoothness is decided EXACTLY: every prime power of every factor-base prime is
divided out (vectorised repeated division), and Q_t is y-smooth iff the
remainder is 1.  No logarithmic thresholds, no floating point.

Only N and the public parameters (y, T) are used.  A factor-base prime dividing
N is reported (it is a factor found by trial division: ACCIDENTAL class).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import gmpy2
import numpy as np


def primes_upto(m: int) -> List[int]:
    if m < 2:
        return []
    s = bytearray([1]) * (m + 1)
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(m) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, m + 1, i)))
    return [i for i in range(m + 1) if s[i]]


def sqrt_mod_prime(a: int, l: int) -> int:
    """A square root of a mod odd prime l (a a nonzero QR). Tonelli-Shanks with
    the deterministic choice of the least non-residue z (exhaustive search)."""
    a %= l
    if l % 4 == 3:
        return pow(a, (l + 1) // 4, l)
    q, s = l - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (l - 1) // 2, l) != l - 1:
        z += 1
    mm, c, t, r = s, pow(z, q, l), pow(a, q, l), pow(a, (q + 1) // 2, l)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % l
            i += 1
        b = pow(c, 1 << (mm - i - 1), l)
        mm, c, t, r = i, b * b % l, t * b * b % l, r * b % l
    return r


@dataclass
class FactorBase:
    primes: List[int]                 # sorted, includes 2
    roots: Dict[int, Tuple[int, ...]]  # l -> residues t0 (mod l) with l | Q_t
    divides_N: List[int]              # factor-base primes dividing N (accidental)


def factor_base(N: int, y: int) -> FactorBase:
    m = int(gmpy2.isqrt(N - 1)) + 1
    ps: List[int] = []
    roots: Dict[int, Tuple[int, ...]] = {}
    hit: List[int] = []
    for l in primes_upto(y):
        if l == 2:
            ps.append(2)
            roots[2] = ()  # handled separately (all t)
            continue
        Nl = N % l
        if Nl == 0:
            hit.append(l)
            ps.append(l)
            roots[l] = ((-m) % l,)
            continue
        if pow(Nl, (l - 1) // 2, l) != 1:
            continue
        s = sqrt_mod_prime(Nl, l)
        ps.append(l)
        roots[l] = tuple(sorted({(s - m) % l, (-s - m) % l}))
    return FactorBase(ps, roots, hit)


def q_values(N: int, T: int) -> Tuple[int, np.ndarray]:
    """m and the array Q[t-1] = (m+t)^2 - N for t = 1..T (uint64, exact)."""
    m = int(gmpy2.isqrt(N - 1)) + 1
    r = m * m - N
    top = r + 2 * m * T + T * T
    if top >= (1 << 63):
        raise ValueError("Q_T too large for exact uint64 arithmetic")
    t = np.arange(1, T + 1, dtype=np.uint64)
    Q = np.uint64(r) + np.uint64(2 * m) * t + t * t
    return m, Q


def _strip(R: np.ndarray, idx: np.ndarray, l: int) -> None:
    """Divide R[idx] by l as many times as possible (exact, in place)."""
    L = np.uint64(l)
    cur = idx
    while cur.size:
        sub = R[cur]
        mask = (sub % L) == 0
        cur = cur[mask]
        if cur.size == 0:
            break
        R[cur] = sub[mask] // L


def smooth_mask(N: int, y: int, T: int, fb: Optional[FactorBase] = None):
    """Boolean mask over t = 1..T: Q_t is y-smooth (exact). Returns (mask, m, Q, fb)."""
    if fb is None:
        fb = factor_base(N, y)
    m, Q = q_values(N, T)
    R = Q.copy()
    # l = 2: remove the full power of two exactly via the lowest set bit
    nz = R != 0
    low = R & (~R + np.uint64(1))
    R[nz] = R[nz] // low[nz]
    for l in fb.primes:
        if l == 2:
            continue
        for t0 in fb.roots[l]:
            first = t0 if t0 >= 1 else t0 + l   # t >= 1
            if first > T:
                continue
            idx = np.arange(first - 1, T, l, dtype=np.int64)
            _strip(R, idx, l)
    return (R == 1), m, Q, fb


def exponents(Qt: int, fb: FactorBase) -> Dict[int, int]:
    """Exact factorisation of a y-smooth Q_t over the factor base (trial division)."""
    e: Dict[int, int] = {}
    x = Qt
    for l in fb.primes:
        if x == 1:
            break
        if x % l == 0:
            k = 0
            while x % l == 0:
                x //= l
                k += 1
            e[l] = k
    if x != 1:
        raise ValueError("not smooth")
    return e


def smooth_relations(N: int, y: int, T: int, fb: Optional[FactorBase] = None):
    """All y-smooth Q_t, t <= T, with exact exponent dicts (vectorised).
    Returns (m, fb, ts, exps) with ts increasing and exps[i] = {l: v_l(Q_{ts[i]})}."""
    mask, m, Q, fb = smooth_mask(N, y, T, fb)
    sel_all = np.nonzero(mask)[0]
    pos = {int(i): j for j, i in enumerate(sel_all)}
    exps: List[Dict[int, int]] = [dict() for _ in range(sel_all.size)]
    Qs = Q[sel_all].copy()
    # l = 2
    low = Qs & (~Qs + np.uint64(1))
    v2 = np.zeros(Qs.size, dtype=np.int64)
    tmp = low.copy()
    while True:
        nz = tmp > 1
        if not nz.any():
            break
        v2[nz] += 1
        tmp[nz] >>= np.uint64(1)
    for j in np.nonzero(v2)[0]:
        exps[int(j)][2] = int(v2[j])
    for l in fb.primes:
        if l == 2:
            continue
        L = np.uint64(l)
        for t0 in fb.roots[l]:
            first = t0 if t0 >= 1 else t0 + l
            if first > T:
                continue
            idx = np.arange(first - 1, T, l, dtype=np.int64)
            idx = idx[mask[idx]]
            if idx.size == 0:
                continue
            js = np.fromiter((pos[int(i)] for i in idx), dtype=np.int64, count=idx.size)
            vals = Q[idx].copy()
            e = np.zeros(idx.size, dtype=np.int64)
            cur = np.arange(idx.size)
            while cur.size:
                mk = (vals[cur] % L) == 0
                cur = cur[mk]
                if cur.size == 0:
                    break
                vals[cur] //= L
                e[cur] += 1
            for j, k in zip(js, e):
                if k:
                    exps[int(j)][l] = exps[int(j)].get(l, 0) + int(k)
    ts = [int(i) + 1 for i in sel_all]
    return m, fb, ts, exps, Q
