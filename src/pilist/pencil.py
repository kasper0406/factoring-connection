"""PTE O2: pencils and linear systems of split polynomials (agent pilist, Task 4).

Exact integer arithmetic (numpy int64 on coefficient vectors that provably fit: |e_i| <= C(M,i) B^i < 2^62
for the cells used; `check_fits` asserts it).

Objects
- A split polynomial P_A = prod_{a in A}(X - a), A a multiset of M integers in [-B, B], is stored as its
  coefficient vector (c_1, ..., c_M) with P = X^M + c_1 X^{M-1} + ... + c_M.
- A *line* (pencil) is {F0 + t D : t in Z} with D primitive; its *members* are the split P_A on it.
- k-projected lines (O2 with remainders of degree < k): lines in the first M - k coordinates.

Lemma PEN (notes/pilist_results.md, PROVED from Binyamini-Cluckers-Novikov): on a pencil with fixed part T
(the multiset intersection of all members) and moving degree e >= 3, distinct members have disjoint moving
roots, at most 2e - 2 members have a repeated moving root, and the number of members is
<= (16/3)(c B^{1/2}(log B)^kappa + 1) + 2e - 2.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations, combinations_with_replacement
from math import comb, gcd
from typing import Dict, List, Sequence, Tuple

import numpy as np


def coeffs_of_roots(roots: Sequence[int]) -> List[int]:
    c = [1]
    for a in roots:
        nc = c + [0]
        for i in range(1, len(nc)):
            nc[i] -= a * c[i - 1]
        c = nc
    return c[1:]


def check_fits(M: int, B: int) -> bool:
    return max(comb(M, i) * B ** i for i in range(M + 1)) < 2 ** 62


def enumerate_split(M: int, B: int) -> Tuple[List[Tuple[int, ...]], np.ndarray]:
    assert check_fits(M, B)
    ms = list(combinations_with_replacement(range(-B, B + 1), M))
    arr = np.array([coeffs_of_roots(A) for A in ms], dtype=np.int64)
    return ms, arr


def _primitive(D: np.ndarray) -> np.ndarray:
    g = np.gcd.reduce(np.abs(D), axis=1)
    g[g == 0] = 1
    D = D // g[:, None]
    # sign: first nonzero coordinate positive
    nz = (D != 0)
    first = np.argmax(nz, axis=1)
    sgn = np.sign(D[np.arange(len(D)), first])
    sgn[sgn == 0] = 1
    return D * sgn[:, None]


def find_lines(arr: np.ndarray, kproj: int = 0, min_members: int = 3) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """All lines in coefficient space (first M - kproj coordinates) carrying >= min_members points.
    Returns [(sorted member indices, primitive direction)], each line once."""
    A = arr[:, : arr.shape[1] - kproj] if kproj else arr
    n = len(A)
    seen = set()
    out = []
    for i in range(n):
        D = A - A[i]
        D[i] = 0
        mask = np.any(D != 0, axis=1)
        idx = np.nonzero(mask)[0]
        if len(idx) == 0:
            continue
        Dp = _primitive(D[idx])
        # group equal directions
        v = np.ascontiguousarray(Dp).view(np.dtype((np.void, Dp.dtype.itemsize * Dp.shape[1]))).ravel()
        order = np.argsort(v, kind="stable")
        vs = v[order]
        breaks = np.nonzero(vs[1:] != vs[:-1])[0] + 1
        groups = np.split(order, breaks)
        for g in groups:
            if len(g) + 1 < min_members:
                continue
            members = tuple(sorted([i] + [int(idx[t]) for t in g]))
            # in kproj > 0 mode, equal projections (D = 0) are excluded above; record each line once
            if members[0] != i:
                continue
            if members in seen:
                continue
            seen.add(members)
            out.append((members, tuple(int(x) for x in Dp[g[0]])))
    return out


def multiset_intersection(ms: Sequence[Tuple[int, ...]]) -> Counter:
    c = Counter(ms[0])
    for m in ms[1:]:
        c &= Counter(m)
    return c


def line_stats(ms_list: List[Tuple[int, ...]], members: Tuple[int, ...]) -> Dict:
    roots = [ms_list[i] for i in members]
    T = multiset_intersection(roots)
    degT = sum(T.values())
    M = len(roots[0])
    e = M - degT
    moving = [Counter(r) - T for r in roots]
    # Claim A: pairwise disjoint moving supports
    disjoint = True
    for x, y in combinations(moving, 2):
        if set(x) & set(y):
            disjoint = False
            break
    nonsqf = sum(1 for mv in moving if any(v > 1 for v in mv.values()))
    return {"m": len(members), "degT": degT, "e": e, "disjoint": disjoint, "nonsqf": nonsqf}


def line_params(arr: np.ndarray, members: Tuple[int, ...], D: Sequence[int], kproj: int = 0) -> List[int]:
    """Integer positions t_j with A_j = A_{m0} + t_j D on the (projected) line."""
    A = arr[:, : arr.shape[1] - kproj] if kproj else arr
    D = np.array(D, dtype=np.int64)
    k = int(np.nonzero(D)[0][0])
    base = A[members[0]]
    ts = []
    for i in members:
        diff = A[i] - base
        t = int(diff[k]) // int(D[k])
        assert np.array_equal(diff, t * D)
        ts.append(t)
    return ts


def small_c_triples(ts: Sequence[int], cmax: int = 2) -> int:
    """Number of member triples whose primitive three-term relation c = (t1-t2, t2-t0, t0-t1)/g has
    max |c_i| <= cmax (c0 P0 + c1 P1 + c2 P2 is then low-degree with small coefficients)."""
    cnt = 0
    for a, b, c in combinations(ts, 3):
        if len({a, b, c}) < 3:      # a repeated position is a two-term (fibre) relation, not a three-term one
            continue
        v = (b - c, c - a, a - b)
        g = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
        if max(abs(x) for x in v) // g <= cmax:
            cnt += 1
    return cnt


def tl_class(ms_list, i0: int, i1: int, i2: int, L: int) -> str:
    """Lemma TL dichotomy for a three-term configuration: 'A' common factor of degree >= M - 2L,
    'B' some point lies in exactly two of the reduced supports (forcing localisation of the third),
    'C' reduced polynomials pairwise coprime."""
    roots = [Counter(ms_list[i]) for i in (i0, i1, i2)]
    g = roots[0] & roots[1] & roots[2]
    M = sum(roots[0].values())
    if M - sum(g.values()) <= 2 * L:
        return "A"
    red = [r - g for r in roots]
    sup = [set(r) for r in red]
    for x in range(3):
        for y in range(x + 1, 3):
            if sup[x] & sup[y]:
                return "B"
    return "C"


def pencil_count_bound_trivial(B: int) -> int:
    """Disjoint moving roots (Claim A) give at most 2B + 1 members on any pencil with e >= 1."""
    return 2 * B + 1
