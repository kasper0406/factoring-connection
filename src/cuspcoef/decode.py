"""OPERATIONAL: factor N = pq from the products A_i = a_N(E_i) = a_p(E_i) a_q(E_i) of a fixed public curve family.

Inputs: N, the public bit size n (instance law p, q in [X, 2X), X = 2^{n/2-1}), the curve list, the list of A_i.
Never sees p, q; never factors N; no PARI.  The integers |A_i| (about n/2 + 2 bits) ARE factored: that is the
"half-size factoring oracle" of the reduction (toy: smallest-prime-factor sieve).

Modes:
  'full' : candidate splits of A_i into x*y with |x|, |y| <= R  (uses the factorisation of A_i);
  'cong' : congruence-only, A_i mod l (no factorisation; Theorem SC says: quadratic-residue content only).
"""
from __future__ import annotations

import heapq
from math import gcd, isqrt, log, prod
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


def legendre(a: int, l: int) -> int:
    a %= l
    if a == 0:
        return 0
    return 1 if pow(a, (l - 1) // 2, l) == 1 else -1


def class_counts(l: int) -> np.ndarray:
    """n[t, d] = #{g in GL2(F_l): tr g = t, det g = d} = l (l + chi(t^2 - 4d)) for d != 0 (row d = 0 is zero)."""
    n = np.zeros((l, l), dtype=np.float64)
    for t in range(l):
        for d in range(1, l):
            n[t, d] = l * (l + legendre(t * t - 4 * d, l))
    return n


def odd_primes_upto(L: int) -> List[int]:
    return [l for l in range(3, L + 1) if all(l % r for r in range(2, isqrt(l) + 1))]


def decoding_primes(X: int) -> List[int]:
    """Least list of consecutive odd primes with product > 2X (then s in [2 sqrt N, 4X] is fixed by s mod 2W)."""
    out: List[int] = []
    l = 3
    while prod(out) <= 2 * X:
        if all(l % r for r in range(2, isqrt(l) + 1)):
            out.append(l)
        l += 2
    return out


# ------------------------------------------------------------------ half-size factoring (toy)

class SPF:
    """Smallest-prime-factor table up to `limit` (toy factoring oracle for |A_i| <= limit)."""

    def __init__(self, limit: int):
        self.limit = limit
        spf = np.zeros(limit + 1, dtype=np.uint32)
        for i in range(2, isqrt(limit) + 1):
            if spf[i] == 0:
                sl = spf[i * i::i]
                sl[sl == 0] = i
        self.spf = spf

    def decompose(self, m: int) -> Dict[int, int]:
        assert 0 < m <= self.limit
        out: Dict[int, int] = {}
        while m > 1:
            r = int(self.spf[m]) or m
            out[r] = out.get(r, 0) + 1
            m //= r
        return out


def divisors(fac: Dict[int, int]) -> List[int]:
    ds = [1]
    for r, e in fac.items():
        ds = [d * r ** k for d in ds for k in range(e + 1)]
    return ds


def candidate_pairs(A: int, R: int, spf: SPF) -> List[Tuple[int, int]]:
    """All (x, y) in Z^2 with x y = A, |x|, |y| <= R (A != 0)."""
    m = abs(A)
    out = []
    for d in divisors(spf.decompose(m)):
        e = m // d
        if d <= R and e <= R:
            sgn = 1 if A > 0 else -1
            out.append((d, sgn * e))
            out.append((-d, -sgn * e))
    return out


# ------------------------------------------------------------------ likelihoods

def loglik_d_fast(l: int, Nl: int, cur_idx: np.ndarray, xs: np.ndarray, ys: np.ndarray, ncur: int) -> np.ndarray:
    """Vectorised: candidates (xs, ys) residues mod l with curve index cur_idx (0..ncur-1)."""
    n = class_counts(l)
    out = np.full(l, -np.inf)
    for d in range(1, l):
        dq = Nl * pow(d, -1, l) % l
        w = n[xs, d] * n[ys, dq]
        per = np.bincount(cur_idx, weights=w, minlength=ncur)
        out[d] = float(np.log(per).sum())
    return out


def loglik_d_cong(l: int, Nl: int, Ares: np.ndarray) -> np.ndarray:
    """Congruence-only: sum_i log sum_{t t' = A_i mod l} n(t, d) n(t', N/d)."""
    n = class_counts(l)
    out = np.full(l, -np.inf)
    counts = np.bincount(Ares % l, minlength=l).astype(float)
    for d in range(1, l):
        dq = Nl * pow(d, -1, l) % l
        lik = np.zeros(l)
        for t in range(l):
            for t2 in range(l):
                lik[t * t2 % l] += n[t, d] * n[t2, dq]
        lik /= lik.sum()
        m = counts > 0
        out[d] = float((counts[m] * np.log(lik[m])).sum()) if np.all(lik[m] > 0) else -np.inf
    return out


def s_posterior(l: int, Nl: int, lld: np.ndarray) -> List[Tuple[float, int]]:
    """Log posterior of s mod l (uniform prior on d = p mod l), sorted descending: [(logpost, s), ...]."""
    acc: Dict[int, List[float]] = {}
    for d in range(1, l):
        s = (d + Nl * pow(d, -1, l)) % l
        acc.setdefault(s, []).append(lld[d])
    out = []
    for s, v in acc.items():
        mx = max(v)
        if mx == -np.inf:
            out.append((-np.inf, s))
        else:
            out.append((mx + log(sum(np.exp(np.array(v) - mx))), s))
    out.sort(reverse=True)
    return out


def crt(res: Sequence[int], mods: Sequence[int]) -> Tuple[int, int]:
    x, m = 0, 1
    for r, mm in zip(res, mods):
        t = ((r - x) * pow(m, -1, mm)) % mm
        x, m = x + m * t, m * mm
    return x % m, m


def list_decode(N: int, X: int, ells: Sequence[int], posts: Sequence[List[Tuple[float, int]]],
                max_combos: int = 3000) -> Tuple[Optional[int], int]:
    """Best-first search over combinations of per-l candidates; returns (p or None, combos tried)."""
    lo, hi = 2 * isqrt(N), 4 * X
    k = len(ells)
    start = tuple([0] * k)

    def score(ix):
        return sum(posts[j][ix[j]][0] for j in range(k))

    heap = [(-score(start), start)]
    seen = {start}
    tried = 0
    while heap and tried < max_combos:
        _, ix = heapq.heappop(heap)
        tried += 1
        s0, W = crt([posts[j][ix[j]][1] for j in range(k)] + [0], list(ells) + [2])
        s = s0 + ((lo - s0 + W - 1) // W) * W if s0 < lo else s0
        while s <= hi:
            D = s * s - 4 * N
            if D >= 0:
                r = isqrt(D)
                if r * r == D and (s - r) % 2 == 0:
                    p = (s - r) // 2
                    if 1 < p < N and N % p == 0:
                        return p, tried
            s += W
        for j in range(k):
            if ix[j] + 1 < len(posts[j]):
                nx = ix[:j] + (ix[j] + 1,) + ix[j + 1:]
                if nx not in seen:
                    seen.add(nx)
                    heapq.heappush(heap, (-score(nx), nx))
    return None, tried


def decode(N: int, nbits: int, curves: Sequence[Tuple[int, int]], A_list: Sequence[int], mode: str = "full",
           spf: Optional[SPF] = None, max_combos: int = 3000) -> dict:
    """Factor N from products A_i (mode 'full' or 'cong').  Returns a dict with 'p' (or None) and diagnostics."""
    X = 1 << (nbits // 2 - 1)
    R = isqrt(4 * 2 * X)                      # floor(2 sqrt(2X))
    ells = decoding_primes(X)
    accidental = []
    keep = []
    for (A, B), a in zip(curves, A_list):
        g = gcd(16 * (4 * A ** 3 + 27 * B ** 2), N)
        if g != 1:
            if g != N:
                accidental.append(g)      # ACCIDENTAL split (bad reduction at p or q); structural decoder ignores it
            continue
        if a != 0:
            keep.append(a)
    posts = []
    ncand = []
    if mode == "full":
        assert spf is not None
        cur, xs, ys = [], [], []
        for i, a in enumerate(keep):
            C = candidate_pairs(a, R, spf)
            ncand.append(len(C))
            for x, y in C:
                cur.append(i)
                xs.append(x)
                ys.append(y)
        cur_a = np.array(cur, dtype=np.int64)
        xs_a = np.array(xs, dtype=np.int64)
        ys_a = np.array(ys, dtype=np.int64)
        for l in ells:
            lld = loglik_d_fast(l, N % l, cur_a, xs_a % l, ys_a % l, len(keep))
            posts.append(s_posterior(l, N % l, lld))
    elif mode == "cong":
        Ares = np.array(keep, dtype=np.int64)
        for l in ells:
            lld = loglik_d_cong(l, N % l, Ares % l)
            posts.append(s_posterior(l, N % l, lld))
    else:
        raise ValueError(mode)
    p, tried = list_decode(N, X, ells, posts, max_combos=max_combos)
    return {"p": p, "tried": tried, "ells": ells, "rank_lists": [[s for _, s in pl] for pl in posts],
            "n_used": len(keep), "mean_cand": float(np.mean(ncand)) if ncand else None,
            "accidental": accidental}
