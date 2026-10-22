"""dsr core (OPERATIONAL): sees only N and public parameters; no factorisation, no group orders, no randomness.

Contents
- factor base and exact smoothness sieves;
- the good sets S_can (least residue B-smooth) and S_bal (short smooth Gauss split), as boolean masks over [0, N);
- deterministic multiplier lists D_box, D_slide and the orbit targets of the DSR orbit reformulation;
- worst-case coverage m*(N, W): least m with S * {w_1..w_m}^{-1} = G (all targets), and first hits for the DSR targets.

Every inversion mod N is preceded by a gcd check; a nonunit would be returned as a factor (never happens after
preprocessing, because all multipliers are B-smooth and N has no prime factor <= B; asserted).
"""
from __future__ import annotations

import math
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "L_half", "factor_base", "smooth_mask", "units_mask", "s_can", "s_bal", "walk_length",
    "d_box", "d_slide", "coverage", "first_hits", "orbit_first_hits", "dsr_output_from_hit",
]


def L_half(N: int, a: float) -> int:
    """floor(L_N[1/2, a]); floating point only to choose the public parameter B."""
    ln = math.log(N)
    return int(math.exp(a * math.sqrt(ln * math.log(ln))))


def _primes_upto(B: int) -> List[int]:
    if B < 2:
        return []
    s = bytearray([1]) * (B + 1)
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(B) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(B + 1) if s[i]]


def factor_base(N: int, a: float) -> Tuple[int, List[int]]:
    B = L_half(N, a)
    P = _primes_upto(B)
    for p in P:  # preprocessing check: a prime <= B dividing N is a factor (would be output)
        if N % p == 0:
            raise ValueError(f"preprocessing: {p} divides N")
    return B, P


def smooth_mask(X: int, P: Sequence[int]) -> np.ndarray:
    """mask[n] = True iff n in [1, X) is P-smooth (exact integer sieve)."""
    arr = np.arange(X, dtype=np.int64)
    for p in P:
        pk = p
        while pk < X:
            arr[pk::pk] //= p
            pk *= p
    m = arr == 1
    m[0] = False
    return m


def units_mask(N: int) -> np.ndarray:
    return np.gcd(np.arange(N, dtype=np.int64), N) == 1


def s_can(N: int, P: Sequence[int]) -> np.ndarray:
    """S_can: residues x in G whose least positive residue is B-smooth (all such x are units here)."""
    return smooth_mask(N, P)


def s_bal(N: int, P: Sequence[int]) -> np.ndarray:
    """S_bal: x = u * v^{-1} mod N with u, v B-smooth, u != 0 (either sign), v > 0, u^2 + v^2 < N, gcd(u, v) = 1."""
    r = math.isqrt(N - 1)
    sm = np.nonzero(smooth_mask(r + 1, P))[0].astype(np.int64)
    mask = np.zeros(N, dtype=bool)
    for v in sm.tolist():
        if v * v >= N:
            break
        assert math.gcd(v, N) == 1
        vinv = pow(v, -1, N)
        lim = math.isqrt(N - 1 - v * v)
        us = sm[: np.searchsorted(sm, lim, side="right")]
        us = us[np.gcd(us, v) == 1]
        x = (us * vinv) % N
        mask[x] = True
        mask[(N - x) % N] = True
    return mask


def walk_length(N: int, B: int) -> int:
    return int(math.ceil(math.log(N) / math.log(B))) + 1


def d_box(P: Sequence[int]) -> Iterator[int]:
    """D_box: the P-smooth integers 1, 2, 3, ... in increasing order (exact, incremental)."""
    n = 0
    while True:
        n += 1
        m = n
        for p in P:
            while m % p == 0:
                m //= p
            if m == 1:
                break
        if m == 1:
            yield n


def _beatty(i: int) -> int:
    """floor(i * (1 + sqrt 5) / 2), exactly."""
    return (i + math.isqrt(5 * i * i)) // 2


def d_slide(P: Sequence[int], t: int) -> Iterator[int]:
    """D_slide: w_s = prod_{i=s}^{s+t-1} p_{l_i}, l_i = floor(i*golden) mod |P|, s = 1, 2, ... (one explicit walk)."""
    k = len(P)
    s = 0
    while True:
        s += 1
        w = 1
        for i in range(s, s + t):
            w *= P[_beatty(i) % k]
        yield w


def coverage(N: int, S: np.ndarray, G: np.ndarray, ws: Iterator[int], cap: int,
             checkpoints: Sequence[int] = ()) -> Dict:
    """Worst-case coverage over ALL targets a in G.

    Returns m* = least m such that every a in G has a*w_s in S for some s <= m (None if > cap), and the number of
    uncovered targets after each checkpoint m.
    """
    S_idx = np.nonzero(S)[0].astype(np.int64)
    covered = ~G  # non-units are not targets
    U = None  # explicit uncovered list once it is small
    ck = set(int(c) for c in checkpoints)
    out_ck: Dict[int, int] = {}
    n_unc = int(G.sum())
    m = 0
    for w in ws:
        m += 1
        if m > cap:
            return {"m_star": None, "uncovered_at": out_ck, "uncovered_final": n_unc}
        wm = w % N
        assert math.gcd(wm, N) == 1, "nonunit multiplier (would be a factor)"
        if U is None:
            covered[(S_idx * pow(wm, -1, N)) % N] = True
            n_unc = int((~covered).sum())
            if n_unc <= len(S_idx) // 4:
                U = np.nonzero(~covered)[0].astype(np.int64)
        else:
            U = U[~S[(U * wm) % N]]
            n_unc = len(U)
        if m in ck:
            out_ck[m] = n_unc
        if n_unc == 0:
            return {"m_star": m, "uncovered_at": out_ck, "uncovered_final": 0}
    return {"m_star": None, "uncovered_at": out_ck, "uncovered_final": n_unc}


def first_hits(N: int, S: np.ndarray, targets: Sequence[int], ws: Iterator[int], cap: int) -> Dict[int, Optional[int]]:
    """For each target a, the first s <= cap with a*w_s mod N in S."""
    todo = {a % N: None for a in targets}
    pending = set(todo)
    s = 0
    for w in ws:
        s += 1
        if s > cap or not pending:
            break
        wm = w % N
        for a in list(pending):
            if S[(a * wm) % N]:
                todo[a] = s
                pending.discard(a)
    return todo


def orbit_first_hits(N: int, S: np.ndarray, P: Sequence[int], H: int, cap: int) -> Dict[int, Optional[int]]:
    """D_orbit: for each q in P, the first j >= 0 (j < cap) with q^(H+j) mod N in S."""
    res = {}
    for q in P:
        x = pow(q, H, N)
        hit = None
        for j in range(cap):
            if S[x]:
                hit = j
                break
            x = (x * q) % N
        res[q] = hit
    return res


def dsr_output_from_hit(N: int, a: int, w: int, P: Sequence[int]) -> Tuple[int, int, int]:
    """Given a*w mod N with B-smooth least residue, return the DSR triple (u, v, w) with v = 1, checked."""
    u = (a * w) % N
    m = u
    for p in P:
        while m % p == 0:
            m //= p
    assert m == 1 and (u - a * w) % N == 0
    return u, 1, w
