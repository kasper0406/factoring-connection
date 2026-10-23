"""Continued-fraction (Euclid-path) structure of a/N.  OPERATIONAL: sees only the integers it is given.

Conventions (notes/classE_results.md section 1).  For 0 < a < N, gcd(a, N) = 1:
  a/N = [0; a_1, ..., a_r],  remainders  r_0 = N, r_1 = a, r_{i+1} = r_{i-1} mod r_i  (r_r = gcd = 1, r_{r+1} = 0),
  convergents h_k/s_k (k = 0..r; h_0 = 0, s_0 = 1; s_r = N),
  identities:  a s_k - N h_k = (-1)^k r_{k+1}  (so a s_k = +-r_{k+1} mod N),     [checked in tests]
               N = s_k r_k + s_{k-1} r_{k+1}   (k >= 1).

Every function here is exact integer arithmetic, O(log N) steps.
"""
from __future__ import annotations

import hashlib
from math import gcd
from typing import Dict, List, Optional, Tuple

import gmpy2


def cf(a: int, N: int) -> List[int]:
    """Partial quotients a_1..a_r of a/N (0 < a < N)."""
    out = []
    x, y = N, a % N
    while y:
        out.append(x // y)
        x, y = y, x % y
    return out


def remainders(a: int, N: int) -> List[int]:
    """r_0 = N, r_1 = a, ..., down to the gcd (last nonzero)."""
    rs = [N, a % N]
    while rs[-1]:
        rs.append(rs[-2] % rs[-1])
    return rs[:-1]


def convergents(qs: List[int]) -> Tuple[List[int], List[int]]:
    """(h_k, s_k) for k = 0..r of [0; a_1..a_r]; h_0 = 0, s_0 = 1."""
    h_prev, h = 1, 0
    s_prev, s = 0, 1
    H, S = [h], [s]
    for x in qs:
        h_prev, h = h, x * h + h_prev
        s_prev, s = s, x * s + s_prev
        H.append(h)
        S.append(s)
    return H, S


def path_gcd_hit(a: int, N: int) -> Optional[int]:
    """Rule 3 / Prop VIS diagnostic that uses only N: the first proper gcd(r_i, N) along the Euclid path
    (equivalently gcd(s_k, N) by a s_k = +-r_{k+1} mod N), else None."""
    for r in remainders(a, N)[1:]:
        g = gcd(r, N)
        if 1 < g < N:
            return g
    return None


def middle_index(S: List[int], N: int) -> int:
    """m with s_m <= sqrt(N) < s_{m+1} (the junction of the two halves)."""
    r = int(gmpy2.isqrt(N))
    m = 0
    while m + 1 < len(S) and S[m + 1] <= r:
        m += 1
    return m


def two_halves(a: int, N: int) -> Dict:
    """Two-halves decomposition (Lemma TH): prefix (a_1..a_m) is the CF of a/N cut at s_m <= sqrt N;
    the reversed suffix (a_r, ..., a_{m+2}) is the corresponding prefix of the CF of a*/N with
    a* = s_{r-1} = (-1)^{r-1} a^{-1} mod N; a_{m+1} is the junction quotient."""
    qs = cf(a, N)
    H, S = convergents(qs)
    m = middle_index(S, N)
    r = len(qs)
    astar = S[r - 1] if r >= 1 else 0
    return {"qs": qs, "m": m, "prefix": qs[:m], "junction": qs[m] if m < r else None,
            "suffix_rev": list(reversed(qs[m + 1:])), "astar": astar, "r": r}


# ------------------------------------------------------------------------------------------------------------
# structured arguments a(N) (registered in notes/classE_prereg.md)
# ------------------------------------------------------------------------------------------------------------
def arguments(N: int) -> Dict[str, int]:
    r = int(gmpy2.isqrt(N))
    return {
        "isqrt": r,
        "isqrt2": int(gmpy2.isqrt(2 * N)),
        "cbrt2": int(gmpy2.iroot(N * N, 3)[0]),
        "div3": N // 3 - 1,
        "fermres": pow(2, N - 1, N),
        "euler3": pow(3, (N - 1) // 2, N),
        "invsqrt": pow(r, -1, N) if gcd(r, N) == 1 else 0,
        "hash2": 1 + int.from_bytes(hashlib.sha256(("classE|" + str(N)).encode()).digest(), "big") % (N - 1),
    }


def sqrt_decomposition(N: int) -> Dict[str, int]:
    """Prop SQ-CF: a = isqrt(N) = r, N = r(r + t) + s', 0 <= s' < r, t in {0,1,2}; CF(r/N) = [0; r+t, CF(s'/r)]."""
    r = int(gmpy2.isqrt(N))
    s = N - r * r
    t, sp = divmod(s, r)
    return {"r": r, "s": s, "t": t, "sprime": sp}
