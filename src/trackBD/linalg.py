"""Track B/D operational linear algebra over Z/N (agent sketchBD).

OPERATIONAL: inputs are N and public data only.  All pivots are gcd-checked.

* elim_profile(M, N)    column-gcd Gaussian elimination.  For each column, after
                        eliminating earlier pivots, g = gcd(N, remaining entries):
                          g == N  -> no pivot (column dependent mod every p | N)
                          1<g<N   -> returns factor g  (pivot profiles differ mod p, q)
                          g == 1  -> choose a unit pivot (an entry, else a seeded random
                                     row combination that is a unit) and continue.
                        Proposition B3 (notes): returns a factor iff the column rank
                        profiles of M mod the primes of N differ (squarefree N).
                        Also records every *accidental* nonunit entry seen (counter).
* Matrix builders (all from N and public parameters):
    mult_matrix(N, f, S)          columns X^i S mod f              (Section B3, M1)
    hankel_proj(N, f, S, ell, L)  [u_{i+j}], u_k = ell(X^k S mod f) (M2, Wiedemann)
    berlekamp_QN(N, f)            columns X^(iN) mod f, and Q_N - I (M3)
    hankel_over_c(seq, L)         Hankel matrix of a scalar sequence    (M4)
* bm_mod_N(seq, N)       Berlekamp-Massey over Z/N, gcd-checked discrepancies.
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Sequence

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import FlintQuotient, power_of_linear


def elim_profile(M: Sequence[Sequence[int]], N: int, seed: int = 0,
                 m: Optional[Metrics] = None) -> Dict:
    m = m if m is not None else NullMetrics()
    rng = random.Random(f"elim|{N}|{seed}")
    A = [[x % N for x in row] for row in M]
    nrows = len(A)
    ncols = len(A[0]) if nrows else 0
    used = [False] * nrows
    profile = []
    acc_nonunits = 0
    for j in range(ncols):
        rows = [i for i in range(nrows) if not used[i]]
        col = [A[i][j] for i in rows]
        g = N
        for v in col:
            g = gcd(g, v)
            m.inc("gcd")
        if g == N:
            continue
        if g != 1:
            return {"factor": g, "column": j, "profile": profile, "acc_nonunits": acc_nonunits}
        piv = None
        for i in rows:
            v = A[i][j]
            if v:
                m.inc("gcd")
                if gcd(v, N) == 1:
                    piv = i
                    break
                acc_nonunits += 1
        if piv is None:
            # combine rows: random combination is a unit (entries generate the unit ideal)
            for _ in range(64):
                coef = {i: rng.randrange(N) for i in rows}
                v = sum(coef[i] * A[i][j] for i in rows) % N
                if gcd(v, N) == 1:
                    newrow = [sum(coef[i] * A[i][k] for i in rows) % N for k in range(ncols)]
                    # replace the first unused row whose coefficient is a unit
                    for i in rows:
                        if gcd(coef[i], N) == 1:
                            A[i] = newrow
                            piv = i
                            break
                    if piv is not None:
                        break
            if piv is None:  # pragma: no cover (probability ~ (2/p)^64)
                return {"factor": None, "column": j, "profile": profile, "stuck": True,
                        "acc_nonunits": acc_nonunits}
        used[piv] = True
        profile.append(j)
        inv = pow(A[piv][j], -1, N)
        m.inc("inv")
        prow = [(x * inv) % N for x in A[piv]]
        A[piv] = prow
        for i in range(nrows):
            if i != piv and A[i][j]:
                t = A[i][j]
                A[i] = [(a - t * b) % N for a, b in zip(A[i], prow)]
                m.inc("rn_mul", ncols)
    return {"factor": None, "profile": profile, "rank": len(profile), "acc_nonunits": acc_nonunits}


# ---------------------------------------------------------------------------
# matrix builders
# ---------------------------------------------------------------------------
def mult_matrix(N: int, f: Sequence[int], S: Sequence[int]) -> List[List[int]]:
    """Matrix (rows = coefficient index, cols = i) of multiplication by S on R_N[X]/(f)."""
    Q = FlintQuotient(N, list(f))
    s = Q.reduce(list(S))
    x = Q.gen()
    cols = []
    cur = s
    for i in range(Q.d):
        cols.append(Q.coeffs(cur))
        cur = Q.mul(cur, x)
    return [[cols[i][k] for i in range(Q.d)] for k in range(Q.d)]


def hankel_proj(N: int, f: Sequence[int], S: Sequence[int], ell: Sequence[int], L: int) -> List[List[int]]:
    """Hankel [u_{i+j}]_{0<=i,j<L}, u_k = <ell, X^k S mod f>."""
    Q = FlintQuotient(N, list(f))
    cur = Q.reduce(list(S))
    x = Q.gen()
    u = []
    for _ in range(2 * L - 1):
        cs = Q.coeffs(cur)
        u.append(sum(a * b for a, b in zip(ell, cs)) % N)
        cur = Q.mul(cur, x)
    return [[u[i + j] for j in range(L)] for i in range(L)]


def berlekamp_QN(N: int, f: Sequence[int], minus_identity: bool = False) -> List[List[int]]:
    """Columns X^(i N) mod (f, N), i = 0..r-1 (the matrix of g(X) -> g(X^N))."""
    Q = FlintQuotient(N, list(f))
    xN = Q.pow(Q.gen(), N)
    cols = []
    cur = Q.one
    for i in range(Q.d):
        cols.append(Q.coeffs(cur))
        cur = Q.mul(cur, xN)
    M = [[cols[i][k] for i in range(Q.d)] for k in range(Q.d)]
    if minus_identity:
        for k in range(Q.d):
            M[k][k] = (M[k][k] - 1) % N
    return M


def hankel_over_c(seq: Sequence[int], L: int) -> List[List[int]]:
    return [[seq[i + j] for j in range(L)] for i in range(L)]


# ---------------------------------------------------------------------------
# Berlekamp-Massey over Z/N
# ---------------------------------------------------------------------------
def bm_mod_N(seq: Sequence[int], N: int, m: Optional[Metrics] = None) -> Dict:
    """Berlekamp-Massey over Z/N.  Every nonzero discrepancy is gcd-checked before it
    is stored for later inversion; a nonunit discrepancy (nonzero mod N) returns
    gcd(d, N).  Returns the connection polynomial C (low->high, C[0]=1) and the
    linear complexity L if no nonunit occurs."""
    m = m if m is not None else NullMetrics()
    s = [x % N for x in seq]
    C = [1]
    B = [1]
    L = 0
    mm = 1
    b = 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            if i < len(C):
                d = (d + C[i] * s[n - i]) % N
        if d == 0:
            mm += 1
            continue
        m.inc("gcd")
        g = gcd(d, N)
        if g != 1:
            return {"factor": g if g < N else None, "step": n, "L": L, "C": C}
        coef = d * pow(b, -1, N) % N
        T = C[:]
        need = len(B) + mm
        if len(C) < need:
            C = C + [0] * (need - len(C))
        for i, bi in enumerate(B):
            C[i + mm] = (C[i + mm] - coef * bi) % N
        if 2 * L <= n:
            L = n + 1 - L
            B = T
            b = d
            mm = 1
        else:
            mm += 1
    while len(C) > 1 and C[-1] == 0:
        C.pop()
    return {"factor": None, "L": L, "C": C}


def bm_field(seq: Sequence[int], p: int) -> Dict:
    """Harness helper: BM over the prime field F_p (p known).  Linear complexity profile."""
    s = [x % p for x in seq]
    C, B, L, mm, b = [1], [1], 0, 1, 1
    prof = []
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            if i < len(C):
                d = (d + C[i] * s[n - i]) % p
        if d == 0:
            mm += 1
        else:
            coef = d * pow(b, -1, p) % p
            T = C[:]
            need = len(B) + mm
            if len(C) < need:
                C = C + [0] * (need - len(C))
            for i, bi in enumerate(B):
                C[i + mm] = (C[i + mm] - coef * bi) % p
            if 2 * L <= n:
                L = n + 1 - L
                B, b, mm = T, d, 1
            else:
                mm += 1
        prof.append(L)
    return {"L": L, "profile": prof}
