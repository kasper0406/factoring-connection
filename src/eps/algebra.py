"""Track A section values, Track B sketches, Berlekamp-Massey lengths, ranks.  OPERATIONAL: input N only.

Section/sketch logic is copied (rule 7) from src/trackA/sections.py and src/trackBD/sketch.py; the
quotient-ring arithmetic uses the shared FlintQuotient (src/common/ring.py; flint audit in its README).
Every inversion mod N is gcd-checked.  Every value whose Jacobi symbol is taken is gcd-checked too, and a
nonunit is recorded through `Nonunits` (rule 3); Jacobi of a nonunit is 0.
"""
from __future__ import annotations

from math import gcd
from typing import List, Optional, Sequence, Tuple

import gmpy2

from src.common.ring import FlintQuotient, power_of_linear


class Nonunits:
    """Collects gcd(value, N) > 1 events: (where, divisor)."""

    def __init__(self, N: int):
        self.N = N
        self.events: List[Tuple[str, int]] = []
        self.checks = 0

    def check(self, where: str, v: int) -> int:
        self.checks += 1
        g = gcd(v % self.N, self.N)
        if g != 1:
            self.events.append((where, g))
        return g


def jac(v: int, N: int, nu: Optional[Nonunits] = None, where: str = "") -> int:
    v %= N
    if nu is not None:
        nu.check(where, v)
    return int(gmpy2.jacobi(v, N))


def section_values(N: int, M: int, c: int) -> List[int]:
    """S_{r,M}(c), r < M, for F_N = (1+X)^N - 1 - X^N (copy of trackA.sections.section_values)."""
    c %= N
    if M == 1:
        return [(pow(1 + c, N, N) - 1 - pow(c, N, N)) % N]
    Q = FlintQuotient(N, [(-c) % N] + [0] * (M - 1) + [1])
    vals = list(Q.coeffs(power_of_linear(Q, N)))
    vals[0] = (vals[0] - 1) % N
    vals[N % M] = (vals[N % M] - pow(c, N // M, N)) % N
    return vals


def sketch(N: int, f: Sequence[int]) -> List[int]:
    """F_N mod (f, N), f monic integer list low->high (copy of trackBD.sketch.sketch)."""
    Q = FlintQuotient(N, list(f))
    S = Q.sub(Q.sub(power_of_linear(Q, N), Q.one), Q.pow(Q.gen(), N))
    return [int(x) for x in Q.coeffs(S)]


def sketch_G(N: int, f: Sequence[int], k: int = 1) -> List[int]:
    """G_k = (1+X)^(N-k) - X^(N-k) mod (f, N) (copy of trackBD.sketch.sketch_G)."""
    Q = FlintQuotient(N, list(f))
    S = Q.sub(power_of_linear(Q, N - k), Q.pow(Q.gen(), N - k))
    return [int(x) for x in Q.coeffs(S)]


def mult_matrix(N: int, f: Sequence[int], S: Sequence[int]) -> List[List[int]]:
    """Matrix (columns = S * X^j mod f) of multiplication by S in (Z/N)[X]/(f), f monic of degree d."""
    d = len(f) - 1
    cols = []
    cur = [int(x) % N for x in S] + [0] * (d - len(S))
    for _ in range(d):
        cols.append(list(cur))
        # cur <- cur * X mod f
        top = cur[-1]
        nxt = [0] + cur[:-1]
        cur = [(nxt[t] - top * f[t]) % N for t in range(d)]
    return [[cols[j][i] for j in range(d)] for i in range(d)]


def rank_mod_N(Mx: List[List[int]], N: int, nu: Nonunits, where: str) -> int:
    """Rank over Z/N by elimination with gcd-checked pivots.  A nonzero nonunit pivot is recorded and the
    elimination returns -1 (rank undefined over Z/N: N split)."""
    A = [list(r) for r in Mx]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i][c] % N:
                piv = i
                break
        if piv is None:
            continue
        if nu.check(where, A[piv][c]) != 1:
            return -1
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], -1, N)
        for i in range(rows):
            if i != r and A[i][c] % N:
                t = A[i][c] * inv % N
                A[i] = [(x - t * y) % N for x, y in zip(A[i], A[r])]
        r += 1
    return r


def det_mod_N(Mx: List[List[int]], N: int, nu: Nonunits, where: str) -> Optional[int]:
    """Determinant over Z/N (fraction-free Bareiss is overkill at d <= 4: gcd-checked Gaussian elimination).
    Returns None if a nonunit pivot appears (recorded)."""
    A = [list(r) for r in Mx]
    n = len(A)
    det = 1
    for c in range(n):
        piv = None
        for i in range(c, n):
            if A[i][c] % N:
                piv = i
                break
        if piv is None:
            return 0
        if nu.check(where, A[piv][c]) != 1:
            return None
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det = det * A[c][c] % N
        inv = pow(A[c][c], -1, N)
        for i in range(c + 1, n):
            t = A[i][c] * inv % N
            if t:
                A[i] = [(x - t * y) % N for x, y in zip(A[i], A[c])]
    return det % N


def bm_mod_N(seq: Sequence[int], N: int, nu: Nonunits, where: str) -> Tuple[int, List[int]]:
    """Berlekamp-Massey over Z/N with gcd-checked inversion of discrepancies.

    Returns (L, C) with C the connection polynomial (C[0] = 1).  If a discrepancy is a nonzero nonunit,
    it is recorded and (-1, []) is returned."""
    s = [x % N for x in seq]
    C, B = [1], [1]
    L, m, b = 0, 1, 1
    for n_ in range(len(s)):
        d = s[n_]
        for i in range(1, L + 1):
            d = (d + C[i] * s[n_ - i]) % N
        if d == 0:
            m += 1
            continue
        if nu.check(where, b) != 1:
            return -1, []
        coef = d * pow(b, -1, N) % N
        T = list(C)
        if len(C) < len(B) + m:
            C = C + [0] * (len(B) + m - len(C))
        for i, bi in enumerate(B):
            C[i + m] = (C[i + m] - coef * bi) % N
        if 2 * L <= n_:
            L = n_ + 1 - L
            B, b, m = T, d, 1
        else:
            m += 1
    return L, C[:L + 1] + [0] * max(0, L + 1 - len(C))
