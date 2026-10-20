"""Bivariate Coppersmith lattice for positive-width WINDOW queries at ladder blocks (E6).

Question (orchestrator, after floorhard): can lattice small-root methods answer the ladder's window
queries {P_j(h)} in [0, rho] on a top block (x ~ sqrt N) cheaply?

Formulation.  P_j(h) = num(h)/D, D = x^{j+1}.  {P_j(h)} in [0, R/D]  <=>  exists r in [0, R] with
f(h, r) := r - g(h) == 0 (mod D), where g = num reduced mod D.  Basic-strategy lattice
(Jochemsz-May style; leading monomial r, f monic in r):
    G_{u,k}(h, r) = h^u f(h, r)^k D^{m-k},   0 <= k <= m,  0 <= u <= j(m-k),
monomials h^a r^b with b <= m, a <= j(m-b); triangular, det = prod X^u Y^k D^{m-k}.
A reduced vector v with sqrt(n) ||v(hX, rY)|| < D^m vanishes over Z at every root with
0 <= h < X, 0 <= r < Y (Howgrave-Graham).  Two such vectors give the roots by a resultant
(heuristic: algebraic independence).  As m -> infinity the determinant condition is X^j * Y < D
(PROVED below in the notes: the exponent sums are j^2 m^3/6, j m^3/6, j m^3/3 against m * j m^2/2).
With Y = rho*D this is X^j < 1/rho: windows can be resolved only on sub-blocks of length
~ rho^{-1/j} ~ H^{1/j}, i.e. ~ H^{1-1/j} lattice calls per block -- worse than the Farey oracle
(H^{2(j-1)/(3j)}) for every j.  HARNESS/DIAGNOSTIC measurement code; nothing here factors.
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

import flint


def _polymul2(a: Dict[Tuple[int, int], int], b: Dict[Tuple[int, int], int]):
    out: Dict[Tuple[int, int], int] = {}
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            k = (i1 + i2, j1 + j2)
            out[k] = out.get(k, 0) + c1 * c2
    return {k: v for k, v in out.items() if v}


def window_lattice(g: List[int], D: int, j: int, m: int, X: int, Y: int):
    """Return (n, monomials, LLL-reduced rows) for f = r - g(h) mod D."""
    f = {(0, 1): 1}
    for i, c in enumerate(g):
        if c % D:
            f[(i, 0)] = f.get((i, 0), 0) - c
    fpow = [{(0, 0): 1}]
    for k in range(m):
        fpow.append(_polymul2(fpow[-1], f))
    monos = [(a, b) for b in range(m + 1) for a in range(j * (m - b) + 1)]
    idx = {mo: t for t, mo in enumerate(monos)}
    rows = []
    for k in range(m + 1):
        for u in range(j * (m - k) + 1):
            row = [0] * len(monos)
            for (a, b), c in fpow[k].items():
                row[idx[(a + u, b)]] += c * D ** (m - k)
            rows.append(row)
    scale = [X ** a * Y ** b for (a, b) in monos]
    M = flint.fmpz_mat([[row[t] * scale[t] for t in range(len(monos))] for row in rows])
    L = M.lll()
    return len(monos), monos, L


def good_vectors(L, n: int, D: int, m: int) -> int:
    """Number of reduced rows v with n * ||v||^2 < D^{2m} (Howgrave-Graham, exact integers)."""
    cnt = 0
    bound = D ** (2 * m)
    for i in range(n):
        s = sum(int(L[i, t]) ** 2 for t in range(n))
        if n * s < bound:
            cnt += 1
    return cnt


def log2_det_condition(j: int, m: int, X: int, Y: int, D: int) -> float:
    """log2 of det^{1/n} * 2^{(n-1)/4} * sqrt(n) / D^m (< 0 means the first vector is guaranteed);
    floats for reporting only."""
    n = sum(j * (m - k) + 1 for k in range(m + 1))
    ld = 0.0
    for k in range(m + 1):
        for u in range(j * (m - k) + 1):
            ld += u * math.log2(X) + k * math.log2(Y) + (m - k) * math.log2(D)
    return ld / n + (n - 1) / 4 + 0.5 * math.log2(n) - m * math.log2(D)
