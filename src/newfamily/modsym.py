"""Plus modular symbols [x]^+ = lambda({oo, x}) for the weight-2 newform of prime level M (default 11, curve 11a).
OPERATIONAL: sees only the rational x = a/N (Manin's continued-fraction trick, O(log N) steps).  Exact (Fractions).

Manin symbols (c:d) in P^1(F_M) <-> g{0, oo} with g in SL2(Z) of bottom row (c, d).  Relations:
  x + xS = 0 (S: (c:d) -> (d:-c)),  x + x tau + x tau^2 = 0 (tau: (c:d) -> (d:-c-d)),  plus quotient x = x eta
  (eta: (c:d) -> (-c:d)).
Hecke T_l (l prime, l != M) acts by Merel's Heilbronn matrices {[[a,b],[c,d]]: ad-bc=l, a>b>=0, d>c>=0} on the right:
  (u:v)[[a,b],[c,d]] = (ua+vc : ub+vd).
lambda is the functional on the plus quotient with lambda o T_2 = a_2 lambda (a_2 = -2 for 11a), normalised to make
[0]^+ = lambda({oo, 0}) a fixed rational.  {oo, a/N} = sum_{k=0}^{r} g_k{0,oo} with
g_k = [[(-1)^{k-1}p_k, p_{k-1}], [(-1)^{k-1}q_k, q_{k-1}]] (convergents p_k/q_k, p_{-1}/q_{-1} = 1/0).
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import List, Tuple

import sympy


def _idx(c: int, d: int, M: int) -> int:
    c %= M
    d %= M
    if d:
        return (c * pow(d, -1, M)) % M
    if c == 0:
        raise ValueError("(0:0) is not in P^1")
    return M


def _rep(i: int, M: int) -> Tuple[int, int]:
    return (1, 0) if i == M else (i, 1)


def heilbronn(l: int) -> List[Tuple[int, int, int, int]]:
    out = []
    for a in range(1, l + 1):
        for d in range(1, l + 1):
            for b in range(0, a):
                for c in range(0, d):
                    if a * d - b * c == l:
                        out.append((a, b, c, d))
    return out


@lru_cache(maxsize=None)
def setup(M: int = 11, a2: int = -2):
    n = M + 1
    rel = []
    for i in range(n):
        c, d = _rep(i, M)
        r = [0] * n
        r[i] += 1
        r[_idx(d, -c, M)] += 1
        rel.append(r)
        r = [0] * n
        r[i] += 1
        r[_idx(d, -c - d, M)] += 1
        r[_idx(-c - d, c, M)] += 1  # tau^2: (c:d) -> (d:-c-d) -> (-c-d : -d+c+d) = (-c-d : c)
        rel.append(r)
        r = [0] * n
        r[i] += 1
        r[_idx(-c, d, M)] -= 1
        rel.append(r)
    R = sympy.Matrix(rel)
    dual = R.nullspace()          # functionals lambda (column vectors) with R lambda = 0
    B = sympy.Matrix.hstack(*dual)  # n x k
    # T_2 on symbols: matrix T (n x n), T[i, j] = coefficient of symbol j in T(i)
    T = sympy.zeros(n, n)
    for i in range(n):
        u, v = _rep(i, M)
        for (a, b, c, d) in heilbronn(2):
            uu, vv = u * a + v * c, u * b + v * d
            if uu % M == 0 and vv % M == 0:
                continue
            T[i, _idx(uu, vv, M)] += 1
    # lambda in span(B) with T lambda = a2 lambda  (lambda(T x) as column action)
    TB = T * B
    # solve TB y = a2 B y, i.e. (TB - a2 B) y = 0
    ns = (TB - a2 * B).nullspace()
    if len(ns) != 1:
        raise RuntimeError(f"eigenspace dim {len(ns)}")
    lam = B * ns[0]
    eigs_all = sorted(set((TB - x * B).nullspace().__len__() and x for x in range(-10, 11)) - {0})
    vals = [Fraction(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1])) for x in lam]
    # normalise: make the value of {oo, 0} = lambda(g_0{0,oo}) with g_0 = [[-0? ...]] -- use x = 0 path below
    tbl = tuple(vals)
    v0 = _value(tbl, M, 0, 1)
    if v0 != 0:
        tbl = tuple(x / v0 for x in tbl)
    return tbl, B.shape[1], eigs_all


def _value(tbl, M: int, a: int, N: int) -> Fraction:
    """lambda({oo, a/N}) for gcd(a, N) = 1, N >= 1."""
    a %= N if N > 1 else 1
    # convergents of a/N = [a0; a1, ..., ar]
    p2, q2 = 0, 1   # p_{-2}/q_{-2}
    p1, q1 = 1, 0   # p_{-1}/q_{-1}
    x, y = a, N
    total = Fraction(0)
    k = 0
    while True:
        t = x // y
        p, q = t * p1 + p2, t * q1 + q2
        sgn = -1 if (k - 1) % 2 else 1   # (-1)^(k-1)
        total += tbl[_idx(sgn * q, q1, M)]
        p2, q2, p1, q1 = p1, q1, p, q
        k += 1
        x, y = y, x - t * y
        if y == 0:
            break
    return total


def modsym_plus(a: int, N: int, M: int = 11) -> Fraction:
    tbl, _, _ = setup(M)
    return _value(tbl, M, a, N)
