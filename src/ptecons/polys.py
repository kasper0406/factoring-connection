"""Exact integer-polynomial helpers for ptecons.  Coefficient lists are LOW-DEGREE FIRST.

Diagnostic/harness-neutral: nothing here sees factors of N.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Iterable, List, Sequence

import flint


def from_roots(roots: Iterable[int]) -> List[int]:
    """prod (X - a) as a coefficient list (low degree first)."""
    c = [1]
    for a in roots:
        nxt = [0] * (len(c) + 1)
        for i, ci in enumerate(c):
            nxt[i + 1] += ci
            nxt[i] -= a * ci
        c = nxt
    return c


def elem_sym(vals: Sequence, k: int):
    """e_k of a sequence (works over any ring supporting + and *)."""
    e = [1] + [0] * k
    for v in vals:
        for j in range(k, 0, -1):
            e[j] = e[j] + e[j - 1] * v
    return e[k]


def all_elem_sym(vals: Sequence) -> list:
    n = len(vals)
    e = [1] + [0] * n
    for v in vals:
        for j in range(n, 0, -1):
            e[j] = e[j] + e[j - 1] * v
    return e


def evaluate(c: Sequence[int], x: int) -> int:
    r = 0
    for ci in reversed(c):
        r = r * x + ci
    return r


def poly_mod(c: Sequence[int], F: Sequence[int]) -> List[int]:
    """Remainder of c modulo a MONIC F over Z (exact)."""
    assert F[-1] == 1
    c = list(c)
    d = len(F) - 1
    for top in range(len(c) - 1, d - 1, -1):
        q = c[top]
        if q:
            for i in range(d + 1):
                c[top - d + i] -= q * F[i]
    out = c[:d] if len(c) >= d else c + [0] * (d - len(c))
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def resultant(F: Sequence[int], G: Sequence[int]) -> int:
    return int(flint.fmpz_poly(list(F)).resultant(flint.fmpz_poly(list(G))))


def norm1(c: Sequence[int]) -> int:
    return sum(abs(x) for x in c)


def centered(x: int, N: int) -> int:
    r = x % N
    return r - N if r > N // 2 else r


def int_roots_split(c: Sequence[int]) -> List[int] | None:
    """If the integer polynomial c splits into linear factors over Z (monic), return the roots
    with multiplicity; else None.  Uses flint factorisation of a POLYNOMIAL (not of N)."""
    p = flint.fmpz_poly(list(c))
    if p.degree() < 1:
        return [] if p.degree() == 0 else None
    cont, facs = p.factor()
    roots: List[int] = []
    for f, e in facs:
        if f.degree() != 1:
            return None
        a, b = int(f[0]), int(f[1])
        if b not in (1, -1) or (-a * b) % 1:
            return None
        roots += [-a * b] * e
    return sorted(roots)


def dickson(n: int, a: int = 1) -> List[int]:
    """D_n(X, a): D_0 = 2, D_1 = X, D_n = X D_{n-1} - a D_{n-2}."""
    D0, D1 = [2], [0, 1]
    if n == 0:
        return D0
    for _ in range(n - 1):
        nxt = [0] + D1
        for i, c in enumerate(D0):
            nxt[i] -= a * c
        D0, D1 = D1, nxt
    return D1


def gap_bound_nonzero(e_vals: Sequence[int]) -> int | None:
    """Lemma GAP-R helper: if e_i = e_{i+1} = 0 for some 1 <= i < len-1 (first such i),
    return i-1 (the max number of nonzero REAL roots allowed); else None."""
    for i in range(1, len(e_vals) - 1):
        if e_vals[i] == 0 and e_vals[i + 1] == 0:
            return i - 1
    return None
