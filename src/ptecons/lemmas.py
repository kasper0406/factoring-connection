"""Diagnostic checkers for the ptecons construction-side lemmas (notes/ptecons_results.md section 3).

None of these is operational; none sees factors of N.

- gap_r_exhaustive: Lemma GAP-R (real multisets with e_i = e_{i+1} = 0 have <= i-1 nonzero entries).
- lin_family_affine_degree: Lemma LIN-R on explicit one/two-parameter families.
- ledger_holds: Lemma LEDGER, |Res(F, P mod F)| = prod |F(a_l)| for monic F.
- cunningham_bound: Lemma CUN, ||P mod (X^n - 1)||_1 >= prod_{|a|>=2} |a| * (1 - eps).
- galois_orbit_collapse: product over a full mu_n orbit of (x - zeta^k t) equals x^n - t^n (one relation).
- root_support_rank: rows <= columns for families whose members have disjoint root supports.
- dickson_integer_split: classify t with D_n(X) - D_n(t) split over Z.
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Sequence, Tuple

import sympy as sp

from .polys import (all_elem_sym, dickson, evaluate, from_roots, gap_bound_nonzero, int_roots_split,
                    norm1, poly_mod, resultant)


def gap_r_exhaustive(M: int, H: int) -> Tuple[int, int]:
    """All integer multisets of size M in [-H, H]; return (#with a double gap, #violations)."""
    tot = viol = 0
    for y in itertools.combinations_with_replacement(range(-H, H + 1), M):
        e = all_elem_sym(y)
        g = gap_bound_nonzero(e[:M + 1])
        if g is not None:
            tot += 1
            if sum(1 for v in y if v) > g:
                viol += 1
    return tot, viol


def lin_family_affine_degree(alphas: Sequence[int], betas: Sequence[Sequence[int]]) -> List[int]:
    """Family P_lam = prod (X - alpha_l - <beta_l, lam>); return, for each coefficient of X^{M-i}
    (i = 1..M), its total degree in lam.  LIN-R predicts deg(e_2) >= 2 or deg(e_3) >= 3 unless
    at most one beta_l is nonzero."""
    s = len(betas[0])
    lam = sp.symbols(f"l0:{s}")
    X = sp.Symbol("X")
    roots = [a + sum(b[j] * lam[j] for j in range(s)) for a, b in zip(alphas, betas)]
    P = sp.Poly(sp.expand(sp.prod([X - r for r in roots])), X)
    M = len(alphas)
    out = []
    for i in range(1, M + 1):
        c = P.coeff_monomial(X ** (M - i))
        out.append(sp.Poly(c, *lam).total_degree() if c != 0 else 0)
    return out


def ledger_holds(F: Sequence[int], roots: Sequence[int]) -> bool:
    P = from_roots(roots)
    r = poly_mod(P, F)
    prod = 1
    for a in roots:
        prod *= evaluate(F, a)
    return abs(resultant(F, r)) == abs(prod)


def cunningham_bound(n: int, roots: Sequence[int]) -> Tuple[int, float]:
    """Return (||P mod X^n - 1||_1, prod_{|a|>=2}|a| * (1 - 2^{-n})^{M/(n-2)}) -- Lemma CUN says
    first >= second whenever n >= 3 (product over the nonreal n-th roots of unity)."""
    F = [-1] + [0] * (n - 1) + [1]
    r = poly_mod(from_roots(roots), F)
    big = 1
    for a in roots:
        if abs(a) >= 2:
            big *= abs(a)
    return norm1(r), big * (1 - 2.0 ** (-n)) ** (len(roots) / (n - 2))


def galois_orbit_collapse(n: int, x: int, t: int) -> bool:
    """prod_k (x - zeta_n^k t) == x^n - t^n (so a full mu_n-orbit gives the single relation x^n - t^n).
    Exact: the product over all roots of Z^n - 1 (monic) of (x - t Z) is Res(Z^n - 1, x - tZ)."""
    return resultant([-1] + [0] * (n - 1) + [1], [x, -t]) == x ** n - t ** n


def root_support_rank(family: Sequence[Sequence[int]]) -> Tuple[int, int]:
    """(#members, #distinct roots used).  A relation matrix over the symbols {m - a} has rank
    <= #distinct roots, so families need #members > #roots to be useful at all."""
    support = set()
    for f in family:
        support.update(f)
    return len(family), len(support)


def dickson_integer_split(n: int, T: int) -> List[int]:
    D = dickson(n)
    good = []
    for t in range(-T, T + 1):
        c = list(D)
        c[0] -= evaluate(D, t)
        if int_roots_split(c) is not None:
            good.append(t)
    return good
