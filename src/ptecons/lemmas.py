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


def cunningham_bound(n: int, roots: Sequence[int]) -> Tuple[int, int]:
    """Lemma CUN (corrected after codex ptecons_2): with r = P mod (X^n - 1) and z = number of nonreal
    n-th roots of unity (n - 2 for even n, n - 1 for odd n),
        prod_{|a|>=2} |a|^{n-2} * (1 - 2^{-n})^M  <=  ||r||_1^z.
    Exact integer form: returns (RHS, LHS) = (||r||_1^z * 2^{nM}, prod|a|^{n-2} * (2^n - 1)^M);
    the lemma says RHS >= LHS for n >= 3."""
    F = [-1] + [0] * (n - 1) + [1]
    r = poly_mod(from_roots(roots), F)
    z = n - 2 if n % 2 == 0 else n - 1
    M = len(roots)
    big = 1
    for a in roots:
        if abs(a) >= 2:
            big *= abs(a)
    return norm1(r) ** z * 2 ** (n * M), big ** (n - 2) * (2 ** n - 1) ** M


def galois_orbit_collapse(n: int, x: int, t: int) -> bool:
    """prod_k (x - zeta_n^k t) == x^n - t^n (so a full mu_n-orbit gives the single relation x^n - t^n).
    Exact: the product over all roots of Z^n - 1 (monic) of (x - t Z) is Res(Z^n - 1, x - tZ)."""
    return resultant([-1] + [0] * (n - 1) + [1], [x, -t]) == x ** n - t ** n


def root_support_rank(family: Sequence[Sequence[int]]) -> Tuple[int, int]:
    """(#members, #distinct roots used).  Over a SYMBOL factor base {m - a}, rank <= #distinct roots,
    so such a family gives no excess unless #members > #roots.  (If pieces are factored into
    primes, columns merge and this count is only an upper bound on columns -- codex ptecons_2.)"""
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


def star_gaussian_check(cs: Sequence[complex]) -> Tuple[int, int, int]:
    """Lemma STAR over Z[i]: for a conjugation-closed multiset of Gaussian integers (given as complex
    ints), return (M', s, bound) where M' = #nonzero, s = degree of the non-leading part of
    prod (Y - c), and bound = 6(s+1) + 2 log2 |g_low| (ceil).  STAR says M' <= bound."""
    import math
    # exact Gaussian-integer arithmetic on (re, im) pairs
    cs = [(int(c.real), int(c.imag)) for c in map(complex, cs)]
    mul = lambda a, b: (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    coeffs = [(1, 0)]
    for c in cs:
        nxt = [(0, 0)] * (len(coeffs) + 1)
        for i, a in enumerate(coeffs):
            nxt[i + 1] = (nxt[i + 1][0] + a[0], nxt[i + 1][1] + a[1])
            ca = mul(c, a)
            nxt[i] = (nxt[i][0] - ca[0], nxt[i][1] - ca[1])
        coeffs = nxt
    M = len(cs)
    g = coeffs[:M]   # Y^0..Y^{M-1}
    nz = [j for j, z in enumerate(g) if z != (0, 0)]
    Mp = sum(1 for c in cs if c != (0, 0))
    if not nz:
        return Mp, -1, 0
    s = max(nz)
    glow2 = g[min(nz)][0] ** 2 + g[min(nz)][1] ** 2       # |g_low|^2, exact
    return Mp, s, 6 * (s + 1) + (glow2.bit_length() if glow2 > 1 else 0)
