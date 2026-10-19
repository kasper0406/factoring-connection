"""Exact lattice facts behind notes/pcest_results.md (Lemma GB, Lemma Q, Lemma Kern).

OPERATIONAL (receives N and public parameters only):
- gbe(N, d, g0, g1): generalized base-(g0/g1) expansion, F(g0, g1) = N exactly.
- field(N, d, g0, g1, R): f = fhat + (g1 x - g0) R(x); returns coefficients f[0..d] (f[i] of x^i).
Homogeneous evaluation: F(a, b) = sum_i f[i] a^i b^(d-i); rational value G(a, b) = g1 a - g0 b.

DIAGNOSTIC (pure lattice algebra, no N):
- rational_pair_lattice: image of the special-q lattice L_l(P) under g -> (G(P)/l, G(P')/l).
- alg_pair_index / kernel checks for c -> (<c, v(P)>, <c, v(P')>).
"""
from __future__ import annotations

from itertools import combinations
from math import gcd
from typing import List, Optional, Sequence, Tuple

import sympy
from sympy.matrices.normalforms import hermite_normal_form


# ------------------------------------------------------------------ Lemma GB
def _centered(x: int, m: int) -> int:
    r = x % m
    return r - m if 2 * r > m else r


def gbe(N: int, d: int, g0: int, g1: int) -> List[int]:
    """Coefficients a[0..d] with sum_i a[i] g0^i g1^(d-i) = N, |a[i]| <= g1/2 for i >= 1.

    Step: a_d = N g0^(-d) mod g1 (centred), then recurse on (N - a_d g0^d)/g1 with d - 1.
    Requires gcd(g0, g1) = 1 and g1 >= 2.
    """
    if gcd(g0, g1) != 1 or g1 < 2:
        raise ValueError("need gcd(g0, g1) = 1 and g1 >= 2")
    a = [0] * (d + 1)
    T = N
    for i in range(d, 0, -1):
        inv = pow(g0, -i, g1)
        a[i] = _centered(T * inv, g1)
        T -= a[i] * g0 ** i
        assert T % g1 == 0
        T //= g1
    a[0] = T
    assert hom_eval(a, g0, g1) == N
    return a


def hom_eval(f: Sequence[int], a: int, b: int) -> int:
    d = len(f) - 1
    return sum(int(f[i]) * a ** i * b ** (d - i) for i in range(d + 1))


def field(N: int, d: int, g0: int, g1: int, R: Sequence[int]) -> List[int]:
    """f = fhat + (g1 x - g0) R(x), R = sum_j R[j] x^j (deg <= d-1). F(g0, g1) = N."""
    assert len(R) == d
    f = gbe(N, d, g0, g1)
    for j, r in enumerate(R):
        f[j + 1] += g1 * r
        f[j] -= g0 * r
    assert hom_eval(f, g0, g1) == N
    return f


def G_value(g0: int, g1: int, a: int, b: int) -> int:
    return g1 * a - g0 * b


# ------------------------------------------------------------------ lattices
def hnf_rows(rows: Sequence[Sequence[int]]) -> sympy.Matrix:
    """Row-style HNF basis of the Z-span of the given rows (returned as rows)."""
    M = sympy.Matrix(rows).T                  # columns = generators
    H = hermite_normal_form(M)                # column-style HNF of the column span
    return H.T


def lattice_index_2d(gens: Sequence[Tuple[int, int]]) -> int:
    """Index in Z^2 of the lattice spanned by gens (assumed full rank)."""
    H = hnf_rows(gens)
    assert H.shape[0] == 2
    return abs(int(H.det()))


def special_q_basis(l: int, a: int, b: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Basis of L_l(P) = {(g0, g1) : g1 a - g0 b = 0 (mod l)}, assuming b != 0 mod l."""
    t = a * pow(b, -1, l) % l
    return (l, 0), (t, 1)


def rational_pair_lattice(l: int, P: Tuple[int, int], Pp: Tuple[int, int]):
    """Return (gens, q, lam): image of L_l(P) under g -> (G(P)/l, G(P')/l).

    Requires l | Delta = a b' - a' b, and b, b' invertible mod l.  Lemma Q: the image is
    Gamma_q = {(n, n') : n' = lam n (mod q)} with q = |Delta|/l and lam a unit mod q
    (for primitive P, P').  Returns lam = None when q = 1.
    """
    (a, b), (ap, bp) = P, Pp
    Delta = a * bp - ap * b
    assert Delta % l == 0 and Delta != 0
    gens = []
    for g0, g1 in special_q_basis(l, a, b):
        n, n2 = G_value(g0, g1, a, b), G_value(g0, g1, ap, bp)
        assert n % l == 0 and n2 % l == 0
        gens.append((n // l, n2 // l))
    q = abs(Delta) // l
    idx = lattice_index_2d(gens)
    assert idx == q, (idx, q)
    if q == 1:
        return gens, 1, None
    lam = None
    for n, n2 in gens + [(gens[0][0] + gens[1][0], gens[0][1] + gens[1][1])]:
        if gcd(n, q) == 1:
            lam = n2 * pow(n, -1, q) % q
            break
    return gens, q, lam


def in_gamma_q(n: int, n2: int, q: int, lam: Optional[int]) -> bool:
    return q == 1 or (n2 - lam * n) % q == 0


def veronese(P: Tuple[int, int], e: int) -> List[int]:
    """v_e(P) = (a^(e-i) b^i)_{0 <= i <= e} (length e + 1)."""
    a, b = P
    return [a ** (e - i) * b ** i for i in range(e + 1)]


def alg_pair_index(P, Pp, d: int) -> int:
    """gcd of the 2x2 minors of [v_{d-1}(P); v_{d-1}(P')] = index of Gamma in Z^2."""
    v, w = veronese(P, d - 1), veronese(Pp, d - 1)
    g = 0
    for i, j in combinations(range(d), 2):
        g = gcd(g, v[i] * w[j] - v[j] * w[i])
    return g


def kernel_basis_Q(P, Pp, d: int) -> List[List[int]]:
    """Q_PP'(x, y) * x^i y^(d-3-i), i = 0..d-3, as coefficient vectors c in the basis of
    v_{d-1} (c_i multiplies a^(d-1-i) b^i), where Q = (b x - a y)(b' x - a' y)."""
    (a, b), (ap, bp) = P, Pp
    # Q(x, y) = b b' x^2 - (a b' + a' b) x y + a a' y^2, coefficients in x^(2-j) y^j
    Qc = [b * bp, -(a * bp + ap * b), a * ap]
    out = []
    for i in range(d - 2):
        c = [0] * d                                # c_j multiplies x^(d-1-j) y^j
        for j, qj in enumerate(Qc):
            c[j + i] += qj
        out.append(c)
    return out


def is_saturated(rows: Sequence[Sequence[int]]) -> bool:
    """True iff the Z-span of rows (r x d, rank r) is saturated in Z^d (gcd of r-minors = 1)."""
    M = sympy.Matrix(rows)
    r, d = M.shape
    g = 0
    for cols in combinations(range(d), r):
        g = gcd(g, int(M[:, list(cols)].det()))
        if g == 1:
            return True
    return g == 1


def fibre_ratio_bound(P, Pp, d: int, H: int) -> float:
    """Diagnostic: prod_i max(1, |k_i|_inf / H) over the kernel basis (crude successive-minima
    proxy). Lemma Kern: fibre counts are uniform within O_d(1) once this is O(1)."""
    prod = 1.0
    for c in kernel_basis_Q(P, Pp, d):
        prod *= max(1.0, max(abs(x) for x in c) / H)
    return prod
