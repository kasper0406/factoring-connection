"""Exact toy checks of the elementary lemmas of notes/mfU_theorem.md (correctness checks only; no asymptotic content).

- base_expansion: Definition 4.1 / Lemma 4.2 with an arbitrary positive integer multiplier w (not necessarily prime).
- forced_residue_ok: Lemma 4.5, F_hat_{g,w}(P) = c_P w (mod |G(P)|) for every integer w, when gcd(|G(P)|, N g1 b) = 1.
- residue_law_deviation: max_r |Pr_w[w = r mod t] - 1/t| for w uniform on the integers of (W, 2W] (Lemma OP uses <= 1/W).
- norm_class_ok: Lemma 8.5's key identity, Norm(delta xi^2) = delta (mod Q^{x2}) in odd degree.
No factor of any test N is used by these functions except as an explicit input of the identity being checked.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import List

import sympy


def base_expansion(N: int, w: int, g0: int, g1: int, d: int) -> List[int]:
    """Coefficients a_0..a_d of F_hat_{g,w}(x, y) = sum a_i x^i y^{d-i} with F_hat(g0, g1) = N w."""
    assert gcd(g0, g1) == 1 and g1 >= 2
    T = N * w
    a = [0] * (d + 1)
    inv = pow(g0, -1, g1)
    for i in range(d, 0, -1):
        r = (T * pow(inv, i, g1)) % g1
        if r > g1 // 2:
            r -= g1
        a[i] = r
        T = (T - r * g0 ** i)
        assert T % g1 == 0
        T //= g1
    a[0] = T
    return a


def form_value(coeffs: List[int], x: int, y: int) -> int:
    d = len(coeffs) - 1
    return sum(c * x ** i * y ** (d - i) for i, c in enumerate(coeffs))


def forced_residue_ok(N: int, w: int, g0: int, g1: int, d: int, P) -> bool:
    a, b = P
    t = abs(g1 * a - g0 * b)
    if t == 0 or gcd(t, N * g1 * b) != 1:
        raise ValueError("hypothesis of Lemma 4.5 fails")
    F = base_expansion(N, w, g0, g1, d)
    cP = (N * pow(b, d, t) * pow(pow(g1, d, t), -1, t)) % t
    return (form_value(F, a, b) - cP * w) % t == 0


def residue_law_deviation(W: int, t: int) -> Fraction:
    counts = [0] * t
    for w in range(W + 1, 2 * W + 1):
        counts[w % t] += 1
    return max(abs(Fraction(c, W) - Fraction(1, t)) for c in counts)


def squarefree_kernel(q: Fraction) -> int:
    """The squarefree integer delta with q/delta a rational square (sign included)."""
    sgn = -1 if q < 0 else 1
    out = sgn
    for n, e in ((abs(q.numerator), 1), (q.denominator, 1)):
        for p, k in sympy.factorint(n).items():
            if k % 2:
                out *= p
    return out


def norm_class_ok(f_coeffs: List[int], xi_coeffs: List[int], delta: int) -> bool:
    """f monic irreducible of odd degree d (coefficients low->high), xi in Z[alpha]; checks that the squarefree kernel of
    Norm(delta * xi^2) equals the squarefree kernel of delta."""
    x = sympy.Symbol("x")
    f = sum(c * x ** i for i, c in enumerate(f_coeffs))
    d = sympy.degree(f, x)
    assert d % 2 == 1
    xi = sum(c * x ** i for i, c in enumerate(xi_coeffs))
    gamma = sympy.rem(sympy.expand(delta * xi ** 2), f, x)
    lc = sympy.LC(sympy.Poly(f, x))
    nrm = Fraction(int(sympy.resultant(f, gamma, x)), int(lc) ** int(sympy.degree(gamma, x)) if gamma.has(x) else 1)
    return squarefree_kernel(nrm) == squarefree_kernel(Fraction(delta))
