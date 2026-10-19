"""Exact arithmetic in Z[beta] = Z[y]/(g), beta = f_d * alpha, for a NON-MONIC f (notes/lv6_results.md §2).

Conventions (as in src/pcest/geometry.py):
- f = [f_0, ..., f_d] (f[i] is the coefficient of x^i), F(a, b) = sum_i f_i a^i b^(d-i), alpha = root of f,
  so a - b*alpha has norm F(a, b)/f_d.
- g(y) = f_d^(d-1) f(y / f_d) = sum_i f_i f_d^(d-1-i) y^i is monic, beta = f_d alpha is a root of g.
- An element of Z[beta] is a list of d integers (coefficient of beta^j at index j).

Everything is exact (python ints). Operational: these functions take only public data (f, N, field).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Sequence, Tuple

import sympy


def monic_g(f: Sequence[int]) -> List[int]:
    d = len(f) - 1
    fd = f[d]
    return [int(f[i]) * fd ** (d - 1 - i) for i in range(d)] + [1]


def reduce_mod_g(c: List[int], g: Sequence[int], m: int | None = None) -> List[int]:
    """Reduce a coefficient list (any length) modulo the monic g (degree d); optional coefficient modulus m."""
    d = len(g) - 1
    c = list(c)
    for i in range(len(c) - 1, d - 1, -1):
        t = c[i]
        if t:
            for j in range(d):
                c[i - d + j] -= t * g[j]
        c[i] = 0
    c = (c + [0] * d)[:d]
    if m is not None:
        c = [x % m for x in c]
    return c


def mul(x: Sequence[int], y: Sequence[int], g: Sequence[int], m: int | None = None) -> List[int]:
    d = len(g) - 1
    prod = [0] * (2 * d - 1)
    for i, xi in enumerate(x):
        if xi:
            for j, yj in enumerate(y):
                prod[i + j] += xi * yj
    if m is not None:
        prod = [p % m for p in prod]
    return reduce_mod_g(prod, g, m)


def one(d: int) -> List[int]:
    return [1] + [0] * (d - 1)


def linear(a: int, b: int, fd: int, d: int) -> List[int]:
    """The element f_d a - b beta = f_d (a - b alpha) of Z[beta]."""
    e = [0] * d
    e[0] = fd * a
    e[1] = -b
    return e


def product(elems: Sequence[Sequence[int]], g: Sequence[int]) -> List[int]:
    """Exact product by a balanced product tree."""
    d = len(g) - 1
    level = [list(e) for e in elems] or [one(d)]
    while len(level) > 1:
        nxt = []
        for i in range(0, len(level) - 1, 2):
            nxt.append(mul(level[i], level[i + 1], g))
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    return level[0]


def gprime_elem(g: Sequence[int]) -> List[int]:
    """g'(beta) as an element of Z[beta]."""
    d = len(g) - 1
    return [(i + 1) * g[i + 1] for i in range(d)]


def disc_poly(coeffs_low_to_high: Sequence[int]) -> int:
    x = sympy.Symbol("x")
    P = sympy.Poly(list(reversed([int(c) for c in coeffs_low_to_high])), x)
    return int(sympy.discriminant(P))


def psi(elem: Sequence[int], beta_image: int, N: int) -> int:
    """Ring homomorphism Z[beta] -> Z/N, beta -> beta_image (requires g(beta_image) = 0 mod N)."""
    r = 0
    for c in reversed(elem):
        r = (r * beta_image + c) % N
    return r


def norm1(x: Sequence[int]) -> int:
    return sum(abs(c) for c in x)


def coeff_bound_sqrt(g: Sequence[int], factors_l1: Sequence[int]) -> int:
    """Rigorous integer bound A on |eta_j| for eta in Z[beta] with eta^2 = g'(beta)^2 * prod(x_t)
    (Lemma SQ(ii)); factors_l1 = [||x_t||_1].  All quantities are integers (no floating point).

    Mh >= Mahler(g) >= max_k max(1, |beta_k|);
    |sigma_i(x)| <= ||x||_1 Mh^(d-1);  |g'(beta_i)| <= d ||g||_1 Mh^(d-1);
    |sigma_i(eta)|^2 <= U := (d||g||_1 Mh^(d-1))^2 prod_t ||x_t||_1 Mh^(d-1);
    |(V^-1)_{j,i}| <= 2^(d-1) Mh / |g'(beta_i)|  and  |g'(beta_i)| >= (d ||g||_1 Mh^(d-1))^-(d-1)
    (because prod_i |g'(beta_i)| = |disc g| >= 1).
    """
    d = len(g) - 1
    Mh = isqrt(sum(c * c for c in g)) + 1
    gl1 = norm1(g)
    gp = d * gl1 * Mh ** (d - 1)
    U = gp * gp
    for n1 in factors_l1:
        U *= n1 * Mh ** (d - 1)
    sig = isqrt(U) + 1
    return d * sig * 2 ** (d - 1) * Mh * gp ** (d - 1)


def class_group_log2_bound(f: Sequence[int]) -> Tuple[int, Dict[str, int]]:
    """Explicit upper bound for log2 |H| (Lemma H of the note):
    |H| <= 2^(d*omega(D)) * 2^d * h_K,  D = disc(f) * f_d,
    h_K <= M_K (d - 1 + log M_K)^(d-1) / (d-1)!,  M_K <= sqrt|Delta_K| <= sqrt|disc g| (Lenstra 1992, Thm 6.5).
    omega(D) is bounded by log2|D| (no factoring needed).  Returns (bound, parts)."""
    d = len(f) - 1
    dg = abs(disc_poly(monic_g(f)))
    D = abs(disc_poly(f) * f[d])
    omega_b = D.bit_length()
    log2_MK = (dg.bit_length() + 1) // 2 + 1
    # (d-1+log M_K)^(d-1)/(d-1)! <= (d + log2_MK)^(d-1)
    log2_h = log2_MK + (d - 1) * (d + log2_MK).bit_length()
    total = d * omega_b + d + log2_h
    return total, {"d_omega": d * omega_b, "units": d, "log2_h": log2_h, "log2_disc_g": dg.bit_length()}
