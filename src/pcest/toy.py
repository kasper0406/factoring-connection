"""Toy special-q counts for designs 1D (g = x - m) and 2D (g = g1 x - g0), experiment E1.

OPERATIONAL part: fields and counts use N and public parameters only.
- 1D field: m uniform in [M/2, M], fhat = centred base-m expansion, f = fhat + (x - m) R.
- 2D field: (g0, g1) uniform coprime in [M/2, M]^2, fhat = gbe(N, d, g0, g1), f = fhat + (g1 x - g0) R.
R has d i.i.d. coefficients uniform in [-H, H].  For a prime l the special-q lattice is
{(a, b) : g1 a = g0 b (mod l)} = {a = b m_l (mod l)}, m_l = g0/g1 mod l, and the box is
hypD's reduced-basis box (one sign per +-pair, primitive, b != 0).
N_rat = #{box points : G/l B-smooth};  N_both = #{... and F(a, b) B-smooth}.
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Tuple

import gmpy2

from ..hypD.specialq import box_points
from ..nfsfruit.smooth import primes_upto
from .geometry import field, hom_eval, G_value


class Smooth:
    """Exact B-smoothness by repeated gcd with the primorial (no factorisation of N)."""

    def __init__(self, B: int):
        P = 1
        for p in primes_upto(B):
            P *= p
        self.P = gmpy2.mpz(P)

    def __call__(self, z: int) -> bool:
        r = gmpy2.mpz(abs(int(z)))
        if r == 0:
            return False
        g = gmpy2.gcd(r, self.P)
        while g > 1:
            r //= g
            g = gmpy2.gcd(r, g)
        return r == 1


def base_m_centred(N: int, d: int, m: int) -> List[int]:
    a, T = [0] * (d + 1), N
    for i in range(d):
        r = T % m
        if 2 * r > m:
            r -= m
        a[i] = r
        T = (T - r) // m
    a[d] = T
    assert hom_eval(a, m, 1) == N
    return a


def field_1d(N: int, d: int, m: int, R) -> List[int]:
    f = base_m_centred(N, d, m)
    for j, r in enumerate(R):
        f[j + 1] += r
        f[j] -= m * r
    assert hom_eval(f, m, 1) == N
    return f


def sample_field(N: int, d: int, M: int, H: int, design: str, rng: random.Random):
    R = [rng.randrange(-H, H + 1) for _ in range(d)]
    if design == "1d":
        m = rng.randrange(M // 2, M + 1)
        return (m, 1), field_1d(N, d, m, R)
    while True:
        g0, g1 = rng.randrange(M // 2, M + 1), rng.randrange(M // 2, M + 1)
        if gcd(g0, g1) == 1:
            return (g0, g1), field(N, d, g0, g1, R)


def counts(N: int, d: int, g: Tuple[int, int], f: List[int], l: int, S: int,
           sm: Smooth) -> Dict:
    g0, g1 = g
    if g1 % l == 0:
        return dict(nrat=0, nboth=0, npts=0, viol=0)
    m_l = g0 * pow(g1, -1, l) % l
    nrat = nboth = npts = viol = 0
    for a, b in box_points(l, m_l, S):
        npts += 1
        G = G_value(g0, g1, a, b)
        if G == 0 or G % l:
            viol += 1
            continue
        if not sm(G // l):
            continue
        nrat += 1
        F = hom_eval(f, a, b)
        if (F * g1 ** d - N * b ** d) % G:
            viol += 1
        if F != 0 and sm(F):
            nboth += 1
    return dict(nrat=nrat, nboth=nboth, npts=npts, viol=viol)
