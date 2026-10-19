"""Toy LV-style special-q counts N_l(A) (experiment E2 of notes/hypD_prereg.md).

OPERATIONAL: field A = (m, r0, r1): monic base-m cubic plus (x - m)(r1 x + r0), so that
    F(a, b) = Fhat(a, b) + b (a - b m)(r1 a + r0 b),   F(a, b) = N b^3 (mod a - b m).
For a prime l, the special-q lattice is {(a, b) : a = b m (mod l)}; its Gauss-reduced
basis (w1, w2) spans the box {i w1 + j w2 : |i| <= S, 0 <= j <= S, (i, j) != 0}.
A box point is a *special-q pair for l* iff gcd(a, b) = 1, b != 0, (a - b m)/l is
B-smooth and F(a, b) is B-smooth.  N_l(A) = number of such pairs (one sign per +-pair).
Smoothness is tested exactly (FactorBase); nothing here touches the factors of N.
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Tuple

import gmpy2

from ..nfsfruit.smooth import FactorBase


def base_m_digits(N: int, m: int) -> List[int]:
    c, t = [], N
    for _ in range(3):
        t, r = divmod(t, m)
        c.append(r)
    assert t == 1, "need m^3 <= N < 2 m^3"
    return c                       # f_hat = x^3 + c2 x^2 + c1 x + c0


def random_m(N: int, rng: random.Random) -> int:
    lo = int(gmpy2.iroot(gmpy2.mpz(N) // 2, 3)[0]) + 1
    hi = int(gmpy2.iroot(gmpy2.mpz(N), 3)[0])
    while True:
        m = rng.randrange(lo, hi + 1)
        if m ** 3 <= N < 2 * m ** 3:
            return m


def reduced_basis(l: int, m: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    a, b = (l, 0), (m % l, 1)
    n2 = lambda w: w[0] * w[0] + w[1] * w[1]
    if n2(a) < n2(b):
        a, b = b, a
    while True:
        nb = n2(b)
        k = (2 * (a[0] * b[0] + a[1] * b[1]) + nb) // (2 * nb)
        a = (a[0] - k * b[0], a[1] - k * b[1])
        if n2(a) >= nb:
            return b, a
        a, b = b, a


def box_points(l: int, m: int, S: int):
    w1, w2 = reduced_basis(l, m)
    for j in range(0, S + 1):
        for i in range(-S, S + 1):
            if j == 0 and i <= 0:
                continue
            a = i * w1[0] + j * w2[0]
            b = i * w1[1] + j * w2[1]
            if b < 0:
                a, b = -a, -b
            if b == 0 or gcd(a, b) != 1:
                continue
            yield a, b


def rational_smooth_points(l: int, m: int, S: int, fb: FactorBase) -> List[Tuple[int, int]]:
    out = []
    for a, b in box_points(l, m, S):
        t = a - b * m
        assert t % l == 0
        if fb.factor(t // l) is not None:
            out.append((a, b))
    return out


def F_value(c: List[int], m: int, r0: int, r1: int, a: int, b: int) -> int:
    c0, c1, c2 = c
    Fh = a ** 3 + c2 * a * a * b + c1 * a * b * b + c0 * b ** 3
    return Fh + b * (a - b * m) * (r1 * a + r0 * b)


def count_Nl(N: int, c: List[int], m: int, r0: int, r1: int,
             pts: List[Tuple[int, int]], fb: FactorBase) -> int:
    n = 0
    for a, b in pts:
        F = F_value(c, m, r0, r1, a, b)
        assert (F - N * b ** 3) % (a - b * m) == 0
        if F != 0 and fb.factor(F) is not None:
            n += 1
    return n
