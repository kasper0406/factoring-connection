"""OPERATIONAL linearised-cell small-square generator (Vallée-type; NOT Vallée's full algorithm).

A cell is x1 + t, |t| <= T.  Since (x1 + t)^2 = x1^2 + 2 x1 t + t^2, the points with
s := (x1^2 + 2 x1 t) mod N (centred) in [-h/2, h/2] and t^2 <= h/4 have |(x1+t)^2 mod N| <= 3h/4.
They are the points of an affine lattice of determinant N in a T x h box, found exactly by
Gauss reduction and enumeration.  No Farey treatment of skewed cells, no quasi-uniformity claim:
this only produces small squares (|z| <= h, h ~ 4 N^{2/3}) for the RK membership tests.
"""
from __future__ import annotations

import random
from math import isqrt
from typing import List, Tuple

import gmpy2


def _gauss(u, v, K):
    """Lagrange-Gauss reduction for the norm (K t)^2 + s^2 on integer vectors (t, s)."""
    def n2(w):
        return (K * w[0]) ** 2 + w[1] ** 2

    def ip(w, z):
        return K * K * w[0] * z[0] + w[1] * z[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        nu = n2(u)
        mu = (2 * ip(u, v) + nu) // (2 * nu)          # round(ip/nu)
        v = (v[0] - mu * u[0], v[1] - mu * u[1])
        if n2(v) >= nu:
            return u, v
        u, v = v, u


def cell_points(N: int, x1: int, T: int, H: int) -> List[int]:
    """All t, |t| <= T, with s = (x1^2 + 2 x1 t) mod N centred and |s| <= H (exact)."""
    c0 = x1 * x1 % N
    K = max(1, H // max(T, 1))
    e1, e2 = _gauss((1, 2 * x1 % N), (0, N), K)
    D = e1[0] * e2[1] - e1[1] * e2[0]                # +-N
    # affine target: points (t, s) = (0, c0) + i e1 + j e2 + k (0, N) collapse into the lattice;
    # move the shift close to the origin by Babai rounding (exact).
    p0 = (0, c0)
    # coordinates of -p0: solve i e1 + j e2 = -p0
    ni = (-p0[0] * e2[1] + p0[1] * e2[0])
    nj = (-e1[0] * p0[1] + e1[1] * p0[0])
    i0 = (2 * ni + D) // (2 * D) if D > 0 else -((2 * ni - D) // (2 * D))
    j0 = (2 * nj + D) // (2 * D) if D > 0 else -((2 * nj - D) // (2 * D))
    q = (p0[0] + i0 * e1[0] + j0 * e2[0], p0[1] + i0 * e1[1] + j0 * e2[1])
    # enumerate: points q + i e1 + j e2 with |t| <= T, |s| <= H.  Brute bound on i, j via
    # |coordinate| <= (box radius + |q|) * |other vector| / |D| in the K-weighted norm.
    Rt, Rs = T + abs(q[0]), H + abs(q[1])
    def bound(w):
        return (Rt * abs(w[1]) + Rs * abs(w[0])) // abs(D) + 2
    Ji, Jj = bound(e2), bound(e1)
    if (2 * Ji + 1) * (2 * Jj + 1) > 200000:
        return []                                  # skewed cell: skipped (disclosed)
    out = []
    for j in range(-Jj, Jj + 1):
        bt, bs = q[0] + j * e2[0], q[1] + j * e2[1]
        lo, hi = -Ji, Ji
        for i in range(lo, hi + 1):
            t, s = bt + i * e1[0], bs + i * e1[1]
            if abs(t) <= T and abs(s) <= H:
                out.append(t)
    return out


def small_squares(N: int, n_cells: int, seed, h: int) -> List[Tuple[int, int]]:
    """Return (x, z) with z = x^2 mod N centred, |z| <= h, one random point per productive cell."""
    rng = random.Random(f"rkrank-sampler|{N}|{seed}")
    T = isqrt(h) // 2
    H = h // 2
    out = []
    for _ in range(n_cells):
        x1 = rng.randrange(1, N)
        ts = cell_points(N, x1, T, H)
        if not ts:
            continue
        t = rng.choice(ts)
        x = (x1 + t) % N
        z = x * x % N
        if z > N // 2:
            z -= N
        assert abs(z) <= h
        out.append((x, z))
    return out


def h_two_thirds(N: int) -> int:
    return 4 * int(gmpy2.iroot(gmpy2.mpz(N) ** 2, 3)[0])
