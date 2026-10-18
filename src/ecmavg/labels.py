"""DIAGNOSTIC / HARNESS (factor-aware): exposure labels modulo known primes.

Never imported by operational code.  For curve j of a family and a prime p, the staged
exposure label is

  k_p(j) = 0                      if p | D_j (degenerate curve or parameters mod p),
           i  (1 <= i <= pi(B))   if Z of [prod_{i' <= i} l_{i'}^{e_{i'}}] P_j vanishes mod p
                                  first at step i  (= index of the largest prime of
                                  ord(P_j mod p), when that order is B-smooth),
           INF (= -1 here)        if the point order is not B-smooth.

Because steps are in increasing order of l, the label for a smaller cutoff B' is obtained by
truncation: k' = k if k <= pi(B') else INF.  So labels are computed once at the largest B.
Vectorised over (curves x primes) with int64 numpy arithmetic (p < 2^31 required).
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

from .curves import curve_params_int, stage_steps

INF = -1


def _init(family: str, J: int, primes: Sequence[int]):
    C, M = J, len(primes)
    X = np.ones((C, M), dtype=np.int64)
    Z = np.ones((C, M), dtype=np.int64)
    A24 = np.zeros((C, M), dtype=np.int64)
    deg = np.zeros((C, M), dtype=bool)
    for j in range(J):
        (x0, z0), (num, den), D = curve_params_int(family, j)
        for i, p in enumerate(primes):
            if D % p == 0:
                deg[j, i] = True
                continue
            X[j, i], Z[j, i] = x0 % p, z0 % p
            A24[j, i] = num * pow(den, -1, p) % p
    return X, Z, A24, deg


def _ladder_vec(s: int, X, Z, A24, P):
    """[s](X:Z) mod P (P broadcast (1, M)), x-only Montgomery ladder, all entries at once."""
    def dbl(x, z):
        t1 = (x + z) % P
        t1 = t1 * t1 % P
        t2 = (x - z) % P
        t2 = t2 * t2 % P
        t3 = (t1 - t2) % P
        return t1 * t2 % P, t3 * ((t2 + A24 * t3 % P) % P) % P

    def add(x1, z1, x2, z2):
        u = (x1 - z1) % P * ((x2 + z2) % P) % P
        v = (x1 + z1) % P * ((x2 - z2) % P) % P
        w = (u + v) % P
        d = (u - v) % P
        return Z * (w * w % P) % P, X * (d * d % P) % P

    r0x, r0z = X, Z
    r1x, r1z = dbl(X, Z)
    for bit in bin(s)[3:]:
        if bit == "1":
            r0x, r0z = add(r0x, r0z, r1x, r1z)
            r1x, r1z = dbl(r1x, r1z)
        else:
            r1x, r1z = add(r0x, r0z, r1x, r1z)
            r0x, r0z = dbl(r0x, r0z)
    return r0x, r0z


def exposure_labels(family: str, J: int, primes: Sequence[int], B: int, Xbound: int,
                    chunk: int = 4096) -> np.ndarray:
    """int16 array (J, len(primes)) of staged labels k_p(j) (INF = -1)."""
    primes = list(primes)
    assert max(primes) < 2 ** 31
    steps = stage_steps(B, Xbound)
    out = np.full((J, len(primes)), INF, dtype=np.int16)
    for c0 in range(0, len(primes), chunk):
        ps = primes[c0:c0 + chunk]
        P = np.array(ps, dtype=np.int64)[None, :]
        X, Z, A24, deg = _init(family, J, ps)
        k = np.full(X.shape, INF, dtype=np.int16)
        k[deg] = 0
        for i, (_l, pe) in enumerate(steps, start=1):
            X, Z = _ladder_vec(pe, X, Z, A24, P)
            hit = (Z % P == 0) & (k == INF)
            k[hit] = i
        out[:, c0:c0 + len(ps)] = k
    return out


def truncate(k: np.ndarray, nsteps: int) -> np.ndarray:
    """Labels for a smaller cutoff with nsteps = pi(B') stages."""
    k2 = k.copy()
    k2[k2 > nsteps] = INF
    return k2


# ---------------------------------------------------------------------------
# Independent slow check (tests only): affine Montgomery group law with y.
# ---------------------------------------------------------------------------
def naive_label(family: str, j: int, p: int, B: int) -> int:
    """Staged label by brute force: the order of P_j mod p from the affine group law on
    By^2 = x^3 + A x^2 + x with B chosen so that P = (x0/z0, 1) lies on it; then the index
    of its largest prime factor among primes <= B (INF if not B-smooth; 0 if degenerate)."""
    from sympy import factorint, primerange
    (x0, z0), (num, den), D = curve_params_int(family, j)
    if D % p == 0:
        return 0
    a24 = num * pow(den, -1, p) % p
    A = (4 * a24 - 2) % p
    x = x0 * pow(z0, -1, p) % p
    Bc = (x ** 3 + A * x * x + x) % p
    if Bc == 0:                        # P is 2-torsion: order 2
        order = 2
    else:
        P0 = (x, 1)

        def add(P1, P2):
            if P1 is None:
                return P2
            if P2 is None:
                return P1
            (x1, y1), (x2, y2) = P1, P2
            if x1 == x2 and (y1 + y2) % p == 0:
                return None
            if P1 == P2:
                lam = (3 * x1 * x1 + 2 * A * x1 + 1) * pow(2 * Bc * y1, -1, p) % p
            else:
                lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
            x3 = (Bc * lam * lam - A - x1 - x2) % p
            y3 = (lam * (x1 - x3) - y1) % p
            return (x3, y3)

        Q, order = P0, 1
        while Q is not None:
            Q = add(Q, P0)
            order += 1
    fac = factorint(order)
    if max(fac) > B:
        return INF
    pl = list(primerange(2, B + 1))
    return pl.index(max(fac)) + 1 if order > 1 else 0
