"""Tests for the explicit Theorem 1 bound (src/sieve2d/bound.py) and for the
case analysis of its proof (notes/sieve2d_results.md §2), checked config by config."""
from itertools import permutations, product
from math import gcd, isqrt, sqrt

import numpy as np
import pytest

from src.dcltheorem.params import primes_upto, weight_set
from src.sieve2d.bound import check_box, err_terms


def gauss_reduce(b1, b2):
    """Lagrange-Gauss reduction (exact integers); returns (shortest, second)."""
    n = lambda v: v[0] * v[0] + v[1] * v[1]
    if n(b1) > n(b2):
        b1, b2 = b2, b1
    while True:
        # round(<b1,b2>/<b1,b1>) exactly
        num = b1[0] * b2[0] + b1[1] * b2[1]
        den = n(b1)
        m = (2 * num + den) // (2 * den)
        b2 = (b2[0] - m * b1[0], b2[1] - m * b1[1])
        if n(b2) >= n(b1):
            return b1, b2
        b1, b2 = b2, b1


def lattice_of(cfg):
    """cfg: tuple of (l, (A, B)) with distinct l.  Lattice {a == r b mod d}."""
    d, r = 1, 0
    for l, (A, B) in cfg:
        ri = A * pow(B, -1, l) % l
        # CRT r mod d, ri mod l
        t = ((ri - r) * pow(d, -1, l)) % l
        r, d = r + d * t, d * l
    return d, r % d


def test_lattice_lemma_random():
    rng = np.random.default_rng(7)
    X = 120
    a = np.arange(1, X + 1)
    A2, B2 = np.meshgrid(a, a, indexing="ij")
    for _ in range(300):
        d = int(rng.integers(2, 5000))
        r = int(rng.integers(0, d))
        cnt = int((((A2 - r * B2) % d) == 0).sum())
        u, _ = gauss_reduce((d, 0), (r, 1))
        rho = sqrt(u[0] ** 2 + u[1] ** 2)
        assert abs(cnt - X * X / d) <= 1 + 3.266 * X / rho + 1e-9


@pytest.mark.parametrize("H,Y,X,jmax", [(2, 64, 400, 4), (3, 144, 300, 3)])
def test_case_analysis_per_config(H, Y, X, jmax):
    V = weight_set(H)
    band = [l for l in primes_upto(2 * Y) if Y < l <= 2 * Y][:5]
    a = np.arange(1, X + 1, dtype=np.int64)
    A2, B2 = np.meshgrid(a, a, indexing="ij")
    xa, xb = A2.ravel(), B2.ravel()
    online = np.zeros(xa.shape, dtype=bool)
    for (A, B) in V:
        online |= (B * xa == A * xb)
    xa, xb = xa[~online], xb[~online]
    Smax = lambda v: max(v)
    for j in range(1, jmax + 1):
        for ls in permutations(band, j):
            for vs in product(V, repeat=j):
                cfg = tuple(zip(ls, vs))
                d, r = lattice_of(cfg)
                cnt = int((((xa - r * xb) % d) == 0).sum())
                u, _ = gauss_reduce((d, 0), (r, 1))
                rho = sqrt(u[0] ** 2 + u[1] ** 2)
                assert rho <= 1.0746 * sqrt(d) + 1e-9
                tv = {}
                for v in V:
                    t = 1
                    for l, w in cfg:
                        if w != v:
                            t *= l
                    tv[v] = t
                par = [v for v in V if u[0] * v[1] - u[1] * v[0] == 0]
                other = 0.0
                for v in V:
                    if v in par:
                        continue
                    rv = sum(1 for _, w in cfg if w != v)
                    assert rv >= -(-j // 2)          # r_{v'} >= ceil(j/2)
                    other += X / (tv[v] * max(v))
                if not par:                          # case (i)
                    bnd = 1 + 3.266 * X / rho + other
                else:                                # case (ii)
                    v = par[0]
                    rv = sum(1 for _, w in cfg if w != v)
                    assert rv <= j // 2
                    g = d // tv[v]
                    if g > H * X:
                        assert cnt == 0
                        continue
                    nv = sqrt(v[0] ** 2 + v[1] ** 2)
                    bnd = 2 * sqrt(2) * X / (tv[v] * nv) + 2 * H * X / g + other
                assert abs(cnt - X * X / d) <= bnd + 1e-9, (cfg, cnt, X * X / d, bnd)


@pytest.mark.parametrize("H,Y,X", [(4, 256, 500), (5, 400, 400)])
def test_theorem1_errors_hold_exhaustively(H, Y, X):
    V = weight_set(H)
    band = [l for l in primes_upto(2 * Y) if Y < l <= 2 * Y]
    r = check_box(H, X, Y, band, V)
    for j in range(5):
        assert r["actual_err"][j] <= r["bound_err"][j]
    if "sum_dev4" in r:
        assert r["sum_dev4"] <= r["dev4_bound"]
        assert r["uncovered"] * r["lam"] ** 4 <= r["sum_dev4"] + 1e-6
