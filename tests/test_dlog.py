"""Tests for agent dlog (notes/dlog_results.md)."""
from fractions import Fraction as Fr
import json
import random

import sympy as sp

from src.dlog import params as dp
from src.dlog import pipeline as pl
from src.thmA import params as tp


def test_constants_exact():
    mu, th = tp.best_mu(9, 13)
    assert th == Fr(16, 39)
    val, a = dp.c_dl(th)
    assert sp.simplify(val - 55 * sp.sqrt(2) / 78) == 0
    assert dp.la_never_binds(th)                  # 2*16/39 <= 1
    # whenever Theorem A's constant is < 1, the cubic linear algebra never binds
    for d in (7, 9, 11, 19):
        k = tp.min_k(d)
        th = tp.best_mu(d, k)[1]
        assert (1 + th) ** 2 < 2 and 2 * th <= 1
    # d = 3: LA binds; constant 3*sqrt(3)/4 < sqrt(2) in the k -> infinity limit
    val, a = dp.c_dl(sp.Rational(3, 4))
    assert sp.simplify(val - 3 * sp.sqrt(3) / 4) == 0 and val < sp.sqrt(2)
    # infimum at d = 9 is 13 sqrt2/20
    assert sp.simplify(dp.c_dl(sp.Rational(3, 10))[0] - 13 * sp.sqrt(2) / 20) == 0


def test_solve_from_bezout():
    rng = random.Random(1)
    n = 2 * 3 * 3 * 5 * 7 * 11
    for _ in range(200):
        X = rng.randrange(n)
        eqs = []
        for _ in range(3):
            B = rng.randrange(n)
            A = rng.randrange(n)
            eqs.append((A, B, (A + X * B) % n))
        Y = pl.solve_from(eqs, n)
        from math import gcd
        g = n
        for _, B, _ in eqs:
            g = gcd(g, B)
        if g == 1:
            assert Y == X
        else:
            assert Y is None


def test_bm_bound_positive():
    for d in (3, 5, 9):
        b = pl.bm_bound(d)
        assert b > 0 and b < 1


def test_end_to_end_small():
    cfg = json.load(open("experiments/configs/dlog_e1.json"))["pipeline"]
    p = 1000000007          # 30 bits
    g = 5                   # primitive root mod 1e9+7
    X = 123456789
    h = pow(g, X, p)
    res = pl.dlog(p, g, h, cfg, seed=7)
    assert res["X"] is not None and pow(g, res["X"], p) == h
    c = res["ctr"]
    for bad in ("ROW_IDENTITY_VIOLATION", "TORSION_VIOLATION_W_not_pm1", "TORSION_VIOLATION_aux",
                "KERNEL_VIOLATION", "NONUNIT_KERNEL_VECTOR", "WRONG_ANSWER_BEFORE_VERIFY"):
        assert bad not in c
    assert c["unit_rank_implied"] <= sum(c["r1r2"]) - 1
