"""Tests for src/ldbound (notes/ldbound_results.md)."""
from fractions import Fraction
from math import gcd
import random

import sympy

from src.ldbound import params
from src.ldbound.checks import (bt_ratio, field_w, lemma12_w, lift_chain, pair_count_direct,
                                pair_count_ramanujan)
from src.pcest.geometry import hom_eval


def test_theta3_formula_and_limits():
    for d in (3, 5, 7, 9, 11, 19):
        s = Fraction(1, d + 1)
        for k in (d, 2 * d, 50, 1000):
            if 4 * k <= d + 1:
                continue
            assert params.theta3(d, k) == 3 * s + Fraction(d, 4 * k)
            if params.theta3(d, k) < 1:
                assert all(v > 0 for v in params.margins(d, k).values()), (d, k)
    assert abs(params.c_const(Fraction(3, 10)) - 0.91924) < 1e-4
    assert params.min_k(9) == 20 and params.min_k(9, design="a2") == 12
    assert params.min_k(11) < params.min_k(9) < params.min_k(7)


def test_mu_optimised_design():
    # codex ldbound_2: d = 9, k = 30, mu = 8/75 gives theta = 197/600
    mu, th = params.best_mu(9, 30)
    assert (mu, th) == (Fraction(8, 75), Fraction(197, 600))
    for d in (5, 7, 9, 11, 13, 19):
        for k in (d, 2 * d, 30, 100, 1000):
            if 4 * k <= d + 1:
                continue
            mu, th = params.best_mu(d, k)
            assert th <= params.theta3(d, k)
            assert th == params.theta_mu(d, k, mu)
            if th < 1:
                assert all(v > 0 for v in params.margins(d, k, mu).values()), (d, k)
            # brute-force check of optimality on a grid
            grid = [Fraction(i, 2000) for i in range(1, 2000) if Fraction(i, 2000) > Fraction(1, 4 * k)]
            assert min(params.theta_mu(d, k, m) for m in grid) >= th
        assert params.c_const(params.theta3opt(d, 10 ** 6)) - params.c_const(Fraction(3, d + 1)) < 1e-4
    assert params.min_k(9, design="a3opt") == 13


def test_field_w_resultant_and_lemma12():
    rng = random.Random(1)
    for _ in range(40):
        d = rng.choice([3, 5, 7])
        N = rng.randrange(10 ** 20, 10 ** 21) | 1
        M = int(round(N ** (1 / (d + 1))))
        while True:
            g0, g1 = rng.randrange(M // 2, M), rng.randrange(M // 2, M)
            if gcd(g0, g1) == 1 and g1 >= 2:
                break
        w = sympy.nextprime(rng.randrange(10 ** 6, 2 * 10 ** 6))
        R = [rng.randrange(-50, 51) for _ in range(d)]
        f = field_w(N, w, d, g0, g1, R)
        assert hom_eval(f, g0, g1) == N * w
        # |a_i| <= g1/2 (i>=1) for the GB part, |a_0| = O(M w)
        for _ in range(20):
            a, b = rng.randrange(-30, 31), rng.randrange(1, 31)
            if gcd(a, b) == 1:
                assert lemma12_w(f, N, w, g0, g1, a, b)


def test_ramanujan_parseval_identity():
    rng = random.Random(2)
    for _ in range(15):
        vals = [rng.randrange(-200, 200) for _ in range(rng.randrange(5, 25))]
        m = rng.randrange(1, 40)
        assert pair_count_direct(vals, m) == pair_count_ramanujan(vals, m)


def test_brun_titchmarsh_ratio():
    for t, W in ((97, 20000), (210, 30000), (1001, 60000)):
        r = bt_ratio(t, W)
        assert r["ok"]
        assert r["ratio_to_uniform"] < 3.0


def test_lift_chain_toy():
    N = 10007 * 10009
    d = 3
    g0, g1 = 71, 89
    w = 1013
    B = 60
    for P, P2 in (((1, 3), (5, 2)), ((-4, 7), (3, 5)), ((2, 9), (-7, 4))):
        out = lift_chain(N, w, d, (g0, g1), P, P2, T=400, B=B)
        assert out["r0_ok"] == 1
        assert out["chain_ok"] == 1
        assert out["Q"] + out["Q2"] > 0
