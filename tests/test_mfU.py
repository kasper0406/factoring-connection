import random
from fractions import Fraction as Fr

import pytest

from src.mfU import params, toy


@pytest.mark.parametrize("a,d0,kappa", [(Fr(3, 5), Fr(1), Fr(1, 8)), (Fr(9, 10), Fr(1), Fr(47, 200)),
                                        (Fr(13, 10), Fr(1, 2), Fr(23, 400)), (Fr(1, 3), Fr(3), Fr(2))])
def test_margins_hold(a, d0, kappa):
    for name, ok, where in params.margins(a, d0, kappa):
        assert ok, (name, where)


def test_cost_bracket_and_terms():
    lo, hi = params.cost_constant(Fr(9, 10), Fr(1), Fr(47, 200), Fr(0))
    assert lo <= hi and hi - lo < Fr(1, 10 ** 3)
    # sigma_r + sigma_F = sigma_bar, and sigma_F > sigma_r
    sr, sF = params.sigma_r(Fr(1), Fr(1)), params.sigma_F(Fr(1), Fr(1), Fr(1, 4))
    assert sr + sF == params.sigma_bar(Fr(1), Fr(1), Fr(1, 4)) and sF > sr
    # c_MF is at least the recursion term with s = 0
    assert lo >= Fr(9, 10) + params.sigma_F(Fr(9, 10), Fr(1), Fr(47, 200))


def test_base_expansion_composite_w():
    rng = random.Random(1)
    for _ in range(300):
        d = rng.choice([3, 5, 7])
        N = rng.randrange(10 ** 20, 10 ** 21) | 1
        M = int(round(N ** (1.0 / d)))
        g0, g1 = rng.randrange(M // 2, M), rng.randrange(M // 2, M)
        if g1 < 2 or toy.gcd(g0, g1) != 1:
            continue
        w = rng.randrange(2, 10 ** 6)          # arbitrary integer, often composite
        a = toy.base_expansion(N, w, g0, g1, d)
        assert toy.form_value(a, g0, g1) == N * w
        assert all(abs(c) <= g1 // 2 + 1 for c in a[1:])


def test_forced_residue_composite_w():
    rng = random.Random(2)
    done = 0
    while done < 300:
        d = rng.choice([3, 5])
        N = rng.randrange(10 ** 15, 10 ** 16) | 1
        g0, g1 = rng.randrange(100, 400), rng.randrange(100, 400)
        if toy.gcd(g0, g1) != 1:
            continue
        P = (rng.randrange(-50, 50), rng.choice([53, 59, 61, 67, 71, 73]))
        t = abs(g1 * P[0] - g0 * P[1])
        if t == 0 or toy.gcd(t, N * g1 * P[1]) != 1 or toy.gcd(P[0], P[1]) != 1:
            continue
        for w in (rng.randrange(1, 10 ** 5) for _ in range(5)):
            assert toy.forced_residue_ok(N, w, g0, g1, d, P)
        done += 1


@pytest.mark.parametrize("W,t", [(1000, 7), (1000, 999), (12345, 101), (5000, 4999)])
def test_residue_law(W, t):
    assert toy.residue_law_deviation(W, t) <= Fr(1, W)


def test_norm_class_odd_degree():
    rng = random.Random(3)
    fs = [[-2, 0, 0, 1], [1, 1, 0, 1], [-3, 1, 0, 0, 0, 1]]   # x^3-2, x^3+x+1, x^5+x-3
    for f in fs:
        for _ in range(10):
            xi = [rng.randrange(-5, 6) for _ in range(len(f) - 1)]
            if all(c == 0 for c in xi):
                continue
            delta = rng.choice([-1, 2, -3, 5, 6, -7, 10, 30, -105])
            assert toy.norm_class_ok(f, xi, delta)


def test_recursion_term_dominates():
    # codex/mfU_2: c_MF = (a + sigma_F + s)/(1 - 2^{-2/5}) always
    flo, _ = params.recursion_factor_bracket()
    for a in (Fr(1, 5), Fr(3, 5), Fr(2), Fr(5)):
        for d0 in (Fr(1, 4), Fr(1), Fr(4)):
            for kappa in (Fr(1, 100), Fr(1, 8), Fr(2)):
                for Cg in (Fr(0), Fr(5)):
                    lo, hi = params.cost_constant(a, d0, kappa, Cg)
                    rec = (a + params.sigma_F(a, d0, kappa) + 2 * kappa * (1 + Cg)) * flo
                    assert lo == rec
