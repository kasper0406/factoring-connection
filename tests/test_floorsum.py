"""Tests for src/floorsum (exact Farey floor-sum oracle, operational ladder, lattice block)."""
import random
from fractions import Fraction
from math import isqrt

import gmpy2
import pytest

from src.analytic.ladder import _Block, iroot_floor
from src.common.instances import balanced_semiprime, prime_instance
from src.floorsum.farey import (Counters, FloorSumConvex, _convergent_le, floor_sum_brute,
                                floor_sum_poly, poly_eval)
from src.floorsum.ladder_ops import _ladder_coef, ladder_factor_op
from src.floorsum.lattice_block import _cond_ok, best_params, lattice_block, lattice_factor


def test_run_lemma_R():
    # floor(P) = floor((p u + r)/q) whenever r <= q P - p u < r + 1
    rng = random.Random(0)
    for _ in range(2000):
        q = rng.randint(1, 50)
        p = rng.randint(-100, 100)
        u = rng.randint(-100, 100)
        r = rng.randint(-100, 100)
        theta = Fraction(rng.randint(0, 999), 1000)
        P = Fraction(p * u + r, q) + theta / q
        assert P.numerator // P.denominator == (p * u + r) // q


def test_convergent():
    assert _convergent_le(355, 113, 10) == (22, 7)
    assert _convergent_le(-7, 3, 5) == (-7, 3)


@pytest.mark.parametrize("deg", [2, 3, 4])
def test_random_convex_exact(deg):
    rng = random.Random(deg)
    for _ in range(300):
        den = rng.randint(1, 10 ** 12)
        H = rng.randint(1, 800)
        coef = [rng.randint(-10 ** 12, 10 ** 12), rng.randint(-10 ** 12, 10 ** 12),
                rng.randint(1, 10 ** 8)] + [0] * (deg - 2)
        if deg >= 3:
            coef[3] = rng.randint(0, 100)
        if deg == 4:
            coef[4] = rng.randint(0, 5)
        a = rng.randint(-50, 50)
        assert floor_sum_poly(coef, den, a, a + H) == floor_sum_brute(coef, den, a, a + H)


def test_concave_exact():
    rng = random.Random(7)
    for _ in range(200):
        den = rng.randint(1, 10 ** 9)
        coef = [rng.randint(-10 ** 9, 10 ** 9), rng.randint(-10 ** 9, 10 ** 9), -rng.randint(1, 10 ** 5)]
        H = rng.randint(1, 500)
        assert floor_sum_poly(coef, den, 0, H) == floor_sum_brute(coef, den, 0, H)


@pytest.mark.parametrize("j", [2, 3, 4])
def test_ladder_block_ranges(j):
    rng = random.Random(10 + j)
    for _ in range(25):
        N = rng.getrandbits(44) | (1 << 43) | 1
        B = iroot_floor(N, j + 2)
        x = rng.randint(2 * B, isqrt(N))
        H = max(1, x // (2 * B))
        blk = _Block(N, x, H, j)
        eta = Fraction(1, 4 * blk.den * blk.rho.denominator)
        for sh in (Fraction(0), -blk.rho - eta, blk.rho):
            c = _ladder_coef(blk, sh.numerator, sh.denominator)
            d = blk.den * sh.denominator
            F = FloorSumConvex(c, d, 0, H)
            for _ in range(3):
                h1 = rng.randint(0, H)
                h2 = rng.randint(h1, H)
                assert F.range_sum(h1, h2) == floor_sum_brute(c, d, h1, h2)


def test_cost_scale_one_block():
    # windows + runs tracks H beta^{1/3} (heuristic cost; generous band)
    N = (1 << 61) + 12345
    j = 2
    B = iroot_floor(N, 4)
    x = isqrt(N) - 10 ** 6
    H = x // (2 * B)
    blk = _Block(N, x, H, j)
    cnt = Counters()
    FloorSumConvex(_ladder_coef(blk, 0, 1), blk.den, 0, H, cnt)
    scale = H * (N / x ** 3) ** (1 / 3)
    assert 0.3 * scale <= cnt["windows"] + cnt["runs"] <= 10 * scale
    assert cnt["windows"] + cnt["runs"] < H / 50


@pytest.mark.parametrize("j", [2, 3])
def test_ladder_op_correct(j):
    for s in range(3):
        N = balanced_semiprime(28, s).N
        f = ladder_factor_op(N, j)
        assert f is not None and 1 < f < N and N % f == 0
    for s in range(2):
        assert ladder_factor_op(prime_instance(26, s).N, j) is None


def test_lattice_block_guarantee():
    inst = balanced_semiprime(40, 3)
    N = inst.N
    p = inst.factors[0]
    m, n, X = best_params(N, p - 5)
    assert X > 5 and _cond_ok(N, p - 5, X, m, n)
    assert p in lattice_block(N, p - 5, m, n, X)


def test_lattice_factor_small():
    inst = balanced_semiprime(32, 1)
    f = lattice_factor(inst.N)
    assert f is not None and inst.N % f == 0 and 1 < f < inst.N


@pytest.mark.parametrize("j", [1, 2, 3, 4])
def test_divisor_summatory_ladder(j):
    from src.floorsum.ladder_ops import divisor_summatory_ladder
    rng = random.Random(100 + j)
    for _ in range(6):
        N = rng.randint(10 ** 5, 3 * 10 ** 6)
        s = isqrt(N)
        D = 2 * sum(N // d for d in range(1, s + 1)) - s * s
        assert divisor_summatory_ladder(N, j) == D


def test_inflection_rejected():
    # P(h) = (h - 50)^3 / 1000 has an inflection inside [0, 100)
    coef = [-125000, 7500, -150, 1]
    with pytest.raises(ValueError):
        floor_sum_poly(coef, 1000, 0, 100)
    assert floor_sum_poly(coef, 1000, 50, 100) == floor_sum_brute(coef, 1000, 50, 100)
