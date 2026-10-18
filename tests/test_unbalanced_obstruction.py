"""Tests for the curvature obstruction (Proposition O) and the regime-II span bound."""
import random
from math import isqrt

from src.common.metrics import Metrics
from src.unbalanced.algorithm import window_bounds
from src.unbalanced.obstruction import (best_window_coverage, covered_count, e_profile,
                                        frac_covered_at, frac_theta_below, prop_o_bound_int)
from src.unbalanced.params import Inst, regime2_span_bound


def test_prop_o_bound_random():
    rng = random.Random(20261018)
    worst = 0.0
    for _ in range(400):
        n = rng.randint(28, 44)
        N = rng.randrange(1 << (n - 1), 1 << n)
        j = rng.randint(0, n // 4)
        R = 1 << j
        x1 = isqrt(N // (2 * R)) + 1
        x2 = isqrt(N // R)
        beta = rng.choice([1, 1, 1, 2, 3, -1, -2])
        alpha = abs(beta) * R + rng.randint(-3, 3) * (1 if rng.random() < .5 else R)
        if beta < 0:
            alpha = -alpha if rng.random() < 0.3 else alpha
        M = rng.choice([1, 3, 7, 11, 31, 101, 1009])
        r = rng.randrange(M)
        c0 = rng.randint(-M * M, M * M)
        W = rng.choice([1, 2, 5, 17, 100, 1000])
        if (x2 - x1) // M > 30000:
            continue
        cov = best_window_coverage(N, alpha, beta, c0, M, r, W, x1, x2)
        b = prop_o_bound_int(N, beta, W, x2)
        assert cov <= b, (N, alpha, beta, M, r, W, cov, b)
        worst = max(worst, cov / b)
    assert worst > 0.2          # the bound is not vacuous: some windows come close


def test_covered_count_matches_window_scan():
    N = 1000003 * 7919
    x1, x2 = isqrt(N // 8) + 1, isqrt(N // 4)
    tot = 0
    for g0 in range(0, 400000, 9973):
        tot += covered_count(N, 4, 1, 0, 13, 5, g0, 50, x1, x2)
    assert tot >= 0
    assert best_window_coverage(N, 4, 1, 0, 13, 5, 50, x1, x2) <= prop_o_bound_int(N, 1, 50, x2)


def test_regime2_span_bound_holds():
    rng = random.Random(7)
    for _ in range(200):
        n = rng.randint(48, 96)
        N = rng.randrange(1 << (n - 1), 1 << n)
        Delta = rng.randint(2, 40)
        j = rng.randint(max(1, (2 * Delta).bit_length()), n // 5)
        for B in [1 << j, (1 << j) + rng.randrange(1 << j), (2 << j) - 1]:
            t = Inst(1, B, max(1, B - Delta), 1, B + Delta, 1, j)
            wb = window_bounds(N, t, Metrics())
            if wb is None:
                continue
            lo, hi = wb
            assert hi - lo <= regime2_span_bound(n, j, Delta) + 2


def test_profile_and_measure():
    assert abs(e_profile(0.2) - 0.2) < 1e-12
    assert e_profile(0.0) == 1 / 6 and abs(e_profile(1 / 3) - 1 / 6) < 1e-12
    assert max(e_profile(k / 3000) for k in range(1001)) <= 0.2 + 1e-12
    assert abs(frac_theta_below(1 / 3) - 1.0) < 1e-12
    f = frac_covered_at(0.005)
    assert 0.1 < f < 0.25          # ~ 30 eta: not almost all
    assert frac_covered_at(1 / 30 + 1e-9) == 1.0


def test_prop_o_prime_arbitrary_sets():
    from src.unbalanced.obstruction import best_set_coverage, prop_o_set_bound_int
    rng = random.Random(4242)
    worst = 0.0
    for _ in range(300):
        n = rng.randint(28, 42)
        N = rng.randrange(1 << (n - 1), 1 << n)
        j = rng.randint(0, n // 3)
        R = 1 << j
        x1, x2 = isqrt(N // (2 * R)) + 1, isqrt(N // R)
        beta = rng.choice([1, 1, 2, -1])
        alpha = abs(beta) * (R + rng.randint(0, R)) * (1 if beta > 0 else rng.choice([1, -1]))
        M = rng.choice([1, 5, 23, 97, 401])
        if (x2 - x1) // M > 30000:
            continue
        r, c0 = rng.randrange(M), rng.randint(-M * M, M * M)
        m = rng.choice([1, 3, 10, 50, 300])
        cov = best_set_coverage(N, alpha, beta, c0, M, r, m, x1, x2)
        b = prop_o_set_bound_int(N, beta, m, x2)
        assert cov <= b, (N, alpha, beta, M, m, cov, b)
        worst = max(worst, cov / b)
    assert worst > 0.1
