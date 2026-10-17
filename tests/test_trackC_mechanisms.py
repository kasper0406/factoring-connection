"""Tests for track C mechanisms (owner factC)."""
import random
from math import gcd

from src.common.metrics import Metrics
from src.trackC.mechanisms import (jacobi, jacobi_prefix_product, power_sum_direct,
                                   power_sum_faulhaber, truncated_FN_test)


def test_faulhaber_matches_direct():
    N = 1000003 * 1000033
    rng = random.Random(3)
    for _ in range(30):
        e = rng.randrange(0, 15)
        A = rng.randrange(1, 10**6)
        B = A + rng.randrange(1, 500)
        r = power_sum_faulhaber(N, A, B, e)
        assert r["S"] == power_sum_direct(N, A, B, e)["S"]


def test_faulhaber_denominator_is_trial_division():
    # e+1 >= 7 forces 1/7 (or 1/5) -> the small prime factor is exposed as a denominator
    r = power_sum_faulhaber(35, 1, 10, 6)
    assert r["S"] is None and r["divisor"] in (5, 7)


def test_odd_power_sum_centered_is_divisible_by_center():
    # Explains the E2 bias: sum_{|t|<=h} (c+t)^e, e odd, is divisible by c.
    for c in range(1, 60):
        for h in range(0, 8):
            for e in (1, 3, 5, 7):
                assert sum((c + t) ** e for t in range(-h, h + 1)) % c == 0


def test_truncated_FN_equals_U():
    for N in range(4, 300):
        lp = next(d for d in range(2, N + 1) if N % d == 0)
        if lp == N:
            continue
        for K in range(2, min(N, lp + 3)):
            assert truncated_FN_test(N, K)["Z"] == (K <= lp)


def test_jacobi_prefix():
    for N in range(3, 400, 2):
        lp = next(d for d in range(2, N + 1) if N % d == 0)
        for K in range(2, min(N, 60)):
            assert (jacobi_prefix_product(N, K)["C"] != 0) == (K <= lp)
