import math
import random

import pytest

from src.sihard.core import (dickman_rho, is_smooth_int, johnson_S, johnson_counts, merge_bohr_counts,
                             merge_interval_average, primes_upto, smooth_in_range, smooth_mask)
from src.sihard.reduction import build_instance, count_exact_covers, count_S_smooth_in_window, random_x3c
from src.sihard import consequences as C


def test_smooth_mask_exact():
    ps = primes_upto(30)
    lo = 10 ** 6
    m = smooth_mask(lo, 3000, ps)
    for i in range(3000):
        assert bool(m[i]) == is_smooth_int(lo + i, 30)


def test_smooth_mask_large_and_powers():
    ps = primes_upto(7)
    lo = 2 ** 60 - 5
    m = smooth_mask(lo, 10, ps)
    assert bool(m[5]) is True  # 2^60
    assert sum(m) == sum(is_smooth_int(lo + i, 7) for i in range(10))


def test_dickman_values():
    assert abs(dickman_rho(2) - (1 - math.log(2))) < 1e-6
    assert abs(dickman_rho(3) - 0.0486083882) < 1e-6
    assert abs(dickman_rho(5) / 3.5472470e-4 - 1) < 1e-3


def test_johnson_S_condition():
    for V, P in ((256, 256), (1024, 2 ** 40), (4096, math.prod(primes_upto(60)))):
        S = johnson_S(V, P)
        assert math.log(S) ** 2 >= math.log(V) * math.log(P) - 1e-9
        assert S == 1 or math.log(S - 1) ** 2 < math.log(V) * math.log(P)


def test_johnson_counts_brute():
    # brute force average over x of certified smooth outputs, tiny scale
    X, M, V, y = 50000, 300, 40, 20
    ps = primes_upto(y)
    sm = smooth_in_range(X, X + M + V + 1, y)
    P = 2 ** 3 * 3 ** 2 * 5 * 7
    S = johnson_S(V, P)
    r = johnson_counts(X, M, V, P, S, sm, ps)
    tot = 0
    for x in range(X, X + M):
        tot += sum(1 for z in range(x, x + V + 1) if math.gcd(z, P) >= S and is_smooth_int(z, y))
    assert r["A"] == pytest.approx(tot / M)
    assert r["A"] <= r["B"] <= r["Bd"]


def test_merge_interval_brute():
    rng = random.Random(1)
    A = [rng.randrange(100, 140) for _ in range(15)]
    B = [rng.randrange(100, 140) for _ in range(15)]
    X, h = 10000, 7
    meas, bound = merge_interval_average(A, B, X, h)
    brute = sum(sum(1 for a in A for b in B if x <= a * b <= x + h) for x in range(X, 2 * X)) / X
    assert meas == pytest.approx(brute)
    assert meas <= bound


def test_merge_bohr_small():
    A, B = [6, 10, 15], [7, 11]
    N, h = 13, 2
    expect = sum(1 for a in A for b in B if (a * b) % N <= h and a * b >= N)
    assert merge_bohr_counts(A, B, [N], h) == [expect]


@pytest.mark.parametrize("seed", range(6))
def test_x3c_reduction_parsimonious(seed):
    rng = random.Random(seed)
    q = 2 + seed % 2
    sets = random_x3c(q, q + 3, rng, plant=seed % 3 != 0)
    inst = build_instance(q, sets, rng)
    assert count_S_smooth_in_window(inst["primes"], inst["x"], inst["h"]) == count_exact_covers(q, inst["sets"])


def test_prop_si_constants():
    k, c = C.prop_si_optimum()
    assert abs(k ** 3 - 2 / 9) < 1e-12 and abs(c ** 3 - 16 / 9) < 1e-12
    # oracle-cost robustness threshold vs MNFS
    assert abs(C.tau_threshold_mnfs() - 0.7573) < 1e-3


def test_johnson_class_collapse_constant():
    # Johnson class with nu forced to the 4/(3 kappa) bound would give (32/9)^(1/3); the theorem gives L[2/3]
    assert abs(C.dl_booting_ecm() - 3 ** (1 / 3)) < 1e-9
    c, cost = C.dl_descent_optimum(kappa=(2 / 9) ** (1 / 3))
    assert 0 < c < 1 and cost < 1.2114


def test_rankin_divisor_count_bound():
    # #{d | P, d <= D} <= #{first-k-primes smooth <= D}
    P = 2 * 7 * 11 * 13 * 17 * 19 * 23
    ps = [2, 7, 11, 13, 17, 19, 23]
    first = primes_upto(30)[:len(ps)]
    for D in (50, 500, 5000, 50000):
        a = sum(1 for d in range(1, D + 1) if P % d == 0)
        b = sum(1 for n in range(1, D + 1) if all(p in first for p in sympy_factors(n)))
        assert a <= b


def sympy_factors(n):
    out = []
    p = 2
    while p * p <= n:
        while n % p == 0:
            out.append(p)
            n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out
