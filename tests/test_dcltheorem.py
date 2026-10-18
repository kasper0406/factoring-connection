"""Tests for src/dcltheorem (DCL-det: deterministic divisor-cover Lehman)."""
from math import gcd, isqrt

import pytest

from src.common.metrics import Metrics
from src.dcltheorem.algorithm import (all_but_one, balanced_range, dcl_det, f_bounds,
                                      metric_E, sqrt_prime)
from src.dcltheorem.harness import (full_closed_form, lemma_checks, ord2, predict,
                                    sample_pair)
from src.dcltheorem.params import Design, schedule_design, theorem_design, weight_set


def test_weight_set_excludes_symmetric_and_noncoprime():
    V = weight_set(6)
    assert (1, 1) not in V and (2, 4) not in V and (2, 3) in V and (3, 2) in V
    assert all(A != B and gcd(A, B) == 1 for A, B in V)


def test_all_but_one_matches_direct():
    N = 1000003 * 999983
    exps = [l * l for l in (7, 11, 19, 23, 31)]
    vals, E = all_but_one(2, exps, N, Metrics())
    for e, v in zip(exps, vals):
        assert v == pow(2, E // e, N)


@pytest.mark.parametrize("l", [43, 47, 59, 61, 73, 97, 101, 113])
def test_sqrt_prime(l):
    for a in range(1, l):
        if pow(a, (l - 1) // 2, l) == 1:
            r = sqrt_prime(a, l, Metrics())
            assert r * r % l == a


def test_f_bounds_cover_balanced_range():
    N, p, q = sample_pair(30, 1)
    a0, b0 = balanced_range(N)
    assert a0 <= p <= b0
    for A, B in [(1, 2), (5, 3), (7, 1)]:
        lo, hi = f_bounds(N, A, B, a0, b0, Metrics())
        for x in (a0, b0, (a0 + b0) // 2, isqrt(A * N // B)):
            if a0 <= x <= b0:
                assert lo * x <= B * x * x + A * N <= hi * x


def test_sample_pair_population():
    for s in range(20):
        N, p, q = sample_pair(36, s)
        assert p < q <= 2 * p and (1 << 35) <= N < (1 << 36)


def test_theorem_design_shapes():
    D = theorem_design(120, 1, 60)
    assert D.K == 4 and all(l % 4 == 3 for b in D.bands for l in b)
    assert D.Y > D.H ** 2                                    # distinct sieve classes


@pytest.mark.parametrize("n", [32, 36])
def test_lemmas_and_certificate(n):
    """Lemmas 1-4, the Lemma 5/7 equivalence on clean blocks, and
    certified => success, on random pairs (harness check)."""
    D = schedule_design(n, 3.0)
    for seed in range(8):
        N, p, q = sample_pair(n, 7000 + seed)
        r = dcl_det(N, D, Metrics(), early_exit=False, record=True)
        assert r["factor"] is None or (1 < r["factor"] < N and N % r["factor"] == 0)
        for b in r["blocks"]:
            if "targets" in b:
                lc = lemma_checks(N, p, q, D, b)
                assert all(v == 0 for v in lc["bad"].values()), lc
        pr = predict(N, p, q, D, [b.get("m", 0) for b in r["blocks"]])
        if pr["certified"]:
            assert r["factor"] in (p, q)


def test_full_closed_form_exact_vs_superset():
    # superset (istar=None) is implied by the exact window
    p, q = 36697, 71249
    d = ord2(q)
    for (A, B, l, m) in [(5, 9, 151, 18), (7, 2, 167, 18)]:
        for istar in range(m):
            if full_closed_form(p, q, A, B, l, m, d, istar):
                assert full_closed_form(p, q, A, B, l, m, d)


def test_metric_E_counts_everything():
    c = {"mulN": 10, "mul_small": 3, "poly_cost": 100, "gcdN": 2, "invN": 1, "bigmul_bits": 129}
    assert metric_E(c, 64) == 10 + 3 + 100 + 6 * 3 + 3


def test_hh_baseline_runs():
    from src.dcltheorem.hh_baseline import hh_params, hh_run
    N, p, q = sample_pair(40, 3)
    r = hh_run(N, hh_params(N))
    assert r["factor"] in (p, q)


def test_exhaustive_outputs_proper_small_N():
    """Regression (GPT-6.1-sol r13): for every 2 <= N <= 10^4 (primes, prime
    powers, even numbers, everything), any output d satisfies 1 < d < N, d | N.
    Two designs per bit length: the C1 schedule and a tiny theorem-shaped one."""
    designs = {}
    for n in range(2, 15):
        designs[n] = [schedule_design(n, 3.0), Design(n=n, H=3, Y=5, K=4)]
    for N in range(2, 10 ** 4 + 1):
        for D in designs[N.bit_length()]:
            r = dcl_det(N, D, Metrics(), early_exit=False)
            d = r["factor"]
            assert d is None or (1 < d < N and N % d == 0), (N, d, r["status"])
            for b in r["blocks"]:
                f = b.get("factor")
                assert f is None or (1 < f < N and N % f == 0), (N, f, b.get("status"))
