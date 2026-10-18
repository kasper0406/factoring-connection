"""Tests for src/unbalanced (DCL-wide)."""
from fractions import Fraction
from math import isqrt

import pytest

from src.common.metrics import Metrics
from src.dcltheorem import algorithm as dalg
from src.dcltheorem.params import Design, schedule_design
from src.unbalanced.algorithm import build_targets, dclw, ratio_range, window_bounds
from src.unbalanced.harness import lemma_checks, predict, sample_pair_ratio
from src.unbalanced.params import (Inst, WideDesign, regime1_insts, regime1_Jmax,
                                   regime1_weights, regime2_design, wide_from_C1)


def test_band0_reproduces_dprime_targets():
    D = schedule_design(40, 3.0)
    W = wide_from_C1(40, [0])
    for s in range(4):
        N, p, q = sample_pair_ratio(40, 0, s)
        for k in range(2):
            t1 = dalg.build_targets(N, D, k, Metrics())
            t2 = build_targets(N, W, k, Metrics())
            key = lambda t: (t["A"], t["B"], t["l"], t["r"], t["c"], t["g_lo"], t["W"])
            assert [key(t) for t in t1] == [key(t) for t in t2]


def test_ratio_range_exact():
    for N in [10007 * 10009, 3 * 2 ** 40 + 7, 99991 * 1000003]:
        for (a, b, c, d) in [(1, 1, 2, 1), (4, 1, 8, 1), (17, 3, 20, 3), (90, 1, 130, 1)]:
            t = Inst(1, 2, a, b, c, d, 0)
            lo, hi = ratio_range(N, t)
            for x in range(max(1, lo - 3), hi + 4):
                inside = Fraction(a, b) < Fraction(N, x * x) <= Fraction(c, d)
                assert inside == (lo <= x <= hi), (N, x)


def test_regime1_span_bound():
    """Every regime-I window span is <= (1+sqrt2) H sqrt N + 2 (Theorem U, Lemma U1)."""
    H = 9
    for N in [2 ** 39 + 15, 2 ** 40 - 3, 3 * 2 ** 38 + 1]:
        for j in range(regime1_Jmax(H) + 1):
            for (A, B) in regime1_weights(H, j):
                wb = window_bounds(N, Inst(A, B, 1 << j, 1, 1 << (j + 1), 1, j), Metrics())
                if wb is None:
                    continue
                lo, hi = wb
                # (1+sqrt2) < 2.4143
                assert (hi - lo) * 10000 <= 24143 * H * (isqrt(N) + 1) + 20000


def test_regime1_weight_bounds():
    H = 12
    for j in range(regime1_Jmax(H) + 1):
        V = regime1_weights(H, j)
        a = max(A for A, _ in V)
        b = max(B for _, B in V)
        assert a * b <= H * H and a * a * (1 << j) <= H * H and b * b <= H * H * (1 << j)


@pytest.mark.parametrize("j", [0, 3, 6])
def test_wide_end_to_end_lemmas(j):
    W = wide_from_C1(40, list(range(7)))
    for s in range(3):
        N, p, q = sample_pair_ratio(40, j, s)
        assert p * q == N and (1 << j) * p < q <= (2 << j) * p
        r = dclw(N, W, early_exit=False, record=True)
        assert r["factor"] is None or (1 < r["factor"] < N and N % r["factor"] == 0)
        for b in r["blocks"]:
            if "targets" in b:
                lc = lemma_checks(N, p, q, W, b)
                assert sum(lc["bad"].values()) == 0
        pr = predict(N, p, q, W, [b.get("m", 0) for b in r["blocks"]])
        if pr["certified"]:
            assert r["factor"] in (p, q)


def test_regime2_end_to_end():
    W = regime2_design(56, 8, 12)
    for s in range(2):
        N, p, q = sample_pair_ratio(56, 8, s)
        r = dclw(N, W, early_exit=False, record=True)
        assert r["factor"] is None or (1 < r["factor"] < N and N % r["factor"] == 0)
        for b in r["blocks"]:
            if "targets" in b:
                assert sum(lemma_checks(N, p, q, W, b)["bad"].values()) == 0


def test_output_validation_small():
    W = wide_from_C1(32, list(range(4)))
    for N in list(range(2, 400)):
        r = dclw(N, W)
        f = r["factor"]
        assert f is None or (1 < f < N and N % f == 0), N
