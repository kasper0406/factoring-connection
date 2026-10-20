"""Tests for src/almostall13 (notes/almostall13_results.md)."""
from __future__ import annotations

import random

from src.almostall13 import rank as R
from src.almostall13 import analyze_e2 as A2
from src.thmA.numfield import hom_eval, make_f


def test_affine_dim_basic():
    assert R.affine_dim([[1, 2], [1, 2]]) == 0
    assert R.affine_dim([[0, 0], [1, 1], [2, 2]]) == 1
    assert R.affine_dim([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]) == 3


def test_left_kernel_annihilates_every_form():
    d = 3
    pts = R.draw_points(d, 2 * d + 2)
    ker = R.left_kernel(pts, d)
    assert len(ker) == 2 * d + 2 - (d + 1)
    rng = random.Random(7)
    for _ in range(20):
        f = [rng.randrange(-10 ** 9, 10 ** 9) for _ in range(d + 1)]
        v = [hom_eval(f, a, b) for (a, b) in pts]
        for lam in ker:
            assert sum(l * y for l, y in zip(lam, v)) == 0


def test_design_values_are_a_form_whatever_N():
    """F(g0, g1) = N*w and F has d+1 coefficients for every N: N only moves F inside Z[x,y]_d."""
    rng = random.Random(3)
    for _ in range(10):
        N = rng.randrange(2 ** 39, 2 ** 40)
        g0, g1 = 1001, 1013
        f = make_f(N, 37, 3, g0, g1, [rng.randrange(-50, 51) for _ in range(3)])
        assert len(f) == 4 and hom_eval(f, g0, g1) == N * 37


def test_run_case_predictions_small():
    r = R.run_case(3, 5)
    assert r["R1"] and r["R2"] and r["R3"] and r["R4"]
    assert r["aF_A"] == 3 and r["aF_C"] == 4


def test_stats_helpers():
    assert abs(A2.mean([1, 2, 3]) - 2) < 1e-12
    assert A2.sd_between([0.5, 0.5, 0.5], 100) == 0.0
    F = A2.anova_F([[1.0, 1.1, 0.9], [5.0, 5.1, 4.9]])
    assert F > 100


# ---------------------------------------------------------------- bounds.py (Lemma M1-AA, Lemma MP, NR-AA)
from fractions import Fraction

from src.almostall13 import bounds as Bd


def test_window_divisors_and_rankin():
    rng = random.Random(11)
    X = 2 ** 12
    for _ in range(30):
        R = rng.randrange(1, 2 ** 40)
        assert Bd.rankin_half_exact_holds(R, X)
    # a primorial has many divisors in the window, still below the Rankin bound
    R = 1
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        R *= p
    cnt = Bd.window_divisors(R, X)
    assert cnt > 10
    assert Bd.rankin_half_exact_holds(R, X)


def test_mp_mixture():
    law = {5: Fraction(1, 2), 7: Fraction(1, 2)}
    out = Bd.mp_mixture(law, Fraction(1, 1000), 6)
    assert out["p0"] == Fraction(1, 1000)
    assert all(r == Fraction(-1, 1000) for r in out["rel_changes"])


def test_nraa_constant_reduces_to_nr():
    from src.nested.barrier import nr_constant
    assert abs(Bd.nraa_constant(4, 1, 0.75, 0.5, 0.0) - nr_constant(4, 1, 0.75, 0.5)) < 1e-12
    assert Bd.nraa_constant(4, 1, 0.75, 0.5, 100.0) >= 4 + 50 - 1e-9
