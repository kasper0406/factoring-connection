import math
import random

import gmpy2
import numpy as np
import pytest

from src.nativeana import chi as C
from src.nativeana import hankel as HK


def test_chi_table_matches_kronecker():
    for D in (-4 * 1009 * 1013, -4 * 77, -4 * 5 * 9):
        tab = C.chi_table(D, 3000)
        for n in range(1, 3001):
            assert tab[n] == gmpy2.kronecker(D, n)


@pytest.mark.parametrize("N", [101 * 113, 1009 * 1021, 30011 * 30029])
def test_afe_exact(N):
    # D = -4N fundamental for N = 1 mod 4
    if N % 4 != 1:
        pytest.skip("need N = 1 mod 4")
    import cypari2
    pari = cypari2.Pari()
    D = -4 * N
    h = int(pari.qfbclassno(D))
    top = int(6.6 * math.sqrt(4 * N / math.pi)) + 2
    tab = C.chi_table(D, top)
    est = C.afe_h(D, tab, [top])[top]
    assert abs(est - h) < 1e-6 * h


def test_L_formula_nonfundamental():
    # h(-4m) = (2 sqrt m / pi) L(1, chi_{-4m}) also for non-fundamental -4m (orders), m > 1
    import cypari2
    pari = cypari2.Pari()
    for m in (45, 49 * 13, 9 * 101, 1001):
        D = -4 * m
        h = int(pari.qfbclassno(D))
        # long smoothed partial sum as L(1)
        T = 400000
        tab = C.chi_table(D, 6 * T)
        n = np.arange(1, 6 * T + 1, dtype=float)
        L = float(np.sum(tab[1:] * np.exp(-(n / T) ** 2) / n))
        assert abs(math.pi * h / (2 * math.sqrt(m)) - L) < 2e-2


def test_hankel_bound_random_and_lsq():
    rng = np.random.default_rng(1)
    for trial in range(200):
        W = int(rng.integers(2, 12))
        T = int(rng.integers(1, W + 1))
        N = W + int(rng.integers(0, 30))
        r = rng.uniform(0.3, 1.0, T)
        nodes = r * np.exp(2j * np.pi * rng.uniform(0, 1, T))
        c = HK.least_squares_coeffs(nodes, N, W) if trial % 2 else rng.normal(size=T) + 1j * rng.normal(size=T)
        assert HK.window_hankel_opnorm(nodes, c, N, W) >= 1 - 1e-8
        assert HK.window_l1_error(nodes, c, N, W) >= 1 - 1e-8


def test_hankel_bound_tight():
    for W in (3, 7, 10):
        N = 40
        nodes, c = HK.tight_example(N, W)
        assert HK.window_l1_error(nodes, c, N, W) < 1e-8


def test_det_counts_formula():
    from src.nativeana.harness import det_counts
    from fractions import Fraction
    for ell, k in ((3, 1), (3, 2), (5, 1)):
        cnt = det_counts(ell, k)
        m = ell ** k
        for a in range(1, m):
            v, aa = 0, a
            while aa % ell == 0:
                aa //= ell
                v += 1
            pred = Fraction(ell ** (3 * k)) * (1 - Fraction(1, ell * ell)) * sum(Fraction(1, ell ** j) for j in range(v + 1))
            assert pred == cnt[a]


def test_kronecker_hurwitz_small():
    from src.nativeana.harness import hurwitz
    for p, q in ((3, 5), (5, 7), (7, 11), (11, 13), (13, 29)):
        N = p * q
        tmax = math.isqrt(4 * N - 1)
        S0 = sum(hurwitz(4 * N - t * t) for t in range(-tmax, tmax + 1))
        assert S0 == 2 * (N + q)
