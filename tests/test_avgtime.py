"""Tests for src/avgtime (notes/avgtime_results.md)."""
from fractions import Fraction as F
from math import isqrt

import pytest
import sympy

from src.avgtime.detector import small_order_detector
from src.avgtime.exponents import (LADDER_NO_DETECTOR, ONE_SIXTH, dprime_admissible,
                                   dprime_exponents, dprime_requirements, ladder,
                                   ladder_bound_terms, ladder_schedule, single_fallback,
                                   single_fallback_closed, sup_eps)
from src.avgtime.ladder import complete_factorisation, ladder_factor, lehman_split
from src.dcltheorem.harness import sample_pair


# ---------------------------------------------------------------- exponents
def test_dprime_requirement_sups():
    sups = {k: sup_eps(a, b) for k, (a, b) in dprime_requirements().items()}
    assert sups["Y>H^2 (Lemma 9)"] == F(1, 7)
    assert sups["Y^4>sqrt2 P (Lemma 6)"] == F(1, 4)
    assert sups["E_a < N^-eps"] is None            # -eps < 0 for every eps > 0
    assert sups["E_c1 < N^-eps"] == F(1, 52)
    assert sups["E_c2 < N^-eps"] == F(1, 16)
    assert sups["useful: cost < 1/5"] == F(1, 25)
    assert dprime_admissible() == (F(1, 52), "E_c1 < N^-eps")


def test_eps_1_55():
    e = dprime_exponents(F(1, 55))
    assert e["cost"] == F(2, 11)
    assert e["E_c1"] == F(-1, 44) and e["E_c1"] < -F(1, 55)
    assert e["E_c2"] < -F(1, 55) and e["E_a"] == -F(2, 55) and e["E_b"] == -F(1, 55)
    # balance with HH: 1/5 - eps = 2/11
    assert F(1, 5) - F(1, 55) == F(2, 11)


@pytest.mark.parametrize("j", [1, 2, 3])
def test_single_fallback_lp_matches_closed_form(j):
    for det in (False, True):
        v, _ = single_fallback(j, det)
        assert v == single_fallback_closed(j)
    assert single_fallback_closed(1) == F(2, 11)
    assert single_fallback_closed(2) == F(3, 17)


def test_ladder_without_detector_tends_to_19_110():
    vals = [ladder(k, 1, False)[0] for k in range(1, 7)]
    assert vals[0] == F(2, 11) and vals[1] == F(35, 198)
    assert all(a > b for a, b in zip(vals, vals[1:]))
    assert all(v > LADDER_NO_DETECTOR for v in vals)
    assert vals[-1] - LADDER_NO_DETECTOR < F(1, 3000)


def test_ladder_with_detector_goes_below_19_110():
    vals = [ladder(k, 1, True)[0] for k in range(1, 7)]
    assert all(a > b for a, b in zip(vals, vals[1:]))
    assert all(v > ONE_SIXTH for v in vals)
    assert vals[-1] < LADDER_NO_DETECTOR


@pytest.mark.parametrize("n", [256, 512, 1024, 2048, 4096])
def test_ladder_schedule(n):
    Hs = ladder_schedule(n)
    assert Hs[0] ** 2 >= n > (Hs[0] - 1) ** 2 and Hs[-1] == (n * n) << (-(-n // 60))
    assert all(a < b for a, b in zip(Hs, Hs[1:]))
    assert all(b ** 5 <= 32 * a ** 6 for a, b in zip(Hs, Hs[1:]))   # H_{i+1}^{5/3} <= 2^{5/3} H_i^2
    assert len(Hs) <= 6 * n.bit_length()                              # O(log n) stages
    terms = ladder_bound_terms(n)
    # every stage >= 2 contributes <= 2^5 n^3 to (term)^3, i.e. <= 2^{5/3} n per stage
    assert all(t <= 32 * n ** 3 for t in terms[1:])


# ---------------------------------------------------------------- detector
def _ord2(r):
    return int(sympy.n_order(2, r))


def test_detector_exhaustive_small_semiprimes():
    ps = list(sympy.primerange(3, 260))
    for T in (1, 7, 30, 100):
        m = isqrt(T - 1) + 1 if T > 1 else 1
        Teff = -(-T // m) * m          # the detector scans exponents 1..Teff (Teff >= T)
        for i, p in enumerate(ps):
            for q in ps[i + 1:]:
                N = p * q
                d = small_order_detector(N, T)
                dp, dq = _ord2(p), _ord2(q)
                pred = min(dp, dq) <= Teff and dp != dq
                assert (d is not None) == pred, (p, q, T)
                if d is not None:
                    assert d in (p, q)
                    assert _ord2(d) == min(dp, dq)


def test_detector_constructed_small_order():
    # 8191 = 2^13 - 1, 131071 = 2^17 - 1 (both orders <= 20)
    assert small_order_detector(8191 * 131071, 20) == 8191
    assert small_order_detector(8191 * 1000003, 20) == 8191
    assert small_order_detector(8191, 100) is None          # prime
    assert small_order_detector(23 * 89, 100) is None       # ord_23 2 = ord_89 2 = 11


# ---------------------------------------------------------------- ladder (operational)
def test_lehman_split_exhaustive():
    for N in range(2, 3000):
        d = lehman_split(N)
        if sympy.isprime(N) or N < 4:
            assert d is None
        else:
            assert d is not None and 1 < d < N and N % d == 0


def _fac(N):
    return sorted(p for p, e in sympy.factorint(N).items() for _ in range(e))


@pytest.mark.parametrize("n", [32, 36, 40])
def test_ladder_factors_sample_pairs(n):
    for s in range(6):
        N, p, q = sample_pair(n, s)
        r = ladder_factor(N, n, [2, 3, 4])
        assert r["factors"] == sorted([p, q])


def test_ladder_complete_on_adversarial_inputs():
    cases = [2, 3, 4, 9, 15, 7919, 7919 ** 2, 3 * 5 * 7 * 11 * 13, 2 ** 20, 8191 * 131071,
             (2 ** 31 - 1), 65537 * 65539 if sympy.isprime(65539) else 65537 * 65543,
             101 * 103 * 107, 999983 * 2]
    for N in cases:
        r = ladder_factor(N, N.bit_length(), [2, 3])
        assert r["factors"] == _fac(N), N
    for N in range(2, 400):
        assert ladder_factor(N, max(4, N.bit_length()), [2])["factors"] == _fac(N)


def test_stage_cost_envelope():
    from src.avgtime.exponents import min_stage_cost, stage_cost_closed
    for e in [F(1, 100), F(1, 50), F(1, 46), F(3, 110), F(1, 20), F(1, 14), F(1, 13), F(1, 12)]:
        assert min_stage_cost(e) == stage_cost_closed(e), e
    # reachability: C(e) - e < z with e >= 1/5 - z  =>  z = 19/110 at e = 3/110
    e = F(3, 110)
    assert stage_cost_closed(e) - e == LADDER_NO_DETECTOR and F(1, 5) - e == LADDER_NO_DETECTOR


def test_dyadic_recount():
    from src.avgtime.exponents import (LADDER_NO_DETECTOR_DYADIC, dprime_dyadic_sup,
                                       min_stage_cost, stage_cost_closed_dyadic)
    assert dprime_dyadic_sup() == F(1, 40)
    for e in [F(1, 60), F(1, 40), F(3, 100), F(1, 25), F(1, 20)]:
        assert min_stage_cost(e, dyadic=True) == stage_cost_closed_dyadic(e)
    e = F(3, 100)
    assert stage_cost_closed_dyadic(e) - e == LADDER_NO_DETECTOR_DYADIC == F(1, 5) - e
    vals = [ladder(k, 1, False, dyadic=True)[0] for k in range(1, 6)]
    assert all(a > b for a, b in zip(vals, vals[1:]))
    assert all(v > LADDER_NO_DETECTOR_DYADIC for v in vals)
    assert vals[-1] < LADDER_NO_DETECTOR          # strictly better than one-threshold accounting
    # with a lambda^-2 tail no detector is needed: balanced cap 1/26 > 1/30
    assert ladder(6, 2, False, dyadic=True)[0] - ONE_SIXTH < F(1, 1000)
