"""Tests for src/mf2L (notes/mf2L_results.md): exact checks of the finite statements and the exponent margins."""
from fractions import Fraction as Fr
from math import gcd
import random

from src.mf2L import checks, params


def test_margins_and_kappa():
    for al in [Fr(3, 8), Fr(19, 50), Fr(7, 18), Fr(39, 100), Fr(399, 1000)]:
        assert params.all_ok(al), (al, params.margins(al))
        assert params.kappa(al) < 1
    assert params.kappa(Fr(3, 8)) == Fr(2, 3)
    # kappa -> 1 as alpha -> 2/5 (parent algebraic family): the residual input gets harder near 2/5
    assert params.kappa(Fr(2, 5) - Fr(1, 10**6)) > Fr(99, 100)
    # the documented lambda-gap: pointwise e^{-u} bounds plus the classical large sieve never suffice
    for v in params.informational(Fr(3, 8)).values():
        assert v.startswith("no")


def test_homogeneity_is_design_invariant():
    rng = random.Random(1)
    n = 0
    while n < 400:
        P1 = (rng.randint(-60, 60), rng.randint(1, 60))
        P2 = (rng.randint(-60, 60), rng.randint(1, 60))
        if gcd(*P1) != 1 or gcd(*P2) != 1 or P1[0] * P2[1] == P1[1] * P2[0]:
            continue
        e = rng.randint(1, 7)
        # an arbitrary integral form: think F = w*Fhat + G*R_c + (any extra randomisation), or G itself (e = 1)
        coeffs = [rng.randint(-10**6, 10**6) for _ in range(e + 1)]
        assert checks.homogeneity(coeffs, P1, P2)
        n += 1


def test_iota_average_is_exact_product():
    rng = random.Random(2)
    for k in [5, 7, 8, 9, 12, 15]:
        a = {rng.randint(1, 200): rng.randint(0, 3) for _ in range(40)}
        b = {rng.randint(1, 200): rng.randint(0, 3) for _ in range(40)}
        tot, prod = checks.iota_average(a, b, k)
        assert tot == prod


def test_dth_power_characters_survive():
    # d odd, k with d | phi(k): characters of order dividing d are trivial on lam^d, so they never cancel over pairs
    assert checks.dth_power_kernel(7, 3)["chars_order_dividing_d"] == 3
    assert checks.dth_power_kernel(31, 5)["chars_order_dividing_d"] == 5
    assert checks.dth_power_kernel(11, 3)["chars_order_dividing_d"] == 1     # gcd(3, 10) = 1
    for p in [3, 5, 7, 11, 13]:
        for d in [1, 3, 5, 7]:
            assert checks.quadratic_odd_power(p, d)


def test_interval_densities_do_not_determine_pair_counts():
    for k in (3, 4, 6):
        out = checks.nid_example(k, 4 * k + 1, Fr(1, 2))
        assert out["max_interval_density_gap"] <= Fr(1, 2) / k * 2      # O(eps/len) on every interval >= k
        assert out["pair_t"] > out["pair_s"]                              # the lattice pair count differs


def test_gauss_transfer_exact():
    rng = random.Random(3)
    for k in (5, 7, 9):
        a = {n: rng.randint(-2, 3) for n in rng.sample(range(1, 60), 12)}
        out = checks.gauss_transfer(k, a)
        assert out["primitive_identity"] and out["orthogonality"]


def test_fixed_iota_phase_is_conjugate():
    # codex mf2L_2's counterexample to v1's chi(iota): k = 5, iota = 2, a = delta_2, b = delta_1
    assert checks.fixed_iota_formula(5, 2, {2: 1}, {1: 1})
    rng = random.Random(4)
    for k in (5, 7, 9):
        for _ in range(3):
            iota = rng.choice([i for i in range(1, k) if gcd(i, k) == 1])
            a = {rng.randint(1, 80): rng.randint(0, 3) for _ in range(10)}
            b = {rng.randint(1, 80): rng.randint(0, 3) for _ in range(10)}
            assert checks.fixed_iota_formula(k, iota, a, b)
