"""Tests for src/repair/verify_r4.py: exhaustive check of the GPT-6.1-sol r4 claims for odd primes below 300."""
import random

import pytest

from src.repair import verify_r4 as V


# --- the class-group machinery is validated independently of the claims ---------------------------------
@pytest.mark.parametrize("D", [-20, -23, -56, -84, -104, -260, -420, -680, -884, -1768, -4004, -8 * 221])
def test_class_number_matches_analytic_formula(D):
    assert len(V.reduced_forms(D)) == V.class_number_analytic(D)


@pytest.mark.parametrize("D", [-884, -680, -8008, -4 * 3 * 7 * 11 * 13 * 17])
def test_composition_is_an_abelian_group_law(D):
    F = V.reduced_forms(D)
    S = set(F)
    one = V.identity_form(D)
    rnd = random.Random(0)
    assert all(V.compose(f, one) == f for f in F)
    assert all(V.compose(f, (f[0], -f[1], f[2])) == one for f in F)
    for _ in range(200):
        a, b, c = (rnd.choice(F) for _ in range(3))
        assert V.compose(a, b) in S
        assert V.compose(a, b) == V.compose(b, a)
        assert V.compose(V.compose(a, b), c) == V.compose(a, V.compose(b, c))


def test_torsion_counts_on_known_groups():
    # Cl(-84) = (Z/2)^2: r2 = 2, r4 = 0.   Cl(-260) = Z/2 x Z/4: r2 = 2, r4 = 1.   Cl(-56) = Z/4: r2 = 1, r4 = 1.
    assert (V.two_four_torsion(-84)["r2"], V.two_four_torsion(-84)["r4"]) == (2, 0)
    assert (V.two_four_torsion(-260)["r2"], V.two_four_torsion(-260)["r4"]) == (2, 1)
    assert (V.two_four_torsion(-56)["r2"], V.two_four_torsion(-56)["r4"]) == (1, 1)


def test_two_is_square_mod_N_by_brute_force():
    P = V.odd_primes_below(60)
    for i, p in enumerate(P):
        for q in P[i + 1:]:
            N = p * q
            brute = any(x * x % N == 2 for x in range(N))
            assert brute == (V.legendre(2, p) == 1 and V.legendre(2, q) == 1)


# --- the r4 claims, all pairs of distinct odd primes below 300 ---------------------------------------------
@pytest.fixture(scope="module")
def result():
    return V.run(300)


def test_scope_counts(result):
    assert result["odd_primes"] == 61 and result["pairs_total"] == 1830
    assert result["pairs_p_eq_q_mod4"] == 902
    assert result["pairs_N_1_mod_8"] == 442 and result["pairs_N_5_mod_8"] == 460


@pytest.mark.parametrize("claim", ["R4_1_identity", "R4_2_table16", "R4_3_mod32_formula", "R4_3_bp_eq_bq",
                                   "R4_3_bit_is_2_square", "R4_4_criterion", "X1_bp_ne_bq", "X1_mod32_no_new_bit",
                                   "X1_mod64_formula", "X2_T_square_iff_2_square", "X3_matrix_minus4N",
                                   "X3_matrix_minus8N", "X4_milovic", "X5_r2_minus4N", "X5_r2_minus8N"])
def test_claim_holds_on_every_pair(result, claim):
    s = result["summary"][claim]
    assert s["checked"] > 0 and s["fail"] == 0, s


def test_discriminant_choice_matters(result):
    # The criterion is false with the other discriminant on most pairs: the eps-dependent choice is essential.
    s = result["summary"]["X6_wrong_disc_agrees"]
    assert s["checked"] == 442 and s["fail"] > 0
