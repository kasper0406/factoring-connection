"""Tests for src/repair/verify_r5.py: two computational claims of the GPT-6.1-sol report review (r5)."""
import pytest

from src.repair import verify_r5 as R
from src.repair.c2r_search import brute_fold


@pytest.fixture(scope="module")
def item3():
    return R.item3(300)


def test_uv_mod4_formula_all_pairs(item3):
    # uv mod 4 = b(1 + 2[N = 1 mod 16]) for every p < q < 300 with N = 1 (mod 8)
    assert item3["pairs_N_1_mod_8"] == 442
    assert item3["uv_claim_fail"] == 0 and item3["s64_formula_fail"] == 0


def test_s_mod_64_determined_without_class_groups(item3):
    assert item3["s64_given_N64_eps_QR2"]["ambiguous_keys"] == 0


def test_s_mod_128_not_determined_by_two_power_ranks(item3):
    # Diagnostic finding: (N mod 128, eps, QR_2, r4/r8 of Cl(-4N), Cl(-8N)) does not determine s mod 128.
    assert item3["s128_given_N128_eps_QR2"]["ambiguous_keys"] > 0
    assert item3["s128_given_N128_eps_QR2_r4r8"]["ambiguous_keys"] > 0


def test_torsion_ranks_on_known_group():
    # Cl(-56) = Z/4: r2 = 1, r4 = 1, r8 = 0.   Cl(-260) = Z/2 x Z/4.
    assert R.two_power_torsion(-56) == {"r2": 1, "r4": 1, "r8": 0}
    assert R.two_power_torsion(-260) == {"r2": 2, "r4": 1, "r8": 0}
    # Cl(-164) has h = 8 and an element of order 8 (cyclic): r8 = 1.
    assert len(R.V.reduced_forms(-164)) == 8
    assert R.two_power_torsion(-164) == {"r2": 1, "r4": 1, "r8": 1}


def test_item10_noncharacter_exact_fold():
    f = R.fold(257, 8, 70, 186)
    assert f["exact"] and f["d"] == 7 and f["one_layer_degree"] == 9
    assert f["gcd_E_pm1"] == 2 and not f["character_exponent"] and f["twist_orders_le_12"] == []
    assert not f["passes_3L8"]  # 7 > 3L/8 = 6.375: compression, but below the pre-registered threshold
    assert brute_fold(257, 70, 186, 8) == (True, 7)  # independent implementation (agent `repair`)
