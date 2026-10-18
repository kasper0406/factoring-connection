"""Tests for src/tower (agent tower128).  All diagnostic: exhaustive checks on small ranges."""
from __future__ import annotations

import pytest

from src.repair import verify_r4 as V
from src.repair import verify_r5 as V5
from src.tower import formula as F
from src.tower.analyze import ambiguity, depth_rule
from src.tower.classgroup import rank, reduced_forms_np, two_sylow
from src.tower.run_tower import ambiguous_prime_form, invariants, pairs
from src.tower.symbols import power_symbol, psi_pm, s2_single


def test_formula_exhaustive_small():
    res = F.run(400)
    assert res["failures"] == {}
    assert res["T1_iii_ambiguous_keys"] == 0          # (N mod 128, eps, psi(p)) determines s mod 128
    assert res["pub_ambiguous_keys"] > 0              # (N mod 128, eps, QR2) does not
    assert all(v["ambiguous_keys"] == 0 for v in res["T2"].values())
    assert res["T1_v_N5mod8_by_omega"]["ambiguous_keys"] == 0


def test_index32_and_characters():
    for x in range(1, 200, 2):
        sgn, a = F.index32(x)
        assert (sgn * pow(5, a, 32) - x) % 32 == 0
        assert F.psi_exp(x) == a % 4
    assert F.psi_exp(5) == 1 and F.psi_exp(9) == 2 and F.psi_exp(7) == 2 and F.psi_exp(15) == 0


def test_reduced_forms_match_bruteforce():
    for D in [-3, -4, -7, -8, -15, -20, -23, -164, -820, -1012, -4 * 17 * 41, -8 * 17 * 41, -3 * 17 * 41]:
        assert set(map(tuple, reduced_forms_np(D).tolist())) == set(V.reduced_forms(D))


@pytest.mark.parametrize("D", [-164, -820, -4 * 17 * 41, -8 * 17 * 41, -12 * 17 * 41, -3 * 17 * 41, -28 * 73 * 89,
                               -7 * 73 * 89, -24 * 73 * 89, -20 * 73 * 89])
def test_two_sylow_matches_bruteforce(D):
    S = two_sylow(D)
    r = V5.two_power_torsion(D)
    assert (rank(S, 1), rank(S, 2), rank(S, 3)) == (r["r2"], r["r4"], r["r8"])
    assert S["h"] == len(V.reduced_forms(D))


def test_cl164_is_cyclic_8():
    assert two_sylow(-164)["type"] == [8]


def test_ambiguous_prime_form():
    N = 17 * 41
    assert ambiguous_prime_form(-4 * N, 2) == V.reduce_form((2, 2, (N + 1) // 2))
    assert ambiguous_prime_form(-8 * N, 2) == V.reduce_form((2, 0, N))
    assert ambiguous_prime_form(-12 * N, 2) is None   # not primitive in the order of conductor 2


def test_T3_T4_small():
    recs = [invariants(pq) for pq in pairs(500)]
    res = depth_rule(recs)
    assert res["T3_checked"] > 50 and res["T3_fail"] == 0 and res["T3_ambiguous_fail"] == 0
    assert res["T4_fail"] == 0 and res["T4p_fail"] == 0
    # (N mod 128, eps, QR2) alone is ambiguous on this range
    from src.tower.analyze import K0
    assert ambiguity(recs, K0)["ambiguous_keys"] > 0


def test_single_prime_prototypes():
    from src.repair.verify_r4 import odd_primes_below
    for p in odd_primes_below(3000):
        if p % 8 == 1:
            assert power_symbol(-4, p, 8) == power_symbol(2, p, 4) * psi_pm(p)
        if p % 8 in (1, 7):
            r = s2_single(p)
            assert all(v for k, v in r.items() if k.startswith("ok"))
