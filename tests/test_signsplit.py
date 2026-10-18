"""Tests for src/signsplit (agent signsplit).  Instances used here (28/36 bits,
seeds >= 100) are disjoint from the pre-registered experiment set."""
from math import gcd

import pytest

from src.common.instances import balanced_semiprime
from src.common.metrics import Metrics
from src.signsplit import accounting as acc
from src.signsplit.hh import (hh_search, prop33, square_classes, v2_counts, v2_search,
                              weight_lattice)
from src.signsplit.identity import (all_eps, carry, carry_range, factor_squarefree,
                                    sigma_of_eps, sign_decomposition, sqrt_mod_prime, tau)
from src.signsplit.run import coverage_shared, coverage_v2
from src.signsplit.verify_identity import run as verify_run


def test_identity_exhaustive_small():
    rec = verify_run(61)
    assert all(v == 0 for v in rec["failures"].values())
    assert rec["counts"]["NEG_varying"] == rec["counts"]["NEG_classes"] > 0


@pytest.mark.parametrize("N", [7909787, 1000003 * 1000033, 4294967311 * 4294967357])
def test_signed_subset_sum_large(N):
    M = 2 * 3 * 5 * 7 * 11 * 13
    if gcd(M, N) != 1:
        pytest.skip("not coprime")
    primes = factor_squarefree(M)
    for gm in list(square_classes(N, M))[:5]:
        rr = []
        for m in primes:
            roots = sqrt_mod_prime((N * pow(gm, -1, m)) % m, m)
            rr.append(min(r for r in roots if r % m) if m > 2 else 1)
        a = 12345
        b = (a * gm) % M + 7 * M
        dec = sign_decomposition(a, b, rr, primes, N)
        lo, hi = carry_range(dec)
        for eps in all_eps(len(dec["Delta"])):
            s = sigma_of_eps(eps, rr, primes)
            S = dec["Delta0"] + sum(e * d for e, d in zip(eps, dec["Delta"]))
            assert (S - tau(a, b, s, N, M)) % (M * M) == 0
            assert lo <= carry(eps, dec) <= hi


def test_anti_lattice_index_and_value():
    N = 2 ** 36 + 2 ** 18 + 1
    while gcd(N, 2310) != 1:
        N += 2
    M = 2310
    primes = factor_squarefree(M)
    for s in (1, 13, 101, 997):
        if gcd(s, M) != 1:
            continue
        A0, g, Q = weight_lattice(N, primes, s, (3, 7))
        assert A0 * Q == M * 21 * 21
        # every lattice vector gives tau == 0 mod 9 and mod 49, for both signs
        for k in range(1, 6):
            a, b = A0 * k, (g * A0 * k) % Q - Q
            for sg in (s, M - s):
                t = tau(a, b, sg, N, M)
                assert t % 9 == 0 and t % 49 == 0


def test_prop33_bound_and_v1_equals_v0_lists():
    inst = balanced_semiprime(36, 101)
    N = inst.N
    for M, m0 in ((30, 3), (210, 1)):
        r0 = hh_search(N, M, m0, "V0", count_only=True)["met"]
        r1 = hh_search(N, M, m0, "V1", count_only=True)["met"]
        assert (r0["giant"], r0["baby"]) == (r1["giant"], r1["baby"])
        t = sum(1 for m in factor_squarefree(M) if m % 2)
        assert r0["groups"] * (1 << t) == r0["giant"]


@pytest.mark.parametrize("seed", [100, 101, 102])
def test_operational_variants_factor_28bit(seed):
    inst = balanced_semiprime(28, seed)
    N, p = inst.N, inst.labels["p"]
    for M in (30, 210):
        m0 = acc.best_shared(N, M, "V0")["m0"]
        for var, A in (("V0", ()), ("V1", ()), ("V3", (3,)), ("V3", (3, 5))):
            res = hh_search(N, M, m0, var, anti=A, seed=0)
            assert res["factor"] in (p, N // p) or res["status"].startswith("fail_")
            assert coverage_shared(N, p, M, m0, A if var == "V3" else ())["covered"]
        res = v2_search(N, M, m0, 1, 1)
        assert res["factor"] in (p, N // p)
        assert coverage_v2(N, p, M, m0, 1, 1)["covered"]


def test_v2_counts_match_operational():
    inst = balanced_semiprime(28, 103)
    N = inst.N
    c = v2_counts(N, 210, 2, 1, 2)
    r = v2_search(N, 210, 2, 1, 2)
    if r["status"] != "gcd_base":
        assert r["met"]["left"] + r["met"]["right"] == c["work"]


def test_merits():
    for row in acc.merit_table():
        assert row["anti_over_lock"] > 1.0          # sign collapse loses to the HH lock
        assert row["pm_offsets"] > row["lock_HH"]
        assert abs(row["anti_over_lock"] - row["anti_over_lock_formula"]) < 1e-12


def test_v0_counts_match_cycloharvey_calibrated_hh():
    from src.cycloharvey.calibration import cost_hh
    for bits in (32, 64, 96):
        N = balanced_semiprime(bits, 1000).N
        for M in (6, 30, 210, 2310):
            for m0 in (1, 3, 17, 250):
                a = cost_hh(N, M, m0)
                b = hh_search(N, M, m0, "V0", count_only=True)["met"]
                assert (a["giant"], a["baby"]) == (b["giant"], b["baby"])


def test_fixedbase_and_table_gray_agree_with_naive():
    inst = balanced_semiprime(36, 104)
    N, p = inst.N, inst.labels["p"]
    for M in (30, 210):
        m0 = acc.best_shared(N, M, "V0")["m0"]
        for var in ("V0", "V1"):
            rn = hh_search(N, M, m0, var, expo="naive")
            rf = hh_search(N, M, m0, var, expo="fixed")
            assert rn["status"] == rf["status"] and rn["factor"] == rf["factor"]
            assert rf["met"]["counts"]["rn_mul"] < rn["met"]["counts"]["rn_mul"]


def test_two_adic_shift_identity():
    from src.signsplit.verify_two_adic import run
    st = run(nmax=60, kmax=5, amax=5)
    assert all(ok for (_, _, ok) in st)
