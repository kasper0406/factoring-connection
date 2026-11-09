"""Tests for agent midk (intermediate-k regime). Exact arithmetic."""
from fractions import Fraction as Fr

import sympy as sp

from src.midk import independence as I
from src.midk import ledger as L


def test_family_A_fallback_quarter_and_Gr_third():
    assert L.family_A_alpha("fallback") == Fr(1, 4)
    assert L.family_A_alpha("Gr") == Fr(1, 3)


def test_amgm_floor_is_half_power_no_lower():
    # min over d of x/d + d*c is 2 sqrt(c x): exponent of x is exactly 1/2 (a=1/2), so alpha>=1/4.
    c = sp.Symbol("c", positive=True)
    lo, dstar = L.min_logX_over_d(1, c)
    assert sp.simplify(lo - 2 * sp.sqrt(c) * sp.sqrt(L.x)) == 0
    # for c >= 1, lo >= 2 sqrt(x): no (d,r) gives a < 1/2, hence no alpha < 1/4.
    assert sp.simplify(lo.subs(c, 1) - 2 * sp.sqrt(L.x)) == 0


def test_aggregate_attaining_scaling_u_is_quarter():
    sc = L.aggregate_optimum_scaling()
    assert sc["alpha"] == Fr(1, 4)
    assert sc["u"][0] == Fr(1, 4)          # u ~ x^{1/4} -> smooth at reasonable cost
    assert sc["logX"] == (Fr(1, 2), Fr(1, 2))  # log X ~ sqrt(x log x)


def test_Gr_and_KL_floors_nondecreasing_in_r():
    for d in (6, 10, 20):
        gr = [v for _, v in L.gr_floor_monotone(d)]
        kl = [v for _, v in L.kl_floor_monotone(d)]
        assert all(gr[i] <= gr[i + 1] for i in range(len(gr) - 1))
        assert all(kl[i] <= kl[i + 1] for i in range(len(kl) - 1))
        assert gr[0] == kl[0]  # r=2 agrees (both the rank-2 floor)


def test_ladder_decreasing_above_third_inf_third():
    tab = [a for _, a in L.ladder_table(12)]
    assert all(tab[i] > tab[i + 1] for i in range(len(tab) - 1))
    assert all(a > Fr(1, 3) for a in tab)
    assert L.ladder_alpha(1) == Fr(2, 5)
    # limit of (j+1)/(3j+2) is 1/3
    j = sp.Symbol("j", positive=True)
    assert sp.limit((j + 1) / (3 * j + 2), j, sp.oo) == Fr(1, 3)


def test_floorsum_deterministic_metric():
    # N^{1/(j+2)}; endpoint j=2 is Coppersmith N^{1/4}; this is not the L[1/3] metric.
    assert L.floorsum_exponent(2) == Fr(1, 4)
    assert L.floorsum_exponent(3) == Fr(1, 5)


# ---- Lemma SPL-INT: dependence that lives only mod N splits N; global does not. ----

def test_rank_discrepant_intertwiner_splits_N():
    # a genuine rank discrepancy across primes (not an entrywise-zero matrix) forces a non-unit pivot
    p, q = 1000003, 1000033
    N = p * q
    for m in (2, 3, 4):
        A, lam = I.rank_discrepant_intertwiner(p, q, m, seed=m)
        fac = I.split_from_intertwiner(N, A, lam)
        assert fac is not None and 1 < fac < N and N % fac == 0
        assert fac in (p, q)


def test_nonscalar_eigenvector_does_not_split():
    # non-scalar is NOT sufficient: [[lam,1],[0,mu]] has a kernel for lam, no split (Opus critic)
    p, q = 99991, 99989
    N = p * q
    T = [[7, 1], [0, 11]]
    assert I.split_from_intertwiner(N, T, lam=7) is None


def test_outcome0_lambda_not_eigenvalue_no_split():
    # codex midk_2: N=21, T=[[0,-1],[1,0]], det(T-lam I)=lam^2+1 is a unit for all lam -> no split
    for lam in range(21):
        assert I.split_from_intertwiner(21, [[0, -1], [1, 0]], lam) is None


def test_global_scalar_symmetry_does_not_split():
    p, q = 99991, 99989
    N = p * q
    for m in (2, 3, 4):
        Ts = I.scalar_intertwiner(N, m, lam=7)
        # a true scalar: no eigen-computation can split (consistent global eigenvector)
        assert I.split_from_intertwiner(N, Ts, lam=7) is None


def test_gauss_kernel_or_factor_dichotomy():
    # unit pivots -> kernel; a non-unit pivot -> factor
    N = 15
    A = [[1, 2], [2, 4]]  # rank 1 mod 15 with unit pivot
    tag, val = I.gauss_kernel_or_factor(A, N)
    assert tag in ("kernel", "injective")
    B = [[3, 0], [0, 1]]  # pivot 3 shares gcd with 15
    tag2, val2 = I.gauss_kernel_or_factor(B, N)
    assert tag2 == "factor" and val2 == 3


def test_automorphism_graph_dependent_but_useless():
    # dependent links at height 1 (favourable failure of the counting half of G_r) but diagonal/useless
    w = I.automorphism_graph_witness(5, 31, 41, r=3)
    assert w["links_dependent_height1"] is True
    assert w["diagonal_useless"] is True


def test_independence_baseline_poisson_band():
    # the r-1 linking conditions are independent: rate ~ p^{-(r-1)}, within 4 sigma
    res = I.linking_independence_count(1009, 2, trials=100000, seed=1)
    assert abs(res["z"]) < 4.0
