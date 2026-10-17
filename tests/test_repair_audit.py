"""Tests for the r3-audit reproductions (agent `repair`).  Small parameters only."""
from fractions import Fraction

from src.repair import audit_trackA as A
from src.repair import audit_trackBD as B
from src.repair import audit_trackC as C


# ---------------- Track A ----------------
def test_selector_counterexample_N35():
    r = A.selector_counterexample()
    assert r["a_5"] == 7 and r["isolates_class_r_on_[0,N]"] and r["value_on_F_N"] == 7
    assert r["observable_rank"] == 6 < r["M"]


def test_corrected_bound_min_M_i_plus_1():
    for q, n, Ms in [(2, 12, (2, 3, 4, 5, 8)), (3, 9, (2, 3, 4))]:
        r = A.corrected_bound_exhaustive(q, n, Ms)
        assert r["violations"] == 0 and r["bound_attained_(delta_r: r+1; class indicator: M)"]


def test_balanced_support_starts_at_p():
    assert A.balanced_scoped_bound(3000)["violations"] == 0


def test_coinflip_rank_one_distinguishes():
    r = A.coinflip_example()
    assert r["distinguishes"] and not r["isolates_class_0_mod_2"] and not r["isolates_class_1_mod_2"]


def test_c2_scope_N35_M7():
    S = A.section(35, 7, 0)
    assert sorted(S) == [1, 2, 3, 4]  # indices 7, 14, 21, 28


def test_c3_gap_and_c4_N6():
    assert A.section(35, 8, 1) == {3: 21}
    S6 = A.section(6, 2, 0)
    assert S6 == {1: 3, 2: 3} and all(A.eval_section(S6, c, 6) == 0 for c in range(6))


def test_c4_repaired_small():
    r = A.c4_repair(nmax_eval=400, nmax_struct=20000)
    assert r["scalar_false_zeros"] == 0 and r["span_bound_violations"] == 0
    assert [3, 5] in [list(x) for x in r["printed_inequality_fails_(p,q)_first"]]


def test_s04_F51_trap():
    S = A.section(51, 4, 0)
    assert S == {3: 34, 6: 34, 9: 17, 12: 17}
    assert all(A.eval_section(S, c, 51) == 0 for c in range(51))


def test_mask_F5_example():
    r = A.mask_trailing_block()
    assert r["u_t=1-2^t, t=0..7"][:5] == [0, 4, 2, 3, 0] and r["linear_complexity"] == 2


def test_fitting_gap_conclusion_holds():
    r = A.fitting_gap_exhaustive(3, 3, 12)
    assert r["violations"] == 0 and r["blocks_starting_before_e"] > 0


def test_dimension_vs_cost():
    assert A.dimension_vs_cost()["small_check_dense_equals_fast"]


# ---------------- Track C ----------------
def test_o1_polynomial_readout_counterexample():
    r = C.o1_poly_readout(primes=(7, 11, 31))
    assert all(x["nonzero_1..p"] and x["zero_at_p+1,p+2"] for x in r["prefix_examples"])
    assert r["interval_example"]["zero_iff_A_mod_p_in_Z_L_for_1<=A<2p"]


def test_o2_strengthened_bound():
    r = C.o2_strengthening(pmax=40)
    assert r["violations_ceil_L_over_2"] == 0 and r["violations_g<=2"] == 0
    assert r["(p,L)_attaining_ceil_L_over_2"] == r["(p,L)_with_exact_fold"]


def test_o2prime_smooth_part():
    r = C.o2prime_consequences()
    assert r["gcd(lcm(1..B),p-1)"] == 8 and r["B_smooth_part"] == 256


def test_o4_point_order():
    r = C.o4_order_example()
    assert r["n*P=O"] and not r["#E divides n"]


# ---------------- Track B/D ----------------
def test_rank_vs_pivot_3_5():
    r = B.rank_vs_pivot()
    assert r["elim_result"]["factor"] == 3 and r["pivots_mod_3"] == [1] and r["pivots_mod_5"] == [0]
    assert r["krylov_violations"] == 0


def test_six_orbit_literal_vs_repaired():
    r = B.sporadic_six(pmax=60)
    assert r["literal_SPORADIC"] == [2, 4] and r["repaired_def_not_0_mod_6"] == 0
    assert r["literal_def_not_0_mod_6"] > 0


def test_n209_degree_dependence():
    assert B.pgap(209, 11, 19, 2) == 0 and B.pgap(209, 11, 19, 6) == Fraction(1, 9)


def test_n247_q_minus_i():
    r = B.n247()
    assert r["mod13"]["Q_N"] == [[1, 0, 0], [0, 12, 0], [0, 0, 1]]
    assert r["mod13"]["rank(Q_N - I)"] == 1 and r["mod19"]["rank(Q_N - I)"] == 2


def test_stirling_nonmonotone_391():
    g = dict(B.stirling_window(391)["gcds_K=1..19"])
    assert (g[12], g[13], g[14]) == (391, 17, 391)


def test_model_M_77():
    assert B.model_M_77()["elim_result"]["factor"] == 7


def test_bm_step_bound():
    r = B.bm_step_bound(bits=(16,), seeds=10)
    assert r["violations"] == 0
