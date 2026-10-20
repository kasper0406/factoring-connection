"""Tests for src/snfsrig (notes/snfsrig_results.md)."""
from fractions import Fraction as Fr

from src.snfsrig import cyclo as C
from src.snfsrig import rigidity as R
from src.snfsrig import sizemodel as S


# ---------------- size model ----------------

def test_one_level_floor_values():
    assert S.one_level_floor(1) == Fr(2, 5)          # any pair with a linear side (RES: theta_coeff >= 1/d)
    assert S.one_level_floor(0) == Fr(1, 2)          # coefficient term not decaying: fixed degree
    assert S.one_level_floor(10 ** 6) > Fr(1, 3)     # -> 1/3 only if the coefficient term vanished
    assert S.one_level_floor(Fr(1, 2)) == Fr(3, 7)


def test_one_level_floor_is_the_lp_minimum():
    # brute force over a rational grid: feasible iff alpha >= floor
    for p in (Fr(1), Fr(2), Fr(1, 2)):
        fl = S.one_level_floor(p)
        al, x, y = S.one_level_optimum(p)
        assert S.one_level_feasible(al, p, x, y)
        grid = [Fr(i, 60) for i in range(0, 61)]
        for alpha in grid:
            feas = any(S.one_level_feasible(alpha, p, x, y) for x in grid for y in grid)
            if alpha < fl:
                assert not feas
            if alpha >= fl and alpha in grid and (1 - 2 * alpha) / p in grid:
                assert feas


def test_theta0_minima():
    for kappa in (Fr(1, 10), Fr(3, 20), Fr(1)):
        # SNFS minimum 1/sqrt(kappa) at d0 = 2 sqrt(kappa); check at a rational kappa = 4 c^2
        pass
    kappa = Fr(1, 4)                                  # sqrt(kappa) = 1/2
    assert S.theta0("SNFS", 1, kappa) == 2            # 1/1 + 1/(4/4) = 2 = 1/sqrt(1/4)
    assert S.theta0("SNFS", 1, kappa) ** 2 == S.theta0_min_squared("SNFS", kappa)
    kappa = Fr(1, 2)                                  # sqrt(8 kappa) = 2
    assert S.theta0("GNFS_heur", 2, kappa) ** 2 == S.theta0_min_squared("GNFS_heur", kappa)
    assert S.theta0("U", 2, kappa) == S.theta0("GNFS_heur", 2, kappa)   # Kern regime: kappa < d0^2/4
    assert S.theta0("U", 1, 1) == 3                   # level+window regime


def test_cost_reproduces_theorem_U_and_snfs_is_lower():
    gk = [Fr(3, 20)]
    gd = [Fr(21, 20)]
    ga = [Fr(11, 20)]
    u = S.optimise("U", 0, grid_a=ga, grid_kappa=gk, grid_d0=gd)
    assert abs(float(u["c"]) - 3.5103) < 1e-3       # lv5unif C_gamma = 0 value (grid)
    sn = S.cost_constant("SNFS", Fr(1, 2), Fr(13, 20), Fr(3, 25), 0)
    assert sn < u["c"]


def test_multilevel_levels():
    assert S.multilevel_levels(Fr(2, 5)) == 1
    assert S.multilevel_levels(Fr(38, 100)) == 2
    assert S.multilevel_levels(Fr(35, 100)) == 6


def test_r2_gnfs_vs_snfs():
    assert S.r2_window_nonempty_exponent(Fr(39, 100))      # GNFS: blocked below 2/5
    assert not S.r2_window_nonempty_exponent(Fr(2, 5))
    assert S.r2s_contraction_exponent(Fr(1, 3)) == 0         # SNFS: constant contraction at 1/3
    assert S.r2s_contraction_exponent(Fr(2, 5)) < 0
    assert S.r2s_contraction_exponent(Fr(3, 10)) > 0         # expansion below 1/3
    ex = S.r2s_alpha_third(b=4, s=Fr(3, 4), phiF=1, phiR=1, b1=1)
    assert ex["contracts"] and ex["heuristic_ok"] and ex["d1_ge_1"]
    assert ex["lam"] == Fr(11, 64)          # (1/2 + 3/16)/4


# ---------------- rigidity ----------------

def test_res_floor_snfs_pair():
    for (r, e, s, d) in [(2, 31, 3, 3), (3, 25, 2, 4), (2, 41, 5, 5)]:
        f0, m, N, j = R.snfs_pair(r, e, s, d)
        out = R.check_res_floor(f0, [-m, 1], m, N)
        assert out["divisible"] and out["nonzero"] and out["hadamard_ok"] and out["abs_ge_N"]


def test_height_rigidity_exhaustive():
    f0, m, N, j = R.snfs_pair(2, 30, 3, 3)            # f0 = x^3 - 3, m = 2^10
    assert j == 0 and f0 == [-3, 0, 0, 1] and R.is_irreducible(f0)
    out = R.exhaustive_small_height(f0, m, N, 3, 3)
    assert out["hr_guarantees_none"]
    assert out["nondivisible"] == []
    assert out["found"] == 2                          # +- f0
    # small N: HR gives no guarantee, non-multiples appear, all consistent with the bound
    f0, m, N, j = R.snfs_pair(2, 9, 3, 3)             # N = 509, m = 8
    out = R.exhaustive_small_height(f0, m, N, 3, 2)
    assert not out["hr_guarantees_none"]
    assert out["consistent_with_HR"]


def test_res_floor_on_nonmultiples():
    f0, m, N, j = R.snfs_pair(2, 9, 3, 3)
    out = R.exhaustive_small_height(f0, m, N, 3, 2)
    for f in out["nondivisible"][:50]:
        chk = R.check_res_floor(list(f) if f[-1] != 0 else list(f), f0, m, N) if f[-1] != 0 else None
        if chk is not None:
            assert chk["divisible"] and chk["abs_ge_N"]


def test_rational_side_rigidity():
    f0, m, N, j = R.snfs_pair(2, 26, 3, 2)
    Hg = N // (2 * m + 2) - 1
    out = R.check_rs(m, N, min(Hg, 20000))
    assert out["hypothesis"] and out["non_multiples"] == []
    # beyond the threshold non-multiples exist: (N - m*g1 ... ) e.g. g = (N - m, -1)? use a small case
    out2 = R.check_rs(5, 37, 10)
    assert not out2["hypothesis"] and out2["non_multiples"]


def test_device_invariance():
    for (r, e, s, d) in [(2, 31, 3, 3), (3, 25, 2, 4)]:
        f0, m, N, j = R.snfs_pair(r, e, s, d)
        ids = R.device_identities(f0, m)
        assert all(v == 0 for v in ids.values())


def test_value_set_bound():
    assert R.value_set_bound(3, 100, 400) == 240000


# ---------------- cyclotomic pieces ----------------

def test_algebraic_pieces_products():
    for (r, e, sg) in [(2, 30, -1), (2, 30, 1), (3, 12, -1), (10, 15, 1)]:
        pcs = C.algebraic_pieces(r, e, sg)
        p = 1
        for v in pcs.values():
            p *= v
        assert p == r ** e + sg


def test_cyclotomic_values():
    assert C.cyclotomic_value(1, 2) == 1
    assert C.cyclotomic_value(6, 2) == 3
    assert C.cyclotomic_value(12, 2) == 13
    assert C.euler_phi_small(2310) == 480
    assert C.mobius_small(30) == -1 and C.mobius_small(12) == 0


def test_aurifeuillian():
    L, M = C.aurifeuillian_2(30)
    assert L * M == 2 ** 30 + 1
    pcs = C.special_pieces(2, 30, 1)
    p = 1
    for v in pcs:
        p *= v
    assert p == 2 ** 30 + 1
    assert len(pcs) > len([v for v in C.algebraic_pieces(2, 30, 1).values() if v > 1])


def test_primorial_largest_piece_shrinks():
    ratios = [C.largest_piece_log_ratio(2, C.primorial(k), -1) for k in (3, 4, 5)]
    assert ratios[0] > ratios[1] > ratios[2]
    assert ratios[2] < Fr(22, 100)                    # phi(2310)/2310 = 0.2078
    lo, hi = C.cost_exponent_ratio(ratios[2])
    assert lo < Fr(54, 100) < hi + Fr(1, 100)
