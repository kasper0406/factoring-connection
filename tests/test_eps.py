"""Tests for the eps/QR_2 battery (agent epsbattery).  Instances here use only small primes or the
'test' seed range; never the bench/main/fresh/supp ranges."""
import ast
import math
import os
import random
from math import comb, gcd

import gmpy2
import numpy as np
import pytest

from src.eps import cf as cfm
from src.eps.algebra import Nonunits, bm_mod_N, det_mod_N, mult_matrix, rank_mod_N, section_values, sketch
from src.eps.analysis import floor_baseline, pc_check, strat_z
from src.eps.forms import small_forms
from src.eps.harness import SEED_BASE, instance, labels
from src.eps.operational import MENU, MENU_NAMES, compute_stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIMES_MID = [x for x in range(257, 1024) if gmpy2.is_prime(x)]


def _pairs(limit_pairs=None, seed=1):
    pr = [(p, q) for i, p in enumerate(PRIMES_MID) for q in PRIMES_MID[i + 1:] if p % 4 == q % 4]
    if limit_pairs:
        random.Random(seed).shuffle(pr)
        pr = pr[:limit_pairs]
    return pr


# ---------------------------------------------------------------- CF recurrence and signals
def test_cf_convergent_identity():
    for N in (1001 * 7 + 6, 65537 * 3, 10403, 2 ** 31 - 1, 99991 * 100003):
        if gmpy2.is_square(N):
            continue
        seq = cfm.cf_sequence(N, 60)
        A_prev, A = 1, seq[0][0]
        B_prev, B = 0, 1
        for i in range(1, 61):
            # A_{i-1}^2 - N B_{i-1}^2 = (-1)^i Q_i
            assert A * A - N * B * B == (-1) ** i * seq[i][2]
            a_i = seq[i][0]
            A_prev, A = A, a_i * A + A_prev
            B_prev, B = B, a_i * B + B_prev


def test_period_detection_matches_pell():
    # N = m^2 + 1 has period 1 (odd) -> pc1 fires at i = 1; N = m^2 + 2 has period 2.
    f = cfm.cf_features(1000 ** 2 + 1, T_max=64, T_list=(64,))
    assert f["first_pc1"] == 1 and f["pc1_T64"] == 1
    f = cfm.cf_features(1000 ** 2 + 2, T_max=64, T_list=(64,))
    assert f["first_pc1"] == 0


@pytest.mark.parametrize("chunk", [0, 1])
def test_positive_controls_never_wrong_and_fire(chunk):
    """Exhaustive-ish: pairs of primes in (256, 1024) with p = q (mod 4); PC => label."""
    pairs = _pairs()[chunk::2][:6000]
    fires = {k: 0 for k in ("pc1", "pc2", "pc3", "pcq2", "pcq3")}
    for p, q in pairs:
        N = p * q
        f = cfm.cf_features(N, T_max=256, T_list=(256,))
        assert f["cf_fb_ok"] == 1 and f["cf_nonunit"] == 0
        eps = 1 if p % 4 == 1 else -1
        two_sq = (p % 8 in (1, 7)) and (q % 8 in (1, 7))
        for k in ("pc1", "pc2", "pc3"):
            if f[f"{k}_T256"]:
                fires[k] += 1
                assert eps == 1, (k, p, q)
        for k in ("pcq2", "pcq3"):
            if f[f"{k}_T256"]:
                fires[k] += 1
                assert two_sq, (k, p, q)
    assert all(v > 10 for v in fires.values()), fires


def test_negative_pell_iff_period_odd_small():
    # brute-force check of PC1 against direct search of x^2 - N y^2 = -1 for tiny N with short periods
    for p in (5, 13, 17, 29, 37, 41):
        for q in (53, 61, 73, 89, 97, 101):
            N = p * q
            f = cfm.cf_features(N, T_max=200, T_list=(200,))
            if f["pc1_T200"]:
                l = f["first_pc1"]
                seq = cfm.cf_sequence(N, l)
                A_prev, A, B_prev, B = 1, seq[0][0], 0, 1
                for i in range(1, l):
                    A_prev, A = A, seq[i][0] * A + A_prev
                    B_prev, B = B, seq[i][0] * B + B_prev
                assert A * A - N * B * B == -1


# ---------------------------------------------------------------- algebra
def test_section_values_naive():
    for N in (7 * 11 * 13, 221, 1003):
        for M, c in ((1, 2), (2, N - 1), (4, 2), (8, 3)):
            got = section_values(N, M, c)
            want = [0] * M
            for k in range(1, N):
                want[k % M] = (want[k % M] + comb(N, k) * pow(c, k // M, N)) % N
            assert got == want


def test_sketch_naive():
    N = 1003
    f = [-2, 0, 0, 1]
    S = sketch(N, f)
    # naive: F_N mod (X^3 - 2)
    want = [0, 0, 0]
    for k in range(1, N):
        want[k % 3] = (want[k % 3] + comb(N, k) * pow(2, k // 3, N)) % N
    assert S == want


def test_bm_and_rank():
    N = 1000003 * 999983
    nu = Nonunits(N)
    fib = [0, 1]
    for _ in range(62):
        fib.append((fib[-1] + fib[-2]) % N)
    L, C = bm_mod_N(fib, N, nu, "t")
    assert L == 2 and C[:3] == [1, N - 1, N - 1]
    rnd = random.Random(5)
    L, _ = bm_mod_N([rnd.randrange(N) for _ in range(64)], N, nu, "t")
    assert L == 32
    assert rank_mod_N([[1, 2], [2, 4]], N, nu, "t") == 1
    M = mult_matrix(N, [-2, 0, 0, 1], [1, 1, 0])     # 1 + X in Z/N[X]/(X^3-2): det = 1 + 2 = 3
    assert det_mod_N(M, N, nu, "t") == 3
    assert nu.events == []


def test_nonunit_recorded():
    N = 101 * 103
    nu = Nonunits(N)
    assert rank_mod_N([[101, 0], [0, 1]], N, nu, "t") == -1
    assert nu.events == [("t", 101)]


# ---------------------------------------------------------------- forms
def _h_analytic(D):
    # h(D) = -(1/|D|) sum_{a=1}^{|D|} (D/a) a  for fundamental D < -4
    s = sum(int(gmpy2.kronecker(D, a)) * a for a in range(1, -D))
    assert s % D == 0
    return -s // (-D)


def test_forms_class_number():
    for p, q in ((13, 17), (11, 19), (29, 37), (7, 23), (43, 47)):
        N = p * q
        F = small_forms(N, a_max=int(math.isqrt(4 * N // 3)) + 1)
        assert len(F) == _h_analytic(-4 * N)
        for a, b, c in F:
            assert b * b - 4 * a * c == -4 * N


# ---------------------------------------------------------------- operational isolation
def test_operational_modules_do_not_see_factors():
    bad_mods = {"src.common.instances", "src.eps.harness", "sympy", "src.trackE.oracles"}
    bad_attrs = {"is_prime", "factor", "factorint", "next_prime", "is_probab_prime", "is_bpsw_prp"}
    for fn in ("cf.py", "algebra.py", "forms.py", "operational.py"):
        tree = ast.parse(open(os.path.join(ROOT, "src", "eps", fn)).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                mod = node.module or ""
                assert mod not in bad_mods and "instances" not in mod and "harness" not in mod, (fn, mod)
            if isinstance(node, ast.Import):
                for al in node.names:
                    assert al.name not in bad_mods, (fn, al.name)
            if isinstance(node, ast.Attribute):
                assert node.attr not in bad_attrs, (fn, node.attr)
            if isinstance(node, ast.arg):
                assert node.arg not in ("p", "q", "phi", "factors", "labels"), (fn, node.arg)


def test_compute_stats_shape_and_determinism():
    for j in range(6):
        g = instance(32, SEED_BASE["test"] + j)
        if g is None:
            continue
        st, aux = compute_stats(g["N"])
        st2, _ = compute_stats(g["N"])
        assert list(st) == MENU_NAMES and st == st2
        assert all(isinstance(v, int) for v in st.values())
        assert st["Bx3m2_rank"] == 3
        lab = labels(g["p"], g["q"])
        if st["pc2_T1024"] or st["pc1_T1024"] or st["pc3_T1024"]:
            assert lab["eps"] == 1
    assert len(MENU) <= 100


def test_harness_labels_qr2():
    lab = labels(5, 13)        # N = 65 = 1 mod 8; (2/5) = (2/13) = -1
    assert lab == {"eps": 1, "qr2": 0}
    lab = labels(7, 23)        # 161 = 1 mod 8, (2/7) = (2/23) = +1
    assert lab == {"eps": -1, "qr2": 1}
    lab = labels(3, 7)         # 21 = 5 mod 8
    assert lab == {"eps": -1, "qr2": None}


# ---------------------------------------------------------------- floor removal and the test statistic
def _synthetic(n=40000, seed=0):
    rng = np.random.default_rng(seed)
    Nmod = (rng.integers(0, 1 << 14, size=n) * 4 + 1) % (1 << 16)    # N = 1 mod 4
    strata = Nmod % 256
    return rng, Nmod, strata


def test_floor_function_is_degenerate():
    rng, Nmod, strata = _synthetic()
    y = rng.integers(0, 2, size=len(Nmod))
    for x in (Nmod % 8, Nmod % 256, ((Nmod >> 3) & 1) * 5 - 2):
        r = strat_z(x.astype(float), y.astype(float), strata)
        assert r["degenerate"] and r["z"] == 0.0


def test_floor_confound_removed():
    """Label depends strongly on N mod 2^8; statistic = floor function + independent noise.
    Raw correlation is huge; the stratified z must be ~N(0,1)."""
    zs = []
    for seed in range(5):
        rng, Nmod, strata = _synthetic(seed=seed)
        prob = np.where((Nmod >> 4) & 1, 0.8, 0.2)
        y = (rng.random(len(Nmod)) < prob).astype(float)
        x = ((Nmod >> 4) & 1) * 3.0 + rng.normal(size=len(Nmod))
        r = strat_z(x, y, strata)
        raw = np.corrcoef(x, y)[0, 1] * math.sqrt(len(x) - 1)
        assert raw > 50
        zs.append(r["z"])
    assert max(abs(z) for z in zs) < 4.0


def test_real_signal_detected():
    rng, Nmod, strata = _synthetic(n=100000, seed=7)
    y = rng.integers(0, 2, size=len(Nmod)).astype(float)
    # statistic agreeing with y with prob 0.52 -> rho = 0.04 -> z ~ 12.6
    agree = rng.random(len(y)) < 0.52
    x = np.where(agree, y, 1 - y)
    r = strat_z(x, y, strata)
    assert 9 < r["z"] < 16


def test_null_calibration():
    rng, Nmod, strata = _synthetic(n=20000, seed=3)
    zs = []
    for _ in range(200):
        y = rng.integers(0, 2, size=len(Nmod)).astype(float)
        x = rng.integers(0, 3, size=len(Nmod)).astype(float)
        zs.append(strat_z(x, y, strata)["z"])
    m2 = float(np.mean(np.square(zs)))
    assert 0.75 < m2 < 1.3


def test_pc_check_one_sided():
    rng, Nmod, strata = _synthetic(n=50000, seed=11)
    y = rng.integers(0, 2, size=len(Nmod))
    x = ((rng.random(len(y)) < 0.02) & (y == 1)).astype(int)
    r = pc_check(x, y, strata)
    assert r["violations"] == 0 and r["consistent"] and r["detected"]


def test_floor_baseline_null():
    rng, Nmod, strata = _synthetic(n=40000, seed=2)
    y = rng.integers(0, 2, size=len(Nmod))
    for row in floor_baseline(Nmod, y, ks=(3, 6, 8)):
        assert abs(row["cv_acc"] - 0.5) < 0.02


def test_mr_ss_controls_never_wrong_and_fire():
    from src.eps.mr import mr_features
    pr = [x for x in range(5, 400) if gmpy2.is_prime(x)]
    fires4 = fires5 = 0
    for i, p in enumerate(pr):
        for q in pr[i + 1:]:
            if p % 4 != q % 4:
                continue
            f = mr_features(p * q)
            if f["pc4_mr"] or f["pc5_ss"]:
                assert p % 4 == 1, (p, q)
            assert f["pc5_ss"] <= f["pc4_mr"]
            fires4 += f["pc4_mr"]
            fires5 += f["pc5_ss"]
    assert fires4 > 10 and fires5 > 5, (fires4, fires5)
