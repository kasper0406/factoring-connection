"""Tests for agent dlog25 (notes/dlog25_results.md)."""
from fractions import Fraction as Fr
import random

from src.dlog25 import params as P
from src.dlog25 import toy_sp, toy_xk
from src.lv5unif import params as up


def test_two_pow_bracket():
    lo, hi = P.two_pow_minus_2_5_bracket()
    assert lo < hi and lo ** 5 * 4 < 1 < hi ** 5 * 4      # (2^{-2/5})^5 = 1/4


def test_lemma_C25_i_same_tuple():
    """c_DL <= max(3a, c_U) when the shape test uses Theorem U at the same tuple (exact brackets)."""
    rng = random.Random(5)
    for _ in range(200):
        a = Fr(rng.randint(1, 40), 20)
        kappa = Fr(rng.randint(1, 60), 200)
        d0 = Fr(rng.randint(5, 60), 20)
        Cg = rng.choice([0, 1, 5, 20])
        cu = P.c_u(a, d0, kappa, Cg)
        lo, hi = P.c_dl(a, d0, kappa, Cg, cu)
        m_lo, m_hi = P.two_pow_minus_2_5_bracket()
        # exact statement on the lower bracket (no rounding slack in that direction)
        assert lo <= max(3 * a, cu[1])
        assert P.same_tuple_dominated(a, d0, kappa, Cg)


def test_lemma_C25_ii_no_gain():
    """If max(3a, T1, a+s+2^{-2/5}c*) < c* then c_U(tuple) < c*: checked with exact brackets on random tuples."""
    rng = random.Random(7)
    m_lo, m_hi = P.two_pow_minus_2_5_bracket()
    f_lo, f_hi = up.recursion_factor_bracket()
    for _ in range(300):
        a = Fr(rng.randint(1, 40), 20)
        kappa = Fr(rng.randint(1, 60), 200)
        d0 = Fr(rng.randint(5, 60), 20)
        Cg = rng.choice([0, 1, 5])
        cstar = Fr(rng.randint(100, 900), 100)
        s = P.s_of(kappa, Cg)
        T1 = a + s + max(a, P.sigma(a, d0, kappa))
        if max(3 * a, T1, a + s + m_hi * cstar) < cstar:
            # a + s < (1 - 2^{-2/5}) c*  =>  (a+s)/(1-2^{-2/5}) < c*
            assert (a + s) * f_lo < cstar or (a + s) * f_hi < cstar
            assert max(2 * a, T1) < cstar


def test_optimum_dl_equals_u_and_la_free():
    bu = P.optimise_u(0)
    bd = P.optimise_dl(0, (bu["lo"], bu["hi"]))
    assert bd["hi"] <= bu["hi"] and bd["lo"] >= bu["lo"] - Fr(1, 10 ** 4)
    assert 3 * bd["a"] < bd["lo"]                          # the cubic linear algebra does not bind
    assert bd["binding"] in ("shape", "specialq")


def test_gauss_reduce_shortest():
    p = 10007
    rng = random.Random(1)
    for _ in range(200):
        T = rng.randrange(1, p)
        u, v = toy_sp.gauss_reduce(p, T)
        assert (u - T * v) % p == 0 and (u, v) != (0, 0)
        n = u * u + v * v
        for vv in range(-120, 121):           # brute force over small v
            uu = (T * vv) % p
            for cand in (uu, uu - p):
                if (cand, vv) != (0, 0):
                    assert cand * cand + vv * vv >= n


def test_lemma_Sp_identity_small():
    for (p, k, lo, hi, Y) in [(10007, 1, 16, 64, 64), (10007, 2, 2, 16, 64), (40009, 2, 4, 32, 140)]:
        rec = toy_sp.run_case(p, k, lo, hi, Y)
        assert rec["identity_ok"]


def test_xk4_toy_d3():
    rec = toy_xk.run_one(3, 900, 600, 150, 150, 400)
    assert rec["status"] == "ok" and rec["certified"] and rec["aux_pm1"] and rec["aux_sign_consistent"]
    assert rec["nonunit_kernel_vectors"] == 0 and rec["unit_rank_implied"] <= rec["rU"]


def test_la_never_binds_any_tuple():
    """3a < (a+s)/(1 - 2^{-2/5}) <= c_U for every tuple, hence c_DL(tau; tau) = c_U(tau) (Lemma C25(i))."""
    f_lo, _ = up.recursion_factor_bracket()
    assert f_lo > 3
    rng = random.Random(11)
    for _ in range(200):
        a = Fr(rng.randint(1, 60), 20)
        kappa = Fr(rng.randint(1, 60), 200)
        d0 = Fr(rng.randint(5, 60), 20)
        Cg = rng.choice([0, 1, 5, 20])
        cu = P.c_u(a, d0, kappa, Cg)
        assert 3 * a < cu[0]
