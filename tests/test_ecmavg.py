"""Tests for src/ecmavg (exponent bookkeeping, operational prefilter, diagnostic labels)."""
from fractions import Fraction as F
import random

import numpy as np
from sympy import primerange

from src.ecmavg import exponents as ex
from src.ecmavg.curves import bound_B, ecm_prefilter, params_for_n, stage_steps
from src.ecmavg.labels import INF, exposure_labels, naive_label, truncate


# ---------------------------------------------------------------- exponents
def test_codex_point_reproduced():
    d = ex.codex_point()
    assert d["ecm"] == d["fallback"] == d["total"] == F(23, 152)


def test_fallback_R_closed_form_matches_search():
    r = random.Random(7)
    for _ in range(40):
        th = F(r.randint(0, 400), 2000)
        if th == 0:
            continue
        assert ex.fallback_R(th) == ex.fallback_R_search(th, 60)


def test_fallback_values():
    assert ex.fallback_R(F(1, 38)) == F(23, 152)
    assert ex.fallback_R(F(2, 25)) == F(3, 25)
    assert ex.fallback_R(F(2, 13)) == F(3, 65)
    assert ex.fallback_AV(F(1, 38)) == F(23, 152)
    assert ex.fallback_AV(F(1, 10)) == F(2, 15)
    assert ex.fallback_AV_floor()[:2] == (F(2, 15), F(1, 60))


def test_optimisations():
    # codex's own cost model, b fixed at 1/8: 23/152 is its optimum
    assert ex.optimise(F(1, 8), 1, False, F(1, 5), "AV") == (F(23, 152), F(1, 38))
    # same model, b free (-> 0): the AV floor 2/15
    assert ex.optimise(F(0), 1, False, F(1, 5), "AV")[0] == F(2, 15)
    # same model with the Theorem R fallback: 1/10
    assert ex.optimise(F(0), 1, False, F(1, 5), "R")[0] == F(1, 10)
    # early stopping, hypothesis only up to J = N^{1/38}: 23/152 for every b <= 23/152
    assert ex.optimise(F(1, 8), 1, True, F(1, 38), "R")[0] == F(23, 152)
    # early stopping, hypothesis up to J = N^{1/5}: exponent b (-> 0)
    assert ex.optimise(F(1, 100), 1, True, F(1, 5), "R")[0] == F(1, 100)


def test_kappa_closed_form():
    for k in [F(1, 10), F(1, 4), F(1, 3), F(1, 2), F(2, 3), F(9, 10)]:
        assert ex.optimise(F(0), k, True, F(1, 5), "R")[0] == ex.kappa_closed_form(k)


# ---------------------------------------------------------------- labels
def test_vectorised_labels_match_naive_group_law():
    primes = list(primerange(101, 700))
    for fam in ("suyama", "mont"):
        B, Xb = 30, 800
        k = exposure_labels(fam, 12, primes, B, Xb)
        for j in range(12):
            for i, p in enumerate(primes):
                assert k[j, i] == naive_label(fam, j, p, B), (fam, j, p)


def test_truncation_equals_recomputation():
    primes = list(primerange(1000, 1400))
    k_big = exposure_labels("suyama", 8, primes, 60, 1500)
    k_small = exposure_labels("suyama", 8, primes, 20, 1500)
    assert np.array_equal(truncate(k_big, len(stage_steps(20, 1500))), k_small)


# ---------------------------------------------------------------- operational
def test_prefilter_matches_labels_on_pairs():
    """Operational run (sees only N) succeeds at the first curve with k_p != k_q and outputs
    the prime with the smaller label (Lemma S)."""
    r = random.Random(3)
    b = 12
    primes = list(primerange(2 ** (b - 1), 2 ** b))
    beta = F(1, 3)
    B = bound_B(b, beta)
    k = exposure_labels("suyama", 16, primes, B, 2 ** b + 2 ** (b // 2 + 1) + 1)
    idx = {p: i for i, p in enumerate(primes)}
    for _ in range(150):
        p, q = sorted(r.sample(primes, 2))
        N = p * q
        assert params_for_n(N.bit_length(), beta)[0] == B
        d, used = ecm_prefilter(N, 16, beta, "suyama")
        kp, kq = k[:, idx[p]], k[:, idx[q]]
        diff = [j for j in range(16) if kp[j] != kq[j]]
        if not diff:
            assert d is None
        else:
            j = diff[0]
            assert used == j + 1
            small = p if (kq[j] == INF or (kp[j] != INF and kp[j] < kq[j])) else q
            assert d == small


def test_pair_stats_brute_force():
    from itertools import combinations
    from src.ecmavg.experiments import pair_stats
    r = np.random.default_rng(5)
    k = r.integers(-1, 4, size=(64, 23)).astype(np.int16)
    ps, npairs = pair_stats(k)
    for J in (8, 64):
        zs = [int((k[:J, i] != k[:J, j]).sum()) for i, j in combinations(range(23), 2)]
        zu = [int(((k[:J, i] != INF) != (k[:J, j] != INF)).sum()) for i, j in combinations(range(23), 2)]
        assert np.isclose(ps[J]["E_Zs_over_J"], np.mean(zs) / J)
        assert np.isclose(ps[J]["Var_Zs_over_J"], np.var(zs) / J)
        assert np.isclose(ps[J]["Var_Zu_over_J"], np.var(zu) / J)
        assert np.isclose(ps[J]["fail_s"], np.mean(np.array(zs) == 0))
