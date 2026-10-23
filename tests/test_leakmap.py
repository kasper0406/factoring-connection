"""Tests for src/leakmap (phase 24)."""
import random
from collections import Counter
from math import gcd

import gmpy2
import numpy as np
import pytest

from src.leakmap import netrees as nt
from src.leakmap import cfgenus as cg
from src.leakmap import jacobi as jc
from src.leakmap.e1 import population, run_cell, primes_top2
from src.leakmap.e2_posthoc import squares_split
from src.eps.cf import cf_features


def test_algebra_pow_matches_repeated_mul():
    alg = nt.LUCAS_ALG
    m = 1009
    e = (5, 7)
    r = alg.one()
    for k in range(40):
        assert alg.pow(e, k, m) == alg.red(r, m)
        r = alg.mul(r, e, m)


def test_lucas_alpha_norm_one():
    # alpha = x in Z[x]/(x^2 - 3x + 1): alpha * (3 - alpha) = 1
    alg = nt.LUCAS_ALG
    assert alg.mul((0, 1), (3, -1), 10007) == (1, 0)


def test_lazy_nonadaptive_no_blowup():
    t = nt.tree_arch(32)
    out = nt.run_on_N(t, 3167 * 3169)
    assert out[0] in ("leaf", "split")


@pytest.mark.parametrize("tname", list(nt.TREES))
def test_diagonal_lemma_small(tname):
    tree = nt.TREES[tname]()
    P = primes_top2(10)
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            p, q = P[i], P[j]
            N = p * q
            out = nt.run_on_N(tree, N)
            tp, tq = nt.diag_transcript(tree, N, p), nt.diag_transcript(tree, N, q)
            if out[0] == "split":
                assert tp != tq and out[1] in (p, q)
            else:
                assert tp == tq == out[1]


def test_formal_transcript_is_prime_only_and_pure_matches_diagonal():
    # for the pure control tree the formal transcript equals the diagonal one (no deviations ever)
    tree = nt.tree_pm1()
    for p, q in population(10, "iid")[:200]:
        assert nt.diag_transcript(tree, p * q, p) == nt.formal_transcript(tree, p)


def test_ne1_bound_and_ptx_hold_small():
    pairs = population(11, "iid")
    for tname in nt.TREES:
        rec = run_cell(11, "iid", tname, pairs)
        assert rec["p1_mismatch"] == 0
        assert rec["p3_violations"] == 0
        assert all(v["ok"] for v in rec["p2"].values())


def test_arch_universal_pin():
    # q = p + k - 1 makes 2^N = 2^k at p for every base (p-1 | N - k)
    p, q = 3167, 3169
    N = p * q
    k = q - p + 1
    assert (N - k) % (p - 1) == 0
    nd = nt.tree_arch()[1][k - 1]
    assert nd.t_target == k
    assert nt.deviation_class(nd, N, p, L=64) == "UNIVERSAL"


def test_lucas_twin_prime_two_sided_universal():
    # q = p + 2, (5|p) = -1, (5|q) = +1: N + 1 = (p+1)^2 is divisible by |T_p| = p+1 and |T_q| = q-1
    p, q = 3167, 3169
    assert gmpy2.legendre(5, p) == -1 and gmpy2.legendre(5, q) == 1
    tree = nt.tree_lucas()
    N = p * q
    out = nt.run_on_N(tree, N)
    assert out[0] == "leaf" and out[1][0] == 1          # alpha^(N+1) = 1 mod N: Lucas pseudoprime, no split
    nd = tree[1][0]
    assert nt.deviation_class(nd, N, p, L=8) == "UNIVERSAL"
    assert nt.deviation_class(nd, N, q, L=8) == "UNIVERSAL"


def test_order_mod_prime():
    for x in (101, 103, 1009):
        for a in (2, 3, 5):
            o = nt.order_mod_prime(a, x)
            assert pow(a, o, x) == 1
            assert all(pow(a, d, x) != 1 for d in range(1, o))


def test_cf_stats_match_eps_and_controls_sound():
    for seed in range(60):
        inst = cg.balanced_semiprime(32, 10_000 + seed)
        N = inst.N
        p, q = inst.labels["factors"]
        s = cg.cf_stats(N)
        f = cf_features(N)
        for k in cg.STATS:
            assert s[k] == f[k]
        for l in cg.ELLS:
            if s[f"pc{l}_even"]:
                assert gmpy2.legendre(l, p) == 1 and gmpy2.legendre(l, q) == 1
            if s[f"pc{l}_odd"]:
                assert gmpy2.legendre(-l, p) == 1 and gmpy2.legendre(-l, q) == 1
        assert squares_split(N)[0] == f["cf64_n_sq_even"]


def test_labels_populations():
    lab = cg.labels(13, 17)          # N = 221 = 1 mod 4, 221 mod 5 = 1, mod 7 = 4, mod 3 = 2
    assert lab["G1"] == 1 and lab["G2"] is None and lab["G4"] is None
    assert lab["G3"] == 0 and lab["NG2"] is None


def test_strat_z_null_and_signal():
    rng = np.random.default_rng(0)
    s = rng.integers(0, 8, 20000)
    y = rng.integers(0, 2, 20000).astype(float)
    x = rng.normal(size=20000)
    assert abs(cg.strat_z(x, y, s)) < 4.5
    assert cg.strat_z(x + 0.2 * y, y, s) > 8
    # cluster z with singleton clusters equals a heteroskedastic z of the same sign
    ids = np.arange(20000)
    assert cg.cluster_z(x + 0.2 * y, y, s, ids, ids + 10 ** 6) > 8


def test_legendre_vec_and_lucas():
    X = jc.primes_in(1000, 3000).astype(np.int64)
    c = np.full_like(X, 12345)
    lv = jc.legendre_vec(c, X)
    for xi, li in zip(X[:50], lv[:50]):
        assert li == gmpy2.legendre(12345, int(xi))
    V = jc.lucas_V(3, 120, X[:20])
    for xi, vi in zip(X[:20], V):
        v0, v1 = 2, 3
        for _ in range(120):
            v0, v1 = v1, (3 * v1 - v0)
        assert vi == v0 % int(xi)


def test_ptj_criterion_toy():
    """For a single Jacobi test on N = 1 mod 4 with iid factors, P_e[(c|N) = 1] = (1 + (beta(c) + e beta(-c))^2)/2
    exactly over a finite prime set; check the identity by brute force."""
    P = [x for x in primes_top2(12)]
    c = 7
    for e in (1, -1):
        S = [x for x in P if (x % 4 == 1) == (e == 1)]
        lam = [int(gmpy2.legendre(c, x)) for x in S]
        m = sum(lam) / len(S)
        tot = sum(1 for a in lam for b in lam if a * b == 1) / len(S) ** 2
        assert abs(tot - (1 + m * m) / 2) < 1e-12
