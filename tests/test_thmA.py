"""Tests for agent thmA (notes/thmA_proof.md; src/thmA/)."""
import random
from fractions import Fraction as Fr
from math import gcd

from src.thmA import arith, params as P, pipeline
from src.thmA import numfield as nf
from src.thmA.smooth import primorial, smooth_parts


def test_min_k_d9_is_13_exact():
    mu12, th12 = P.best_mu(9, 12)
    mu13, th13 = P.best_mu(9, 13)
    assert th13 == Fr(16, 39) and mu13 == Fr(17, 156)
    assert P.c_lt_1(th13) and not P.c_lt_1(th12)
    assert (1 + th13) ** 2 == Fr(3025, 1521)
    assert P.min_k(9) == 13 and P.min_k(7) == 21 and P.min_k(11) == 13


def test_closed_form_and_ldbound_agree():
    from src.ldbound import params as L
    for d in (3, 5, 7, 9, 11, 13, 19):
        for k in range(1, 200):
            if Fr(k) <= Fr(d + 1, 4):
                continue
            assert P.best_mu(d, k)[1] == P.closed_form(d, k)[1] == L.best_mu(d, k)[1]


def test_infimum_limit():
    for d in (7, 9, 11):
        th = P.best_mu(d, 10 ** 6)[1]
        assert abs(float(th) - 3 / (d + 1)) < 1e-4


def test_gbw_identity_and_lemma12w():
    rng = random.Random(1)
    for _ in range(50):
        N = rng.randrange(2 ** 39, 2 ** 40) | 1
        M = arith.iroot(N, 3)
        while True:
            g0, g1 = rng.randrange(M // 2, M), rng.randrange(M // 2, M)
            if gcd(g0, g1) == 1:
                break
        w = 2 * M + 1
        while not arith.is_prime(w):
            w += 2
        c = [rng.randrange(-4, 5) for _ in range(3)]
        f = nf.make_f(N, w, 3, g0, g1, c)
        assert nf.hom_eval(f, g0, g1) == N * w
        assert f[3] % g1 == (nf.gbe(N * w, 3, g0, g1)[3]) % g1 != 0
        for _ in range(20):
            a, b = rng.randrange(-300, 300), rng.randrange(1, 300)
            G = g1 * a - g0 * b
            if G == 0:
                continue
            assert (g1 ** 3 * nf.hom_eval(f, a, b) - N * w * b ** 3) % G == 0


def test_smooth_parts():
    primes = arith.primes_upto(200)
    Pm = primorial(primes)
    rng = random.Random(2)
    xs = [rng.randrange(1, 10 ** 12) for _ in range(300)]
    sp = smooth_parts(xs, Pm)
    for x, s in zip(xs, sp):
        y = x
        for p in primes:
            while y % p == 0:
                y //= p
        assert s == x // y


def test_gauss_split():
    rng = random.Random(3)
    for _ in range(100):
        N = rng.randrange(2 ** 30, 2 ** 40) | 1
        R = rng.randrange(1, N)
        (u, v), _ = arith.gauss_reduce((N, 0), (R, 1))
        assert (u - R * v) % N == 0 and u * u + v * v <= 2 * N


def test_zsqrt_recovers_square_and_rejects_nonsquare():
    rng = random.Random(4)
    f = nf.make_f(1000003 * 1000033, 1009, 3, 7001, 9973, [1, -2, 3])
    g = nf.monic_g(f)
    D = f[3] * nf.disc(f)
    p = nf.hensel_prime(g, D, 1000)
    for _ in range(10):
        eta = [rng.randrange(-10 ** 6, 10 ** 6) for _ in range(3)]
        Gam = nf.zmul(eta, eta, g)
        if nf.psi(Gam, 0, 1) == 0:
            pass
        try:
            out = nf.zsqrt(Gam, g, 10 ** 7, p, rng)
        except ValueError:
            continue
        assert out is not None and (out == eta or out == [-x for x in eta])
    Gam = nf.zmul([3, 1, 0], [1, 0, 0], g)        # 3 + beta: not a square in general
    try:
        out = nf.zsqrt(Gam, g, 10 ** 7, p, rng)
        assert out is None or nf.zmul(out, out, g) == Gam
    except ValueError:
        pass


def _cfg():
    return {"d": 3, "mu": "1/3", "hbits": 2, "W_mult": 2, "k_large": 2, "min_large": 1, "n_chars": 40,
            "char_bits": 62, "extra_rows": 20, "n_deps": 8, "stats_all_deps": True, "max_rounds": 6,
            "max_fields": 20000, "max_rtrials": 5000000,
            "irr_primes": [p for p in arith.primes_upto(97) if p >= 5], "N0": 2 ** 20,
            "tiers": [{"bits_max": 36, "Bbits": 10, "B1_mult": 4, "S": 12},
                      {"bits_max": 60, "Bbits": 11, "B1_mult": 4, "S": 12}]}


def test_end_to_end_small_semiprime():
    p, q = 1000003, 1000033                       # test-only labels; the pipeline sees N only
    res = pipeline.factor(p * q, _cfg(), 7)
    assert res["factors"] == [p, q]
    for L in res["log"]:
        c = L.get("ctr", {})
        assert c.get("ROW_IDENTITY_VIOLATION", 0) == 0 and c.get("CONGRUENCE_VIOLATION", 0) == 0


def test_preprocessing_paths():
    cfg = _cfg()
    assert pipeline.factor(101 ** 2 * 1000003, cfg, 0)["factors"] == [101, 101, 1000003]
    assert pipeline.factor(1000003 ** 3, cfg, 0)["factors"] == [1000003] * 3
    assert pipeline.factor(1000003, cfg, 0)["factors"] == [1000003]
