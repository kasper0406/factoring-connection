import math
import random

from src.crtdecode.core import (certified, is_smooth, prefix_incidences, prime_powers_upto, primes_upto,
                                smooth_in_range, window_profile)
from src.crtdecode.decoders import (_int_pow_ceil, divisor_trick, heavy_hitter, information_decode_and_certify,
                                    random_tests)
from src.crtdecode.theory import MNFS, a_star, c_of_a, c_of_a_grid, dictionary


def _brute_smooth(z, y):
    for p in primes_upto(y):
        while z % p == 0:
            z //= p
    return z == 1


def test_window_profile_matches_trial_division():
    rng = random.Random(3)
    for _ in range(5):
        X = rng.randrange(1 << 29, 1 << 30)
        h, y = 3000, 97
        prof = window_profile(X, h, y)
        truth = [_brute_smooth(X + j, y) for j in range(h + 1)]
        assert [bool(b) for b in prof["smooth"]] == truth
        assert all(is_smooth(X + j, y) == truth[j] for j in range(h + 1))
        inc = sum(1 for (_q, _e, m) in prime_powers_upto(y, X + h) for j in range(h + 1) if (X + j) % m == 0)
        assert prof["incidences"] == inc
        if prof["first"] >= 0:
            jl = prof["first"]
            assert prefix_incidences(X, jl, y) == sum(
                1 for (_q, _e, m) in prime_powers_upto(y, X + jl) for j in range(jl + 1) if (X + j) % m == 0)


def test_smooth_in_range():
    assert smooth_in_range(100, 200, 7) == [z for z in range(100, 201) if _brute_smooth(z, 7)]


def test_certified_exact_rules():
    assert certified(64, 64, 64) and not certified(32, 64, 64)
    assert certified(1000, 1000, 999) and not certified(1000, 1000, 1001)
    # q^e with q = 3, e = 5, region 100: b certified iff 3^(b*b) >= 100^5
    for b in range(1, 6):
        assert certified(3 ** b, 3 ** 5, 100) == (3 ** (b * b) >= 100 ** 5)


def test_int_pow_ceil():
    for h in (1 << 10, 1 << 14, 1 << 16):
        for th in (1.0, 1.25, 1.5):
            v = _int_pow_ceil(h, th)
            k = h.bit_length() - 1
            a = int(round(k * th * 4))
            assert v ** 4 >= (1 << a) and (v - 1) ** 4 < (1 << a)


def test_decoders_output_smooth_and_respect_budget():
    rng = random.Random(5)
    h, y = 1 << 10, 64
    rep = smooth_in_range(h // 2, h, y)
    rep2 = smooth_in_range(h // 2, 8 * h, y)
    for _ in range(10):
        X = rng.randrange(1 << 29, 1 << 30)
        for r in (divisor_trick(X, h, y, rep, 500), heavy_hitter(X, h, y, rep2, 500, 1.25),
                  random_tests(X, h, y, list(range(h + 1)), 500)):
            assert r.cost <= 500
            assert r.cost == r.tests + r.pointers
            if r.found:
                assert _brute_smooth(r.z, y) and X <= r.z <= X + h


def test_information_decoder_exact_and_shallow():
    rng = random.Random(7)
    bits, h, y = 34, 1 << 10, 64
    pi = len(primes_upto(y))
    for _ in range(6):
        X = rng.randrange(1 << (bits - 1), 1 << bits)
        prof = window_profile(X, h, y)
        truth = [int(j) for j in prof["smooth"].nonzero()[0]]
        r = information_decode_and_certify(X, h, y, bits)
        assert r["smooth_offsets"] == truth
        assert r["output"] == (truth[0] if truth else None)
        assert r["depth"] <= 3 * pi + 1


def test_lemma_gc_zoom_chain_closure():
    """Lemma GC: with a monotone threshold T on gcd(z, P1), every z in the progression of a certified divisor
    d = gcd(z1, P1) >= T is itself certified by call 1."""
    rng = random.Random(11)
    P1 = 2 ** 6 * 3 ** 4 * 5 ** 3 * 7 ** 2 * 11 * 13 * 17 * 19 * 23
    h = 5000
    for _ in range(50):
        X = rng.randrange(10 ** 9, 2 * 10 ** 9)
        T = rng.randrange(50, 5000)
        C1 = {j for j in range(h + 1) if math.gcd(X + j, P1) >= T}
        for j1 in C1:
            d = math.gcd(X + j1, P1)
            prog = range(j1 % d, h + 1, d)
            assert set(prog) <= C1


def test_theory_constants():
    assert abs(a_star() - 1.1415601395) < 1e-9
    assert abs(c_of_a(1.0) - 6 ** (1 / 3)) < 1e-12
    assert abs(c_of_a_grid(a_star()) - MNFS) < 1e-5
    d = dictionary(1024)
    assert abs(d["t_over_k2"] - 0.5) < 1e-9   # t = k^2/2 at kappa*, nu*
    assert d["t"] > d["k"] and d["t"] < d["johnson_t"]


def test_no_leakage_in_operational_decoders():
    import inspect
    import src.crtdecode.decoders as D
    src = inspect.getsource(D)
    for bad in ("window_profile", "factorint", "sympy", "flint"):
        assert bad not in src


def test_crt_query_tree_shallow_and_exact():
    from src.crtdecode.decoders import crt_query_tree
    rng = random.Random(13)
    bits, h, y = 34, 1 << 10, 64
    for _ in range(6):
        X = rng.randrange(1 << (bits - 1), 1 << bits)
        prof = window_profile(X, h, y)
        r = crt_query_tree(X, h, y, bits)
        assert r["reconstructed"]
        assert r["output"] == (prof["first"] if prof["first"] >= 0 else None)
        assert r["depth"] <= 6


def test_oracle_rejects_unsupported_general_modulus():
    import pytest
    from src.crtdecode.decoders import JohnsonOracle
    with pytest.raises(NotImplementedError):
        JohnsonOracle(15).call(30, 0, 2)
