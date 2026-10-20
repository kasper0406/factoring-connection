import math, random
import numpy as np
import pytest
from sympy import isprime

from src.prepread import slp as S, spectra as SP, compsearch as CSR, dlred
from src.prepread.run_experiments import tonelli, preimages, haar_wreath


def test_compsearch_matches_bruteforce():
    for k, ps in [(4, [13]), (5, [29]), (5, [11, 13]), (6, [37, 41]), (6, [17, 19, 23])]:
        r = CSR.search(k, ps)
        b, arg = CSR.brute_best(k, ps)
        assert r['best'] == b and r['nties'] == len(arg)
        assert all(tuple(a) in set(arg) for a in r['argmax'])


def test_fibre_matches_slp_eval():
    rng = random.Random(1)
    for _ in range(20):
        cs = [rng.choice(CSR.CS) for _ in range(rng.randrange(1, 7))]
        p = rng.choice([31, 37, 101, 103])
        v = S.eval_mod(S.composition_slp(cs), np.arange(p, dtype=np.int64), p)
        assert int(np.count_nonzero(v == 0)) == CSR.fibre(cs, p)


def test_reduced_poly_agrees_with_eval():
    p = 101
    xs = np.arange(1, p, dtype=np.int64)
    for i in range(30):
        prog = S.random_slp(12, random.Random(i))
        co = S.reduced_poly(prog, p)
        direct = S.eval_mod(prog, xs, p)
        via = np.array([sum(int(co[e]) * pow(int(x), e, p) for e in range(len(co))) % p for x in xs])
        assert np.array_equal(direct, via)


def test_weil_width_and_bound():
    p = 211
    g = SP.primitive_root(p)
    orb = SP.orbit(p, g)
    co = np.zeros(p - 1, dtype=np.int64); co[3] = 1; co[p - 1 - 2] = 5   # y^3 + 5 y^{-2}
    assert S.weil_width(co, p) == 5
    co2 = np.zeros(p - 1, dtype=np.int64); co2[7] = 2; co2[9] = 1
    assert S.weil_width(co2, p) == 9
    # Weil: |A| <= w sqrt(p)/(p-1) for nontrivial characters
    for prog in [S.power_slp(3, plus_x=True), S.composition_slp([1, 2])]:
        v = S.eval_mod(prog, orb, p)
        w = S.weil_width(S.reduced_poly(prog, p), p)
        M = SP.spectrum_max(v, p, range(1, p), np.arange(p - 1) != 0)[0]
        assert M <= w * math.sqrt(p) / (p - 1) + 1e-12


def test_spectrum_max_matches_direct():
    p = 103
    g = SP.primitive_root(p)
    v = S.eval_mod(S.random_slp(8, random.Random(3)), SP.orbit(p, g), p)
    mask = SP.large_order_mask(p)
    M, c, j = SP.spectrum_max(v, p, range(1, p), mask)
    assert abs(abs(SP.A_value(v, p, c, j)) - M) < 1e-12
    assert mask[j]


def test_dlog_reduction_planted():
    p = 1009
    g = SP.primitive_root(p)
    orb = SP.orbit(p, g)
    logs = np.empty(p, dtype=np.int64); logs[orb] = np.arange(p - 1)
    vp = (p * logs) // (p - 1)
    f = lambda ys: np.exp(2j * np.pi * vp[ys] / p)
    Fj = SP.A_value((p * np.arange(p - 1)) // (p - 1), p, 1, 1)
    rng = np.random.default_rng(0)
    for x in [0, 1, 5, 500, 1007]:
        xr, q = dlred.recover_dlog(p, g, int(orb[x]), f, 1, Fj, rng, gpow=orb)
        assert xr == x


def test_dlog_reduction_noncoprime_frequency():
    # character of order (p-1)/gcd(j, p-1); cofactor brute force
    p = 1009
    g = SP.primitive_root(p)
    orb = SP.orbit(p, g)
    logs = np.empty(p, dtype=np.int64); logs[orb] = np.arange(p - 1)
    j = 6
    f = lambda ys: np.exp(2j * np.pi * (j * logs[ys] % (p - 1)) / (p - 1))
    rng = np.random.default_rng(1)
    for x in [3, 77, 1000]:
        xr, _ = dlred.recover_dlog(p, g, int(orb[x]), f, j, 1.0, rng, gpow=orb)
        assert xr == x


def test_pohlig_hellman_slp_correlation():
    p = 12289                      # 3 * 2^12 + 1
    assert isprime(p)
    e, m = 12, 3 * 5 * 7 * 9       # odd part of p-1 divides m
    w = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)
    ys = np.arange(1, p, dtype=np.int64)
    L, ops = S.ph_phi(ys, w, p, e, m)
    zeta = pow(w, m, p)
    assert all(pow(zeta, int(l), p) == pow(int(y), m, p) for y, l in zip(ys[:300], L[:300]))
    # e(c L / p) correlates with the order-2^e character e(L/2^e) for c = round(p / 2^e)
    c = round(p / 2 ** e)
    corr = abs(np.mean(np.exp(2j * np.pi * (c * L % p) / p) * np.exp(-2j * np.pi * L / 2 ** e)))
    assert corr > 0.8
    assert ops < 400


def test_tonelli_and_preimages():
    for p in [17, 97, 193, 7681, 12289]:
        for a in [a for a in range(1, 60) if a % p]:
            r = tonelli(a, p)
            if pow(a, (p - 1) // 2, p) == 1:
                assert r * r % p == a % p
            else:
                assert r is None
        T = preimages([0, 3, 5], 1, p)
        assert sorted(T) == sorted(x for x in range(p) if (x * x + 1) % p in (0, 3, 5))


def test_haar_wreath_recursion():
    pk = haar_wreath(6)
    assert pk[1] == 0.5 and abs(pk[2] - 0.375) < 1e-15
    # brute force over [C2]^2 acting on 4 leaves: fraction with a fixed leaf = 3/8
    import itertools
    cnt = 0
    for s, a, b in itertools.product([0, 1], repeat=3):
        fixed = (s == 0) and (a == 0 or b == 0)
        cnt += fixed
    assert cnt / 8 == pk[2]
