"""Tests for src/abcsunit (agent abcsunit, phase 19)."""
from math import gcd

import numpy as np
import pytest

from src.abcsunit import core


def test_primes_and_smooth_count():
    assert core.primes_upto(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    P = core.primes_upto(13)
    for X in (1, 10, 1000, 123457):
        vals = core.smooth_numbers(X, P)
        brute = [n for n in range(1, X + 1) if core.is_smooth(n, P)] if X < 200000 else None
        assert vals.size == core.smooth_count(X, P)
        if brute is not None:
            assert list(map(int, vals)) == brute


def test_cap():
    with pytest.raises(core.TooMany):
        core.smooth_numbers(10 ** 9, core.primes_upto(50), cap=1000)


def test_congruent_pairs_bruteforce():
    P = core.primes_upto(7)
    N = 91
    vals = core.smooth_numbers(5000, P)
    pairs, _, _ = core.congruent_pair_count(vals, N)
    L = list(map(int, vals))
    brute = sum(1 for i in range(len(L)) for j in range(i + 1, len(L)) if (L[j] - L[i]) % N == 0)
    assert pairs == brute
    rel = core.relations_from_buckets(core.congruent_buckets(vals, N), N, coprime_only=False)
    assert len(rel) == brute
    for u, v, k in rel:
        assert v - u == k * N and core.is_smooth(u, P) and core.is_smooth(v, P)


def test_prop_AB_lower_bound_with_H():
    # sum_r c_r^2 >= S^2 / (#occupied classes) >= S^2/|H|
    P = core.primes_upto(11)
    p, q = 101, 113
    N = p * q
    vals = core.smooth_numbers(10 ** 7, P)
    c = core.residue_class_sizes(vals, N).astype(object)
    S = int(vals.size)
    H = core.subgroup_order_diag(P, p, q)
    assert H == core.subgroup_order_bruteforce_diag(P, N)
    assert len(c) <= H
    assert sum(int(x) ** 2 for x in c) * H >= S * S


def test_subgroup_order_small_cases():
    for p, q, P in [(7, 11, [2]), (7, 11, [2, 3]), (31, 43, [2, 5]), (101, 103, [2, 3, 5, 7])]:
        assert core.subgroup_order_diag(P, p, q) == core.subgroup_order_bruteforce_diag(P, p * q)


def test_f2_dependencies():
    rows = [0b011, 0b110, 0b101, 0b111]
    rank, deps = core.f2_dependencies(rows, 3)
    assert rank == 3 and len(deps) == 1
    for m in deps:
        acc = 0
        for i, r in enumerate(rows):
            if m >> i & 1:
                acc ^= r
        assert acc == 0


def test_square_and_factor_toy():
    p, q = 1009, 1013
    N = p * q
    P = core.primes_upto(40)
    vals = core.smooth_numbers(N * 30, P)
    rel = core.relations_from_buckets(core.congruent_buckets(vals, N), N, coprime_only=False,
                                      consecutive_only=True)
    rows, exps = [], []
    for u, v, k in rel:
        b, e = core.relation_row(u, v, P)
        rows.append(b)
        exps.append(e)
    rank, deps = core.f2_dependencies(rows, len(P))
    assert rank <= len(P)
    gs = {core.square_from_dependency(m, rel, exps, P, N) for m in deps}
    for g in gs:
        assert g in (1, p, q, N)
    assert p in gs or q in gs  # at least one fruitful dependency at this size


def test_primitive_triples_match_von_kaenel_matschke():
    # VKM 2016 Theorem A: #solutions up to symmetry for S = first n primes: 1, 4, 17, 63.
    # A solution x + y = 1 up to the 6-element symmetry group <-> a primitive a + b = c, a <= b.
    for n, expected, H in [(1, 1, 10 ** 6), (2, 4, 10 ** 9), (3, 17, 10 ** 12), (4, 63, 10 ** 13)]:
        P = core.primes_upto(30)[:n]
        vals = core.smooth_numbers(H, P)
        tr = core.primitive_smooth_triples(vals, H)
        assert len(tr) == expected, (n, len(tr))
        for a, b, c in tr:
            assert a + b == c and gcd(a, b) == 1


def test_radicals():
    P = core.primes_upto(13)
    assert core.rad_smooth(2 ** 5 * 3 * 13 ** 2, P) == 2 * 3 * 13
    assert core.rad_diag(2 ** 3 * 97 ** 2) == 2 * 97


def test_sigma_haar_mean_is_one():
    from src.abcsunit import model
    for l in (2, 3, 5, 7):
        for a in (0.3, 0.5, 0.8):
            mean = sum(l ** (-m) * (1 - 1 / l) * model.sigma_small(l, m, a) for m in range(0, 60))
            assert abs(mean - 1.0) < 1e-6, (l, a, mean)


def test_pair_model_fast_matches_exact():
    from src.abcsunit import model
    p, q = 2003, 2011
    N = p * q
    P = core.primes_upto(40)
    X = 6 * N
    vals = core.smooth_numbers(X, P)
    H = core.subgroup_order_diag(P, p, q)
    a = model.pair_model(vals, N, X, H, P, vals)
    b = model.pair_model_fast(vals, N, X, H, P, vals, usample=10 ** 9)
    assert abs(a - b) / a < 1e-9


def test_lemma_FR_critic_example_zrank_deficit_still_fruitful():
    # Opus proof critic: N = 1147 = 31*37, S = {2,3,5,7}: Z-rank 2 < s = 4, r = 0, rho = 2 -> fruitful.
    N, P = 1147, [2, 3, 5, 7]
    rels = [(225 ** 2, 256 ** 2), (175 ** 2, 972 ** 2)]
    rels = [(min(a, b), max(a, b), (max(a, b) - min(a, b)) // N) for a, b in rels]
    rows, exps = [], []
    for u, v, k in rels:
        assert (v - u) % N == 0
        b, e = core.relation_row(u, v, P)
        rows.append(b)
        exps.append(e)
    rank, deps = core.f2_dependencies(rows, len(P))
    assert rank == 0 and len(deps) == 2
    gs = {core.square_from_dependency(m, rels, exps, P, N) for m in (deps[0], deps[1], deps[0] ^ deps[1])}
    assert {31, 37} <= gs


def test_closure_lemma_identity():
    # codex closure: c = a(b+N)/(b-a), c+N = b(a+N)/(b-a); exponent vector of (c, c+N) is a Z-combination.
    N = 1147
    for a in range(1, 300):
        for b in range(a + 1, 300):
            num = a * (b + N)
            if num % (b - a) == 0:
                c = num // (b - a)
                assert c + N == b * (a + N) // (b - a) and (b * (a + N)) % (b - a) == 0
