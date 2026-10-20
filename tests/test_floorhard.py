import random
from math import gcd

import gmpy2

from src.floorhard import core


def test_window_identity():
    rng = random.Random(1)
    for _ in range(200):
        f = [rng.randint(-50, 50) for _ in range(rng.randint(1, 4))]
        D = rng.randint(1, 40)
        h1 = rng.randint(-20, 20)
        h2 = h1 + rng.randint(0, 30)
        u = rng.randint(0, D - 1)
        v = rng.randint(u + 1, D)
        want = sum(1 for h in range(h1, h2) if u <= core.poly_eval(f, h) % D < v)
        assert core.window_count(core.floor_sum_brute, f, D, h1, h2, u, v) == want


def test_floor_sum_np_matches_brute():
    rng = random.Random(2)
    for _ in range(50):
        f = [rng.randint(-1000, 1000) for _ in range(3)]
        D = rng.randint(1, 500)
        h1 = rng.randint(-100, 100)
        h2 = h1 + rng.randint(0, 200)
        assert core.floor_sum_np(f, D, h1, h2) == core.floor_sum_brute(f, D, h1, h2)


def test_full_period_formula_small():
    for m in range(3, 400, 2):
        if not core.squarefree_odd(m):
            continue
        for a in (1, 2, 7):
            if gcd(a, m) == 1:
                assert core.full_period_brute(m, a) == core.full_period_rhs(m, a)


def test_class_number_decode_and_factor():
    rng = random.Random(5)
    for m in (7 * 11, 13 * 19, 5 * 23, 7 * 19, 11 * 31, 13 * 29):   # mix of m = 3 mod 4 and p = q = 3 mod 4
        M, hM, calls = core.decode_class_number_semiprime(m, rng)
        assert hM == core.forms.class_number(-M)
        g = core.factor_from_class_number(m, M, hM, rng)
        assert g is not None and m % g == 0 and 1 < g < m


def test_subset_sum_reduction_parsimonious():
    rng = random.Random(3)
    for _ in range(30):
        k = rng.randint(2, 6)
        w = [rng.randint(1, 10) for _ in range(k)]
        T = sum(x for x in w if rng.random() < 0.5)
        inst = core.subset_sum_to_qc(w, T)
        assert core.qc_count_crt(inst) == core.subset_sum_count(w, T)


def test_subset_sum_via_floor_sums_tiny():
    inst = core.subset_sum_to_qc([1, 2], 2)
    assert core.qc_count_floorsums(inst, oracle=core.floor_sum_np) == core.subset_sum_count([1, 2], 2) == 1


def test_taylor_embedding_windows():
    rng = random.Random(4)
    for j in (2, 3, 4, 5):
        for _ in range(4):
            b = rng.randint(2, 40)
            A = rng.randint(1, 4)
            cp = rng.randint(1, 2 * b)
            emb = core.TaylorEmbedding(A, b, cp, j)
            assert emb.H >= cp and emb.x * emb.x <= emb.N
            orc = core.PrefixOracle(cp)
            for r in range(b):
                want = sum(1 for h in range(cp) if (A * h * h - r) % b == 0)
                assert emb.window_count(orc, r, 0, cp) == want


def test_prime_embedding():
    rng = random.Random(6)
    emb = core.find_prime_embedding(2, 21, 30, 2, rng)
    assert emb is not None and gmpy2.is_prime(emb.N)
    orc = core.PrefixOracle(30)
    for r in range(21):
        want = sum(1 for h in range(30) if (2 * h * h - r) % 21 == 0)
        assert emb.window_count(orc, r, 0, 30) == want


def test_factor_via_taylor_class():
    for m in (15, 77, 221, 3 * 5 * 7, 45, 1001):
        g, calls, N = core.factor_via_taylor_class(m)
        assert g is not None and 1 < g < m and m % g == 0
        assert calls <= 2 * (m - 1).bit_length() + 4


def test_subset_sum_deterministic_and_range():
    rng = random.Random(8)
    for _ in range(40):
        k = rng.randint(1, 6)
        w = [rng.randint(1, 12) for _ in range(k)]
        T = rng.randint(-3, sum(w) + 3)
        inst = core.subset_sum_to_qc(w, T, deterministic=True)
        assert core.qc_count_crt(inst) == (core.subset_sum_count(w, T) if T >= 0 else 0)
    assert core.subset_sum_to_qc([1, 2], 9) is None      # referee counterexample (aliasing)
