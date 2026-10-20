import random
from fractions import Fraction

import flint

from src.common.instances import balanced_semiprime
from src.topblock import core
from src.topblock import lattices as L


def test_identities_T1_T2():
    rng = random.Random(1)
    for _ in range(300):
        x = rng.randrange(2, 10 ** 6)
        h = rng.randrange(0, 10 ** 5)
        N = rng.randrange(1, 10 ** 12)
        q = rng.randrange(-10 ** 6, 10 ** 6)
        for j in range(0, 7):
            assert core.check_T1_T2(x, j, h, N, q)


def test_window_equals_near_divisor_set():
    inst = balanced_semiprime(22, 3)
    N, p = inst.N, inst.factors[0]
    rng = random.Random(2)
    for j in (2, 3, 4):
        x, H = core.top_block_containing(N, p, j)
        for _ in range(15):
            a = Fraction(rng.randrange(0, 900), 1000)
            b = a + Fraction(rng.randrange(1, 100), 1000)
            assert core.window_set_direct(N, x, j, H, a, b) == core.window_set_neardiv(N, x, j, H, a, b)


def test_divisors_land_in_divisor_window():
    inst = balanced_semiprime(30, 5)
    N, p = inst.N, inst.factors[0]
    for j in (2, 3, 4, 5):
        x, H = core.top_block_containing(N, p, j)
        h = p - x
        if h >= 1:
            a, b = core.divisor_window(N, x, j, h, h + 1)
            assert a <= core.frac_P(N, x, j, h) <= b


def test_locate_via_oracle_and_bound():
    for seed in range(4):
        inst = balanced_semiprime(32, seed)
        N, p = inst.N, inst.factors[0]
        for j in (2, 3):
            x, H = core.top_block_containing(N, p, j)
            orc = core.WindowOracle(N, x, j, H)
            divs, calls = core.locate_divisor_via_window_oracle(N, x, j, H, orc.query)
            assert p in divs
            assert calls <= core.weighted_neardiv_bound(N, x, j, H) + 1e-9


def test_window_count_via_divisor_finder():
    inst = balanced_semiprime(26, 7)
    N, p = inst.N, inst.factors[0]
    rng = random.Random(3)
    for j in (2, 3):
        x, H = core.top_block_containing(N, p, j)
        orc = core.WindowOracle(N, x, j, H)
        for _ in range(10):
            h1 = rng.randrange(0, H - 1)
            h2 = rng.randrange(h1 + 1, H + 1)
            a = Fraction(rng.randrange(0, 950), 1000)
            b = a + Fraction(rng.randrange(1, 40), 1000)
            assert orc.query(h1, h2, a, b) == core.window_count_via_divisor_finder(N, x, j, h1, h2, a, b)[0]


def test_F_relation_substitutes_to_zero_and_w_to_constant():
    N, x = 1000003 * 999983, 999000
    for j in (2, 3):
        win = [(0, j, 1)]
        # F = (x+h) r + N(-h)^{j+1}: monomials (e, t)
        monos = [(0, 1), (1, 1), (j + 1, 0)]
        coeffs = [x, 1, N * (-1) ** (j + 1)]
        assert L.substitute_univariate(N, x, win, monos, coeffs) == 0
        # (x+h) * w(h, r(h)) = -N D  with w = r - num(h) (unreduced)
        numc = L._num_shifted_coeffs(N, x, j, 0)
        monos = [(0, 1)] + [(a, 0) for a in range(j + 1)]
        coeffs = [1] + [-c for c in numc]
        G = L.substitute_univariate(N, x, win, monos, coeffs)
        assert G == flint.fmpz_poly([-N * x ** (j + 1)])


def test_lattices_small():
    inst = balanced_semiprime(48, 999)
    N, p = inst.N, inst.factors[0]
    X = 2 ** 8
    x = p - 200
    assert L.univariate_attempt(N, x, 3, 7, X)["success"]
    r = L.mixed_attempt(N, x, [(0, 2, 1)], 1, 6, X)
    assert r["success"] and r["good_nontrivial"] >= 1
    Hs = 16
    assert L.hm_attempt(N, x, Hs, X // Hs, 2, 4)["success"]
    s = L.univariate_attempt(N, x, 3, 8, X, use_svp=True)
    assert s["svp_success_hg"]
