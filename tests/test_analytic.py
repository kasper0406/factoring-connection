"""Tests for src/analytic (phase 17, agent `analytic`)."""
import math
import random
from fractions import Fraction
from math import isqrt

import pytest

from src.analytic.hiary import (hiary_factor, convergents, best_convergent, nearest_int,
                                integer_roots_quadratic, icbrt_ceil, icbrt_floor)
from src.analytic.hull import hull_vertices
from src.analytic.cubic import (weyl_oracle, cubic_decode, zw_gcd, zw_norm, qw_mul, qw_conj,
                                qw_norm, round_to_norm)
from src.common.instances import balanced_semiprime, prime_instance


def _is_prime(n):
    if n < 2:
        return False
    return all(n % k for k in range(2, isqrt(n) + 1))


def test_convergents_last_is_value():
    for num, den in [(355, 113), (10 ** 9 + 7, 12345), (17, 1)]:
        b, q = list(convergents(num, den))[-1]
        assert Fraction(b, q) == Fraction(num, den)


def test_best_convergent_dirichlet():
    # |num/den - b/q| < 1/(q*qmax) style bound (Lemma 3.1 with q' > qmax)
    for num, den, qmax in [(10 ** 12 + 39, 999983 ** 2, 40), (7 ** 20, 3 ** 25, 1000)]:
        b, q, _ = best_convergent(num, den, qmax)
        assert q <= qmax
        assert abs(Fraction(num, den) - Fraction(b, q)) < Fraction(1, q * qmax)


def test_nearest_int_ties_down():
    assert nearest_int(5, 2) == 2 and nearest_int(7, 2) == 3 and nearest_int(-5, 2) == -3
    assert nearest_int(10, 3) == 3 and nearest_int(11, 3) == 4


def test_integer_roots():
    assert integer_roots_quadratic(-6, 1, 1) == [-3, 2]
    assert integer_roots_quadratic(6, -3, 0) == [2]
    assert integer_roots_quadratic(1, 0, 1) == []


def test_icbrt():
    for m in [0, 1, 7, 8, 9, 26, 27, 28, 10 ** 30 + 1]:
        assert icbrt_floor(m) ** 3 <= m < (icbrt_floor(m) + 1) ** 3
        assert (icbrt_ceil(m) - 1) ** 3 < m <= icbrt_ceil(m) ** 3


def test_hiary_exhaustive_small():
    for N in range(4, 3000):
        f = hiary_factor(N)
        if _is_prime(N):
            assert f is None
        else:
            assert f is not None and 1 < f < N and N % f == 0


def test_hiary_semiprimes_and_primes():
    for bits in (32, 40):
        for seed in range(3):
            inst = balanced_semiprime(bits, seed)
            f = hiary_factor(inst.N)
            assert f in inst.labels["factors"]
            assert hiary_factor(prime_instance(bits, seed).N) is None


def test_hiary_block_scan_finds_all_divisors_in_blocks():
    # skip trial division: every divisor in [N^{1/3}, sqrt N] must be found by the block phase
    rng = random.Random(5)
    for _ in range(200):
        p = rng.randrange(200, 2000)
        q = rng.randrange(p, 4 * p)
        N = p * q
        lo = icbrt_ceil(17 * N) + 2
        if lo > p:
            continue
        st = {}
        f = hiary_factor(N, trial_division=False, stats=st, lo=lo, hi=isqrt(N))
        assert f is not None and 1 < f < N and N % f == 0 and st.get("found_in") == "block"


def test_divisors_are_hull_vertices():
    for N in [1009 * 1021, 3 * 5 * 7 * 11 * 13 * 17, 2 ** 20, 999983 * 7]:
        lo, hi = 2, isqrt(N)
        V = set(hull_vertices(N, lo, hi))
        for d in range(lo, hi + 1):
            if N % d == 0:
                assert (d, N // d) in V


def test_hull_vertices_convex_and_minimal():
    N = 123456791
    V = hull_vertices(N, 500, 11000)
    for (x, y) in V:
        assert y == -((-N) // x)
    for a, b, c in zip(V, V[1:], V[2:]):
        cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        assert cross > 0


def test_zw_arithmetic():
    x, y = (3, 5), (2, -7)
    assert zw_norm(x) * zw_norm(y) == qw_norm(qw_mul(tuple(map(Fraction, x)), tuple(map(Fraction, y))))
    assert qw_norm(qw_conj((Fraction(3), Fraction(5)))) == zw_norm(x)
    # gcd of pi_p pi_q and pi_p conj(pi_q) has norm p (7 = N(3+w)? check via search)
    p_el = next((a, b) for a in range(-5, 6) for b in range(-5, 6) if zw_norm((a, b)) == 7)
    q_el = next((a, b) for a in range(-5, 6) for b in range(-5, 6) if zw_norm((a, b)) == 13)
    P = tuple(map(Fraction, p_el))
    Q = tuple(map(Fraction, q_el))
    e1 = tuple(int(c) for c in qw_mul(P, Q))
    e2 = tuple(int(c) for c in qw_mul(P, qw_conj(Q)))
    assert zw_norm(zw_gcd(e1, e2)) == 7


def test_weyl_oracle_crt_and_vanishing():
    # S(a, pq) = S(a qbar, p) S(a pbar, q); S vanishes if a factor is 2 mod 3
    p, q = 7, 13
    N = p * q
    for a in (1, 2, 5, 11):
        lhs = weyl_oracle(a, N)
        rhs = weyl_oracle(a * pow(q, -1, p) % p, p) * weyl_oracle(a * pow(p, -1, q) % q, q)
        assert abs(lhs - rhs) < 1e-9
    for a in (1, 2, 3):
        assert abs(weyl_oracle(a, 11 * 13)) < 1e-9
        assert abs(weyl_oracle(a, 5 * 11)) < 1e-9


def test_weyl_prime_value_formula():
    # S(b, p) = 2 Re(chibar(b) G(chi)), |S| <= 2 sqrt(p), and G^3 = p J  => |S|^2 relation
    for p in (7, 13, 19, 31, 37):
        for b in (1, 2, 3):
            assert abs(weyl_oracle(b, p)) <= 2 * math.sqrt(p) + 1e-9


@pytest.mark.parametrize("p,q", [(7, 13), (1009, 1021), (601, 607), (4003, 4027)])
def test_cubic_decode_factors(p, q):
    N = p * q
    cache = {}

    def orc(a):
        if a not in cache:
            cache[a] = Fraction(weyl_oracle(a, N))
        return cache[a]
    f = cubic_decode(N, orc, random.Random(0), attempts=8)
    assert f in (p, q)


def test_cubic_decode_no_info_when_factor_is_2_mod_3():
    N = 1013 * 1019  # 1013 = 2 mod 3
    f = cubic_decode(N, lambda a: Fraction(weyl_oracle(a, N)), random.Random(0), attempts=2)
    assert f is None


# ---------------------------------------------------------------- ladder (reduction L_j)
from src.analytic.ladder import floor_sum_linear, ladder_factor


def test_floor_sum_linear_bruteforce():
    rng = random.Random(3)
    for _ in range(3000):
        n, m = rng.randrange(0, 60), rng.randrange(1, 80)
        a, b = rng.randrange(-300, 300), rng.randrange(-300, 300)
        assert floor_sum_linear(n, m, a, b) == sum((a * i + b) // m for i in range(n))


@pytest.mark.parametrize("j", [1, 2, 3])
def test_ladder_exhaustive_small(j):
    for N in range(4, 1500):
        f = ladder_factor(N, j)
        if _is_prime(N):
            assert f is None
        else:
            assert f is not None and 1 < f < N and N % f == 0


@pytest.mark.parametrize("j", [1, 2, 3])
def test_ladder_semiprimes(j):
    for seed in range(3):
        inst = balanced_semiprime(32, seed)
        assert ladder_factor(inst.N, j) in inst.labels["factors"]
        assert ladder_factor(prime_instance(32, seed).N, j) is None
