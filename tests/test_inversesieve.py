"""Tests for src/inversesieve (exactness of forbidden sets, U(N), statistics, ledger, Theorem S/L lemmas)."""
import itertools
import random
from fractions import Fraction as Fr
from math import gcd, isqrt

import gmpy2
import numpy as np
import pytest

from src.inversesieve.design import (band, weights, farey_set, forbidden_set, candidate_interval,
                                     sqrt_mod_3mod4)
from src.inversesieve.usets import U_of_N, U_bruteforce, is_pair_uncovered, uncovered_partners
from src.inversesieve import stats as st
from src.inversesieve import theory as th


def _primes(lo, hi):
    return [x for x in range(lo, hi) if gmpy2.is_prime(x)]


def test_sqrt_mod():
    for l in (7, 11, 19, 43, 103):
        for a in range(l):
            r = sqrt_mod_3mod4(a, l)
            sq = {x for x in range(l) if x * x % l == a}
            assert set(r) == sq


@pytest.mark.parametrize("kind", ["signed", "signed_dw", "positive"])
def test_forbidden_set_bruteforce(kind):
    rng = random.Random(1)
    V = weights(4, kind)
    for l in band(40):
        for _ in range(5):
            N = rng.randrange(10 ** 6, 10 ** 7)
            if N % l == 0:
                continue
            brute = sorted({x for x in range(l) for (A, B) in V if (B * x * x - A * N) % l == 0})
            assert forbidden_set(N, l, V).tolist() == brute


def test_theorem_S_exact_size_signed():
    """Theorem S(b): for sign-symmetric V (incl. (1,1)) and l = 3 mod 4, l > 2H^2, l !| N: |R_l(N)| = |V|."""
    rng = random.Random(2)
    for H in (3, 5, 7):
        V = weights(H, "signed")
        for l in band(2 * H * H + 1, 6 * H * H):
            if l <= 2 * H * H:
                continue
            assert len(farey_set(l, V)) == len(V)
            for _ in range(10):
                N = rng.randrange(10 ** 8, 10 ** 9)
                if N % l:
                    assert len(forbidden_set(N, l, V)) == len(V)


def test_theorem_S_mult_energy_identity():
    """Theorem S(d): E_x(R_l(N)) = E_x(F_l) for sign-symmetric V."""
    rng = random.Random(3)
    V = weights(6, "signed")
    for l in band(80, 400):
        if l <= 72:
            continue
        N = rng.randrange(10 ** 8, 10 ** 9)
        if N % l == 0:
            continue
        assert st.mult_energy_mod(forbidden_set(N, l, V), l) == st.mult_energy_mod(farey_set(l, V), l)


def test_U_matches_bruteforce_and_true_factor():
    rng = random.Random(4)
    V = weights(3, "signed")
    pr = band(30)
    for _ in range(6):
        p = random.choice(_primes(300, 400))
        q = random.choice([x for x in _primes(p + 1, 2 * p)])
        N = p * q
        U = U_of_N(N, pr, V)
        assert U.tolist() == U_bruteforce(N, pr, V)
        assert (p in set(U.tolist())) == is_pair_uncovered(p, q, pr, V)


def test_partners_match_pair_test():
    V = weights(3, "signed")
    pr = band(30)
    p = 331
    Q = set(uncovered_partners(p, pr, V, p + 1, 2 * p).tolist())
    for q in range(p + 1, 2 * p + 1):
        assert (q in Q) == is_pair_uncovered(p, q, pr, V)


def test_interval_cover_optimal():
    rng = random.Random(5)
    for _ in range(200):
        xs = np.array(sorted(rng.sample(range(60), rng.randint(1, 12))))
        L = rng.randint(1, 10)
        # DP optimum
        best = {}
        def f(i):
            if i >= len(xs):
                return 0
            if i not in best:
                j = i
                while j < len(xs) and xs[j] <= xs[i] + L - 1:
                    j += 1
                best[i] = 1 + f(j)
            return best[i]
        assert st.interval_cover_number(xs, L) == f(0)


def test_ap_triples_bruteforce():
    rng = random.Random(6)
    for _ in range(30):
        xs = np.array(sorted(rng.sample(range(2000), 12)), dtype=np.int64)
        Lm = rng.randint(5, 60)
        brute = sum(1 for a, b, c in itertools.combinations(xs.tolist(), 3) if gcd(b - a, c - a) >= Lm)
        assert st.ap_triples(xs, Lm) == brute


def test_energies_bruteforce():
    rng = random.Random(7)
    l = 31
    A = np.array(sorted(rng.sample(range(l), 7)))
    brute = sum(1 for a, b, c, d in itertools.product(A.tolist(), repeat=4) if (a + b - c - d) % l == 0)
    assert st.energy_mod(A, l) == brute
    xs = np.array([1, 4, 5, 9, 13])
    assert st.energy_int(xs) == sum(1 for a, b, c, d in itertools.product(xs.tolist(), repeat=4) if a + b == c + d)


def test_local_pair_density_exact_over_period():
    """g[d] is the exact density of {x : x, x+d both admissible} over a full CRT period."""
    V = weights(2, "signed")
    pr = [19, 23]
    N = 10 ** 6 + 3
    g = st.local_pair_density(N, pr, V, 6)
    M = 19 * 23
    adm = [all(x % l not in set(forbidden_set(N, l, V).tolist()) for l in pr) for x in range(M)]
    for d in range(7):
        cnt = sum(1 for x in range(M) if adm[x] and adm[(x + d) % M])
        assert abs(g[d] - cnt / M) < 1e-12


def test_lehman_containers():
    N = 1000003 * 1999993
    cen = st.lehman_centres(N, 5)
    a0, b0 = candidate_interval(N)
    assert all(a0 - 1 <= c <= b0 for c in cen)
    half = 7
    xs = np.arange(a0, b0 + 1)
    assert st.container_measure(cen, half, a0, b0) == int(st.in_containers(xs, cen, half).sum())


def test_ledger_exact():
    L = th.ledger(Fr(1, 5))
    assert L["X_N"] == Fr(6, 31) and L["Y_P"] == Fr(10, 31) and L["lambda_P"] == Fr(2, 31)
    assert L["container_mass_fraction_P"] == Fr(-7, 62)
    assert th.moment_order_needed() == Fr(7, 2)
    assert th.max_forgetful_modulus_exponent() == Fr(69, 62)
    assert th.exponent(Fr(0)) == Fr(1, 6) and th.exponent(Fr(1, 4)) == Fr(1, 5)


def test_theorem_L_partner_bound():
    """Lemma L1: for x coprime to M, #{q : x q in [k Delta, (k+1) Delta), x q = d (mod M)} <= Delta/(x M) + 1."""
    rng = random.Random(8)
    for _ in range(300):
        M = rng.randint(1, 60)
        x = rng.randint(50, 400)
        if gcd(x, M) != 1:
            continue
        Delta = rng.randint(1, 5000)
        k = rng.randint(0, 50)
        d = rng.randrange(M)
        cnt = sum(1 for q in range(0, (k + 1) * Delta // x + 2)
                  if k * Delta <= x * q < (k + 1) * Delta and (x * q - d) % M == 0)
        assert cnt <= Delta / (x * M) + 1


def test_theorem_L_toy_rule_mass():
    """Theorem L on a toy population (all prime pairs of S_20): a rule that sees only N mod M and N to
    precision Delta >= P*M and catches every pair must have average mass >= |S_n| / (4 P) per n-bit input.
    The extremal rule C(cell) = {p : (p,q) in cell} is computed exactly."""
    n = 20
    P = 2 ** (n / 2)
    pairs = []
    for p in _primes(int(isqrt(1 << (n - 2))) + 1, isqrt((1 << n) - 1) + 1):
        for q in _primes(p + 1, 2 * p + 1):
            if (1 << (n - 1)) <= p * q < (1 << n):
                pairs.append((p, q))
    for M in (1, 7, 7 * 11):
        Delta = int(P * M) * 2
        cells = {}
        for p, q in pairs:
            N = p * q
            cells.setdefault((N % M, N // Delta), set()).add(p)
        mass = 0
        for (d, k), S in cells.items():
            nN = sum(1 for N in range(max(k * Delta, 1 << (n - 1)), min((k + 1) * Delta, 1 << n)) if N % M == d)
            mass += len(S) * nN
        avg = mass / (1 << (n - 1))
        assert avg >= len(pairs) / (8 * P)
