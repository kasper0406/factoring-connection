"""Tests for track C baselines (owner factC).  Harness: uses reference factors."""
import random
from math import gcd, prod

import pytest

from src.common.metrics import Metrics
from src.trackC.baselines import (bgs_interval, find_nonunit, interval_predicate,
                                  least_factor_by_search, naive_interval, ps_interval,
                                  U_predicate)
from src.trackC.polyback import PolyRing, multipoint_eval, product_tree_linear


def ref_interval(N, A, B):
    v = prod(range(A, B)) % N
    return v, gcd(v, N) > 1


def least_prime(n):
    d = 2
    while d * d <= n:
        if n % d == 0:
            return d
        d += 1
    return n


@pytest.mark.parametrize("N", [15, 561, 2**61 - 1, (2**61 - 1) * (2**31 - 1), 10007 * 10009])
def test_poly_backends_agree(N):
    rng = random.Random(N)
    for _ in range(10):
        a = [rng.randrange(N) for _ in range(rng.randrange(1, 90))]
        b = [rng.randrange(N) for _ in range(rng.randrange(1, 90))]
        Rf, Rp = PolyRing(N, "flint"), PolyRing(N, "py")
        assert Rf.coeffs(Rf.mul(Rf.make(a), Rf.make(b))) == Rp.coeffs(Rp.mul(Rp.make(a), Rp.make(b)))
        n = rng.randrange(1, 60)
        assert Rf.coeffs(Rf.mullow(Rf.make(a), Rf.make(b), n)) == Rp.coeffs(Rp.mullow(Rp.make(a), Rp.make(b), n))


@pytest.mark.parametrize("backend", ["flint", "py"])
def test_divrem_and_multipoint(backend):
    N = 1000003 * 999983
    rng = random.Random(1)
    R = PolyRing(N, backend)
    for _ in range(10):
        a = R.make([rng.randrange(N) for _ in range(rng.randrange(1, 70))])
        bl = [rng.randrange(N) for _ in range(rng.randrange(1, 30))] + [1]
        b = R.make(bl)
        q, r = R.divrem_monic(a, b)
        assert R.coeffs(R.add(R.mul(q, b), r)) == R.coeffs(a)
        assert R.length(r) < len(bl)
        pts = [rng.randrange(N) for _ in range(rng.randrange(1, 40))]
        assert multipoint_eval(R, a, pts) == [R.evaluate(a, x) for x in pts]


def _check(N, A, B, backend="flint"):
    v, nu = ref_interval(N, A, B)
    for fn in (naive_interval, ps_interval, bgs_interval):
        r = fn(N, A, B, Metrics()) if fn is naive_interval else fn(N, A, B, Metrics(), backend=backend)
        if r["aux_hit"]:
            assert 1 < r["divisor"] and N % r["divisor"] == 0
            assert gcd(r["witness"], N) == r["divisor"]
            continue
        assert r["value"] == v, (fn.__name__, N, A, B)
        assert r["nonunit"] == nu
        if nu:
            assert A <= r["witness"] < B and gcd(r["witness"], N) == r["divisor"] > 1
    ip = interval_predicate(N, A, B, "bgs", Metrics(), backend)
    assert ip["nonunit"] == nu and ip["value"] == v


def test_exhaustive_small():
    rng = random.Random(7)
    for N in range(4, 400):
        for _ in range(3):
            A = rng.randrange(1, N)
            B = rng.randrange(A + 1, min(N, A + 300) + 1)
            _check(N, A, B)


@pytest.mark.parametrize("backend", ["flint", "py"])
def test_random_composites(backend):
    rng = random.Random(11)
    Ns = [101 * 103, 1009 * 1013, 65537 * 65539 * 3 // 3, 7**9, 3 * 11 * 17 * 1009,
          561, 41041, 2**31 - 1, (2**31 - 1) * (2**61 - 1), 1000003 * 1000033]
    for N in Ns:
        for _ in range(6):
            A = rng.randrange(1, 5000)
            B = A + rng.randrange(65, 3000)
            _check(N, A, B, backend)


def test_prefix_predicate_and_search():
    from src.common.instances import balanced_semiprime, multi_prime, prime_power
    for inst in [balanced_semiprime(24, s) for s in range(4)] + [prime_power(24, 2, 0), multi_prime(32, 3, 1)]:
        N, p = inst.N, inst.labels["factors"][0]
        for K in [2, p - 1, p, p + 1, p + 50]:
            assert U_predicate(N, K, "bgs", Metrics())["U"] == (K <= p)
            assert U_predicate(N, K, "ps", Metrics())["U"] == (K <= p)
        r = least_factor_by_search(N, "bgs", Metrics())
        assert r["p"] == p


def test_find_nonunit():
    N = 35
    assert find_nonunit([1, 2, 3, 4, 6, 8], N, Metrics()) is None
    i = find_nonunit([1, 2, 3, 4, 14, 8], N, Metrics())
    assert i == 4
