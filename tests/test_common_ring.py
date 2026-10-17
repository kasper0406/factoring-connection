import random
from math import comb

import pytest

from src.common.metrics import Metrics
from src.common.ring import (FlintQuotient, NonUnit, PolyQuotient, Zmod, binomial_quotient,
                             ext_ring, flint_binomial_power, inv_mod, power_of_linear,
                             truncated)


def naive_binom_mod(N, e, r, c):
    out = [0] * r
    for k in range(e + 1):
        out[k % r] = (out[k % r] + comb(e, k) * pow(c, k // r, N)) % N
    return out


def test_inv_nonunit():
    assert inv_mod(2, 15) * 2 % 15 == 1
    with pytest.raises(NonUnit) as ei:
        inv_mod(6, 15)
    assert ei.value.divisor == 3 and ei.value.is_proper
    with pytest.raises(NonUnit) as ei:
        Zmod(15).inv(0)
    assert ei.value.divisor == 15 and not ei.value.is_proper


@pytest.mark.parametrize("N", [15, 561, 1001, 3 ** 5])
@pytest.mark.parametrize("r,c", [(1, 1), (3, 1), (4, 7), (5, 2)])
def test_binomial_quotient_vs_naive(N, r, c):
    for e in (1, 2, 5, N - 1, N):
        R = binomial_quotient(Zmod(N), r, c)
        assert list(power_of_linear(R, e)) == naive_binom_mod(N, e, r, c)


def test_truncated_and_F():
    for N in (9, 15, 35, 45, 49, 143):
        K = 16
        R = truncated(Zmod(N), K)
        got = list(power_of_linear(R, N))
        assert got == [comb(N, k) % N for k in range(K)]


def test_tower_equals_flat():
    # R_N[s]/(s^2-c)[X]/(X^M - s)  ~=  R_N[X]/(X^{2M} - c) via s = X^M
    N, M, c = 1001, 3, 5
    Z = Zmod(N)
    E = ext_ring(Z, [(-c) % N, 0, 1])
    s = E.gen()
    T = PolyQuotient(E, [E.neg(s)] + [E.zero] * (M - 1) + [E.one])
    a = power_of_linear(T, N)
    flat = naive_binom_mod(N, N, 2 * M, c)
    for i in range(M):
        s0, s1 = a[i]
        assert s0 == flat[i] and s1 == flat[i + M]


def test_counters():
    m = Metrics()
    R = binomial_quotient(Zmod(10007 * 10009, m), 8, 3)
    power_of_linear(R, 1000)
    assert m.counts["q_sqr"] == 9 and m.counts["q_mul"] == 5
    assert m.maxima["quot_dim"] == 8 and m.counts["rn_mul"] > 0


@pytest.mark.parametrize("N", [561, 10007 * 10009, (1 << 61) - 1, (1 << 89) - 1, 2 ** 64 + 13])
def test_flint_agrees(N):
    rng = random.Random(N)
    for f in ([3, 0, 0, 1], [0] * 9 + [1], [1, 2, 3, 4, 1], [N - 5] + [0] * 16 + [1]):
        P = PolyQuotient(Zmod(N), f)
        Q = FlintQuotient(N, f)
        for e in (1, 2, 7, 1000, N):
            assert list(power_of_linear(P, e)) == Q.coeffs(power_of_linear(Q, e))
        x = [rng.randrange(N) for _ in range(len(f) - 1)]
        y = [rng.randrange(N) for _ in range(len(f) - 1)]
        assert list(P.mul(P.from_coeffs(x), P.from_coeffs(y))) == \
            Q.coeffs(Q.mul(Q.from_coeffs(x), Q.from_coeffs(y)))
    assert flint_binomial_power(N, 4, 7, 99) == naive_binom_mod(N, 99, 4, 7)
