from fractions import Fraction as F

import gmpy2

from src.polyimpl import primefind, residual


def test_phi_d_matches_definition():
    for p in (5, 7, 101, 65537):
        assert residual.phi_d(1, p) * residual.phi_d(2, p) == p * p - 1
        assert residual.phi_d(4, p) * residual.phi_d(2, p) * residual.phi_d(1, p) == p ** 4 - 1
        assert residual.phi_d(3, p) * residual.phi_d(1, p) == p ** 3 - 1
        assert residual.phi_d(6, p) * residual.phi_d(3, p) * residual.phi_d(2, p) \
            * residual.phi_d(1, p) == p ** 6 - 1


def test_b2_bound_is_floor_cuberoot():
    for k in (20, 26, 32):
        b = residual.b2_bound(k)
        assert b ** 3 <= 2 ** (2 * k) < (b + 1) ** 3


def test_largest_prime_factor_small():
    assert residual.largest_prime_factor(1) == 1
    assert residual.largest_prime_factor(2 * 3 * 3 * 101) == 101
    assert residual.largest_prime_factor(65536) == 2


def test_exposure_profile_example():
    # p = 65537: p - 1 = 2^16 is smooth; p + 1 = 2 * 3 * 10923 (10923 = 3 * 11 * 331).
    prof = residual.exposure_profile(65537, 1000)
    assert prof[1] is True and prof[2] is True


def test_sample_primes_are_prime_and_in_block():
    ps = residual.sample_primes(16, 50)
    assert len(set(ps)) == 50
    assert all(gmpy2.is_prime(p) and 2 ** 15 <= p < 2 ** 16 for p in ps)


def test_primefind_merikoski_and_threshold():
    e = primefind.exponents(F(1, 2), F(18, 19), F(0))
    assert e["best"] == F(19, 36) and e["best"] > F(1, 2)
    # beta >= 1/2 with theta <= 1 can never beat 1/2, whatever alpha >= 0
    for name, (b, t) in primefind.KNOWN.items():
        assert primefind.exponents(b, t, F(0))["best"] >= F(1, 2), name
    # factoring can matter only when 2 beta < theta and the strip route fails
    assert primefind.factoring_matters(F(1, 10), F(3, 10))
    assert not primefind.factoring_matters(F(1, 10), F(9, 10))  # strip suffices
    assert not primefind.factoring_matters(F(1, 2), F(18, 19))
