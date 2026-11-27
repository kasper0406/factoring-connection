"""Check second moments of the remaining terms after quartic Gauss completion."""

from math import fsum, sqrt

from stationary_cubic import direct_quartic, nonunit_gauss_values, quartic_gauss_correction


def predicted_absolute_mean(p, q):
    """The average absolute summand including the zero terms of the correction."""
    def one_side(prime, other):
        local_magnitude = sqrt(prime + (prime % 4 == 3))
        return (1 - 1 / other) * (prime - 1) * (local_magnitude + 1) / prime**2

    return p * q * (one_side(p, q) + one_side(q, p))


def check(p, q):
    n = p * q
    modulus = n * n
    gauss_values = nonunit_gauss_values(n)
    active = {s: value for s, value in gauss_values.items() if s % n}
    assert len(active) == (p + q - 2) * n
    magnitudes = [abs(value - n) for value in active.values()]
    second = fsum(value**2 for value in magnitudes) / modulus
    predicted_second = 4 * modulus * (1 - 1 / p) * (1 - 1 / q)
    first = fsum(magnitudes) / modulus
    assert abs(second - predicted_second) <= 1e-11 * predicted_second
    assert abs(first - predicted_absolute_mean(p, q)) <= 1e-11 * first

    # The same second moments apply at every unit scale; the summand phase
    # only rotates. Check the residual still has constant size in two sums.
    for scale in (1, n + 1):
        public, residual = quartic_gauss_correction(n, scale, gauss_values)
        assert abs(residual) <= 4 + 1e-10
        assert abs(1 + public + residual - direct_quartic(n, scale)) < 1e-9
    print(f'N={n}: second moment {second:.6g}, '
          f'best one-term second moment {first**2:.6g}, active fraction {len(active) / modulus:.6g}')


if __name__ == '__main__':
    for p, q in ((3, 5), (5, 7), (7, 11), (11, 13), (13, 17)):
        check(p, q)
