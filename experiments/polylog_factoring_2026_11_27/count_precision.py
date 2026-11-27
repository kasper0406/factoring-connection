"""Check additive amplification and a relative estimate for idempotent counts."""

from fractions import Fraction
from math import ceil, isqrt, prod

from sympy import integer_nthroot

from idempotent_count import decode_radix, decode_two_components, radix_bits


def power_sum(n, total, exponent):
    """Evaluate the power-sum recurrence for a pair with product n."""
    previous, current = 2, total
    if exponent == 0:
        return previous
    for _ in range(1, exponent):
        previous, current = current, total * current - n * previous
    return current


def decode_power_count(n, exponent, observation):
    """Invert the monotone count for exactly two prime-power components."""
    if exponent < 1:
        raise ValueError('Expected a positive power')
    excess = observation - n**exponent - 1
    minimum = isqrt(4 * n)
    minimum += minimum * minimum < 4 * n
    left, right = minimum, n + 1
    while left < right:
        middle = (left + right) // 2
        if power_sum(n, middle, exponent) < excess:
            left = middle + 1
        else:
            right = middle
    choices = (left, max(minimum, left - 1))
    nearest = min(choices, key=lambda s: abs(power_sum(n, s, exponent) - excess))
    return decode_two_components(n, n + 1 + nearest)


def padded_modulus(n, exponent_bound):
    """Choose a known factor so absolute error M**exponent_bound is harmless."""
    exponent_bound = Fraction(exponent_bound)
    if n % 2 == 0 or not 0 <= exponent_bound < 1:
        raise ValueError('Expected odd n and an error exponent in [0,1)')
    bits = ceil((exponent_bound * n.bit_length() + 2) / (1 - exponent_bound))
    multiplier = 1 << bits
    return n * multiplier, multiplier


def relative_interval(n, tolerance):
    """Bound the unweighted-selector count using only small trial divisions."""
    tolerance = Fraction(tolerance)
    if n < 2 or not 0 < tolerance <= 1:
        raise ValueError('Expected n >= 2 and relative tolerance in (0,1]')
    bound = ceil(2 * n.bit_length() / tolerance)
    remainder = n
    known = 1
    divisions = 0
    for candidate in range(2, bound + 1):
        if candidate * candidate > remainder:
            if remainder > 1:
                known *= remainder + 1
                remainder = 1
            break
        divisions += 1
        if remainder % candidate:
            continue
        component = 1
        while remainder % candidate == 0:
            remainder //= candidate
            component *= candidate
        known *= component + 1
    low = known * remainder
    high = low if remainder == 1 else low * (1 + tolerance)
    return low, high, divisions


def integer_error(modulus, exponent):
    return int(integer_nthroot(modulus**exponent.numerator, exponent.denominator)[0])


def check():
    pairs = ((5, 7), (9, 5), (101, 109), (503, 937), (101**2, 109),
             (1000000007, 1000000009), (3**51, 5**32))
    power_checks = padding_checks = interval_checks = 0
    largest_power_noise = largest_padding_noise = max_divisions = 0
    for components in pairs:
        n = prod(components)
        for exponent_bound in (Fraction(1, 4), Fraction(2, 5), Fraction(9, 20)):
            exponent = ceil(1 / (1 - 2 * exponent_bound)) + 1
            count = prod(q**exponent + 1 for q in components)
            noise = integer_error(n**exponent, exponent_bound)
            for sign in (-1, 1):
                assert decode_power_count(n, exponent, count + sign * noise) == tuple(sorted(components))
                power_checks += 1
            largest_power_noise = max(largest_power_noise, noise.bit_length())

    for components in ((101, 1000000007), (101**2, 1000000007)):
        n = prod(components)
        exponent = 5 * n.bit_length()
        modulus = n**exponent
        count = prod(q**exponent + 1 for q in components)
        assert decode_power_count(n, exponent, count + isqrt(modulus)) == components

    close_pair = (1000000007, 1000000009)
    n = prod(close_pair)
    for exponent in (1, 2, 8, 32, 128):
        modulus = n**exponent
        public_estimate = modulus + 1 + isqrt(4 * modulus)
        assert abs(prod(q**exponent + 1 for q in close_pair) - public_estimate) <= isqrt(modulus)

    fixtures = (*pairs, (9, 25, 7), (3**25, 5**7, 7**8, 11**14, 13**5, 17**4, 19, 23))
    for components in fixtures:
        n = prod(components)
        radix = 1 << radix_bits(n)
        count = prod(q + radix for q in components)
        for exponent_bound in (Fraction(1, 2), Fraction(9, 10)):
            modulus, multiplier = padded_modulus(n, exponent_bound)
            noise = integer_error(modulus, exponent_bound)
            for sign in (-1, 1):
                padded_count = (multiplier + radix) * count + sign * noise
                recovered = Fraction(padded_count, multiplier + radix)
                assert abs(recovered - count) < Fraction(1, 2)
                assert decode_radix(n, recovered) == tuple(sorted(components))
                padding_checks += 1
            largest_padding_noise = max(largest_padding_noise, noise.bit_length())

        true_count = prod(q + 1 for q in components)
        for tolerance in (Fraction(1, 2), Fraction(1, 10), Fraction(1, 100)):
            low, high, divisions = relative_interval(n, tolerance)
            assert low <= true_count <= high <= low * (1 + tolerance)
            interval_checks += 1
            max_divisions = max(max_divisions, divisions)

    print(f'{power_checks} power-lift checks, maximum noise {largest_power_noise} bits')
    print('Two imbalanced and five near-balanced square-root-error controls passed')
    print(f'{padding_checks} padding checks, maximum noise {largest_padding_noise} bits')
    print(f'{interval_checks} certified relative intervals, maximum {max_divisions} trial divisors')


if __name__ == '__main__':
    check()
