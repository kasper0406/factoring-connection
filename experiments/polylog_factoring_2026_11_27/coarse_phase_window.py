"""Average cubic or quartic phase magnitudes on a long scale window."""

import sys
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path

from mpmath import mp

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.lattice_window import interval_square_roots
from experiments.polylog_factoring_2026_11_27.stationary_dual import prime_power_base


LOCAL = {'cubic': {-4: 1, 0: 2, 4: 1}, 'quartic': {-1: 2, 0: 5, 1: 2}}


def find_small_component(modulus, cutoff):
    """Trial-divide only through the precision-dependent cutoff."""
    for prime in range(3, min(cutoff, isqrt(modulus)) + 1, 2):
        if modulus % prime == 0:
            component = prime
            while modulus % (component * prime) == 0:
                component *= prime
            return component
    return None


def find_low_frequency_component(modulus, labels, cutoff):
    """Search every low mixed frequency through its integer quadratic quotient."""
    square = modulus**2
    for left, right in combinations(labels, 2):
        last = (cutoff + max(abs(left), abs(right)))**2 // square
        for quotient in range(1, last + 1):
            discriminant = (left - right)**2 + 4 * quotient * square
            root = isqrt(discriminant)
            if root**2 != discriminant:
                continue
            for signed_root in (root, -root):
                frequency = (left + right + signed_root) // 2
                if abs(frequency) <= cutoff:
                    divisor = gcd(frequency - left, modulus)
                    if 1 < divisor < modulus:
                        return divisor
    return None


def find_band_component(modulus, period, labels, center, cutoff):
    """Use a public band to locate a mixed frequency and verify a component."""
    half = pow(2, -1, period)
    for left, right in combinations(labels, 2):
        shift = (left + right) * half % period
        target = (shift**2 - left * right) % period
        low = center - cutoff - shift
        for root in interval_square_roots(period, low, low + 2 * cutoff, target):
            divisor = gcd(root + shift - left, modulus)
            if 1 < divisor < modulus:
                return divisor
    return None


def sum_on_multiples(frequency, modulus, start, stop, step, *, complex_sum=False):
    """Sum a geometric progression over step-multiples in [start, stop)."""
    first = -(-start // step)
    count = (stop - 1) // step - first + 1
    if count <= 0:
        return 0, 0
    phase = frequency * step % modulus
    if not phase:
        return count, count
    period = 2 * modulus
    center = phase * (2 * first + count - 1) % period
    width = phase * count % period
    angle = mp.mpf(center) / modulus
    rotate = mp.expjpi(angle) if complex_sum else mp.cospi(angle)
    total = rotate * mp.sinpi(mp.mpf(width) / modulus)
    return total / mp.sinpi(mp.mpf(phase) / modulus), count


def factored_lines(modulus, local, component):
    complement = modulus // component
    idempotent = complement**2 * pow(complement**2, -1, component**2)
    for left, right in product(local, repeat=2):
        yield left * idempotent + right * (1 - idempotent), local[left] * local[right]


def factored_average(modulus, local, start, stop, component, *, center=0):
    square = modulus**2
    complement = modulus // component
    prime = prime_power_base(component)
    other = prime_power_base(complement)
    divisors = ((1, 1), (prime, -1), (other, -1), (prime * other, 1))
    total = 0
    for frequency, weight in factored_lines(modulus, local, component):
        total += weight * mp.fsum(
            sign * sum_on_multiples(frequency - center, square, start, stop, step, complex_sum=center != 0)[0]
            for step, sign in divisors
        )
    units = sum(sign * sum_on_multiples(0, square, start, stop, step)[1]
                for step, sign in divisors)
    if not units:
        raise ValueError('The interval contains no unit scales')
    return total / units


def validate_window(modulus, family, length, precision):
    if modulus < 15 or modulus % 2 == 0 or family not in LOCAL:
        raise ValueError('Expected an odd modulus with two components and a supported phase')
    if (family == 'cubic' and modulus % 3 == 0) or length < 1 or precision < 1:
        raise ValueError('Expected nondegenerate phase, positive length and positive precision')


def window_average(modulus, local, mass, start, length, precision, center):
    trial_bound = max(32 * mass * precision, max(local) - min(local))
    square = modulus**2
    center %= square
    component = find_small_component(modulus, trial_bound)
    if component is None:
        cutoff = min(square // 2, (mass * precision * square + length - 1) // length)
        if length < trial_bound:
            cutoff = square // 2
        # Quotient enumeration stays within a polynomial factor of block search here.
        if not center and cutoff <= modulus * modulus.bit_length()**2:
            component = find_low_frequency_component(modulus, local, cutoff)
        else:
            component = find_band_component(modulus, square, local, center, cutoff)
    stop = start + length
    if component is not None:
        return factored_average(modulus, local, start, stop, component, center=center), component
    public = mp.fsum(weight**2 * sum_on_multiples(frequency - center, square, start, stop, 1,
                                                 complex_sum=center != 0)[0]
                     for frequency, weight in local.items())
    return public / length, None


def average_squared(modulus, family, start, length, precision, *, center=0):
    """Approximate the modulated mean squared magnitude to error 1/precision."""
    validate_window(modulus, family, length, precision)
    local = LOCAL[family]
    working_bits = 8 * modulus.bit_length() + 2 * length.bit_length() + precision.bit_length() + 48
    with mp.workprec(working_bits):
        return window_average(modulus, local, sum(local.values())**2, start, length, precision, center)


def cubic_magnitude_weights(precision):
    scale = 4 / mp.pi
    local = {0: scale}
    for harmonic in range(1, 12 * precision + 1):
        weight = scale * (-1)**(harmonic + 1) / (4 * harmonic**2 - 1)
        local[4 * harmonic] = local[-4 * harmonic] = weight
    return local


def quartic_magnitude_weights(precision):
    coefficients = [mp.mpf(1)]
    for degree in range(1, (512 * precision).bit_length() + 1):
        coefficients.append(coefficients[-1] * (3 - 2 * degree) / (4 * degree))
    local = {}
    for offset in range(len(coefficients)):
        weight = 2 * mp.fsum(left * right for left, right in
                            zip(coefficients[offset:], coefficients))
        local[offset] = local[-offset] = weight
    return local


def average_magnitude(modulus, family, start, length, precision, *, center=0):
    """Approximate the modulated mean magnitude to error 1/precision."""
    validate_window(modulus, family, length, precision)
    working_bits = 8 * modulus.bit_length() + 2 * length.bit_length() + 8 * precision.bit_length() + 48
    with mp.workprec(working_bits):
        if family == 'cubic':
            local, mass = cubic_magnitude_weights(precision), 9
        else:
            local, mass = quartic_magnitude_weights(precision), 64
        return window_average(modulus, local, mass, start, length, 2 * precision, center)


def average_on_stride(modulus, family, start, length, precision, *, center=0, squared=False):
    """Average at scales start + N*j after modulation in the reduced frequency circle."""
    validate_window(modulus, family, length, precision)
    if gcd(start, modulus) != 1:
        raise ValueError('Expected a unit starting scale')
    center %= modulus
    square = modulus**2
    working_bits = 8 * modulus.bit_length() + 2 * length.bit_length() + 8 * precision.bit_length() + 48
    with mp.workprec(working_bits):
        if squared:
            local = LOCAL[family]
            mass, accuracy = sum(local.values())**2, precision
        elif family == 'cubic':
            local, mass, accuracy = cubic_magnitude_weights(precision), 9, 2 * precision
        else:
            local, mass, accuracy = quartic_magnitude_weights(precision), 64, 2 * precision
        cutoff = min(modulus // 2, (mass * accuracy * modulus + length - 1) // length)
        component = find_small_component(modulus, max(local) - min(local))
        if component is None:
            component = find_band_component(modulus, modulus, local, center, cutoff)
        lines = (factored_lines(modulus, local, component) if component is not None
                 else ((frequency, weight**2) for frequency, weight in local.items()))
        total = mp.fsum(
            weight * mp.expjpi(mp.mpf(2 * (frequency * start % square)) / square)
            * sum_on_multiples(frequency - center, modulus, 0, length, 1, complex_sum=True)[0]
            for frequency, weight in lines
        )
        return total / length, component


def check():
    checks = 0
    for left, right in ((5, 7), (25, 11), (11, 13)):
        modulus = left * right
        square = modulus**2
        left_idempotent = right**2 * pow(right**2, -1, left**2)
        right_idempotent = (1 - left_idempotent) % square
        for family in LOCAL:
            for start, length in ((1, 70), (square - 20, 65), (-20, 150)):
                actual, factor = average_squared(modulus, family, start, length, 64)
                assert factor in (left, right)
                values = []
                for scale in range(start, start + length):
                    if gcd(scale, modulus) != 1:
                        continue
                    x, y = ((scale * frequency % square) / square
                            for frequency in (left_idempotent, right_idempotent))
                    if family == 'cubic':
                        values.append((4 * mp.cospi(4 * x) * mp.cospi(4 * y))**2)
                    else:
                        values.append((5 + 4 * mp.cospi(2 * x)) * (5 + 4 * mp.cospi(2 * y)))
                assert abs(actual - mp.fsum(values) / len(values)) < mp.mpf('1e-10')
                magnitude, found = average_magnitude(modulus, family, start, length, 2)
                assert found in (left, right)
                assert abs(magnitude - mp.fsum(mp.sqrt(value) for value in values) / len(values)) < mp.mpf('0.5')
                checks += 1

    # Catch a mixed low frequency by its quotient without using trial division.
    assert find_low_frequency_component(35, LOCAL['quartic'], 71) in (5, 7)

    coarse_checks = 0
    for left, right in ((7001, 7103), (10007, 10009)):
        modulus = left * right
        for family, local in LOCAL.items():
            for start, length in ((1, modulus), (modulus // 2, modulus)):
                actual, factor = average_squared(modulus, family, start, length, 2)
                assert factor is None
                with mp.workprec(128):
                    control = factored_average(modulus, local, start, start + length, left)
                    assert abs(actual - control) < mp.mpf('0.5')
                coarse_checks += 1
    modulus = 7001 * 7103
    for family in LOCAL:
        actual, factor = average_magnitude(modulus, family, 1, modulus, 1)
        assert factor is None
        with mp.workprec(128):
            control_weights = cubic_magnitude_weights(4) if family == 'cubic' else quartic_magnitude_weights(4)
            control = factored_average(modulus, control_weights, 1, 1 + modulus, 7001)
            assert abs(actual - control) < 1
    print(f'{2 * checks} direct unit averages and {coarse_checks + 2} unfactored long windows matched')


def check_modulation():
    modulus = 25 * 11
    square = modulus**2
    left_idempotent = 11**2 * pow(11**2, -1, 25**2)
    right_idempotent = (1 - left_idempotent) % square
    for center in (-3, 24_580, square - 4):
        for family in LOCAL:
            with mp.workprec(128):
                start, length = square - 10, 55
                terms = []
                for scale in range(start, start + length):
                    if gcd(scale, modulus) != 1:
                        continue
                    x, y = (mp.mpf(scale * idempotent % square) / square
                            for idempotent in (left_idempotent, right_idempotent))
                    value = ((4 * mp.cospi(4 * x) * mp.cospi(4 * y))**2 if family == 'cubic'
                             else (5 + 4 * mp.cospi(2 * x)) * (5 + 4 * mp.cospi(2 * y)))
                    terms.append(value * mp.expjpi(mp.mpf(-2 * center * scale % (2 * square)) / square))
                observed, _ = average_squared(modulus, family, start, length, 64, center=center)
                assert abs(observed - mp.fsum(terms) / len(terms)) < mp.mpf('1e-25')

    left, right = 7001, 7103
    modulus = left * right
    square = modulus**2
    planted = right**2 * pow(right**2, -1, left**2)
    for family, label in (('cubic', 4), ('quartic', 1)):
        center = (label * planted + modulus // 5) % square
        observed, factor = average_squared(modulus, family, 1, 2 * modulus, 1, center=center)
        assert factor in (left, right)
        with mp.workprec(128):
            exact = factored_average(modulus, LOCAL[family], 1, 1 + 2 * modulus, left, center=center)
            assert abs(observed - exact) < mp.mpf('1e-25')

    center = square // 3
    observed, factor = average_squared(modulus, 'cubic', 1, 2 * modulus, 1, center=center)
    assert factor is None
    with mp.workprec(128):
        exact = factored_average(modulus, LOCAL['cubic'], 1, 1 + 2 * modulus, left, center=center)
        assert abs(observed - exact) < 1
    print('6 direct modulated averages, 2 planted bands and 1 empty band matched')


def check_stride():
    modulus = 25 * 11
    square = modulus**2
    left_idempotent = 11**2 * pow(11**2, -1, 25**2)
    right_idempotent = (1 - left_idempotent) % square
    for family, squared in product(LOCAL, (False, True)):
        with mp.workprec(128):
            terms = []
            for index in range(31):
                scale = 3 + modulus * index
                x, y = (mp.mpf(scale * idempotent % square) / square
                        for idempotent in (left_idempotent, right_idempotent))
                value = ((4 * mp.cospi(4 * x) * mp.cospi(4 * y))**2 if family == 'cubic'
                         else (5 + 4 * mp.cospi(2 * x)) * (5 + 4 * mp.cospi(2 * y)))
                terms.append((value if squared else mp.sqrt(value)) * mp.expjpi(mp.mpf(-42 * index) / modulus))
            observed, factor = average_on_stride(modulus, family, 3, 31, 4, center=21, squared=squared)
            assert factor in (25, 11)
            assert abs(observed - mp.fsum(terms) / len(terms)) < mp.mpf('0.25')

    left, right = 7001, 7103
    modulus = left * right
    folded = right * pow(right, -1, left)
    length = isqrt(modulus)
    center = (folded + 10) % modulus
    _, factor = average_on_stride(modulus, 'quartic', 1, length, 1, center=center, squared=True)
    assert factor in (left, right)
    center = modulus // 3
    estimate, factor = average_on_stride(modulus, 'quartic', 1, length, 1, center=center, squared=True)
    assert factor is None
    with mp.workprec(128):
        exact = mp.fsum(
            weight * mp.expjpi(mp.mpf(2 * (frequency % modulus**2)) / modulus**2)
            * sum_on_multiples(frequency - center, modulus, 0, length, 1, complex_sum=True)[0]
            for frequency, weight in factored_lines(modulus, LOCAL['quartic'], left)
        ) / length
        assert abs(exact - estimate) < 1
    print('4 direct stride averages, 1 planted reduced band and 1 empty reduced band matched')


if __name__ == '__main__':
    check()
    check_modulation()
    check_stride()
