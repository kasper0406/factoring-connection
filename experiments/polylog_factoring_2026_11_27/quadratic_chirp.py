"""Projected quadratic chirps and stable recovery from unit chirps."""

import sys
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path

from mpmath import mp
from sympy import I, Matrix, Rational, eye, jacobi_symbol

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.coarse_phase_window import LOCAL, find_small_component


RAW = {'cubic': ({-2: 1, 2: 1}, 3), 'quartic': ({-1: 2, 0: 1}, 4)}


def phase(argument, modulus):
    return mp.expjpi(mp.mpf(2 * (argument % modulus)) / modulus)


def gauss_coefficient(modulus, quadratic):
    return int(jacobi_symbol(quadratic, modulus)) * (1 if modulus % 4 == 1 else 1j)


def projected_actual_raw(modulus, family, quadratic, offset, precision):
    """Return a 1/precision estimate or a split on a squarefree semiprime."""
    if modulus < 15 or modulus % 2 == 0 or gcd(quadratic, modulus) != 1 or precision < 1:
        raise ValueError('Expected an odd semiprime, a unit reduced coefficient and positive precision')
    if family == 'cubic' and modulus % 3 == 0:
        return None, 3
    weights, degree = RAW[family]
    component = find_small_component(modulus, (2 * degree**2 * precision)**2)
    if component is not None:
        return None, component
    for label, weight in weights.items():
        common = gcd(offset + label, modulus)
        if 1 < common < modulus:
            return None, common
        if common == modulus:
            linear = (offset + label) // modulus
            argument = -linear**2 * pow(4 * quadratic, -1, modulus)
            with mp.workprec(4 * modulus.bit_length() + precision.bit_length() + 40):
                return weight**2 * gauss_coefficient(modulus, quadratic) * phase(argument, modulus), None
    return mp.mpc(0), None


def clean_unit_chirps(modulus, target, label, rank):
    """Choose 2*rank unit chirps with no deleted labels, or return a split."""
    common = gcd(target, modulus)
    if common != 1:
        if common < modulus:
            return [], common
        raise ValueError('Expected a unit target scale')
    count = 2 * rank
    if count * (count + 2) >= modulus:
        return [], find_small_component(modulus, isqrt(modulus))

    square = modulus**2
    for base in (1, count + 1, count + 2):
        queries = []
        for index in range(1, count + 1):
            d = base * index
            common = gcd(d, modulus)
            if common != 1:
                return [], common
            quadratic = pow(4 * d, -1, square)
            offset = -target * pow(2 * d, -1, square) % square
            divisors = tuple(gcd(offset + frequency, modulus) for frequency in (-label, 0, label))
            for divisor in divisors:
                if 1 < divisor < modulus:
                    return [], divisor
            if modulus in divisors:
                break
            queries.append((quadratic, offset, d))
        if len(queries) == count:
            return queries, None
    raise ArithmeticError('Two excluded coefficients cannot spoil three disjoint query sets')


def extrapolate_zero(samples, error):
    """Recover a missing initial value from 2R rational samples of R unit modes."""
    rank = len(samples) // 2
    if rank < 1 or len(samples) != 2 * rank or error <= 0:
        raise ValueError('Expected 2R samples and a positive error bound')
    matrix = Matrix(rank, rank, lambda row, col: samples[row + col + 1])
    target = Matrix(samples[:rank])
    adjoint = matrix.conjugate().T
    normal = adjoint * matrix + error**2 * eye(rank)
    coefficients = -normal.inv(method='DM') * adjoint * target
    return -(target.T * coefficients)[0]


def idempotents(components):
    square = prod(components)**2
    return [(square // q**2) * pow(square // q**2, -1, q**2) for q in components]


def collision_count(prime, exponent):
    """Number of congruent pairs for x^4-2x^2 modulo one odd prime power."""
    sign = 1 if prime % 4 == 1 else -1
    return prime**(exponent - 1) * ((5 * exponent - 2) * prime - 5 * exponent - sign)


def literal_quartic_centered_reference(prime_powers, quadratic):
    """Factor-supplied value of the literal all-scale, squared-magnitude chirp."""
    components = [p**exponent for p, exponent in prime_powers]
    square = prod(components)**2
    slope = -pow(4 * quadratic, -1, square)
    return mp.fprod(collision_count(*power) + 4 * phase(slope * projector, square)
                    for power, projector in zip(prime_powers, idempotents(components)))


def literal_quartic_queries(modulus, target):
    """Select five centered literal chirps, or split on a nonunit coefficient."""
    square = modulus**2
    queries = []
    for index in range(1, 6):
        common = gcd(index * target, modulus)
        if common == modulus:
            raise ValueError('Expected a unit target scale')
        if common != 1:
            return [], common
        queries.append(pow(4 * index * target, -1, square))
    return queries, None


def quartic_from_literal(modulus, target, samples, error):
    """Recover magnitude from five literal chirps already divided by sqrt(N)."""
    if len(samples) != 5 or not 0 < error < 1:
        raise ValueError('Expected five samples and a positive error below one')
    z = phase(-target, modulus**2)
    public = 16 / mp.sqrt(modulus)
    values = [value - public * z**index for index, value in enumerate(samples, 1)]
    h1, h2, h3, h4 = (right - left for left, right in zip(values, values[1:]))
    threshold = mp.root(error**2, 3)
    if max(abs(h2), abs(h3)) < threshold:
        squared = 9 * (5 + 4 * z.real)
    else:
        slope = (h3 + z * h1) / h2 if abs(h2) >= abs(h3) else (h4 + z * h2) / h3
        squared = 25 + 20 * slope.real + 4 * (2 * z.real + (slope**2 / z).real - 2)
    return mp.sqrt(min(81, max(1, squared)))


def polynomial(value, family):
    if family == 'cubic':
        return value * (3 - value**2)
    return value**2 * (value**2 - 2)


def projected_lift_reference(components, family, quadratic, offset):
    """Independent control using all quotient roots and the local quadratic means."""
    left, right = components
    modulus = left * right
    roots = [[value for value in range(p) if (polynomial(value, family) + offset) % p == 0]
             for p in components]
    parts = []
    for u, v in product(*roots):
        value = u + left * ((v - u) * pow(left, -1, right) % right)
        linear = (polynomial(value, family) + offset) // modulus
        derivative = 3 * (1 - value**2) if family == 'cubic' else 4 * value * (value**2 - 1)
        factors = []
        for p in components:
            coefficient = -pow(4 * quadratic * (modulus // p), -1, p)
            if derivative % p:
                factors.append(gauss_coefficient(p, coefficient) / mp.sqrt(p))
            else:
                factors.append(phase(coefficient * linear**2, p))
        parts.append(mp.fprod(factors))
    return gauss_coefficient(modulus, quadratic) * mp.fsum(parts)


def check_projection():
    for family, offset in product(RAW, (0, 1, 2, 5, -2, 1 + 5003 * 5101)):
        modulus = 5003 * 5101
        observed, factor = projected_actual_raw(modulus, family, 1, offset, 2)
        assert factor is None
        expected = projected_lift_reference((5003, 5101), family, 1, offset)
        assert abs(observed - expected) < mp.mpf(1) / 2
    assert projected_actual_raw(35, 'cubic', 1, 0, 8)[1] in (5, 7)


def check_extrapolation():
    modulus = 101 * 109
    i, j = idempotents((101, 109))
    assert clean_unit_chirps(modulus, 658, 1, 4)[1] == 109
    for family, target in product(LOCAL, (1, 17, 19)):
        local = LOCAL[family]
        label = max(local)
        square = modulus**2
        weights = [(a * i + b * j, local[a] * local[b]) for a, b in product(local, repeat=2)]
        queries, factor = clean_unit_chirps(modulus, target, label, 4)
        assert factor is None and len(queries) == 8
        values = []
        scale = 1 << 28
        for _, _, d in queries:
            sample = mp.fsum(weight * phase(target * frequency - d * frequency**2, square)
                             for frequency, weight in weights)
            real, imag = (int(mp.nint(part * scale)) for part in (sample.real, sample.imag))
            values.append(Rational(real, scale) + I * Rational(imag, scale))
        predicted = complex(extrapolate_zero(values, Rational(1, scale)))
        expected = mp.fsum(weight * phase(target * frequency, square) for frequency, weight in weights)
        assert abs(predicted - complex(expected)) < 1e-4


def check_literal():
    import numpy as np

    def complete(modulus, family):
        square = modulus**2
        residues = np.fromiter((polynomial(value, family) % square for value in range(square)),
                               dtype=np.int64, count=square)
        return modulus * np.fft.ifft(np.bincount(residues, minlength=square))

    raw_checks = 0
    for components in ((5, 7), (7, 11), (11, 13)):
        modulus = prod(components)
        square = modulus**2
        grid = np.arange(square, dtype=np.int64)
        for family, offset in product(RAW, (0, 1, 2, 17)):
            weight = np.exp(2j * np.pi * ((modulus * grid**2 + offset * grid) % square) / square)
            literal = np.sum(complete(modulus, family) * weight) / (modulus * np.sqrt(modulus))
            assert abs(literal - complex(projected_lift_reference(components, family, 1, offset))) < 1e-9
            raw_checks += 1

    squared_checks = 0
    for prime_powers in (((5, 1), (7, 1)), ((5, 2), (7, 1)), ((3, 3), (5, 2))):
        modulus = prod(p**exponent for p, exponent in prime_powers)
        square = modulus**2
        grid = np.arange(square, dtype=np.int64)
        magnitude = abs(complete(modulus, 'quartic'))**2
        for quadratic in (1, 2, 4):
            literal = np.sum(magnitude * np.exp(2j * np.pi * ((quadratic * grid**2) % square) / square)) / modulus
            expected = literal_quartic_centered_reference(prime_powers, quadratic)
            assert abs(literal - complex(expected)) < 1e-7
            squared_checks += 1
    print(f'{raw_checks} literal projected chirps and {squared_checks} literal squared chirps matched')


def check_five_queries():
    checks = 0
    for prime_powers in (((101, 1), (109, 1)), ((11, 2), (10007, 1)), ((29, 3), (31, 2))):
        left, right = (p**exponent for p, exponent in prime_powers)
        modulus = left * right
        for target in (1, 17, left**2 + right**2):
            queries, factor = literal_quartic_queries(modulus, target)
            assert factor is None
            error = mp.mpf(2)**-28
            samples = [literal_quartic_centered_reference(prime_powers, quadratic) / mp.sqrt(modulus)
                       + error * phase(index, 7) for index, quadratic in enumerate(queries)]
            observed = quartic_from_literal(modulus, target, samples, error)
            expected = mp.fprod(abs(1 + 2 * phase(-target * component, modulus**2))
                                for component in idempotents((left, right)))
            assert abs(observed - expected) < mp.mpf('0.01')
            checks += 1
    print(f'{checks} noisy five-chirp point evaluations matched')


if __name__ == '__main__':
    with mp.workprec(100):
        check_projection()
        check_extrapolation()
        check_literal()
        check_five_queries()
    print('12 unfactored projections and 6 noisy extrapolations matched')
