"""Evaluate and check complete cubic and quartic sums at a square modulus."""

import argparse
from collections import Counter
from cmath import exp
from math import cos, fsum, gcd, log, prod, tau


CUBIC_DIFFERENCES = (0, 2, -6, -6)
QUARTIC_DIFFERENCES = (0, -1, 10, 36, 24)


def phase(modulus, residue):
    return exp(1j * tau * (residue % modulus) / modulus)


def polynomial_phases(modulus, scale, initial_differences, length):
    differences = [scale * value % modulus for value in initial_differences]
    for _ in range(length):
        yield differences[0]
        for index in range(len(differences) - 1):
            differences[index] = (differences[index] + differences[index + 1]) % modulus


def direct_cubic(n, scale):
    """Return the normalized sum of e(c(3x-x^3)/n²) over all x modulo n²."""
    if n < 3 or n % 2 == 0 or gcd(n, scale) != 1:
        raise ValueError("use an odd modulus greater than one and a unit scale")
    modulus = n * n
    residues = polynomial_phases(modulus, scale, CUBIC_DIFFERENCES, (modulus + 1) // 2)
    # Opposite residues give conjugate phases; the term at zero occurs once.
    real_sum = fsum(cos(tau * residue / modulus) for residue in residues)
    return (2 * real_sum - 1) / n


def direct_quartic(n, scale):
    """Return the normalized sum of e(c(x^4-2x²)/n²) over all x modulo n²."""
    if n < 3 or n % 2 == 0 or gcd(n, scale) != 1:
        raise ValueError("use an odd modulus greater than one and a unit scale")
    modulus = n * n
    residues = polynomial_phases(modulus, scale, QUARTIC_DIFFERENCES, (modulus + 1) // 2)
    # Opposite residues give equal phases; the term at zero occurs once.
    total = sum(phase(modulus, residue) for residue in residues)
    return (2 * total - 1) / n


def lift_roots(n, polynomial, derivative):
    """Scan modulo n and lift each simple root to n² for the checks."""
    roots = []
    for lower in range(n):
        value = polynomial(lower)
        if value % n == 0:
            correction = -(value // n) * pow(derivative(lower), -1, n) % n
            lift = lower + n * correction
            assert polynomial(lift) % (n * n) == 0
            roots.append(lift)
    return roots


def direct_inverse(n, scale, inverses):
    modulus = n * n
    return sum(phase(modulus, x + scale * scale * inverse) for x, inverse in inverses) / n


def quartic_idempotents(n, components):
    modulus = n * n
    idempotents = []
    for component in components:
        local = component * component
        cofactor = modulus // local
        idempotents.append(cofactor * pow(cofactor, -1, local))
    return idempotents


def quartic_product(n, scale, idempotents):
    return prod(1 + 2 * phase(n * n, -scale * value) for value in idempotents)


def check_quartic_spectrum(n, components, idempotents):
    values = [log(abs(quartic_product(n, 1 + n * digit, idempotents))) for digit in range(n)]
    for component, idempotent in zip(components, idempotents):
        coefficient = sum(value * phase(n, -idempotent * digit)
                          for digit, value in enumerate(values)) / n
        assert gcd(idempotent, n) == n // component
        lower_bound = 0.24 if component >= 5 else 0.16
        assert abs(coefficient) > lower_bound


def nonunit_gauss_values(n):
    modulus = n * n
    square_counts = {}
    values = {}
    for frequency in range(modulus):
        common = gcd(frequency, modulus)
        if common == 1:
            continue
        reduced = modulus // common
        if reduced not in square_counts:
            square_counts[reduced] = Counter(x * x % reduced for x in range(reduced))
        primitive = frequency // common
        values[frequency] = common * sum(count * phase(reduced, primitive * square)
                                         for square, count in square_counts[reduced].items())
    return values


def quartic_gauss_correction(n, scale, nonunit_values):
    modulus = n * n
    inverse = pow(4 * scale, -1, modulus)
    multiples = 0j
    partial = 0j
    for frequency, value in nonunit_values.items():
        kernel = phase(modulus, -inverse * frequency * frequency - frequency - scale)
        term = kernel * (value - n) / modulus
        if frequency % n == 0:
            multiples += term
        else:
            assert 1 < gcd(frequency, n) < n
            partial += term
    return multiples, partial


def check(components):
    n = prod(components)
    modulus = n * n
    scales = (1, 2, 1 + n, 1 + n * (n // 2), modulus - 1)
    quartic_roots = lift_roots(n, lambda a: a * a * a - a, lambda a: 3 * a * a - 1)
    assert len(quartic_roots) == 3 ** len(components)
    idempotents = quartic_idempotents(n, components)
    nonunit_values = nonunit_gauss_values(n) if n <= 105 else {}
    cubic_roots = []
    inverses = []
    if n % 3:
        cubic_roots = lift_roots(n, lambda a: a * a - 1, lambda a: 2 * a)
        assert len(cubic_roots) == 2 ** len(components)
        inverses = [(x, pow(x, -1, modulus)) for x in range(1, modulus) if gcd(x, n) == 1]

    for scale in scales:
        quartic = direct_quartic(n, scale)
        stationary_quartic = sum(phase(modulus, scale * (a ** 4 - 2 * a * a))
                                 for a in quartic_roots)
        assert abs(quartic - stationary_quartic) < 5e-9
        assert abs(quartic - quartic_product(n, scale, idempotents)) < 5e-9
        assert abs(quartic) >= 1 - 5e-9
        if nonunit_values:
            multiples, partial = quartic_gauss_correction(n, scale, nonunit_values)
            assert abs(multiples - 2 ** len(components) * phase(modulus, -scale)) < 5e-9
            assert abs(quartic - (1 + multiples + partial)) < 5e-9
        if cubic_roots:
            cubic = direct_cubic(n, scale)
            stationary_cubic = sum(phase(modulus, 2 * scale * a) for a in cubic_roots)
            assert abs(cubic - stationary_cubic) < 5e-9
            assert abs(cubic - direct_inverse(n, scale, inverses)) < 5e-9

    check_quartic_spectrum(n, components, idempotents)
    cubic_checks = len(scales) if cubic_roots else 0
    gauss_checks = len(scales) if nonunit_values else 0
    print(f"N={n} cubic_checks={cubic_checks} quartic_checks={len(scales)} "
          f"gauss_checks={gauss_checks} "
          f"phase_evaluations_per_direct_sum={(modulus + 1) // 2}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modulus', type=int, help='evaluate directly without factoring')
    parser.add_argument('--scale', type=int, default=1)
    parser.add_argument('--phase', choices=('cubic', 'quartic'), default='cubic')
    args = parser.parse_args()
    if args.modulus is None:
        for components in ((5,), (3, 5), (5, 7), (9, 5), (5, 11), (3, 5, 7),
                           (25, 7), (5, 49), (5, 7, 11)):
            check(components)
        return

    evaluator = direct_cubic if args.phase == 'cubic' else direct_quartic
    try:
        value = evaluator(args.modulus, args.scale)
    except ValueError as error:
        parser.error(str(error))
    print(f"normalized_sum={value} magnitude={abs(value)} "
          f"phase_evaluations={(args.modulus * args.modulus + 1) // 2}")


if __name__ == '__main__':
    main()
