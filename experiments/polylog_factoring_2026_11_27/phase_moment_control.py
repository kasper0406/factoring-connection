"""Check joint low moments of the cubic and quartic stationary traces."""

from collections import Counter
from fractions import Fraction
from math import cos, fsum, gcd, prod, tau


def laurent_coefficients(stencil):
    """Build a local Laurent polynomial, using squared magnitudes for quartics."""
    coefficients = Counter({0: 1})
    for phase, scale, count in stencil:
        term = {-scale: 1, scale: 1} if phase == 'cubic' else {-scale: 2, 0: 5, scale: 2}
        for _ in range(count):
            next_coefficients = Counter()
            for exponent, coefficient in coefficients.items():
                for shift, weight in term.items():
                    next_coefficients[exponent + shift] += coefficient * weight
            coefficients = next_coefficients
    return coefficients


def unit_grid_moment(prime, component, coefficients):
    modulus = component * component
    total = Fraction(0)
    for exponent, coefficient in coefficients.items():
        if exponent % modulus == 0:
            total += coefficient
        elif exponent % (modulus // prime) == 0:
            total -= Fraction(coefficient, prime - 1)
    return total


def literal_local(prime, component, stencil):
    modulus = component * component
    values = []
    for residue in range(modulus):
        if gcd(residue, prime) != 1:
            continue
        term = 1.0
        for phase, scale, count in stencil:
            cosine = cos(tau * ((scale * residue) % modulus) / modulus)
            term *= (2 * cosine if phase == 'cubic' else 5 + 4 * cosine) ** count
        values.append(term)
    return fsum(values) / len(values)


if __name__ == '__main__':
    stencils = ((('cubic', 1, 4),), (('quartic', 1, 3),),
                (('cubic', 1, 3), ('cubic', 2, 2), ('quartic', 3, 1)),
                (('quartic', 1, 1), ('quartic', 2, 2), ('cubic', 3, 1)),
                (('quartic', 2, 2), ('quartic', 3, 2), ('cubic', 1, 1)))
    fixtures = (((13, 13), (17, 17)), ((17, 17), (19, 19), (23, 23)),
                ((5, 25), (19, 19)), ((23, 23), (29, 29)),
                ((7, 49), (13, 13), (31, 31)))
    checks = 0
    for components in fixtures:
        for stencil in stencils:
            coefficients = laurent_coefficients(stencil)
            degree = max(map(abs, coefficients))
            assert all(degree < q * q // p for p, q in components)
            local_values = []
            for p, q in components:
                observed = literal_local(p, q, stencil)
                exact = unit_grid_moment(p, q, coefficients)
                assert abs(observed - exact) <= 1e-10 * max(1, abs(exact))
                local_values.append(exact)
            assert prod(local_values) == coefficients[0]**len(components)
            checks += 1

    alias = (('quartic', 1, 5),)
    coefficients = laurent_coefficients(alias)
    assert unit_grid_moment(5, 5, coefficients) != coefficients[0]
    assert unit_grid_moment(7, 7, coefficients) == coefficients[0]
    print(f'{checks} joint moments and two alias boundary cases checked')
