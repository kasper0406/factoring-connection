"""Decode an odd integer using an exact quartic collision count oracle."""

from math import factorial

from sympy import Poly, QQ, Rational, symbols


def decode(modulus, count):
    if modulus < 3 or modulus % 2 == 0:
        raise ValueError("Expected an odd integer at least three")
    variable = symbols("j")
    degree_bound = modulus.bit_length() - 1
    differences = []
    polynomial = Poly(0, variable, domain=QQ)
    basis = Poly(1, variable, domain=QQ)
    for order in range(degree_bound + 1):
        power = modulus ** (order + 1)
        value, remainder = divmod(count(power), power // modulus)
        if remainder:
            raise ValueError("Inconsistent collision count")
        # Update the right edge of the forward difference table.
        for index, previous in enumerate(differences):
            differences[index], value = value, value - previous
        differences.append(value)
        polynomial += basis * Rational(value, factorial(order))
        basis *= Poly(variable - order - 1, variable, domain=QQ)

    return factors_from_roots(modulus, polynomial.ground_roots())


def divisors_from_roots(modulus, roots):
    degree_bound = modulus.bit_length() - 1
    candidates = set()
    for root in roots:
        for exponent in range(1, degree_bound + 1):
            scaled = 5 * exponent * root
            if scaled == 2:
                continue
            for character in (-1, 1):
                candidate = (scaled + character) / (scaled - 2)
                if candidate.q != 1 or candidate < 3:
                    continue
                prime = int(candidate)
                if modulus % prime == 0:
                    candidates.add(prime)
    return sorted(candidates)


def factors_from_roots(modulus, roots):
    # Every true prime occurs. Its composite multiples become redundant.
    remaining = modulus
    factors = []
    for candidate in divisors_from_roots(modulus, roots):
        while remaining % candidate == 0:
            remaining //= candidate
            factors.append(candidate)
    if remaining != 1:
        raise ValueError("Collision values failed to decode")
    return factors
