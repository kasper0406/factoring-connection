"""Relative precision controls for the normalized quartic collision mean on odd inputs."""

from fractions import Fraction

from sympy import primerange


def distinct_count(n, approximate):
    """Recover omega(n) from two rational estimates with relative error at most 10%."""
    if n < 1 or n % 2 == 0:
        raise ValueError("Expected a positive odd integer")
    if n == 1:
        return 0
    j = 8 * n.bit_length()
    power = n**j
    ratio = Fraction(approximate(power*power)) / Fraction(approximate(power))
    return nearest_power_of_two_exponent(ratio)


def nearest_power_of_two_exponent(ratio):
    # Round log2(ratio) using exact geometric-midpoint comparisons.
    square = ratio**2
    count = 0
    while square > 2 * 4**count:
        count += 1
    return count


def squarefree_power_mean(n, components, j, inverse_error):
    """Approximate A(n**j) for squarefree odd n of known component count.

    The estimate is an upper bound within relative error 1/inverse_error.
    """
    if n < 1 or n % 2 == 0 or j < 1 or inverse_error < 1 or components < 0:
        raise ValueError("Invalid squarefree-power parameters")
    return power_mean_from_profile(n, {1: components}, j, inverse_error)


def power_mean_from_profile(n, profile, j, inverse_error):
    """Approximate A(n**j) given the number of prime divisors of each multiplicity."""
    if n < 1 or n % 2 == 0 or j < 1 or inverse_error < 1:
        raise ValueError("Expected positive odd n, exponent and inverse error")
    if any(e < 1 or count < 0 for e, count in profile.items()):
        raise ValueError("Invalid exponent profile")
    remaining = dict(profile)
    exact = Fraction(1)
    cutoff = 2 * sum(remaining.values()) * (inverse_error + 1)
    rest = n
    for p in primerange(3, min(n, cutoff) + 1):
        if rest % p:
            continue
        e = 0
        while rest % p == 0:
            rest //= p
            e += 1
        if remaining.get(e, 0) == 0:
            raise ValueError("Small prime valuation inconsistent with profile")
        remaining[e] -= 1
        base = 5*e*j - 2
        character = 1 if p % 4 == 1 else -1
        exact *= Fraction(base*p - 5*e*j - character, p)
    for e, count in remaining.items():
        exact *= (5*e*j - 2)**count
    return exact


def multiplicity_summary(n, approximate):
    """Return omega(n) and a multiplicative estimate of the product of exponents.

    Queries use odd moduli. A 10% relative oracle puts the returned estimate
    between 63/80 and 11/10 times the true exponent product.
    """
    if n < 1:
        raise ValueError("Expected a positive integer")
    bits = n.bit_length()
    part, known_count, known_product = n, 0, 1
    for p in primerange(2, min(n, 16 * bits) + 1):
        if part % p:
            continue
        exponent = 0
        while part % p == 0:
            part //= p
            exponent += 1
        known_count += 1
        known_product *= exponent
    if part == 1:
        return known_count, Fraction(known_product)

    j = 16 * bits
    power = part**j
    first = Fraction(approximate(power))
    second = Fraction(approximate(power*power))
    residual_count = nearest_power_of_two_exponent(second / first)
    estimate = known_product * first / (5*j)**residual_count
    return known_count + residual_count, estimate


def mobius_from_relative_mean(n, approximate):
    count, estimate = multiplicity_summary(n, approximate)
    if estimate > Fraction(7, 5):
        return 0
    return (-1)**count
