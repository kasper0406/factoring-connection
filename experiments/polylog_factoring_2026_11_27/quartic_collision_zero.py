"""Decode from yes/no divisibility of quartic collision counts."""

from math import gcd
from sympy import Poly, Rational, nextprime, symbols

if __package__:
    from .quartic_collision_powers import factors_from_roots
else:
    from quartic_collision_powers import factors_from_roots


def rational_lift(residue, period, bound):
    """Bounded rational reconstruction by the extended Euclidean algorithm."""
    last, current = period, residue
    previous_denominator, denominator = 0, 1
    while current > bound:
        quotient, remainder = divmod(last, current)
        last, current = current, remainder
        previous_denominator, denominator = denominator, previous_denominator - quotient * denominator
    if abs(denominator) > bound or gcd(current, denominator) != 1:
        raise ValueError("No bounded rational reconstruction")
    return Rational(current, denominator)


def decode_zero(modulus, divides):
    """divides(n, j, ell) answers whether ell divides H(n**j)."""
    if modulus < 3 or modulus % 2 == 0:
        raise ValueError("Expected odd modulus >= 3")
    bits = modulus.bit_length()
    # The ith distinct odd prime is at least 2i+1 (with i starting at 1).
    degree_bound, minimum_product = 0, 1
    while minimum_product * (2 * degree_bound + 3) <= modulus:
        minimum_product *= 2 * degree_bound + 3
        degree_bound += 1
    height_product = 5**degree_bound * modulus
    coefficient_bound = 2**degree_bound * height_product
    pairs = degree_bound * (degree_bound - 1) // 2
    bad_product = modulus * 2**pairs * height_product**degree_bound
    threshold = max(bad_product, 2 * coefficient_bound**2)
    degree, period, residues = -1, 1, []
    prime = max(5, bits)
    while period <= threshold:
        prime = int(nextprime(prime))
        if modulus % prime in (0, 1):
            continue
        roots = [j % prime for j in range(1, prime + 1) if divides(modulus, j, prime)]
        if len(roots) == prime or len(roots) < degree:
            continue
        if len(roots) > degree:
            degree, period, residues = len(roots), 1, [0] * (len(roots) + 1)
        coefficients = [1]
        for root in roots:
            extended = [0] + coefficients
            for index, coefficient in enumerate(coefficients):
                extended[index] = (extended[index] - root * coefficient) % prime
            coefficients = extended
        inverse = pow(period, -1, prime)
        for index, coefficient in enumerate(coefficients):
            correction = (coefficient - residues[index]) * inverse % prime
            residues[index] += period * correction
        period *= prime

    lifted = [rational_lift(value, period, coefficient_bound) for value in residues]
    polynomial = Poly.from_list(list(reversed(lifted)), symbols("j"))
    return factors_from_roots(modulus, polynomial.ground_roots())
