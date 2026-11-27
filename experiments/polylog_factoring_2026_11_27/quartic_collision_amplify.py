"""Recover factors from a unit-additive approximation to the normalized quartic mean."""

from fractions import Fraction

from sympy import nextprime

if __package__:
    from .quartic_collision_powers import decode
else:
    from quartic_collision_powers import decode


def collision_padding(n):
    multiplier, collisions = 1, 1
    prime = 5
    while collisions <= 2*n*multiplier:
        prime = int(nextprime(prime))
        if n % prime:
            multiplier *= prime
            character = 1 if prime % 4 == 1 else -1
            collisions *= 3*prime - 5 - character
    return multiplier, collisions


def decode_coarse_mean(n, approximate):
    """approximate(m) returns a rational within one of H(m)/m for each odd m."""
    if n < 1:
        raise ValueError("Expected a positive integer")
    twos = []
    while n % 2 == 0:
        twos.append(2)
        n //= 2
    if n == 1:
        return twos
    multiplier, collisions = collision_padding(n)
    scale = Fraction(n*multiplier, collisions)

    def count(power):
        estimate = Fraction(approximate(multiplier*power))
        return round(scale*estimate) * (power//n)

    return twos + decode(n, count)
