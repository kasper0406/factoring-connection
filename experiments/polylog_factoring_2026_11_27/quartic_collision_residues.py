"""Recover exact quartic collision counts using small auxiliary prime residues."""

from itertools import count
from math import isqrt


def recover_count(modulus, residue):
    """residue(m, prime) returns the collision count at m modulo prime."""
    if modulus < 1 or modulus % 2 == 0:
        raise ValueError("Expected a positive odd modulus")
    value, period = 0, 1
    for prime in count(7, 2):
        if modulus % prime == 0 or any(
            prime % d == 0 for d in range(3, isqrt(prime) + 1, 2)
        ):
            continue
        character = 1 if prime % 4 == 1 else -1
        local_count = 3 * prime - 5 - character
        observed = residue(modulus * prime, prime)
        target = observed * pow(local_count, -1, prime) % prime
        digit = (target - value) * pow(period, -1, prime) % prime
        value += period * digit
        period *= prime
        if period > modulus**2:
            return value
