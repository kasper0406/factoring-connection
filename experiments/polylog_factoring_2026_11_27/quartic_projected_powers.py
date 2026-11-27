"""Projected quartic chirp at regular offsets and the exact offset-one power decoder."""

from math import gcd
from sympy import primerange
from mpmath import mp


def regular_quartic_projection(n, offset, inverse_error):
    """Approximate J at a unit outer parameter, gcd(n,offset*(offset-1))=1."""
    if n < 3 or n % 2 == 0 or inverse_error < 1 or gcd(n, offset*(offset-1)) != 1:
        raise ValueError("Expected an odd modulus, positive inverse error and a regular offset")
    rest, count = n, 1
    for p in primerange(3, 65):
        if rest % p:
            continue
        count *= sum((x**4 - 2*x*x + offset) % p == 0 for x in range(p))
        if not count:
            return mp.mpf(0)
        while rest % p == 0:
            rest //= p
    # Every remaining prime contributes at most four roots, hence at most p^(1/3).
    if (2*inverse_error*count)**6 <= (n//rest)**3 * rest:
        return mp.mpf(0)
    offset %= rest
    count *= sum((x**4 - 2*x*x + offset) % rest == 0 for x in range(rest))
    with mp.workprec(n.bit_length() + inverse_error.bit_length() + 32):
        return count / mp.sqrt(n)


def square_part_from_projection(n, approximate):
    """Return (number of distinct primes, largest s with s² dividing n), for odd n.

    approximate(m) estimates the real J at outer parameter and offset one
    with absolute error strictly below 1/2.
    """
    if n < 3 or n % 2 == 0:
        raise ValueError("Expected odd n >= 3")
    diagonal = round(approximate(n**4))
    multiplicity, rem = divmod(diagonal, n)
    if rem or multiplicity <= 0 or multiplicity & (multiplicity-1):
        raise ValueError("Invalid fourth-power projection")
    square_part, rem = divmod(round(approximate(n*n)), multiplicity)
    if rem or square_part <= 0 or n % square_part**2:
        raise ValueError("Invalid square projection")
    return multiplicity.bit_length()-1, square_part
