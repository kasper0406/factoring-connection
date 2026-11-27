"""Split an odd composite using approximate quartic collision counts by difference."""

from fractions import Fraction
from math import gcd
from secrets import randbelow


def sample_interval(n, interval_count, minimum_width, draw=randbelow):
    """Count ordered collisions with (x-y) mod n in each requested half-open interval.

    The callback takes (n, low, high, inverse_error) and returns a nonnegative
    rational estimate, with relative error at most 1/inverse_error (zero at zero).
    """
    if not 1 <= minimum_width <= n:
        raise ValueError("Invalid minimum interval width")
    if n < 3 or n % 2 == 0:
        raise ValueError("Expected odd n >= 3")
    low, high = 0, n
    inverse_error = 16 * n.bit_length()
    while high - low >= 2 * minimum_width:
        middle = (low + high) // 2
        left = Fraction(interval_count(n, low, middle, inverse_error))
        right = Fraction(interval_count(n, middle, high, inverse_error))
        if min(left, right) < 0 or left + right == 0:
            raise ValueError("Invalid estimates at a positive-mass interval")
        prob = left / (left + right)
        if draw(prob.denominator) < prob.numerator:
            high = middle
        else:
            low = middle
    return low, high


def sample_difference(n, interval_count, draw=randbelow):
    return sample_interval(n, interval_count, 1, draw)[0]


def try_split(n, interval_count, draw=randbelow):
    difference = sample_difference(n, interval_count, draw)
    divisor = gcd(difference, n)
    return divisor if 1 < divisor < n else None
