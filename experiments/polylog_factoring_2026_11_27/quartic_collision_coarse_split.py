"""Extract a large gcd from a short terminal interval of collision differences."""

from math import gcd
from secrets import randbelow

from flint import fmpz_mat
from sympy import Poly, integer_nthroot, symbols

if __package__:
    from .quartic_collision_sampling import sample_interval
else:
    from quartic_collision_sampling import sample_interval


def interval_divisor(n, low, high, dimension):
    if not 0 <= low < high <= n or dimension not in (4, 8):
        raise ValueError("Invalid interval or lattice dimension")
    if high - low <= 2:
        for value in range(low, high):
            divisor = gcd(value, n)
            if 1 < divisor < n:
                return divisor
        return None

    radius = high - low - 1
    rows = [[n] + [0] * (dimension - 1)]
    for index in range(1, dimension):
        row = [0] * dimension
        row[index - 1] = low * radius**(index - 1)
        row[index] = radius**index
        rows.append(row)
    short = [int(entry) for entry in fmpz_mat(rows).lll().tolist()[0]]
    if sum(map(abs, short))**(dimension // 2) >= n:
        raise ValueError("Interval exceeds the certified large-gcd radius")
    coefficients = [entry // radius**i for i, entry in enumerate(short)]
    polynomial = Poly.from_list(coefficients[::-1], symbols("x"))
    for root in polynomial.ground_roots():
        if root.is_Integer and 0 <= root <= radius:
            divisor = gcd(low + int(root), n)
            if 1 < divisor < n:
                return divisor
    return None


def try_coarse_split(n, interval_count, draw=randbelow, *, semiprime=False):
    """The larger radius assumes n has exactly two distinct prime factors, both simple."""
    dimension, exponent = (4, 6) if semiprime else (8, 28)
    width = max(1, int(integer_nthroot(n, exponent)[0]) // 8)
    low, high = sample_interval(n, interval_count, width, draw)
    return interval_divisor(n, low, high, dimension)
