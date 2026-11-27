"""Count roots of one in a window given the odd prime-power components."""

from bisect import bisect_left
from math import prod


def root_count(prime_powers, start, stop):
    """Count in [start, stop) given distinct, pairwise coprime odd prime powers."""
    components = tuple(prime_powers)
    n = prod(components)
    if not components or not 0 <= start <= stop <= n:
        raise ValueError("expected a nonempty factorization and a canonical window")

    def signs(parts):
        sums = [0]
        for q in parts:
            quotient = n // q
            weight = quotient * pow(quotient, -1, q)
            sums = [(value + sign*weight) % n for value in sums for sign in (-1, 1)]
        return sums

    middle = len(components)//2
    right = sorted(signs(components[middle:]))
    total = 0
    for left in signs(components[:middle]):
        for wrap in (0, n):
            total += (bisect_left(right, stop + wrap - left)
                      - bisect_left(right, start + wrap - left))
    return total
