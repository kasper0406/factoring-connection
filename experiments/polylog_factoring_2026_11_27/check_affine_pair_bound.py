"""Verify the uniform constant and the pair-count identity on small moduli."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod


def check_constant():
    bound = Fraction(1)
    for i in range(1, 13):
        numerator = 10_000 + Fraction(isqrt((2*i+1) * 10**8) + 1, i)
        bound *= numerator / 14_142
    assert bound < 9
    assert (Fraction(14_142, 10_000) - 1)**2 > Fraction(27, 169)
    return bound


def check(factors):
    n = prod(factors)
    roots = [r for r in range(n) if r*r % n == 1]
    units = [r for r in range(n) if gcd(r, n) == 1]
    widths = [1 << j for j in range(n.bit_length()) if 1 << j <= n//4]
    for width in widths:
        observed = 0
        for a in units:
            images = [a*r % n for r in roots]
            for x, y in combinations(images, 2):
                distance = (x-y) % n
                observed += 2*(max(0, width-distance) + max(0, width-(n-distance)))
        predicted = Fraction(0)
        for mask in range(1, 1 << len(factors)):
            m = prod(factor for i, factor in enumerate(factors) if mask & (1 << i))
            d = n//m
            reduced_units = [x for x in range(1, m) if gcd(x, m) == 1]
            overlap = sum(
                max(0, width - d*x) + max(0, width - (n - d*x)) for x in reduced_units
            )
            predicted += len(roots)*Fraction(overlap, len(reduced_units))
        assert Fraction(observed, len(units)) == predicted
    return len(widths)


if __name__ == "__main__":
    bound = check_constant()
    widths = sum(check(factors) for factors in [(3, 5), (9, 5), (3, 5, 7), (9, 5, 7), (3, 5, 7, 11)])
    print(f"Uniform constant below {float(bound):.3f}; {widths} exact pair-count checks.")
