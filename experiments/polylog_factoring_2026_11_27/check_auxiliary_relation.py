"""Extract a verified root or factor from an explicit integer polynomial."""

from math import gcd
from sympy import Poly, symbols


def extract(n, center, y, multiplicity, auxiliary):
    """Use the lowest nonstandard coefficient layer in the monic expansion."""
    if n <= 1 or multiplicity < 1:
        raise ValueError("Expected n > 1 and positive multiplicity")
    t, = auxiliary.gens
    quadratic = Poly((t + center)**2 - y, t, domain="ZZ")
    remaining = auxiliary.set_domain("ZZ")
    best = multiplicity
    relation = None
    index = 0
    while not remaining.is_zero and index < best:
        remaining, linear = remaining.div(quadratic)
        slope, intercept = map(int, (linear.nth(1), linear.nth(0)))
        weight = index
        while weight < best and slope % n == 0 and intercept % n == 0:
            if slope == intercept == 0:
                break
            slope //= n
            intercept //= n
            weight += 1
        if weight < best and (slope % n or intercept % n):
            best = weight
            relation = slope, intercept
        index += 1
    if relation is None:
        return None
    slope, intercept = relation
    divisor = gcd(slope, n)
    if 1 < divisor < n:
        return "factor", divisor
    if divisor == 1:
        root = (center - intercept * pow(slope, -1, n)) % n
        if root**2 % n == y % n:
            return "root", root
    raise ValueError("Coefficient fails the promised root relation")


if __name__ == "__main__":
    t = symbols("t")
    counts = dict(root=0, factor=0, standard=0)
    for n, prime in ((35, 5), (77, 7), (125, 5), (175, 7)):
        for center in (-3, 0, 8):
            root = 2
            y = root**2 % n
            f = (t + center)**2 - y
            for m in range(1, 5):
                standard = (f + n)**m
                assert extract(n, center, y, m, Poly(standard, t)) is None
                counts["standard"] += 1
                for i in range(m):
                    for slope, kind in ((1, "root"), (prime, "factor")):
                        h = standard + n**(m-1-i)*f**i*slope*(t+center-root)
                        assert extract(n, center, y, m, Poly(h, t)) == (
                            kind, root if kind == "root" else prime)
                        counts[kind] += 1
    print(counts)
