"""Exact finite controls for the affine prime-window reader. Run directly or as a module."""

from itertools import product
from math import gcd, isqrt, prod
from random import Random

from sympy import Poly, Symbol

if __package__:
    from .affine_prime_window import (
        _integer_roots, affine_divisors, prime_window_divisor, window_idempotent, window_radius,
    )
else:
    from affine_prime_window import (
        _integer_roots, affine_divisors, prime_window_divisor, window_idempotent, window_radius,
    )


def compare(p, q, x, radius):
    n = p*q
    selected = [r for r in (p, q) if min(x % r, (-x) % r) <= radius]
    stats = {}
    value = prime_window_divisor(n, x, radius=radius, stats=stats)
    assert value == prod(selected), (p, q, x, radius, value, selected)
    e = window_idempotent(n, value)
    assert e*e % n == e and gcd(e-1, n) == value and gcd(e, n) == n//value
    assert all(e % r == int(r in selected) for r in (p, q))
    return value, stats


def exhaustive_checks():
    count = 0
    for p, q in ((2, 3), (5, 7), (67, 71), (73, 83)):
        n = p*q
        for radius in range(window_radius(n)+1):
            hits = set()
            for x in range(n):
                hits.add(compare(p, q, x, radius)[0])
                count += 1
            assert hits == {1, p, q, n}
    return count


def lattice_checks():
    p, q = 1009, 1321
    n = p*q
    radius = window_radius(n)
    assert radius == 4
    count = 0
    # Distinct hits, the same hit, signed bounds, and just outside both boundaries.
    for u, v in product(range(-radius-1, radius+2), repeat=2):
        x = u+p*((v-u)*pow(p, -1, q) % q)
        _, stats = compare(p, q, x-n, radius)
        assert stats['method'] == ('central' if u == v and abs(u) <= radius else 'lattice')
        count += 1
    for radius in (0, 1, 2, 3):
        compare(p, q, -3*n+p+radius, radius)
        count += 1
    return count


def arbitrary_divisor_checks():
    # These moduli have repeated factors. Trial enumeration is an independent oracle.
    rng = Random(413)
    count = 0
    for n in (2**12*3**4*5*7, 2**9*3**6*11, 10007*12401):
        lower = isqrt(n)
        radius = window_radius(n)
        for center in [n-3, 3, 0, 2**12+radius, 10007+radius, *(rng.randrange(-n,n) for _ in range(12))]:
            direct = {gcd(n, center+y) for y in range(-radius, radius+1)}
            expected = tuple(sorted(g for g in direct if g >= lower and g > 1))
            assert affine_divisors(n, center, radius, lower) == expected
            count += 1
    # At the inclusive lower bound the fixture's composite gcd must be returned.
    n, lower, radius = 2**12*3**4*5*7, 2**12, 5
    assert affine_divisors(n, lower+radius, radius, lower) == (lower,)
    assert affine_divisors(n, lower+radius, radius, lower+1) == ()
    return count+2


def root_checks():
    t = Symbol('t')
    rng = Random(229)
    radius = 7
    count = 0
    for _ in range(80):
        poly = Poly.from_list([1]+[rng.randrange(-9,10) for _ in range(5)], t)
        # Repeated roots and zeros exactly on a half-integer bisection boundary.
        poly *= Poly((t-rng.randrange(-8,9))**2*(2*t-rng.choice([-3,1,3,5])), t)
        assert sorted(_integer_roots(poly, radius)) == [x for x in range(-radius, radius+1) if poly.eval(x) == 0]
        count += 1
    radius = 1 << 100
    poly = Poly((t-radius)*(t+radius)**2*(2*t-1), t)
    assert sorted(_integer_roots(poly, radius)) == [-radius, radius]
    return count+1


def large_windows():
    p, q = 16777213, 20000003
    n = p*q
    radius = window_radius(n)
    assert radius == 534
    both = radius+p*((-2*radius*pow(p, -1, q)) % q)
    centers = ((2**24+512, p), (isqrt(n), 1), (both, n), (-radius, n))
    records = []
    for center, expected in centers:
        result, stats = compare(p, q, center, radius)
        assert result == expected
        if min(center % n, (-center) % n) <= radius:
            assert stats['method'] == 'central'
        else:
            assert stats['method'] == 'lattice'
            assert stats['offsets_tested'] < stats['dimension'] < 2*radius+1
        records.append((center, result, stats['dimension'], stats['offsets_tested']))
    return n, radius, records


if __name__ == '__main__':
    print('Exhaustive residue/window pairs:', exhaustive_checks())
    print('Planted local boundary pairs and narrower windows:', lattice_checks())
    print('Arbitrary-composite affine gcd controls:', arbitrary_divisor_checks())
    print('Integer-root isolation controls:', root_checks())
    print('Large modulus, radius, (center, divisor, dimension, gcd offsets):', large_windows())
