"""Exact-domain, recovery, and work comparisons for the pooled norm search."""

from itertools import combinations
from math import gcd, isqrt

from sympy import jacobi_symbol, primerange

if __package__:
    from .check_large_conic_batch import direct
    from .large_conic_batch import NormHit, batch_conic_factor
    from .pooled_conic_batch import FormHit, _form_tiles, pooled_conic_factor, pooled_reduced_form_factor
else:
    from check_large_conic_batch import direct
    from large_conic_batch import NormHit, batch_conic_factor
    from pooled_conic_batch import FormHit, _form_tiles, pooled_conic_factor, pooled_reduced_form_factor


def check_tiles():
    for h in range(1, 65):
        for form in [(1, 0, 2), (3, -2, 6), (3, 2, 6), (2, 2, 3), (1, 0, 1)]:
            u, v, w = form
            tiles = _form_tiles(h, form)
            points = [(a,b) for (width, start), rows in tiles.items()
                      for b in rows for a in range(start, start+width)]
            bound = isqrt(2*h)
            expected = {(a,b) for a in range(-bound, bound+1) for b in range(1, bound+1)
                        if u*a*a+v*a*b+w*b*b <= h}
            assert len(points) == len(expected) and set(points) == expected


def check_cohort():
    eligible = mixed = 0
    for p, q in combinations(primerange(11, 85), 2):
        if q >= 2*p:
            continue
        n = p*q
        for x in range(10):
            a = x*x+1
            result = pooled_conic_factor(n, x)
            assert bool(result) == bool(direct(n, x))
            if result:
                d, r, s = result
                value = r*r+a*s*s
                assert 0 < value <= isqrt(8*n) and gcd(value, n) == d in (p, q)
            eligible += 1
            if jacobi_symbol(-a, n) == -1:
                form, (d, i, j) = pooled_reduced_form_factor(n, x)
                u, v, w = form
                assert d in (p, q) and gcd(n, u*i*i+v*i*j+w*j*j) == d
                mixed += 1
    return eligible, mixed


def check_examples():
    n, x = 19*23, 4
    assert pooled_conic_factor(n, x) is None
    assert pooled_reduced_form_factor(n, x) == FormHit((3, -2, 6), NormHit(23, 1, 2))
    stats = {}
    assert pooled_conic_factor(91, 4, stats=stats) == NormHit(7, 2, 1)
    assert stats['saturated_blocks'] == 1
    n, x = 10412720415486989, 8192
    assert pooled_conic_factor(n, x, stats=stats) == NormHit(83317541, 4026, 1)
    print(f'54-bit planted case: {stats["block_values"]} block values, '
          f'{stats["scalar_gcds"]} scalar gcds, maximum polynomial degree {stats["max_degree"]}')
    # Both primes are inert in Q(sqrt(-2)), so the full fixed-form search is a miss.
    n, x = 131101*163901, 1
    old, new = {}, {}
    assert pooled_conic_factor(n, x, stats=new) is None
    assert batch_conic_factor(n, x, stats=old) is None
    assert new['block_values'] < old['block_values'] and new['scalar_gcds'] == 0
    print(f'Complete 35-bit miss: old {old["block_values"]} block values and {old["scalar_gcds"]} '
          f'scalar gcds; pooled {new["block_values"]} and {new["scalar_gcds"]}; '
          f'polynomial degree {old["max_degree"]} and {new["max_degree"]}, respectively')


if __name__ == '__main__':
    check_tiles()
    eligible, mixed = check_cohort()
    print(f'{eligible} fixed-form scans agree with direct enumeration; all {mixed} mixed inputs split')
    check_examples()
