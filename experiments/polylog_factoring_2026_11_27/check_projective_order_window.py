"""Exhaust the small boxes and recover large power-image cofactors from public inputs."""

from math import gcd

from sympy import Poly, Symbol

if __package__:
    from .projective_order_window import _rational_roots, projective_order_window, small_power_image_factor
else:
    from projective_order_window import _rational_roots, projective_order_window, small_power_image_factor


def brute(n, e, h_bound, k_bound, lower, offset):
    factors = set()
    for h in range(1, h_bound+1):
        for k in range(1, k_bound+1):
            if gcd(h, k) != 1:
                continue
            common = gcd(k-offset*e*h, n)
            if lower <= common < n:
                factors.add(common)
            quotient, tail = divmod(e*h, k)
            candidate = quotient-offset
            if not tail and lower <= candidate < n and n % candidate == 0:
                factors.add(candidate)
    return tuple(sorted(factors))


def check():
    t = Symbol('T')
    poly = Poly(t*(5*t-3)**2*(17*t-10)*(2*t-3)*(t-5)*(t*t-2)*(t+1), t)
    assert _rational_roots(poly, 17, 10) == [(1, 5), (2, 3), (5, 3), (17, 10)]
    assert _rational_roots(poly, 5, 3) == [(2, 3), (5, 3)]

    p, q = 1000000007, 1000000009
    n = p*q
    small_boxes = 0
    for h_bound, k_bound in ((1, 1), (4, 4), (3, 8), (8, 3), (9, 13)):
        for offset in (-1, 1):
            # Place two different primitive slopes on the two components.
            left = offset * k_bound * pow(h_bound, -1, p) % p
            right = offset * (k_bound+1) * pow(h_bound+1, -1, q) % q
            both = left + ((right-left)*pow(p, -1, q) % q)*p
            for e in (1, 53, n+1, max(1, both), offset*(p-1) % n, 19*n+73):
                expected = brute(n, e, h_bound, k_bound, p, offset)
                stats = {}
                actual = projective_order_window(n, e, h_bound, k_bound, lower=p,
                                                 order_offset=offset, stats=stats)
                assert actual == expected
                if h_bound*k_bound > 16:
                    assert stats['method'] == 'lattice'
                small_boxes += 1
    left, right = -11*pow(5, -1, p) % p, -7*pow(4, -1, q) % q
    two_slopes = left + ((right-left)*pow(p, -1, q) % q)*p
    assert projective_order_window(n, two_slopes, 5, 11) == (p, q)

    # These factors construct the challenge exponent and check the answer.
    # The public-envelope call in each case receives only (n,e) and the sign.
    p127 = (1 << 127)-1
    q128 = 204169420152563078078024764459060926911
    p521 = (1 << 521)-1
    q522 = int('8237757192156731657978280958897671860723322360171966491273356151022651820077187'
               '262547071568793745465972755573669777029644546385599659972575088833949338068639')
    fixtures = ((p127, q128, (1 << 21)-1, (1 << 17)+1, -1),
                (p127, q128, 1 << 27, (1 << 15)+1, 1),
                (p521, q522, (1 << 52)-1, (1 << 47)+3, -1))
    reports = []
    for target, other, h, k, offset in fixtures:
        n = target*other
        assert target < other < 2*target and (target+offset) % h == 0 and gcd(h, k) == 1
        e = (target+offset)//h*k
        assert (target+offset)//gcd(target+offset, e) == h
        stats = {}
        assert projective_order_window(n, e, h, k, order_offset=offset, stats=stats) == (target,)
        assert stats['method'] == 'lattice' and stats['candidates'] < h
        stats = {}
        assert small_power_image_factor(n, e, order_offset=offset, stats=stats) == target
        assert stats['method'] == 'lattice' and stats['max_h'] >= h and stats['k_bound'] >= k
        reports.append((n.bit_length(), offset, h, k, stats['dimension'], stats['boxes']))

    # e > n^(19/25) and h*k > n^(1/4): this control requires the modulus n*e.
    n, h, k = p127*q128, 3, (1 << 69)+3
    e = (p127-1)//h*k
    assert (p127-1)//gcd(p127-1, e) == h and e**25 > n**19 and (h*k)**4 > n
    stats = {}
    assert projective_order_window(n, e, h, k, order_certificate=True, stats=stats) == (p127,)
    assert stats['certificate'] == 'n_times_e'
    assert small_power_image_factor(n, e, stats=stats) == p127
    assert stats['certificate'] == 'n_times_e'
    reports.append((n.bit_length(), -1, h, k, stats['dimension'], stats['boxes']))

    # Generic gcd n can hide the factor on an exact order witness.
    assert projective_order_window(35, 30, 1, 5, lower=5) == (7,)
    # A fixed enormous exponent has an empty public envelope.
    assert small_power_image_factor(n, n*n) is None
    assert small_power_image_factor(3, 1) is None
    return small_boxes, reports


if __name__ == '__main__':
    boxes, reports = check()
    print(f'Matched exhaustive factors in {boxes} boxes and recovered both prime components in one box.')
    for bits, offset, h, k, dimension, boxes_used in reports:
        print(f'{bits} bits, order r{offset:+d}, h={h}, k={k}: public lattice dimension {dimension}, '
              f'{boxes_used} boxes')
