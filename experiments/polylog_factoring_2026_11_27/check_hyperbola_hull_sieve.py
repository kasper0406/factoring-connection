"""Independent finite checks for the clipped hyperbola hull and its divisor stream."""

from itertools import product
from math import gcd, isqrt
from random import Random

if __package__:
    from .hyperbola_hull_sieve import corner_chain, divisor_events, initial_hull
else:
    from hyperbola_hull_sieve import corner_chain, divisor_events, initial_hull


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def lower_chain(points):
    """Andrew's monotone chain for the given traversal order."""
    hull = []
    for p in points:
        while len(hull) >= 2:
            a, b = hull[-2:]
            if det((b[0]-a[0], b[1]-a[1]), (p[0]-a[0], p[1]-a[1])) > 0:
                break
            hull.pop()
        hull.append(p)
    return hull


def brute_hull(n, lo, hi):
    return lower_chain((x, (n+x-1)//x) for x in range(lo, hi+1))


def brute_corner(u, v):
    # The sail is the side facing the origin of the hull of the nonzero
    # lattice points in conv(0,u,v). Its other side is the segment [u,v].
    signed_area = det(u, v)
    sign = 1 if signed_area > 0 else -1
    xs, ys = zip((0, 0), u, v)
    pts = []
    for p in product(range(min(xs), max(xs)+1), range(min(ys), max(ys)+1)):
        a, b = sign*det(p, v), sign*det(u, p)
        if a >= 0 and b >= 0 and 0 < a+b <= abs(signed_area):
            pts.append(p)
    pts.sort()
    polygon = lower_chain(pts)[:-1] + lower_chain(reversed(pts))[:-1]
    at = polygon.index(u)
    chain = [u]
    while chain[-1] != v:
        at = (at-sign) % len(polygon)
        chain.append(polygon[at])
    return chain


def check_initial_hulls():
    rng = Random(2491)
    cases = []
    for n in range(1, 17):
        cases.extend((n, lo, hi) for lo in range(1, n+2) for hi in range(lo, n+3))
    for n in range(17, 131):
        cases.append((n, 1, n))
        for _ in range(5):
            lo = rng.randint(1, n)
            cases.append((n, lo, rng.randint(lo, n+3)))
    for n, lo, hi in cases:
        assert initial_hull(n, lo, hi) == brute_hull(n, lo, hi), (n, lo, hi)
    return len(cases)


def check_corners():
    vectors = [v for v in product(range(-3, 4), repeat=2) if gcd(*v) == 1]
    cases = 0
    for u, v in product(vectors, repeat=2):
        if det(u, v):
            assert corner_chain(u, v) == brute_corner(u, v), (u, v)
            cases += 1
    # A long collinear run in the sail's lattice points has just two extremes.
    huge = 2**333+11
    assert corner_chain((1, 0), (1, huge)) == [(1, 0), (1, huge)]
    # Unimodular shearing preserves the entire hull, including corner order.
    def shear(p):
        return p[0]+huge*p[1], p[1]
    u, v = (-3, 2), (1, 5)
    assert corner_chain(shear(u), shear(v)) == [shear(p) for p in brute_corner(u, v)]
    return cases+2


def brute_events(start, stop, lo, hi):
    return [(n, x, n//x) for n in range(start, stop+1)
            for x in range(lo, hi+1) if n % x == 0]


def check_events():
    rng = Random(31512)
    cases = [(1, 180, 1, None), (6, 6, 2, 3), (36, 37, 1, 36),
             (37, 39, 17, 18), (97, 101, 10, 11), (2, 18, 20, 23)]
    cases.extend((70, 120, x, x) for x in range(1, 46))
    for _ in range(250):
        start = rng.randint(1, 200)
        stop = start+rng.randrange(16)
        lo = rng.randint(1, 90)
        cases.append((start, stop, lo, rng.randint(lo, 110)))
    total = 0
    for start, stop, lo, hi in cases:
        stats = {}
        got = list(divisor_events(start, stop, lo, hi, stats=stats))
        assert got == brute_events(start, stop, lo, stop if hi is None else hi), (start, stop, lo, hi)
        assert stats['events'] == len(got)
        total += len(got)
    return len(cases), total


def check_large_window():
    # For |t|<=W, (X+t)(X-t)=X²-t², while (X+t)(X-t+1)>X²+1.
    x, w = 2**160, 2**78
    n = x*x
    lo, hi = x-w, x+w
    ends = [(lo, hi+1), (hi, lo+1)]
    assert initial_hull(n-1, lo, hi) == [ends[0], (x-1, x+1), (x+1, x-1), ends[1]]
    assert initial_hull(n, lo, hi) == [ends[0], (x, x), ends[1]]
    assert initial_hull(n+1, lo, hi) == ends
    got = list(divisor_events(n-1, n+1, lo, hi))
    assert got == [(n-1, x-1, x+1), (n-1, x+1, x-1), (n, x, x)]
    # Many consecutive levels are empty on this entire column band.
    assert not list(divisor_events(n+1, n+x//2, lo, hi))
    return n.bit_length(), hi-lo+1


def check_consecutive_batch():
    first, stop = 10000123, 10000634
    expected = []
    for n in range(first, stop+1):
        divisors = set()
        for d in range(1, isqrt(n)+1):
            if n % d == 0:
                divisors.update((d, n//d))
        expected.extend((n, d, n//d) for d in sorted(divisors))
    for lo, hi in ((1, stop), (2200, 4300)):
        stats = {}
        assert list(divisor_events(first, stop, lo, hi, stats=stats)) == [
            row for row in expected if lo <= row[1] <= hi]
        print(f"512 consecutive products, columns {lo}..{hi}, {stats['initial_vertices']} "
              f"initial vertices, {stats['events']} events")


if __name__ == '__main__':
    hulls, corners = check_initial_hulls(), check_corners()
    windows, events = check_events()
    bits, columns = check_large_window()
    print(f"{hulls} clipped hulls, {corners} cone checks, {events} events in {windows} windows")
    print(f"{bits}-bit square, {columns} columns, three events and a certified empty level interval")
    check_consecutive_batch()
