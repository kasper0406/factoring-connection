"""Check signed primitive sphere counts, exact Hopf fibers, and perpendicular forms."""

from collections import defaultdict
from math import gcd, isqrt, prod
from random import Random

if __package__:
    from .check_four_square_incidence import prime_powers, quadruples
    from .ternary_sphere import hopf_point, perpendicular_form, primitive_ternary_sphere
else:
    from check_four_square_incidence import prime_powers, quadruples
    from ternary_sphere import hopf_point, perpendicular_form, primitive_ternary_sphere


def primitive_points(n):
    for x in range(-n, n + 1):
        y_bound = isqrt(n*n - x*x)
        for y in range(-y_bound, y_bound + 1):
            rest = n*n - x*x - y*y
            z = isqrt(rest)
            if z*z == rest and gcd(x, y, z) == 1:
                yield x, y, z
                if z:
                    yield x, y, -z


def gaussian_gcd_norm(a, b, c, d):
    while c or d:
        norm = c*c + d*d
        re = (2*(a*c + b*d) + norm) // (2*norm)
        im = (2*(b*c - a*d) + norm) // (2*norm)
        a, b, c, d = c, d, a - re*c + im*d, b - im*c - re*d
    return a*a + b*b


def check_radius(n):
    factors = prime_powers(n)
    points = set(primitive_points(n))
    counted = 6*prod(r**(e-1)*(r - (1 if r % 4 == 1 else -1)) for r, e in factors)
    assert len(points) == counted
    fibers = defaultdict(set)
    examined = primitive_quaternions = 0
    for row in quadruples(n):
        examined += 1
        point = hopf_point(row)
        content = gcd(*point)
        assert sum(v*v for v in point) == n*n and point[0] % 2 == 1
        assert content == gaussian_gcd_norm(*row)
        if gcd(*row) == 1:
            primitive_quaternions += 1
        if content == 1:
            assert gcd(*row) == 1
            fibers[point].add(row)
    assert primitive_quaternions == 8*prod(r**(e-1)*(r+1) for r, e in factors)
    assert set(fibers) == {point for point in points if point[0] % 2}
    for lifts in fibers.values():
        a, b, c, d = next(iter(lifts))
        assert lifts == {(a, b, c, d), (-b, a, -d, c),
                         (-a, -b, -c, -d), (b, -a, d, -c)}
    for point in points:
        a, middle, c = perpendicular_form(point)
        assert a > 0 and c > 0 and gcd(a, middle, c) == 1
        assert middle*middle - 4*a*c == -4*n*n
    return examined, len(fibers), len(points)


if __name__ == "__main__":
    totals = [0, 0, 0]
    for radius in range(1, 100, 2):
        for i, count in enumerate(check_radius(radius)):
            totals[i] += count
    rng = Random(28022026)
    radii = (1, 3, 5, 25, 45, 65, 3**4*5**2, 1009*1013)
    for radius in radii:
        for _ in range(8):
            point = primitive_ternary_sphere(radius, rng)
            assert gcd(*point) == 1 and sum(t*t for t in point) == radius*radius
    print("50 odd radii through 99:", totals[0], "quaternion lifts,",
          totals[1], "four-element fibers,", totals[2], "primitive points and forms.")
    print("Public sampler produced verified points on", 8*len(radii), "seeded calls.")
