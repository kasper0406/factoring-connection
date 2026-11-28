"""Uniform primitive integer points on the sphere of a positive odd radius."""

from math import gcd

if __package__:
    from .four_square_sampler import norm_minus_one, primitive_four_squares
else:
    from four_square_sampler import norm_minus_one, primitive_four_squares


def hopf_point(quaternion):
    """Map (a+bi,c+di) to a point whose radius is a²+b²+c²+d²."""
    a, b, c, d = quaternion
    return a*a + b*b - c*c - d*d, 2*(a*c + b*d), 2*(b*c - a*d)


def primitive_ternary_sphere(n, rng, *, split=None):
    """Sample every signed, ordered primitive solution x²+y²+z²=n² uniformly.

    The optional split is a reusable solution of a²+b² = −1 modulo n.
    The algorithm needs only the odd radius n and uniform random integers.
    """
    if n < 1 or n % 2 == 0:
        raise ValueError("the radius must be positive and odd")
    if split is None:
        split = norm_minus_one(n, rng)
    while True:
        point = hopf_point(primitive_four_squares(n, rng, split=split))
        if gcd(*point) == 1:
            shift = rng.randrange(3)
            return point[shift:] + point[:shift]


def perpendicular_form(point):
    """Return the integral binary norm form on the full perpendicular lattice."""
    x, y, z = point
    if gcd(x, y, z) != 1:
        raise ValueError("the input point must be primitive")
    g = gcd(x, y)
    if g == 0:
        return 1, 0, 1
    if y:
        s = pow(x // g, -1, abs(y) // g)
        t = (g - s*x) // y
    else:
        s, t = x // g, 0
    # The cross product of these two vectors is exactly (x,y,z).
    e = (-y // g, x // g, 0)
    f = (-s*z, -t*z, g)
    return sum(a*a for a in e), 2*sum(a*b for a, b in zip(e, f)), sum(b*b for b in f)
