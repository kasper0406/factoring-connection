"""Sample primitive or unrestricted points on the ternary sphere of norm (x²+1)T²."""

from itertools import permutations, product
from math import gcd

if __package__:
    from .fixed_field_packing import fixed_field_parameters, packing_params
    from .four_square_incidence import sample_content
    from .four_square_sampler import _multiply, norm_minus_one, primitive_four_squares
    from .short_jacobi_grams import first_vector_orbits
else:
    from fixed_field_packing import fixed_field_parameters, packing_params
    from four_square_incidence import sample_content
    from four_square_sampler import _multiply, norm_minus_one, primitive_four_squares
    from short_jacobi_grams import first_vector_orbits


def base_points(x, *, primitive=True):
    """List signed, ordered points of squared norm x²+1, optionally including imprimitive ones."""
    points = set()
    for shape in first_vector_orbits(x):
        if not primitive or gcd(*shape) == 1:
            for order in set(permutations(shape)):
                for signs in product((-1, 1), repeat=3):
                    points.add(tuple(s*a for s, a in zip(signs, order)))
    return sorted(points)


def primitive_base_points(x):
    return base_points(x, primitive=True)


def conjugate_vector(q, vector):
    """Return the three imaginary coordinates of q vector conjugate(q)."""
    a, b, c, d = q
    return _multiply(_multiply(q, (0, *vector)), (a, -b, -c, -d))[1:]


def primitive_scaled_sphere(x, conductor, rng, *, split=None):
    """Sample exactly uniformly for any x≥0, odd T=conductor≥1; no factors needed.

    The work to enumerate the base shell depends on x. A supplied split is
    any pair of residues with squared norm −1 modulo T.
    """
    if conductor < 1 or conductor % 2 == 0:
        raise ValueError("expected a positive odd conductor")
    base = primitive_base_points(x)
    if split is None:
        split = norm_minus_one(conductor, rng)
    while True:
        s = base[rng.randrange(len(base))]
        q = primitive_four_squares(conductor, rng, split=split)
        w = conjugate_vector(q, s)
        if gcd(*w) == 1:
            return w


def scaled_sphere(x, conductor, rng, *, split=None):
    """Sample all points uniformly for an odd conductor coprime to x²+1."""
    if conductor < 1 or conductor % 2 == 0 or gcd(conductor, x*x+1) != 1:
        raise ValueError("expected a positive odd conductor coprime to the base norm")
    base = base_points(x, primitive=False)
    if split is None:
        split = norm_minus_one(conductor, rng)
    while True:
        s = base[rng.randrange(len(base))]
        g = sample_content(conductor, 1, rng)
        m = conductor//g
        q = primitive_four_squares(m, rng, split=split)
        w = conjugate_vector(q, s)
        if gcd(*w, m) == 1:
            return tuple(g*a for a in w)


def packed_sphere_point(n, x, rng, *, gaussian=None, primitive=False):
    """Sample the whole sphere or its primitive points for a public packing."""
    if gaussian is None:
        params = fixed_field_parameters(n, x)
        t = 1
    else:
        b, t = gaussian
        params = packing_params(n, x, b, t)
    root = params.B * pow(t*n, -1, params.T) % params.T
    sampler = primitive_scaled_sphere if primitive else scaled_sphere
    return sampler(x, params.T, rng, split=(root, 0))


def packed_primitive_point(n, x, rng, *, gaussian=None):
    return packed_sphere_point(n, x, rng, gaussian=gaussian, primitive=True)
