"""Probe for a proper denominator in the left order of a sampled level-11 ideal."""

from math import gcd
from random import SystemRandom

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

if __package__:
    from . import level11_sampler as sampler
else:
    import level11_sampler as sampler


def _field_product(a, b):
    x, y = a
    z, w = b
    return x*z - 3*y*w, x*w + y*z + y*w


def _bar(a):
    x, y = a
    return x+y, -y


def _product(left, right):
    """Product of a+b*j and c+d*j in the order basis (1, omega, j, omega*j)."""
    a, b = left[:2], left[2:]
    c, d = right[:2], right[2:]
    ac = _field_product(a, c)
    bd = _field_product(b, _bar(d))
    ad = _field_product(a, d)
    bc = _field_product(b, _bar(c))
    return ac[0]-bd[0], ac[1]-bd[1], ad[0]+bc[0], ad[1]+bc[1]


def _conjugate(point):
    x, y, z, w = point
    return x+y, -y, -z, -w


def _unit_numerators(n, inverse, line):
    """Enumerate n times every unit in the left order of the given line ideal."""
    basis = sampler._line_basis(n, inverse, line)
    columns = [tuple(map(int, basis.col(j))) for j in range(4)]
    products = [_product(a, _conjugate(b)) for a in columns for b in columns]
    # I*bar(I) = n*O_left(I), since I is locally principal and has norm n.
    scaled_order = hermite_normal_form(Matrix(products).T)
    units = sampler._shell_points(scaled_order, n*n)
    if len(units) not in (4, 6):
        raise ArithmeticError("unexpected unit count in a maximal order of discriminant 11")
    return units


def _integrality_divisor(n, numerators):
    # All noncentral units in either class agree up to sign and scalar addition.
    unit = next(x for x in numerators if any(x[1:]))
    return gcd(n, *unit)


def automorphism_probe(n, rng=None):
    """One factor-blind trial. Return a verified proper divisor, or None."""
    if n <= 1:
        raise ValueError("the input must exceed one")
    if n % 11 == 0:
        return 11 if n > 11 else None
    if rng is None:
        rng = SystemRandom()
    inverse = sampler._inverse_splitting(n, rng)
    line = sampler._primitive_pair(n, rng)
    divisor = _integrality_divisor(n, _unit_numerators(n, inverse, line))
    return divisor if 1 < divisor < n else None
