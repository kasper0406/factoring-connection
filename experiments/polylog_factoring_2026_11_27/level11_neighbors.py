"""Find a level-11 norm generator from a line and its fixed norm-two translate."""

from math import gcd

from sympy import Matrix

if __package__:
    from .level11_sampler import _line_points
else:
    from level11_sampler import _line_points


def norm_two_translate(n, inverse, line):
    """Apply 1+j to a projective line in the supplied matrix splitting."""
    # inverse maps row-major matrix entries to coordinates on 1, omega, j, omega*j.
    coordinates = Matrix((1, 0, 1, 0))
    forward = inverse.adjugate() * pow(int(inverse.det()), -1, n)
    action = (forward * coordinates).reshape(2, 2)
    return tuple(int(z) % n for z in action * Matrix(line))


def repair_line(n, inverse, line):
    """Return (generator, used_twist) for L or (1+j)L, with at most two shell tests.

    Supply an odd positive n coprime to 11, its inverse level-11 matrix
    splitting, and any primitive pair representing L. The generator choice
    is deterministic once these inputs are fixed.
    """
    if n < 1 or gcd(n, 22) != 1:
        raise ValueError("the norm must be positive and coprime to 22")
    u, v = line
    if gcd(u, v, n) != 1:
        raise ValueError("the input coordinates must generate the unit ideal modulo n")
    line = (u % n, v % n)
    points = _line_points(n, inverse, line)
    if points:
        return points[0], False
    points = _line_points(n, inverse, norm_two_translate(n, inverse, line))
    if not points:
        raise ArithmeticError("the two neighboring level-11 ideals cannot both be nonprincipal")
    return points[0], True
