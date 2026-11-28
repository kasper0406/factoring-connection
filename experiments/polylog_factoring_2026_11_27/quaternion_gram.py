"""Primitive four-dimensional embeddings of a binary Gram of odd square determinant."""

from math import isqrt

if __package__:
    from .four_square_sampler import _doubled_units, _left_gcd, _multiply
    from .ternary_sphere import primitive_ternary_sphere
else:
    from four_square_sampler import _doubled_units, _left_gcd, _multiply
    from ternary_sphere import primitive_ternary_sphere


def _lifts(odd_norm, product):
    """The eight integral pairs with bar(u)*v=product and Norm(u)=odd_norm."""
    generator = _left_gcd((2*odd_norm, 0, 0, 0), tuple(2*x for x in product))
    pairs = []
    for unit in _doubled_units():
        associate = _multiply(generator, unit)
        if any(c % 4 for c in associate):
            continue
        a, b, c, d = (value // 4 for value in associate)
        u = (a, -b, -c, -d)
        numerator = _multiply(u, product)
        if any(value % odd_norm for value in numerator):
            raise ArithmeticError("quaternion quotient must be integral")
        pairs.append((u, tuple(value // odd_norm for value in numerator)))
    if len(pairs) != 8:
        raise ArithmeticError("expected eight distinct integral associates")
    return pairs


def embed_four_quaternion(gram, rng, *, split=None):
    """Sample (u,v) with the supplied Gram (A,B,C), where u·v=B.

    Require A>0, AC-B²=n² for a positive odd n. The distribution is exactly
    uniform on primitive frames when gcd(A,B,C)=1. In general it is uniform
    on frames for which the imaginary coordinates of bar(u)*v are coprime.
    A split is an optional reusable solution of x²+y²=−1 modulo n.
    """
    a, b, c = gram
    determinant = a*c - b*b
    if a <= 0 or determinant <= 0:
        raise ValueError("expected a positive definite binary Gram matrix")
    n = isqrt(determinant)
    if n*n != determinant or n % 2 == 0:
        raise ValueError("the determinant must be a positive odd square")
    point = primitive_ternary_sphere(n, rng, split=split)
    # If a is even, c is odd since the valuation of n²+b² at 2 is at most 1.
    swapped = a % 2 == 0
    pairs = _lifts(c if swapped else a, (b, *point))
    u, v = pairs[rng.randrange(len(pairs))]
    return (v, u) if swapped else (u, v)
