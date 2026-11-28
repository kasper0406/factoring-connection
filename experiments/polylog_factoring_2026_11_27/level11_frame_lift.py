"""Decode a supplied intermediate level-11 integral frame on a balanced semiprime."""

from itertools import product
from math import gcd, isqrt

from sympy import Matrix

if __package__:
    from .level11_sampler import _inverse_splitting, norm
else:
    from level11_sampler import _inverse_splitting, norm


REFERENCE_NORM = {"principal": 1, "nonprincipal": 3}


def frame_lift_exponent(n):
    """Least k >= 2 with (2^k)^4 >= n."""
    if n < 1:
        raise ValueError("expected a positive integer")
    return max(2, ((n - 1).bit_length() + 3) // 4)


def _coordinate_lifts(residue, modulus, bound):
    first = -bound + ((residue + bound) % modulus)
    return range(first, bound + 1, modulus)


def factor_from_frame(n, k, frame, *, base="principal"):
    """Return the smaller factor, or None if no lift yields a proper divisor.

    The guarantee assumes n=p*q, p<q<=2p, and frame comes from a relative
    norm-p or norm-q intermediate ideal in framed_line_classifier's basis.
    If the level 2^k is at least n^(1/4), at most 1,600 candidate
    quadruples are tested. Lower levels use the same finite search.
    """
    if n < 2 or gcd(n, 22) != 1:
        raise ValueError("expected a modulus greater than one and coprime to 22")
    if not isinstance(k, int) or k < 2:
        raise ValueError("expected a frame exponent of at least two")
    target, entries = frame
    if target not in REFERENCE_NORM or base not in REFERENCE_NORM:
        raise ValueError("unknown ideal class")
    if len(entries) != 4:
        raise ValueError("expected four row-major matrix entries")
    a, b, c, d = entries
    if (a*d - b*c) % 2 == 0:
        raise ValueError("the integral frame must be invertible modulo 2")

    modulus = 1 << k
    target_norm = REFERENCE_NORM[target]
    norm_scale = REFERENCE_NORM[base] * target_norm
    # At a pure power of two the splitting is deterministic.
    inverse = _inverse_splitting(modulus, None)
    beta = tuple(int(x) * target_norm % modulus for x in inverse * Matrix(entries))
    max_norm = norm_scale * isqrt(2*n)
    x_bound = isqrt(12 * max_norm // 11)
    y_bound = isqrt(4 * max_norm // 11)
    bounds = (x_bound, y_bound, x_bound, y_bound)
    lifts = [_coordinate_lifts(x, modulus, bound) for x, bound in zip(beta, bounds)]
    for candidate in product(*lifts):
        divisor, remainder = divmod(norm(candidate), norm_scale)
        if not remainder and 1 < divisor < n and n % divisor == 0:
            return min(divisor, n // divisor)
    return None
