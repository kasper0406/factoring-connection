"""Common Jacobi tests on a candidate smaller prime and its floor quotient."""

from math import isqrt
from random import SystemRandom

if __package__:
    from .public_factor_jacobi_sieve import jacobi
else:
    from public_factor_jacobi_sieve import jacobi


def paired_moduli(n, failure_bits=20, *, rng=None):
    """Common denominators with error at most 2^-failure_bits on balanced prime pairs."""
    if n < 15 or n % 2 == 0 or failure_bits < 0:
        raise ValueError("expected an odd balanced semiprime and a nonnegative error exponent")
    rounds = (isqrt(n).bit_length()+failure_bits+3)//2
    height = n*(rounds+1)
    source = SystemRandom() if rng is None else rng
    return [1+4*n*source.randrange(height) for _ in range(rounds)]


def quotient_candidates(n, lower=None, upper=None):
    low = isqrt((n-1)//2)+1
    root = isqrt(n)
    high = root
    if lower is not None:
        low = max(low, lower)
    if upper is not None:
        high = min(high, upper)
    square_below = root*root
    for x in range(low | 1, high+1, 2):
        y = n//x
        if y % 2 and isqrt(x)**2 != x and isqrt(y)**2 != y and x*y != square_below:
            yield x, y


def paired_scan(n, moduli):
    """Explicit candidate scan, returning survivors and their count after every row."""
    pairs = list(quotient_candidates(n))
    counts = [len(pairs)]
    for d in moduli:
        pairs = [(x, y) for x, y in pairs if jacobi(d, x) == 1 and jacobi(d, y) == 1]
        counts.append(len(pairs))
    return pairs, counts


def composite_single_coordinate_witness(lower, upper, moduli):
    """Find a product of distinct sign-matched seeds, or None in this seed interval."""
    if lower < 4 or upper < lower or any(d < 1 or d % 4 != 1 for d in moduli):
        raise ValueError("expected a positive interval and positive denominators congruent to 1 mod 4")
    earlier = {}
    for u in range((isqrt(lower-1)+1) | 1, isqrt(upper)+1, 2):
        signs = tuple(jacobi(u, d) for d in moduli)
        if 0 in signs:
            continue
        if signs in earlier:
            v = earlier[signs]
            if isqrt(u*v)**2 != u*v:
                return u*v, v, u
        else:
            earlier[signs] = u
    return None
