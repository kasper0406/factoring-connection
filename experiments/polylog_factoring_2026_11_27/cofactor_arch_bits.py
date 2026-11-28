"""Bit-mask zero tests on balanced Archimedean factor completions."""

from fractions import Fraction
from math import isqrt


def _first_avoiding(start, mask):
    """Least integer >= start whose word has no selected bit."""
    value = start
    while value & mask:
        period = 1 << (value & mask).bit_length()
        value = (value // period + 1)*period
    return value


def arch_outcomes(n, mask, *, lower=None, upper=None, stats=None):
    """All answers on u*v=n, sqrt(n/2)<=u<=sqrt(n), using floor(u), floor(v).

    An answer is true when each of the two words contains a selected bit.
    Optional nonnegative integer endpoints lower and upper restrict u further.
    Floor labels on real coordinates need not be primes or multiply to n.
    """
    if n < 1 or mask < 0:
        raise ValueError("expected a positive integer and a nonnegative mask")
    if (lower is not None and lower < 0) or (upper is not None and upper < 0):
        raise ValueError("expected nonnegative interval endpoints")
    low, high = Fraction(n, 2), Fraction(n)
    if lower is not None:
        low = max(low, Fraction(lower**2))
    if upper is not None:
        high = min(high, Fraction(upper**2))
    if stats is not None:
        stats.update(cut_points=0, full_period=False)
    if low > high:
        return frozenset()

    def floor_root(z):
        return isqrt(z.numerator // z.denominator)

    left = floor_root(low), floor_root(high)
    right = floor_root(n*n/high), floor_root(n*n/low)
    effective = mask & ((1 << right[1].bit_length())-1)
    if not effective:
        return frozenset({False})
    if all(_first_avoiding(first, effective) > last for first, last in (left, right)):
        return frozenset({True})

    # A full bit-one interval for u has length P/2. On its interior v=n/u
    # varies by strictly more than P/2, so v also attains a one label.
    # Check sqrt(high)-sqrt(low) >= 2P using rational comparisons only.
    twice_period = 4*(effective & -effective)
    residual = high-low-twice_period**2
    if residual >= 0 and residual**2 >= 4*twice_period**2*low:
        if stats is not None:
            stats["full_period"] = True
        return frozenset({False, True})

    # z=u² makes endpoints and bit-change boundaries rational. Generate edges
    # in the two coordinate intervals separately: on a restricted interval
    # there can be a large gap between the possible left and right factors.
    cuts = {low, high}
    bits = effective
    while bits:
        weight = bits & -bits
        for side, (first, last) in enumerate((left, right)):
            start = max(1, (first+weight-1)//weight)*weight
            for edge in range(start, last+1, weight):
                boundary = Fraction(edge**2) if side == 0 else Fraction(n*n, edge**2)
                if low <= boundary <= high:
                    cuts.add(boundary)
        bits -= weight
    boundaries = sorted(cuts)
    if stats is not None:
        stats["cut_points"] = len(boundaries)
    def occupied(z):
        return bool((floor_root(z) & effective) and (floor_root(n*n/z) & effective))

    # A coincident boundary can give an isolated yes because v is decreasing.
    yes = any(occupied(z) for z in boundaries) or any(
        occupied((a+b)/2) for a, b in zip(boundaries, boundaries[1:]))
    return frozenset({False, True}) if yes else frozenset({False})


def floor_lift(n, modulus, a, b):
    """Lift a*b=n mod modulus to balanced x,y with y=floor(n/x), or return None.

    The lift is guaranteed when n >= 4096*modulus**6. A returned pair
    satisfies its defining constraints at any height.
    """
    if n < 1 or modulus < 1 or (a*b-n) % modulus:
        raise ValueError("expected a positive input and a compatible residue pair")
    # Select the first congruent sum h above 2*sqrt(n).
    h = isqrt(4*n-1)+1
    h += (a+b-h) % modulus
    discriminant = h*h-4*n
    root = isqrt(discriminant)
    if root*root < discriminant:
        root += 1
    upper = (h-root)//2
    x = upper-(upper-a) % modulus
    y = h-x
    if x > 0 and y <= 2*x and n//x == y:
        return x, y
    return None
