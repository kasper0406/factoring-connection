"""Construct balanced integer quotient pairs satisfying the specified low product bits."""

from math import gcd, isqrt
from typing import NamedTuple


class QuotientPair(NamedTuple):
    u: int
    v: int
    line: int
    remainder: int


def inverse_line(n, line, period):
    """Unique odd x,y with xy=n and 3x+2y=line modulo a power of two."""
    if n % 2 == 0 or line % 2 == 0 or period < 2 or period & (period-1):
        raise ValueError("expected an odd input and line, and a power-of-two period")
    x, precision = 1, 2
    while precision < period:
        precision = min(precision*precision, period)
        value = 3*x*x-line*x+2*n
        x = (x-value*pow(6*x-line, -1, precision)) % precision
    return x, n*pow(x, -1, period) % period


def _lift(n, line, x, y, period):
    # u=x (mod period) also needs (line-3u)/2=y (mod period).
    residue = x + period*(((line-3*x-2*y)//period) % 2)
    disc = line*line-24*n
    ceiling_root = 0 if disc == 0 else 1+isqrt(disc-1)
    below = (line-ceiling_root)//6
    u = below-(below-residue) % (2*period)
    if u < 1:
        return None
    v = (line-3*u)//2
    remainder = n-u*v
    if not (u < v < 2*u and 0 <= remainder < u):
        return None
    return QuotientPair(u, v, line, remainder)


def quotient_witnesses(n, mask, budget, *, wheel=1):
    """Find both predicate outcomes in the floor-quotient relaxation.

    The predicate is (u & mask != 0 and v & mask != 0). A witness has odd
    u<v<2u, v=floor(n/u), and n-u*v divisible by both the bit period and
    the optional odd wheel. The method tests at most budget tangent lines.
    """
    if n < 3 or n % 2 == 0 or mask < 2 or mask % 2 or budget < 1:
        raise ValueError("expected an odd input, positive non-units bit mask and positive budget")
    if wheel < 1 or wheel % 2 == 0 or gcd(n, wheel) != 1:
        raise ValueError("the extra wheel must be odd, positive and coprime to the input")
    period = 1 << mask.bit_length()
    line = 1+isqrt(24*n-1)
    # u=1, v=n mod wheel also presieves both candidate factors.
    line += (3+2*n-line) % wheel
    if line % 2 == 0:
        line += wheel
    inv = pow(period, -1, wheel) if wheel > 1 else 0
    found = {}
    for step in range(budget):
        x, y = inverse_line(n, line, period)
        outcome = bool(x & mask and y & mask)
        if outcome not in found:
            x += period*((1-x)*inv % wheel)
            y += period*((n-y)*inv % wheel)
            point = _lift(n, line, x, y, period*wheel)
            if point is not None:
                found[outcome] = point
                if len(found) == 2:
                    return found, step+1
        line += 2*wheel
    return found, budget
