"""Batch a²+(x²+1)b² up to floor(sqrt(8n)) over Z/nZ, with no factors supplied."""

from math import gcd, isqrt
from typing import NamedTuple

if __package__:
    from .sqrt_prime_marks import _DecimalRing
else:
    from sqrt_prime_marks import _DecimalRing

from src.trackC.polyback import multipoint_eval


class NormHit(NamedTuple):
    factor: int
    a: int
    b: int


def _quadratic_product(ring, constant, width):
    level = [ring.make([constant + j*j, 2*j, 1]) for j in range(width)]
    while len(level) > 1:
        next_level = [ring.mul(level[j], level[j + 1]) for j in range(0, len(level) - 1, 2)]
        if len(level) % 2:
            next_level.append(level[-1])
        level = next_level
    return level[0]


def batch_conic_factor(n, x, *, stats=None):
    """Return the first NormHit in (b,a) order, or None after the full exact scan.

    If n=pq with p<q<2p, a representation a²+(x²+1)b² in {p,2p,q,2q}
    guarantees a split. The scan can also find a multiple 3p.
    """
    if n < 9 or n % 2 == 0 or x < 1:
        raise ValueError("expected an odd n >= 9 and a positive x")
    limit = isqrt(8*n)
    norm = x*x + 1
    ring = _DecimalRing(n)
    work = {} if stats is None else stats
    work.update(rows=0, block_values=0, scalar_gcds=0, saturated_blocks=0, max_degree=0)
    for b in range(1, isqrt(limit // norm) + 1):
        constant = norm*b*b
        length = 1 + isqrt(limit - constant)
        width = isqrt(length)
        blocks = length // width
        starts = range(0, blocks*width, width)
        polynomial = _quadratic_product(ring, constant, width)
        values = multipoint_eval(ring, polynomial, starts)
        work["rows"] += 1
        work["block_values"] += blocks
        work["max_degree"] = max(work["max_degree"], 2*width)
        # The tail is shorter than a block. Inspect it after all full blocks.
        start = blocks * width
        stop = length
        for base, value in zip(starts, values):
            divisor = gcd(value, n)
            if divisor != 1:
                work["saturated_blocks"] += divisor == n
                start, stop = base, base + width
                break
        for a in range(start, stop):
            work["scalar_gcds"] += 1
            divisor = gcd(constant + a*a, n)
            if 1 < divisor < n:
                return NormHit(divisor, a, b)
    return None
