"""Search the exact norm ellipse by sharing block products across several rows."""

from collections import defaultdict
from math import gcd, isqrt
from typing import NamedTuple

from sympy import jacobi_symbol

if __package__:
    from .large_conic_batch import NormHit, _quadratic_product
    from .sqrt_prime_marks import _DecimalRing
else:
    from large_conic_batch import NormHit, _quadratic_product
    from sqrt_prime_marks import _DecimalRing

from src.trackC.polyback import multipoint_eval


def ellipse_tiles(limit, norm):
    """Return dyadic column intervals and rows, partitioning a²+norm*b²<=limit."""
    return _form_tiles(limit, (1, 0, norm), nonnegative_a=True)


def _form_tiles(limit, form, *, nonnegative_a=False):
    leading, middle, last = form
    delta = 4*leading*last-middle*middle
    intervals = defaultdict(list)
    for b in range(1, isqrt((4*leading*limit)//delta)+1):
        root = isqrt(4*leading*limit-delta*b*b)
        start = -((root+middle*b)//(2*leading))
        end = 1+(root-middle*b)//(2*leading)
        if nonnegative_a:
            start = max(0, start)
        while start < end:
            width = 1 << ((end-start).bit_length()-1)
            if start:
                width = min(width, start & -start)
            intervals[width, start].append(b)
            start += width
    return intervals


def _product(ring, level):
    while len(level) > 1:
        following = [ring.mul(level[i], level[i+1]) for i in range(0, len(level)-1, 2)]
        if len(level) % 2:
            following.append(level[-1])
        level = following
    return level[0]


def _pooled_poly(ring, form, rows, stride):
    leading, middle, last = form
    if leading == 1 and middle == 0:
        return _product(ring, [_quadratic_product(ring, last*b*b, stride) for b in rows])
    factors = [ring.make([leading*j*j+middle*b*j+last*b*b, 2*leading*j+middle*b, leading])
               for b in rows for j in range(stride)]
    return _product(ring, factors)


def _statistics(stats):
    work = {} if stats is None else stats
    work.update(tiles=0, quadratic_factors=0, block_values=0, scalar_gcds=0,
                saturated_blocks=0, max_degree=0)
    return work


def pooled_conic_factor(n, x, *, stats=None):
    """Return a proper gcd with a²+(x²+1)b² in [1,floor(sqrt(8n))], or None.

    Both coordinates are nonnegative and b>=1. The order of visited points
    differs from the older row-wise scan. All returned hits obey the same ellipse.
    """
    if n < 9 or n % 2 == 0 or x < 0:
        raise ValueError('expected odd n>=9 and nonnegative x')
    norm, limit = x*x+1, isqrt(8*n)
    return _search(_DecimalRing(n), (1, 0, norm), ellipse_tiles(limit, norm), _statistics(stats))


def _search(ring, form, intervals, work):
    leading, middle, last = form
    n = ring.N
    for width, start in sorted(intervals, key=lambda k: (-k[0], k[1])):
        rows = intervals[width, start]
        # A power-of-two stride divides the full tile and approximates sqrt(width/#rows).
        ideal = max(1, isqrt(width//len(rows)))
        stride = 1 << (ideal.bit_length()-1)
        polynomial = _pooled_poly(ring, form, rows, stride)
        starts = range(start, start+width, stride)
        values = multipoint_eval(ring, polynomial, starts)
        work['tiles'] += 1
        work['quadratic_factors'] += stride*len(rows)
        work['block_values'] += len(starts)
        work['max_degree'] = max(work['max_degree'], 2*stride*len(rows))
        for base, value in zip(starts, values):
            divisor = gcd(value, n)
            if divisor == 1:
                continue
            work['saturated_blocks'] += divisor == n
            for b in rows:
                for a in range(base, base+stride):
                    work['scalar_gcds'] += 1
                    divisor = gcd(leading*a*a+middle*a*b+last*b*b, n)
                    if 1 < divisor < n:
                        return NormHit(divisor, a, b)
            # Every examined norm lies strictly between zero and n for n>=9.
            raise ArithmeticError('a nonunit block contained no proper factor')
    return None


def reduced_forms(norm):
    """Primitive positive reduced forms of discriminant -4*norm."""
    for a in range(1, isqrt(4*norm//3)+1):
        for b in range(-a, a+1):
            c, remainder = divmod(b*b+4*norm, 4*a)
            if remainder or c < a or gcd(a, b, c) != 1:
                continue
            if b < 0 and (-b == a or a == c):
                continue
            yield a, b, c


class FormHit(NamedTuple):
    form: tuple
    hit: NormHit


def pooled_reduced_form_factor(n, x, *, stats=None):
    """Split a balanced odd semiprime with a mixed Jacobi symbol at -(x²+1)."""
    norm = x*x+1
    if n < 9 or n % 2 == 0 or x < 0 or jacobi_symbol(-norm, n) != -1:
        raise ValueError('require odd n>=9, nonnegative x, and Jacobi(-A,n)=-1')
    ring, limit = _DecimalRing(n), isqrt(2*n)
    work = _statistics(stats)
    work['forms'] = 0
    for form in reduced_forms(norm):
        work['forms'] += 1
        # A prime r represented with second coordinate zero equals the leading coefficient.
        divisor = gcd(form[0], n)
        if 1 < divisor < n:
            return FormHit(form, NormHit(divisor, 1, 0))
        result = _search(ring, form, _form_tiles(limit, form), work)
        if result:
            return FormHit(form, result)
    return None
