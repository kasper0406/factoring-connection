"""Find a factor supported by a marked projective slope on the binary conic."""

from math import gcd, isqrt
from pathlib import Path
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.genusleak.forms import reduce_form
from experiments.polylog_factoring_2026_11_27.pooled_conic_batch import (
    _form_tiles, _search, _statistics,
)
from experiments.polylog_factoring_2026_11_27.sqrt_prime_marks import _DecimalRing


def slope_form(A, ell, line):
    """Reduced norm form on {(k,j): k*right-j*left=0 mod ell}."""
    left, right = line
    if A < 1 or ell < 1 or gcd(ell, left, right) != 1:
        raise ValueError('require A,ell positive and a unimodular projective generator')
    g = gcd(left, ell)
    step = ell//g
    offset = right*pow(left//g, -1, step) % step if step > 1 else 0
    # The initial lattice basis is (g,offset),(0,step).
    form = g*g+A*offset*offset, 2*A*offset*step, A*step*step
    return reduce_form(form)


def pooled_slope_factor(n, A, D, ell, line, *, stats=None):
    """Return a proper divisor of n or None, given the transport bound D and slope line.

    If a coprime pair (k,j) on the line has k²+A*j² <= 2D*floor(sqrt(2n))
    and 1<gcd(k²+A*j²,n)<n, a factor is guaranteed. Incidental proper factors
    are also returned. The caller must still verify its original point query.
    """
    if n < 3 or n % 2 == 0 or D <= 0:
        raise ValueError('require odd n>=3 and positive transport denominator D')
    form = slope_form(A, ell, line)
    if gcd(ell, n) != 1:
        raise ValueError('the fine slope modulus must be coprime to n')
    bound = 2*D*isqrt(2*n)
    work = _statistics(stats)
    work.update(bound=bound, form=form, axis_gcds=0, trial_divisions=0)
    common = gcd(n, D)
    if 1 < common < n:
        work['method'] = 'transport gcd'
        return common
    if bound >= n:
        # In this range sqrt(n) <= 2 sqrt(2) D, even if D shares all of n.
        work['method'] = 'small input'
        for divisor in range(3, isqrt(n)+1, 2):
            work['trial_divisions'] += 1
            if n % divisor == 0:
                return divisor
        return None
    a, _, c = form
    if a > bound:
        work['method'] = 'empty ellipse'
        return None
    work['axis_gcds'] = 1
    common = gcd(a, n)
    if common > 1:
        # a<=bound<n, hence this divisor is proper.
        work['method'] = 'shortest vector'
        return common
    if c > bound:
        # Collinear vectors are multiples of the shortest. Only its two signs
        # can be primitive in the ambient integer coordinates.
        work['method'] = 'single direction'
        return None
    work['method'] = 'pooled ellipse'
    result = _search(_DecimalRing(n), form, _form_tiles(bound, form), work)
    return result.factor if result else None
