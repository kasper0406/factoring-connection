"""Read a centered residue window at every prime of a balanced squarefree semiprime.

The public promise is n=p*q for primes p<q<=2*p. Only n, x and the requested
radius enter the algorithm. The lattice is the affine-hint basis of
src/trackE/coppersmith.py, with f(T)=T+x. A certified row vanishes over Z at
*every* in-window residue, even when the two primes have different residues.
"""

from functools import lru_cache
from math import gcd, isqrt, lcm
from pathlib import Path
import sys

from sympy import Matrix, Poly, Symbol

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.trackE.coppersmith import build_basis, rigorous_ok, lll_reduce as _native_lll


def lll_reduce(rows):
    try:
        return _native_lll(rows)
    except ModuleNotFoundError as error:
        if error.name != 'flint':
            raise
    basis = Matrix(rows)
    try:
        return basis.lll().tolist()
    except (AssertionError, OverflowError):
        # Some large coefficients exceed SymPy's floating-point rounding step.
        from experiments.polylog_factoring_2026_11_27.level11_sampler import _reduced_basis
        return _reduced_basis(basis).tolist()


def window_radius(n):
    """Return floor(n^(1/4)/8), using only integer square roots."""
    if n < 2:
        raise ValueError('n must be at least 2')
    return isqrt(isqrt(n)) // 8


@lru_cache(maxsize=128)
def _parameters(n, radius, lower):
    # For fixed d the bound is strongest at m=floor(log_n(lower^d)):
    # increasing m by one multiplies the two sides' ratio by (n^(m+1)/lower^d)^4.
    lower_power, n_power, m = lower, 1, 0
    for d in range(2, max(8, n.bit_length()+2) + 1):
        lower_power *= lower
        while n_power*n <= lower_power:
            m += 1
            n_power *= n
        best = min(d-1, max(1, m))
        if rigorous_ok(n, lower, radius, d, best):
            return d, best
    raise ValueError('no affine lattice certificate for this radius and divisor lower bound')


def _integer_roots(polynomial, radius):
    """Isolate possible integers with Sturm bisection on half-integer boundaries."""
    if polynomial.degree() < 1:
        return []
    first = polynomial.sqf_part()
    second = first.diff()
    sturm = [first, second]
    while second.degree() > 0:
        remainder = first.prem(second)
        # prem multiplies the ordinary remainder by LC(second)^(deg gap + 1).
        scale_is_negative = second.LC() < 0 and (first.degree()-second.degree()+1) % 2
        remainder = remainder.primitive()[1]
        first, second = second, remainder if scale_is_negative else -remainder
        sturm.append(second)
    rows = [list(map(int, member.all_coeffs())) for member in sturm]

    def variations(twice_x):
        previous, changes = 0, 0
        for row in rows:
            value = 0
            for i, coeff in enumerate(row):
                value = value*twice_x + (coeff << i)
            if value:
                changes += previous*value < 0
                previous = 1 if value > 0 else -1
        return changes

    roots = []
    stack = [(-radius, radius, variations(-2*radius-1), variations(2*radius+1))]
    while stack:
        lo, hi, left, right = stack.pop()
        if left == right:
            continue
        if lo == hi:
            if polynomial.eval(lo) == 0:
                roots.append(lo)
            continue
        mid = (lo+hi)//2
        middle = variations(2*mid+1)
        stack.extend(((lo, mid, left, middle), (mid+1, hi, middle, right)))
    return roots


def affine_divisors(n, center, radius, lower, *, stats=None):
    """All distinct gcd(n,center+y)>1 at least lower with |y|<=radius, sorted.

    Accepts any n>=2. Radii beyond the exact public lattice certificate (through
    dimension max(8,bit_length(n)+2)) raise ValueError. Radii <=2 use direct gcds.
    """
    if n < 2 or not isinstance(radius, int) or radius < 0 or not 1 <= lower <= n:
        raise ValueError('require n>=2, integer radius>=0 and 1<=lower<=n')
    if stats is None:
        stats = {}
    stats.update(method='tiny', radius=radius, dimension=0, offsets_tested=0)
    if radius <= 2:
        offsets = range(-radius, radius+1)
    else:
        d, m = _parameters(n, radius, lower)
        stats.update(method='lattice', dimension=d)
        c = center % n
        if c > n//2:
            c -= n
        row = list(map(int, lll_reduce(build_basis(n, c, radius, d, m))[0]))
        if not any(row) or d*sum(v*v for v in row) >= lower**(2*m):
            raise ArithmeticError('reduced row failed the exact short-vector certificate')
        coeffs = [v // radius**j for j, v in enumerate(row)]
        polynomial = Poly.from_list(coeffs[::-1], Symbol('T'), domain='ZZ')
        offsets = _integer_roots(polynomial, radius)
    divisors = set()
    for offset in offsets:
        divisor = gcd(center+offset, n)
        if divisor > 1 and divisor >= lower:
            divisors.add(divisor)
        stats['offsets_tested'] += 1
    return tuple(sorted(divisors))


def prime_window_divisor(n, x, *, radius=None, stats=None):
    """Return the product of the r | n with x mod r in [−radius,radius].

    The input promises n=p*q, p<q<=2*p. The radius defaults to window_radius(n)
    and may be any smaller nonnegative integer. An optional dict receives the
    method, radius, lattice dimension, and number of candidate offsets tested.
    """
    maximum = window_radius(n)
    radius = maximum if radius is None else radius
    if not isinstance(radius, int) or not 0 <= radius <= maximum:
        raise ValueError('radius must be an integer in the certified range')
    x %= n
    if min(x, n-x) <= radius:
        if stats is not None:
            stats.update(method='central', radius=radius, dimension=0, offsets_tested=0)
        return n
    lower = isqrt((n-1)//2) + 1  # ceil(sqrt(n/2))
    return lcm(*affine_divisors(n, x, radius, lower, stats=stats))


def window_idempotent(n, divisor):
    """Encode a unitary divisor by the scalar e=1 on it and e=0 on n/divisor."""
    if n < 2 or divisor < 1 or n % divisor or gcd(divisor, n//divisor) != 1:
        raise ValueError('a positive unitary divisor of n is required')
    if divisor == 1:
        return 0
    complement = n//divisor
    return complement * pow(complement, -1, divisor) % n
