"""Recover balanced divisors with an unknown large common divisor of e and r±1.

The public inputs are N, e, a lower bound G on that common divisor, and an
optional interval. No factorization of N, e, or the common divisor is used.
"""

from math import gcd, isqrt
from pathlib import Path
import sys

from sympy import Poly, Symbol

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.affine_prime_window import _integer_roots, lll_reduce
from src.trackE.coppersmith import build_basis


def _box(modulus, exponent, center, radius):
    if not all(isinstance(v, int) for v in (modulus, exponent, center)):
        raise ValueError('modulus, exponent and center must be integers')
    if modulus < 2 or exponent < 2 or gcd(modulus, exponent) != 1:
        raise ValueError('require coprime modulus and exponent, both at least two')
    upper = isqrt(2*modulus)
    radius = upper if radius is None else radius
    if not isinstance(radius, int) or radius < 0:
        raise ValueError('radius must be a nonnegative integer')
    lower = max(2, isqrt((modulus-1)//2)+1, center-radius)
    return lower, min(modulus-1, upper, center+radius), radius


def _certifies(modulus, exponent, lower, bound, radius, dimension, m, n):
    # Exact form of sqrt(D) 2^((D-1)/4) det^(1/D) < lower^m bound^n.
    d = dimension
    det4 = (modulus**(2*m*(m+1)) * exponent**(2*n*(n+1))
            * radius**(2*d*(d-1)))
    return (1 << (d*(d-1))) * d**(2*d) * det4 < (lower**m * bound**n)**(4*d)


def _multiplicity(base, lower, dimension):
    """Maximize the determinant margin in one prime-independent modulus."""
    target = lower**dimension
    lo, hi = 0, dimension
    while hi-lo > 1:
        mid = (lo+hi)//2
        if base**mid <= target:
            lo = mid
        else:
            hi = mid
    return lo


def _parameters(modulus, exponent, lower, bound, radius, max_dimension):
    for d in range(2, max_dimension+1):
        # Raising m to m+1 improves the margin exactly when N^(m+1) < lower^d.
        # The e term has the identical independent optimization.
        m = _multiplicity(modulus, lower, d)
        n = _multiplicity(exponent, bound, d)
        if _certifies(modulus, exponent, lower, bound, radius, d, m, n):
            return d, m, n
    return None


def minimum_hint_bound(modulus, exponent, dimension, *, center=0, radius=None):
    """Smallest integer G certified in this exact dimension, or None.

    The determinant test is monotone in G. This binary search uses public
    integer arithmetic only. Passing its result to weighted_order_divisors
    takes a single lattice reduction for the chosen sign.
    """
    lower, upper, radius = _box(modulus, exponent, center, radius)
    if not isinstance(dimension, int) or dimension < 2:
        raise ValueError('dimension must be an integer at least two')
    if lower > upper:
        return None
    if lower == upper:
        return 1
    m = _multiplicity(modulus, lower, dimension)

    def works(bound):
        n = _multiplicity(exponent, bound, dimension)
        return _certifies(modulus, exponent, lower, bound, radius, dimension, m, n)

    if not works(exponent):
        return None
    lo, hi = 0, exponent
    while hi-lo > 1:
        mid = (lo+hi)//2
        if works(mid):
            hi = mid
        else:
            lo = mid
    return hi


def _weighted_basis(modulus, exponent, sign, center, radius, d, m, n):
    # With t=r-center and C=0 mod N, C=sign mod e, f(t)=t+center-C.
    product = modulus*exponent
    constant = (center-sign*modulus*pow(modulus, -1, exponent)) % product
    if 2*constant > product:
        constant -= product
    # Above M=max(m,n), T^(i-M) f^M is triangularly equivalent to f^i.
    rows = build_basis(1, constant, radius, d, max(m, n))
    for i, row in enumerate(rows):
        weight = modulus**max(m-i, 0) * exponent**max(n-i, 0)
        for j in range(i+1):
            row[j] *= weight
        # Unimodular size reduction preserves divisibilities and keeps the
        # large CRT representative out of nonleading coefficients where possible.
        for j in reversed(range(i)):
            pivot = rows[j][j]
            quotient = (row[j]+pivot//2)//pivot
            if quotient:
                for k in range(j+1):
                    row[k] -= quotient*rows[j][k]
    return rows


def weighted_order_divisors(modulus, exponent, bound, *, sign=1, center=0, radius=None,
                            max_dimension=None, stats=None):
    """Return all proper balanced divisors r in the box with gcd(e,r-sign)>=G.

    Balanced means sqrt(N/2) <= r <= sqrt(2N). The search box is
    |r-center| <= radius. The defaults cover the whole balanced interval.
    Raises ValueError if the public bound lacks a determinant certificate
    through max_dimension. The default cap is linear in the input bit length.
    """
    lower, upper, radius = _box(modulus, exponent, center, radius)
    if sign not in (-1, 1) or not isinstance(bound, int) or not 1 <= bound <= exponent:
        raise ValueError('require sign in {-1,1} and integer bound in [1,exponent]')
    if max_dimension is None:
        max_dimension = max(8, modulus.bit_length()+exponent.bit_length()+2)
    if not isinstance(max_dimension, int) or max_dimension < 2:
        raise ValueError('max_dimension must be an integer at least two')
    if stats is not None:
        stats.update(dimension=0, modulus_multiplicity=0, exponent_multiplicity=0)
    if lower > upper:
        return ()
    if lower == upper:
        candidates = [lower]
    else:
        params = _parameters(modulus, exponent, lower, bound, radius, max_dimension)
        if params is None:
            raise ValueError('no weighted lattice certificate through the dimension cap')
        d, m, n = params
        rows = _weighted_basis(modulus, exponent, sign, center, radius, d, m, n)
        row = list(map(int, lll_reduce(rows)[0]))
        if not any(row) or d*sum(v*v for v in row) >= (lower**m * bound**n)**2:
            raise ArithmeticError('LLL row failed the exact short-vector certificate')
        powers = [radius**j for j in range(d)]
        if any(v % power for v, power in zip(row, powers)):
            raise ArithmeticError('LLL row failed the exact coordinate scaling check')
        polynomial = Poly.from_list([v//power for v, power in zip(row, powers)][::-1],
                                    Symbol('T'), domain='ZZ')
        candidates = [center+t for t in _integer_roots(polynomial, radius)]
        if stats is not None:
            stats.update(dimension=d, modulus_multiplicity=m, exponent_multiplicity=n)
    return tuple(sorted(r for r in candidates if lower <= r <= upper
                        and modulus % r == 0 and gcd(exponent, r-sign) >= bound))
