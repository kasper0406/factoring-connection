"""Recover balanced divisors from partial congruence agreement.

A signed order hint is the single-coordinate case gcd(e,r±1)>=G. Multiple
coordinates may carry separate lower bounds or preselected multiplicities.
No factorization of any input modulus is used.
"""

from math import gcd, isqrt, prod
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
    return _profile_certifies((modulus, exponent), (lower, bound), radius, dimension, (m, n))


def _profile_certifies(moduli, bounds, radius, dimension, weights):
    # Exact form of sqrt(D) 2^((D-1)/4) det^(1/D) < lower^m bound^n.
    d = dimension
    det4 = radius**(2*d*(d-1)) * prod(e**(2*m*(m+1)) for e, m in zip(moduli, weights))
    divisor_bound = prod(g**m for g, m in zip(bounds, weights))
    return (1 << (d*(d-1))) * d**(2*d) * det4 < divisor_bound**(4*d)


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
    return _profile_parameters((modulus, exponent), (lower, bound), radius, max_dimension)


def _profile_parameters(moduli, bounds, radius, max_dimension):
    for d in range(2, max_dimension+1):
        weights = tuple(_multiplicity(e, g, d) for e, g in zip(moduli, bounds))
        if _profile_certifies(moduli, bounds, radius, d, weights):
            return (d, *weights)
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


def combine_congruences(congruences):
    """Return (product, CRT residue) for nonempty pairwise-coprime (e,a) pairs.

    gcd(product,r-residue) = product_i gcd(e_i,r-a_i). An aggregate lower
    bound G can thus be supplied as one hint (*combine_congruences(pairs), G).
    """
    product, residue = 1, 0
    for e, a in congruences:
        if not isinstance(e, int) or e < 2 or not isinstance(a, int) or gcd(e, product) != 1:
            raise ValueError('require pairwise-coprime integer moduli >=2 and integer residues')
        residue += product*((a-residue)*pow(product, -1, e) % e)
        product *= e
    if product == 1:
        raise ValueError('at least one congruence is required')
    return product, residue


def _weighted_basis(modulus, exponent, sign, center, radius, d, m, n):
    crt = combine_congruences(((modulus, 0), (exponent, sign)))
    return _profile_basis((modulus, exponent), crt, center, radius, d, (m, n))


def _profile_basis(moduli, crt, center, radius, d, weights):
    product, residue = crt
    constant = (center-residue) % product
    if 2*constant > product:
        constant -= product
    # Above the maximum multiplicity M, T^(i-M) f^M is triangularly equivalent to f^i.
    rows = build_basis(1, constant, radius, d, max(weights))
    for i, row in enumerate(rows):
        weight = prod(e**max(m-i, 0) for e, m in zip(moduli, weights))
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
    if sign not in (-1, 1):
        raise ValueError('require sign in {-1,1}')
    profile_stats = None if stats is None else {}
    try:
        return weighted_congruence_divisors(modulus, [(exponent, sign, bound)], center=center,
                                            radius=radius, max_dimension=max_dimension,
                                            stats=profile_stats)
    finally:
        if profile_stats:
            m, n = profile_stats['multiplicities']
            stats.update(dimension=profile_stats['dimension'], modulus_multiplicity=m,
                         exponent_multiplicity=n)


def _problem(modulus, pairs, center, radius):
    if not pairs:
        raise ValueError('at least one hint is required')
    moduli = (modulus, *(e for e, _ in pairs))
    crt = combine_congruences(((modulus, 0), *pairs))
    box = _box(modulus, crt[0]//modulus, center, radius)
    return moduli, crt, box


def _divisor_roots(moduli, crt, center, box, dimension, weights, stats, divisor_bound=None):
    lower, upper, radius = box
    if stats is not None:
        stats.update(dimension=0, multiplicities=(0,)*len(moduli), evaluation_bound_squared=0)
    if lower > upper:
        return ()
    if lower == upper:
        return (lower,) if moduli[0] % lower == 0 else ()
    rows = _profile_basis(moduli, crt, center, radius, dimension, weights)
    row = list(map(int, lll_reduce(rows)[0]))
    evaluation_bound_squared = dimension*sum(v*v for v in row)
    if not any(row) or (divisor_bound is not None and evaluation_bound_squared >= divisor_bound**2):
        raise ArithmeticError('LLL row failed the exact short-vector certificate')
    powers = [radius**j for j in range(dimension)]
    if any(v % power for v, power in zip(row, powers)):
        raise ArithmeticError('LLL row failed the exact coordinate scaling check')
    polynomial = Poly.from_list([v//power for v, power in zip(row, powers)][::-1],
                                Symbol('T'), domain='ZZ')
    if stats is not None:
        stats.update(dimension=dimension, multiplicities=tuple(weights),
                     evaluation_bound_squared=evaluation_bound_squared)
    candidates = (center+t for t in _integer_roots(polynomial, radius))
    return tuple(sorted(r for r in candidates if lower <= r <= upper and moduli[0] % r == 0))


def weighted_congruence_divisors(modulus, hints, *, center=0, radius=None,
                                max_dimension=None, stats=None):
    """Recover each balanced r in the box satisfying gcd(e_i,r-a_i)>=G_i.

    Supply triples (e_i,a_i,G_i). All e_i and N must be pairwise coprime.
    A product lower bound alone can be used with combine_congruences. It
    generally gives a different certificate from supplying separate G_i.
    """
    hints = tuple(hints)
    pairs = tuple((e, a) for e, a, _ in hints)
    moduli, crt, box = _problem(modulus, pairs, center, radius)
    lower, upper, radius = box
    bounds = (lower, *(g for _, _, g in hints))
    if any(not isinstance(g, int) or not 1 <= g <= e for e, _, g in hints):
        raise ValueError('each hint bound must be an integer between one and its modulus')
    if max_dimension is None:
        max_dimension = max(8, 2+sum(e.bit_length() for e in moduli))
    if not isinstance(max_dimension, int) or max_dimension < 2:
        raise ValueError('max_dimension must be an integer at least two')
    if stats is not None:
        stats.update(dimension=0, multiplicities=(0,)*len(moduli), evaluation_bound_squared=0)
    params = (0, *(0 for _ in moduli))
    if lower < upper:
        params = _profile_parameters(moduli, bounds, radius, max_dimension)
        if params is None:
            raise ValueError('no weighted lattice certificate through the dimension cap')
    d, *weights = params
    divisor_bound = prod(g**m for g, m in zip(bounds, weights))
    roots = _divisor_roots(moduli, crt, center, box, d, weights, stats, divisor_bound)
    return tuple(r for r in roots if all(gcd(e, r-a) >= g for e, a, g in hints))


def fixed_weight_divisors(modulus, congruences, dimension, multiplicities, *, center=0,
                          radius=None, stats=None):
    """Read verified divisors at public multiplicities without gcd lower bounds.

    Supply (e_i,a_i) pairs and one multiplicity for N followed by one per e_i.
    A balanced divisor in the box is guaranteed to appear if
    (r^m_0 product_i gcd(e_i,r-a_i)^m_i)^2 exceeds evaluation_bound_squared
    in stats. Other roots that are verified divisors are also returned.
    """
    pairs, weights = tuple(congruences), tuple(multiplicities)
    moduli, crt, box = _problem(modulus, pairs, center, radius)
    if not isinstance(dimension, int) or dimension < 2:
        raise ValueError('dimension must be an integer at least two')
    if len(weights) != len(moduli) or any(not isinstance(m, int) or not 0 <= m < dimension
                                          for m in weights):
        raise ValueError('supply one integer multiplicity in [0,dimension) per modulus')
    return _divisor_roots(moduli, crt, center, box, dimension, weights, stats)
