"""Recover a divisor from a rational root k/h of T+e modulo an unknown large divisor."""

from fractions import Fraction
from math import gcd, isqrt

from sympy import Poly, Rational, Symbol

if __package__:
    from .affine_prime_window import _parameters, build_basis, lll_reduce, primitive_sturm_rows, window_radius
else:
    from affine_prime_window import _parameters, build_basis, lll_reduce, primitive_sturm_rows, window_radius


def _rational_roots(polynomial, h_bound, k_bound):
    """Positive roots k/h in lowest terms, with k<=k_bound and h<=h_bound."""
    if polynomial.degree() < 1:
        return []
    # Sturm variations count roots in (left,right], including a root at a
    # dyadic subdivision point exactly once. Positive row scaling keeps signs.
    chain = primitive_sturm_rows(polynomial)

    def variations(numerator, denominator):
        previous = changes = 0
        for row in chain:
            value, scale = 0, 1
            for coefficient in row:
                value = value*numerator + scale*coefficient
                scale *= denominator
            if value:
                sign = 1 if value > 0 else -1
                changes += previous*sign < 0
                previous = sign
        return changes

    roots = set()
    stack = [(0, k_bound, 1, variations(0, 1), variations(k_bound, 1))]
    while stack:
        left, right, denominator, at_left, at_right = stack.pop()
        if at_left == at_right:
            continue
        if 2*h_bound*h_bound*(right-left) <= denominator:
            # Two different rationals with denominator <= H differ by >=1/H².
            candidate = Fraction(left+right, 2*denominator).limit_denominator(h_bound)
            k, h = candidate.numerator, candidate.denominator
            if 0 < k <= k_bound and polynomial.eval(Rational(k, h)) == 0:
                roots.add((h, k))
            continue
        mid, doubled = left + right, 2*denominator
        at_mid = variations(mid, doubled)
        stack.append((2*left, mid, doubled, at_left, at_mid))
        stack.append((mid, 2*right, doubled, at_mid, at_right))
    return sorted(roots)


def _order_certificate(n, e):
    modulus = n*e
    lower = (e+1)//2  # k*r = e*h*r/(r+offset) >= e/2 on an exact order pair
    # The dyadic divisor-window bound from affine-divisor-scales.md works
    # for every lower scale. Its proof chooses dimension bits, multiplicity k.
    bits = (modulus-1).bit_length()
    k = lower.bit_length()-1
    exponent = k*(k-1)//max(1, bits-1) - 2
    capacity = 1 << exponent if exponent >= 0 else 0
    return modulus, lower, capacity


def projective_order_window(n, e, h_bound, k_bound, *, lower=None, order_offset=-1,
                            order_certificate=False, stats=None):
    """Return sorted proper factors detected from 1<=h<=H, 1<=k<=K, gcd(h,k)=1.

    By default the certificate covers every pair for which k-order_offset*e*h
    vanishes modulo a divisor of n >= lower. With order_certificate=True,
    it instead covers exact order pairs k*r=abs(k-order_offset*e*h), k | e,
    for a divisor r >= lower of n. It uses the larger modulus n*e in this case.
    The default lower assumes n=p*q with p<q<=2p.
    The default order_offset=-1 targets r-1; use +1 for r+1. Both the
    congruence gcd and the exact factor candidate e*h/k-order_offset are checked.
    An uncertifiable large box raises ValueError.
    """
    if n < 2 or e < 1 or h_bound < 1 or k_bound < 1:
        raise ValueError('n>=2 and positive exponent and bounds required')
    if order_offset not in (-1, 1):
        raise ValueError('order_offset must be -1 or 1')
    lower = isqrt((n-1)//2) + 1 if lower is None else lower
    if not 1 <= lower <= n:
        raise ValueError('lower must be between 1 and n')
    stats = {} if stats is None else stats
    modulus, divisor_lower = n, lower
    if order_certificate:
        modulus, divisor_lower, _ = _order_certificate(n, e)
    stats['certificate'] = 'n_times_e' if order_certificate else 'n'
    if h_bound*k_bound <= 16:
        stats.update(method='tiny', dimension=0)
        pairs = [(h, k) for h in range(1, h_bound+1) for k in range(1, k_bound+1) if gcd(h, k) == 1]
    else:
        dimension, multiplicity = _parameters(modulus, h_bound*k_bound, divisor_lower)
        center = (-order_offset*e) % modulus
        if center > modulus//2:
            center -= modulus
        # k^j h^(d-1-j) gets column weight K^j H^(d-1-j).
        degree = dimension - 1
        h_weights = [h_bound**(degree-j) for j in range(dimension)]
        weights = [k_bound**j * h_weights[j] for j in range(dimension)]
        affine = build_basis(modulus, center, k_bound, dimension, multiplicity)
        basis = [[v*w for v, w in zip(row, h_weights)] for row in affine]
        row = list(map(int, lll_reduce(basis)[0]))
        if not any(row) or dimension*sum(v*v for v in row) >= divisor_lower**(2*multiplicity):
            raise ArithmeticError('reduced row failed the exact projective certificate')
        coefficients = [v//w for v, w in zip(row, weights)]
        polynomial = Poly.from_list(coefficients[::-1], Symbol('T'), domain='ZZ')
        pairs = _rational_roots(polynomial, h_bound, k_bound)
        stats.update(method='lattice', dimension=dimension, multiplicity=multiplicity)
    factors = set()
    for h, k in pairs:
        common = gcd(k - order_offset*e*h, n)
        if common < lower:
            continue
        if 1 < common < n:
            factors.add(common)
        quotient, remainder = divmod(e*h, k)
        candidate = quotient - order_offset
        if not remainder and lower <= candidate < n and candidate > 1 and n % candidate == 0:
            factors.add(candidate)
    stats.update(h_bound=h_bound, k_bound=k_bound, candidates=len(pairs))
    return tuple(sorted(factors))


def small_power_image_factor(n, e, *, order_offset=-1, stats=None):
    """Try the public balanced-semiprime envelope with no image-order bound supplied.

    Its final H is floor(sqrt((lower-1)*capacity/e)), using the larger
    certified product bound from the n and n*e lattices, and
    K=floor(e*H/(lower-1)). If (r+order_offset)/gcd(r+order_offset,e)<=H at either prime r,
    the corresponding proper factor is detected by the final box or earlier.
    """
    if n < 3 or e < 1 or order_offset not in (-1, 1):
        raise ValueError('n>=3, a positive exponent and order_offset -1 or 1 are required')
    stats = {} if stats is None else stats
    lower = isqrt((n-1)//2) + 1
    ordinary = window_radius(n)
    _, _, lifted = _order_certificate(n, e)
    use_lifted = lifted > ordinary
    maximum = isqrt((lower-1)*max(ordinary, lifted)//e)
    stats.update(max_h=maximum, boxes=0)
    h_bound = min(1, maximum)
    while h_bound:
        k_bound = e*h_bound // (lower-1)
        if k_bound:
            stats['boxes'] += 1
            factors = projective_order_window(n, e, h_bound, k_bound, lower=lower,
                                              order_offset=order_offset, order_certificate=use_lifted,
                                              stats=stats)
            if factors:
                return factors[0]
        if h_bound == maximum:
            break
        h_bound = min(2*h_bound, maximum)
    return None
