"""Audit full-input predictors constrained by a public integer polynomial."""

from math import isqrt
from sympy import Poly, symbols

if __package__:
    from .quartic_collision_powers import divisors_from_roots
else:
    from quartic_collision_powers import divisors_from_roots


def split_or_weight_bound(n, template):
    """Return (verified split, None) or (None, upper bound on product of hitting primes).

    n is an odd composite. Every prediction must be a modular root of template.
    """
    if n < 9 or n % 2 == 0 or template.is_zero or not template.domain.is_ZZ:
        raise ValueError("Expected odd composite and nonzero integral polynomial")
    for divisor in divisors_from_roots(n, template.ground_roots()):
        if divisor < n:
            return divisor, None
    k, minimum = 0, 1
    while minimum * (2*k+3) <= n:
        minimum *= 2*k+3
        k += 1
    # Integer 1-norm bounds the Euclidean norm in Hadamard's inequality.
    height = sum(abs(int(c)) for c in template.all_coeffs())
    return None, (10**k * n)**template.degree() * height**k


def square_root_template(n):
    t = (isqrt(n) - 1) // 2 * 2 + 1
    sign = 1 if t % 4 == 1 else -1
    j = symbols("j")
    # The factor j covers the zero fallback at noninvertible denominators.
    return t, Poly(j*(5*(t-1)*j - (2*t+sign)), j)


if __name__ == "__main__":
    for n in (10403, 1893893, 150235091, 1500011500021):
        estimate, template = square_root_template(n)
        divisor, bound = split_or_weight_bound(n, template)
        print(n, estimate, divisor, bound.bit_length() if bound is not None else None)
