"""Factor reduction on a sphere of prescribed squareclass x²+1."""

from math import gcd
from typing import NamedTuple

from sympy import jacobi_symbol

if __package__:
    from .packed_ternary_decision import _primes_until, _residue_candidates
    from .short_jacobi_grams import first_vector_orbits
else:
    from packed_ternary_decision import _primes_until, _residue_candidates
    from short_jacobi_grams import first_vector_orbits


class Parameters(NamedTuple):
    B: int
    P: int
    Q: int
    T: int
    M: int


def packing_params(n, x, B, t):
    """Use the Gaussian integer B − i*t*n. No factorization is required."""
    a = x*x + 1
    if n < 3 or n % 2 == 0 or x < 0 or jacobi_symbol(-a, n) != -1:
        raise ValueError("expected an odd input and a nonnegative mixed Jacobi parameter")
    if B <= 0 or t <= 0:
        raise ValueError("both Gaussian parameters must be positive")
    c = n*t
    p = B*B - c*c - 2*x*B*c
    q = 2*a*B*t
    radius_factor = B*B + c*c
    if p <= 0 or gcd(p, q) != 1 or gcd(radius_factor, 2*a*n) != 1:
        raise ValueError("the Gaussian parameters must give positive coprime scales and a unit conductor")
    return Parameters(B, p, q, radius_factor, a*radius_factor**2)


def fixed_field_parameters(n, x):
    """Choose an input-only packing with t=1 and B a power of two times A."""
    a = x*x + 1
    b = 2*a
    threshold = (2*x + 2)*n
    while b <= threshold:
        b *= 2
    return packing_params(n, x, b, 1)


def recover_factor(n, x, decide, gaussian=None):
    """Recover a prime from exact existence of points in specified residue classes.

    decide(M, modulus, residue) may count only primitive points. The promise
    is n=p*q, distinct odd primes, and Jacobi(−x²−1,n)=−1. An optional pair
    (B,t) supplies the Gaussian parameters in place of the default.
    """
    params = fixed_field_parameters(n, x) if gaussian is None else packing_params(n, x, *gaussian)
    first_vectors = tuple(first_vector_orbits(x))
    excluded = 210*n*params.Q*params.M
    e2 = e3 = 0
    modulus = 1
    for prime in _primes_until(excluded, n**6):
        query_modulus = params.Q*prime
        pairs = set()
        for u in first_vectors:
            for z in _residue_candidates(u, -n*x, n*n, prime):
                residue = tuple((params.P*i - params.Q*j) % query_modulus for i, j in zip(u, z))
                if decide(params.M, query_modulus, residue):
                    a, b, c = (i*i % prime for i in z)
                    pairs.add(((a*b + b*c + a*c) % prime, a*b*c % prime))
        pairs.discard((0, 0))
        if len(pairs) > 1:
            raise ValueError("inconsistent symmetric coefficients at the test prime")
        r2, r3 = next(iter(pairs), (0, 0))
        inverse = pow(modulus, -1, prime)
        e2 += modulus*((r2 - e2)*inverse % prime)
        e3 += modulus*((r3 - e3)*inverse % prime)
        modulus *= prime
    divisor = gcd(n, e2, e3)
    if not (0 < e2 < n**4 and 0 <= e3 < n**6 and 1 < divisor < n):
        raise ValueError("the reconstructed coefficients did not give a factor")
    return divisor
