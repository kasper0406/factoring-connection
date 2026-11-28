"""Certify both the squared radius and modulus in a packed ternary decision."""

from dataclasses import dataclass
from math import gcd

from sympy import jacobi_symbol

if __package__:
    from .factored_gaussian_packing import certified_prime
    from .fixed_field_packing import Parameters, packing_params
else:
    from factored_gaussian_packing import certified_prime
    from fixed_field_packing import Parameters, packing_params


@dataclass(frozen=True)
class FullyFactoredPacking:
    parameters: Parameters
    auxiliary_prime: int
    norm_factorization: dict
    scale_factorization: dict
    box: int
    prime_pairs: int


def _box_prime(lo, rng, prime_test):
    # Uniform on primes in [lo,2lo]. The algorithm only uses lo >= 4.
    while True:
        p = rng.randrange(lo | 1, 2*lo+1, 2)
        if prime_test(p):
            return p


def small_factorization(a):
    """Return the exact prime factorization of a positive short base by trial division."""
    if a < 1:
        raise ValueError("the base must be positive")
    base = {}
    remaining = a
    divisor = 2
    while divisor*divisor <= remaining:
        while remaining % divisor == 0:
            base[divisor] = base.get(divisor, 0)+1
            remaining //= divisor
        divisor += 1 if divisor == 2 else 2
    if remaining > 1:
        base[remaining] = 1
    return base


def fully_factored_packing(n, x, rng, *, prime_test=certified_prime):
    """Find prime B,q,T=B²+4n²q² by a widening search.

    Existence and positive density in all sufficiently large boxes hold for
    each fixed n. The theorem does not bound the initial successful box
    uniformly in n. prime_test must certify primes.
    """
    a = x*x+1
    if n < 3 or n % 2 == 0 or x < 0 or jacobi_symbol(-a, n) != -1:
        raise ValueError("expected an odd input and a nonnegative mixed Jacobi parameter")
    base = small_factorization(a)
    box = 1
    prime_pairs = 0
    while True:
        size = (2*n*box).bit_length()
        for _ in range(8*size**2):
            b = _box_prime(2*n*box, rng, prime_test)
            q = _box_prime(16*box, rng, prime_test)
            prime_pairs += 1
            if b == q or gcd(b*q, 2*a*n) != 1:
                continue
            norm = b*b + 4*n*n*q*q
            if gcd(norm, 2*a*n) != 1 or not prime_test(norm):
                continue
            params = packing_params(n, x, b, 2*q)
            radius_factors = dict(base)
            radius_factors[norm] = 2
            scale_factors = dict(base)
            scale_factors[2] = scale_factors.get(2, 0)+2
            scale_factors[b] = scale_factors[q] = 1
            return FullyFactoredPacking(params, q, radius_factors, scale_factors, box, prime_pairs)
        box *= 2
