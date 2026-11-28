"""Generate factored coordinates for the prime-conductor three-square reduction."""

from collections import Counter
from dataclasses import dataclass
from math import gcd

from sympy import jacobi_symbol

if __package__:
    from .factored_gaussian_packing import certified_prime
    from .fixed_field_packing import Parameters, packing_params
    from .fully_factored_packing import small_factorization
else:
    from factored_gaussian_packing import certified_prime
    from fixed_field_packing import Parameters, packing_params
    from fully_factored_packing import small_factorization


def uniform_factored(bound, rng, prime_test=certified_prime):
    """Kalai's exact uniform draw from 1..bound, with prime multiplicities."""
    if bound < 1:
        raise ValueError("expected a positive upper bound")
    while True:
        value = 1
        factors = Counter()
        item = rng.randrange(1, bound+1)
        while item > 1:
            if prime_test(item):
                value *= item
                if value > bound:
                    break
                factors[item] += 1
            item = rng.randrange(1, item+1)
        if value <= bound and rng.randrange(bound) < value:
            return value, factors


@dataclass(frozen=True)
class CoordinatePacking:
    parameters: Parameters
    t: int
    norm_factorization: dict
    scale_factorization: dict
    height: int
    draws: int

    def query_factorization(self, ell):
        """Include the extra marking prime chosen by recover_factor."""
        return dict(Counter(self.scale_factorization) + Counter([ell]))


def factored_coordinate_packing(n, x, rng, *, prime_test=certified_prime):
    """Sample a prime conductor T and give the complete factors of M and Q.

    The primality predicate must be exact. The bit-complexity bound in the
    accompanying note uses a polynomial-time predicate.
    """
    a = x*x+1
    if n < 3 or n % 2 == 0 or x < 0 or jacobi_symbol(-a, n) != -1:
        raise ValueError("expected a mixed Jacobi parameter on an odd input")
    base = Counter(small_factorization(a))
    height = n*n
    draws = 0
    while True:
        for _ in range(height.bit_length()**2):
            draws += 1
            b, b_factors = uniform_factored(height, rng, prime_test)
            t, t_factors = uniform_factored(height // n, rng, prime_test)
            if b <= n or t < 2:
                continue
            conductor = b*b + (n*t)**2
            if conductor > height*height or gcd(conductor, 2*a*n) != 1:
                continue
            if not prime_test(conductor):
                continue
            params = packing_params(n, x, b, t)
            radius = base + Counter({conductor: 2})
            scale = base + b_factors + t_factors + Counter([2])
            return CoordinatePacking(params, t, dict(radius), dict(scale), height, draws)
        height *= height
