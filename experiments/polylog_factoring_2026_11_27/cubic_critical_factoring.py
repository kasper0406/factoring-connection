"""Factor from additive-unit estimates of the projected cubic at a=1, B=2."""

from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from math import prod
from sympy import Poly, integer_nthroot, mobius, symbols


def padding(n):
    value = 1
    while value <= 4*n*n:
        value *= 5
    return value


def inverse_log_coefficient(m):
    numerator, power = 0, 1
    while m % power == 0:
        numerator += power * int(mobius(m//power))
        power *= 2
    return Fraction(numerator, m)


def factor_with_projection(n, approximate):
    """approximate(modulus) returns a rational real estimate with absolute error <=1."""
    if n < 1:
        raise ValueError("Expected a positive integer")
    factors = []
    for prime in (2, 3, 5):
        while n % prime == 0:
            factors.append(prime)
            n //= prime
    if n == 1:
        return sorted(factors)
    bits = n.bit_length()
    pad = padding(n)
    limit = bits + 5

    @lru_cache(None)
    def exact_product(t):
        auxiliary = pad**t
        multiplier = Fraction(n**(2*t)*auxiliary**2, auxiliary**3+1)
        return round(multiplier * Fraction(approximate((n*pad)**(4*t))))

    weights = [inverse_log_coefficient(m) for m in range(1, limit+1)]
    powers = []
    for j in range(1, bits+1):
        with localcontext() as context:
            context.prec = 3*bits*j + 3*bits + 30
            log_n = Fraction(Decimal(n).ln())
            power_sum = Fraction(0)
            for m, weight in enumerate(weights, 1):
                if weight:
                    t = j*m
                    log_product = Fraction(Decimal(exact_product(t)).ln()) - 3*t*log_n
                    power_sum += weight * log_product
        powers.append(round(power_sum * n**(3*j)))

    coefficients = [1]
    for order in range(1, bits+1):
        numerator = -sum(coefficients[order-i] * powers[i-1] for i in range(1, order+1))
        coefficient, remainder = divmod(numerator, order)
        if remainder:
            raise ValueError("Inconsistent projected power sums")
        coefficients.append(coefficient)
    roots = Poly.from_list(coefficients, symbols("x")).ground_roots()
    components = []
    for root in roots:
        if root == 0:
            continue
        ratio, exact = integer_nthroot(int(root), 3)
        if not exact or ratio <= 0 or n % ratio:
            raise ValueError("Invalid component root")
        component = n//ratio
        components.append(component)
        for exponent in range(component.bit_length(), 0, -1):
            prime, exact = integer_nthroot(component, exponent)
            if exact:
                factors.extend([int(prime)]*exponent)
                break
    if prod(components) != n:
        raise ValueError("Recovered components do not reconstruct input")
    return sorted(factors)
