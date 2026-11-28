"""Evaluate gcds of both sides of a stretch of the central-binomial recurrence."""

from math import gcd, isqrt

if __package__:
    from .adaptive_binomial_lift import _fold
    from .sqrt_prime_marks import _DecimalRing
else:
    from adaptive_binomial_lift import _fold
    from sqrt_prime_marks import _DecimalRing
from src.trackC.polyback import multipoint_eval, product_tree_linear


def transition_gcds(n, h, d, length):
    """Numerator and denominator gcds for length steps after t=floor(h*n/d).

    A single product polynomial covers both progressions. No factor of n
    or central binomial coefficient is needed.
    """
    if n < 3 or n % 2 == 0 or d < 8 or d & (d-1) or h < 1 or h % 2 == 0 or length < 0:
        raise ValueError("expected odd n and h, dyadic d>=8 and a nonnegative length")
    if length == 0:
        return 1, 1
    remainder = h*n % d
    block = isqrt(length)
    rows = length//block
    ring = _DecimalRing(n)
    polynomial = product_tree_linear(ring, [d*j for j in range(1, block+1)])[-1][0]
    starts = [d*block*i-remainder for i in range(rows)]
    half = d//2
    values = multipoint_eval(ring, polynomial, [v-half for v in starts] + starts)
    numerator, denominator = 1, 1
    for j in range(rows):
        numerator = numerator*values[j] % n
        denominator = denominator*values[rows+j] % n
    for j in range(rows*block+1, length+1):
        term = d*j-remainder
        numerator = numerator*(term-half) % n
        denominator = denominator*term % n
    return gcd(numerator, n), gcd(denominator, n)


def opposite_window(n, a, d):
    """Two gcds for an opposite-lift state. A true state gives one proper factor.

    The public promise is n=pq with odd primes p<q<2p. An inconsistent
    state can also expose a factor. The input guard guarantees d/2<p.
    """
    if (n < 15 or n % 2 == 0 or d < 16 or d & (d-1) or not 0 < a < d//2
            or (a*a+n) % d or d//2 >= isqrt((n-1)//2)+1):
        raise ValueError("expected a balanced opposite state at d>=16 below the factor bound")
    h = _fold((d//2-1)*pow(a, -1, d), d//2)
    length = (isqrt(2*n)+1)//d+2
    return transition_gcds(n, h, d, length)
