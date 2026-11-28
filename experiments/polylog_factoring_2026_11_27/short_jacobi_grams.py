"""Deterministic small-coordinate search for the mixed ternary Gram."""

from math import gcd, isqrt

from sympy import jacobi_symbol


def select_x(n, limit):
    """Return (first negative-Jacobi x, None), (None, proper gcd), or (None, None).

    Test x=0,...,limit in order. The large-sieve guarantee concerns products
    of two distinct balanced odd primes. The tests themselves need only n.
    """
    if n <= 1 or n % 2 == 0 or limit < 0:
        raise ValueError("expected odd n>1 and a nonnegative limit")
    for x in range(limit+1):
        a = x*x+1
        symbol = jacobi_symbol(-a, n)
        if symbol == -1:
            return x, None
        if symbol == 0:
            divisor = gcd(a, n)
            if divisor < n:
                return None, divisor
    return None, None


def first_vector_orbits(x):
    """Yield each sorted nonnegative ternary first vector of squared norm x²+1."""
    if x < 0:
        raise ValueError("expected a nonnegative x")
    a = x*x+1
    for i in range(isqrt(a//3)+1):
        rest = a-i*i
        for j in range(i, isqrt(rest//2)+1):
            k = isqrt(rest-j*j)
            if j*j+k*k == rest:
                yield i, j, k
