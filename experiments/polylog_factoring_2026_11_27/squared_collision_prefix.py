"""Locate the first prime using coarse multiplicative estimates of second moments."""

from fractions import Fraction
from math import isqrt


def least_factor_from_second_moment(n, estimate, contraction=Fraction(1, 2), *, quartic=False):
    """At t, require two-sided factor max(1,contraction*sqrt(t)/local_bound).

    The quartic option assumes a product of two distinct odd primes. Otherwise
    the callback estimates the gcd second moment and n can be any integer >=2.
    """
    c = Fraction(contraction)
    if n < 2 or not 0 < c < 1:
        raise ValueError("Expected n >= 2 and a contraction strictly between zero and one")
    if n % 2 == 0:
        return 2
    local_bound = 9 if quartic else 1
    limit = isqrt(n)
    low = min(limit, 2 + int((local_bound/c)**2) + int(1/(1-c)))
    for d in range(3, low+1, 2):
        if n % d == 0:
            return d

    def detected(t):
        # Compare against (local_bound**2 * t) times the allowed factor.
        value = Fraction(estimate(n, t))
        return value**2 > (local_bound*c)**2 * t**3

    while low < limit:
        high = min(limit, int(low/c))
        if detected(high):
            while high-low > 1:
                middle = (low+high)//2
                if detected(middle):
                    high = middle
                else:
                    low = middle
            return high
        low = high
    return n


def critical_semiprime_estimate(n, t):
    """Native factor-sqrt(t) estimate on distinct semiprimes for 1 <= t <= sqrt(n)."""
    from math import gcd

    if not 1 <= t <= isqrt(n):
        raise ValueError("Cutoff outside the short-prefix range")
    p = gcd(t, n)
    if p > 1:
        # The other prime exceeds sqrt(n).
        return t + (p*p - 1)*(t//p)
    return isqrt(t**3)
