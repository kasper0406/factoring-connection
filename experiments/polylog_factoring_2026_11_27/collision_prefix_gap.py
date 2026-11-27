"""Deterministic splits from coarse collision counts."""

from math import isqrt


def factor_from_gcd_prefix(n, estimate):
    """Factor n with 1/8-relative estimates of sum(gcd(d,n), 1<=d<=t)."""
    if n < 1:
        raise ValueError("Expected positive n")
    factors = []
    while n % 2 == 0:
        factors.append(2)
        n //= 2
    while n > 1:
        low, high = 2, isqrt(n) + 1
        sentinel = high
        while low < high:
            middle = (low + high) // 2
            if 5 * estimate(n, middle) > 6 * middle:
                high = middle
            else:
                low = middle + 1
        if low == sentinel:
            factors.append(n)
            break
        if n % low:
            raise ValueError("Prefix estimates failed the promised accuracy")
        factors.append(low)
        n //= low
    return factors


def split_balanced_quartic(n, interval_count):
    """Split p*q, p<q<2p, from 1/2-relative quartic difference interval counts."""
    from sympy import integer_nthroot
    if __package__:
        from .quartic_collision_coarse_split import interval_divisor
    else:
        from quartic_collision_coarse_split import interval_divisor
    root = isqrt(n)
    width = max(1, int(integer_nthroot(n, 6)[0]) // 8)
    low, high = root + 1, isqrt(2*n) + 1
    # Divide evenly so the last interval cannot be an arbitrarily short tail.
    block_size = max(1, root // 64)
    blocks = (high - low + block_size - 1) // block_size
    for i in range(blocks):
        start = low + i * (high - low) // blocks
        end = low + (i+1) * (high - low) // blocks
        if root < 64:
            if n % start == 0:
                return start
            continue
        if interval_count(n, start, end, 2) <= 16*(end-start):
            continue
        while end-start >= 2*width:
            middle = (start+end)//2
            if interval_count(n, start, middle, 2) > 16*(middle-start):
                end = middle
            else:
                start = middle
        divisor = interval_divisor(n, start, end, 4)
        if divisor is not None:
            return divisor
        raise ValueError("Positive interval yielded no verified divisor")
    raise ValueError("No balanced semiprime split found")


def least_factor_from_noisy_prefix(n, estimate, relative_error):
    """Return the least prime factor (n itself for prime n), for relative error < 1/3."""
    from fractions import Fraction

    error = Fraction(relative_error)
    if n < 2 or not 0 <= error < Fraction(1, 3):
        raise ValueError("Expected n >= 2 and a relative error below 1/3")
    slack = 2 - (1 + error) / (1 - error)
    limit = isqrt(n)
    start = min(limit, 1 + int(4 / slack))
    for divisor in range(2, start + 1):
        if n % divisor == 0:
            return divisor
    def detected(t):
        return estimate(n, t) > (1 + error) * t

    low = start
    while low < limit:
        high = min(limit, low + max(1, int(slack * low / 4)))
        if detected(high):
            while high - low > 1:
                middle = (low + high) // 2
                if detected(middle):
                    high = middle
                else:
                    low = middle
            return high
        low = high
    return n
