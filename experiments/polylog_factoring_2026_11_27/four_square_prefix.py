"""Deterministic relative approximations to prefixes of one scalar annihilator."""

from bisect import bisect_right
from collections.abc import Mapping
from fractions import Fraction
from math import ceil, gcd, isqrt


def _precision(epsilon):
    precision = Fraction(epsilon)
    if not 0 < precision <= 1:
        raise ValueError("epsilon must be positive and at most one")
    return precision


def _lower_step_table(maximum, value_at, delta):
    """Keep each new integer value exceeding (1+delta) times the last step."""
    starts = [0]
    values = [value_at(0)]
    scale = delta.denominator
    growth = scale + delta.numerator
    last_value = value_at(maximum)
    while last_value * scale > values[-1] * growth:
        low, high = starts[-1] + 1, maximum
        threshold = values[-1] * growth
        while low < high:
            middle = (low + high) // 2
            if value_at(middle) * scale > threshold:
                high = middle
            else:
                low = middle + 1
        starts.append(low)
        values.append(value_at(low))

    def lookup(b):
        return values[bisect_right(starts, b) - 1]

    return lookup


def approximate_factored_prefix(factors, B, epsilon, include_zero):
    """Lower estimate for the prefix with the odd prime factorization supplied.

    `factors` is a mapping or an iterable of (prime, exponent) pairs. The
    exact count is at most (1+epsilon) times the result, also at exact zero.
    Scalar representatives are integers from 0 or 1 through B, inclusive.
    """
    if B < 0:
        raise ValueError("the cutoff must be nonnegative")
    precision = _precision(epsilon)
    pairs = list(factors.items() if isinstance(factors, Mapping) else factors)
    delta = precision / (4 * max(1, len(pairs)))
    include = int(bool(include_zero))

    def previous(b):
        return b + include

    modulus = 1
    for prime, exponent in pairs:
        if prime < 3 or prime % 2 == 0 or exponent < 1 or gcd(modulus, prime) != 1:
            raise ValueError("expected distinct odd prime bases and positive exponents")
        # Multiplicative weights psi((prime^e)^2) = p^(2e) + p^(2e-1).
        weights = [(1, 1)]
        power = 1
        for _ in range(exponent):
            power *= prime
            weights.append((power, power * power + power * (power // prime)))
        modulus *= power

        def next_value(b):
            return sum(weight * previous(b // divisor) for divisor, weight in weights)

        previous = _lower_step_table(B, next_value, delta)
    if B >= modulus:
        raise ValueError("the cutoff must be smaller than the modulus")
    return 8 * previous(B)


def strip_small_prime_powers(n, limit):
    """Trial-divide an odd positive n. Return known powers, cofactor and trial count."""
    residual = n
    factors = []
    candidate = 3
    trials = 0
    while candidate <= limit and candidate * candidate <= residual:
        trials += 1
        exponent = 0
        while residual % candidate == 0:
            residual //= candidate
            exponent += 1
        if exponent:
            factors.append((candidate, exponent))
        candidate += 2
    if residual > 1 and candidate * candidate > residual:
        factors.append((residual, 1))
        residual = 1
    return factors, residual, trials


def approximate_closed_prefix(n, B, epsilon):
    """Lower estimate for the zero-inclusive prefix on an unfactored odd n.

    The exact count is at most (1+epsilon) times the result. Trial division
    removes small prime factors; the remaining cofactor stays intact.
    """
    if n < 1 or n % 2 == 0 or not 0 <= B < n:
        raise ValueError("expected an odd positive n and 0 <= B < n")
    precision = _precision(epsilon)
    limit = ceil(12 * n.bit_length() / precision)
    factors, residual, _ = strip_small_prime_powers(n, limit)

    # Linear term B times sum_{d|n} d eta(d), with each large-prime
    # multiplier rounded down to its leading prime power.
    mass = residual
    for prime, exponent in factors:
        power = prime ** exponent
        mass *= power + 2 * ((power - 1) // (prime - 1))

    # For the correction, enumerate small complementary divisors e=n/d.
    # The omitted mass is at most n² eta(n) times sum_{e>K} 1/e².
    correction = 0
    for complement in range(1, min(n, ceil(4 / precision)) + 1):
        if n % complement:
            continue
        divisor = n // complement
        weighted = divisor
        for prime, _ in factors:
            if divisor % prime == 0:
                weighted = weighted // prime * (prime + 1)
        correction += weighted * (divisor - B % divisor)
    return 8 * (B * mass + correction)


def factor_from_open_prefix(n, oracle):
    """Find the least prime factor with an oracle(n, B, epsilon) for 1 <= t <= B.

    For an odd composite n, estimates within ±25% of the positive-scalar
    count suffice. Up to the first prime factor the baseline is 8B. At and
    beyond that factor the additional terms raise the count to at least 24B.
    """
    if n <= 1 or n % 2 == 0:
        raise ValueError("expected an odd composite")
    low, high = 1, isqrt(n)
    while low < high:
        middle = (low + high) // 2
        if oracle(n, middle, Fraction(1, 4)) > 12 * middle:
            high = middle
        else:
            low = middle + 1
    if low <= 1 or n % low:
        raise ValueError("no divisor found under the composite and accuracy promises")
    return low
