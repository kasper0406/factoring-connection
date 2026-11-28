"""Test one public center at each dyadic scale of a possible large divisor."""

from random import SystemRandom

if __package__:
    from .affine_prime_window import affine_divisors
else:
    from affine_prime_window import affine_divisors


def divisor_scales(n):
    """Yield (lower divisor bound, symmetric radius) for large divisors of n."""
    if n < 4:
        raise ValueError("expected n >= 4")
    bits = (n - 1).bit_length()
    for k in range((n.bit_length() - 1) // 2, bits):
        exponent = k * (k-1) // (bits-1) - 2
        radius = 1 << exponent if exponent >= 0 else 0
        yield 1 << k, radius


def try_affine_scales(n, center=None, rng=None):
    """Return a verified proper divisor or None. Each center is one trial.

    With a uniform center and composite n the success probability is at least
    n^(-1/4)/256. The scan needs no balance or squarefreeness promise.
    """
    if n < 4:
        raise ValueError("expected n >= 4")
    if center is None:
        center = (SystemRandom() if rng is None else rng).randrange(n)
    center %= n
    distance_to_zero = min(center, n - center)
    for lower, radius in divisor_scales(n):
        # If the window contains a multiple of n, every other offset has
        # gcd at most 2*radius, strictly below this band's divisor bound.
        if distance_to_zero <= radius:
            continue
        for divisor in affine_divisors(n, center, radius, lower):
            if 1 < divisor < n:
                return divisor
    return None
