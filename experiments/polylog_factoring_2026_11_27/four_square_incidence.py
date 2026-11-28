"""Sample integer four-square points with modular scalar annihilators."""

from math import gcd
from random import SystemRandom

if __package__:
    from .four_square_sampler import primitive_four_squares
else:
    from four_square_sampler import primitive_four_squares


def _other_prime_powers(n, m):
    """Return the full prime powers of n supported outside m, multiplied together."""
    return n // gcd(n, pow(m, n.bit_length(), n))


def _gcd_of_residues(modulus, count, rng):
    result = modulus
    for _ in range(count):
        result = gcd(result, rng.randrange(modulus))
        if result == 1:
            break
    return result


def sample_content(n, count, rng):
    """Sample g with weight g^count times the number of points with content g."""
    if n < 1 or n % 2 == 0:
        raise ValueError("n must be positive and odd")
    if count < 0 or count == 2:
        raise ValueError("the factor-free sampler covers counts 0, 1, and 3 or more")
    h = abs(count - 2)
    while True:
        v = _gcd_of_residues(n, h, rng)
        g = v if count < 2 else n // v
        m = n // g
        # Correct the prime factors missing from the gcd proposal's weight.
        if _gcd_of_residues(_other_prime_powers(n, n // v), h, rng) != 1:
            continue
        # A uniform projective line over each remaining prime r has r+1
        # directions. Requiring a unit first coordinate accepts r of them.
        other = _other_prime_powers(n, m)
        while True:
            x, y = rng.randrange(other), rng.randrange(other)
            if gcd(x, y, other) == 1:
                break
        if gcd(x, other) == 1:
            return g


def _point_and_scalars(n, g, count, rng):
    m = n // g
    point = tuple(g * value for value in primitive_four_squares(m * m, rng))
    scalars = tuple(m * rng.randrange(g) for _ in range(count))
    return point, scalars


def sample_incidence(n, count, rng=None):
    """Uniform (point, scalars) with sum(point[i]²)=n² and each t*point=0 mod n."""
    if rng is None:
        rng = SystemRandom()
    return _point_and_scalars(n, sample_content(n, count, rng), count, rng)


def sample_with_fixed_scalar(n, fixed, rng=None):
    """Uniform (point, t) conditional on both fixed*point=t*point=0 mod n."""
    if n < 1 or n % 2 == 0:
        raise ValueError("n must be positive and odd")
    if rng is None:
        rng = SystemRandom()
    modulus = gcd(n, fixed)
    point, (residue,) = sample_incidence(modulus, 1, rng)
    return tuple(n // modulus * x for x in point), residue + modulus * rng.randrange(n // modulus)


def sample_two_from_factors(prime_powers, rng=None):
    """Two-scalar reference sampler, given distinct odd primes and their exponents."""
    if rng is None:
        rng = SystemRandom()
    n = g = 1
    for prime, exponent in prime_powers:
        if prime < 3 or prime % 2 == 0 or exponent < 1 or gcd(prime, n) != 1:
            raise ValueError("expected distinct odd prime bases and positive exponents")
        # Exponents 0,...,a-1 have weight p+1 each. Exponent a has weight p.
        ticket = rng.randrange(exponent * (prime + 1) + prime)
        g *= prime ** (ticket // (prime + 1))
        n *= prime ** exponent
    return _point_and_scalars(n, g, 2, rng)
