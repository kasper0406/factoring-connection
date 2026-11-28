"""Cubic Jacobi information on the Eisenstein norm-one torus, using the public modulus."""

from math import gcd
from random import SystemRandom


OMEGA = (0, 1)
ONE = (1, 0)


def norm(z):
    a, b = z
    return a*a-a*b+b*b


def _product(x, y):
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c-b*d


def torus_product(n, x, y):
    a, b = _product(x, y)
    return a % n, b % n


def torus_power(n, z, exponent):
    """Nonnegative integer powering in (Z/n)[omega]."""
    if n < 2 or exponent < 0:
        raise ValueError("invalid ring modulus or exponent")
    result = ONE
    while exponent:
        if exponent & 1:
            result = torus_product(n, result, z)
        z = torus_product(n, z, z)
        exponent //= 2
    return result


def _remainder(z, den):
    c, d = den
    den_norm = norm(den)
    v, w = _product(z, (c-d, -d))
    # Round each coefficient of z/den to a nearest integer.
    quotient = ((2*v+den_norm)//(2*den_norm), (2*w+den_norm)//(2*den_norm))
    e, f = _product(quotient, den)
    return z[0]-e, z[1]-f


def cubic_jacobi(n, z):
    """Return log_omega of the cubic Jacobi symbol (z/(n)), or None for a nonunit."""
    if n <= 1 or gcd(n, 3) != 1:
        raise ValueError("expected a positive rational denominator >1 coprime to 3")
    z = z[0] % n, z[1] % n
    if gcd(n, norm(z)) != 1:
        return None
    den = (n, 0)
    answer = 0
    while norm(den) > 1:
        rem = _remainder(z, den)
        c, d = den
        ramified = 0
        a, b = rem
        while (a+b) % 3 == 0:
            a, b = (2*a-b)//3, (a+b)//3
            ramified += 1
        rotations = 0
        while b % 3:
            a, b = -b, a-b
            rotations += 1
        # Supplements for 1-omega and omega at the primary denominator c+d*omega.
        answer -= ramified*((c*c-1)//3) + rotations*((c*c-c*d-1)//3)
        z, den = den, (a, b)
    return answer % 3


def torus_sample(n, rng=None):
    """Return (uniform norm-one element, None), or (None, an incidental proper factor)."""
    if n < 5 or gcd(n, 6) != 1:
        raise ValueError("expected an odd modulus >=5 coprime to 3")
    generator = SystemRandom() if rng is None else rng
    while True:
        x = generator.randrange(n), generator.randrange(n)
        value = norm(x)
        common = gcd(value, n)
        if common == n:
            continue
        if common > 1:
            return None, common
        # Hilbert 90: x/conj(x) = x²/Norm(x).
        scale = pow(value, -1, n)
        a, b = _product(x, x)
        return (scale*a % n, scale*b % n), None


def rotation_labels(n, z):
    """Cubic Jacobi exponents for z, omega*z, omega²*z."""
    label = cubic_jacobi(n, z)
    if label is None:
        raise ValueError("rotations require a unit")
    step = ((n*n-1)//3) % 3
    return tuple((label+i*step) % 3 for i in range(3))


def mixed_probe(n, rng=None):
    """For n=pq=-1 mod 3, certify different local endpoint labels or give a proper factor."""
    if n < 5 or n % 6 != 5:
        raise ValueError("expected an odd semiprime congruent to -1 modulo 3")
    if n % 9 != 8:
        return OMEGA, None
    generator = SystemRandom() if rng is None else rng
    while True:
        z, divisor = torus_sample(n, generator)
        if divisor is not None or cubic_jacobi(n, z) != 0:
            return z, divisor


def torus_mark(r, z):
    """Verification with a known prime r: log_omega of z^((r-epsilon_r)/3)."""
    if r < 2 or r == 3 or norm(z) % r != 1:
        raise ValueError("expected a norm-one element at a prime different from 3")
    epsilon = 1 if r % 3 == 1 else -1
    endpoint = torus_power(r, z, (r-epsilon)//3)
    return (ONE, OMEGA, (r-1, r-1)).index(endpoint)
