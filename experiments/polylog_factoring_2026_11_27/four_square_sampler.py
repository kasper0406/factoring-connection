"""Sample signed, ordered primitive four-square representations of an odd norm."""

from itertools import product
from math import gcd, isqrt


def _multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )


def _left_gcd(dividend, divisor):
    """Greatest common left divisor in the Hurwitz order, in doubled coordinates."""
    while any(divisor):
        a, b, c, d = divisor
        denominator = a * a + b * b + c * c + d * d
        numerator = _multiply((a, -b, -c, -d), dividend)
        # The quotient divisor^{-1} dividend is numerator / denominator.
        # Compare the nearest integral and nearest all-half-integral points.
        candidates = []
        for parity in (0, 1):
            doubled = tuple(
                2 * ((2 * x + (1 - parity) * denominator) // (2 * denominator)) + parity
                for x in numerator
            )
            error = sum((q * denominator - 2 * x) ** 2 for q, x in zip(doubled, numerator))
            candidates.append((error, doubled))
        _, quotient = min(candidates)
        multiple = _multiply(divisor, quotient)
        remainder = tuple(x - y // 2 for x, y in zip(dividend, multiple))
        dividend, divisor = divisor, remainder
    return dividend


def gaussian_norm_pair(modulus, target, rng):
    """Find a²+b² = target modulo an odd modulus for any unit target."""
    if modulus < 1 or modulus % 2 == 0:
        raise ValueError("the modulus must be positive and odd")
    if modulus == 1:
        return (0, 0)
    inverse_target = pow(target, -1, modulus)
    while True:
        x, y = rng.randrange(modulus), rng.randrange(modulus)
        t = (x * x + y * y) * inverse_target % modulus
        if t == 1:
            return x, y
        if t % 4 != 1 or gcd(t, modulus) != 1:
            continue
        # A prime t = 1 mod 4 passes this root check with probability 1/2.
        # For every accepted t the exact sum-of-squares check suffices.
        root = pow(rng.randrange(1, t), (t - 1) // 4, t)
        if (root * root + 1) % t:
            continue
        previous, z = t, root
        while z * z > t:
            previous, z = z, previous % z
        w = isqrt(t - z * z)
        if z * z + w * w != t:
            continue
        inverse = pow(t, -1, modulus)
        return (x * z + y * w) * inverse % modulus, (y * z - x * w) * inverse % modulus


def norm_minus_one(modulus, rng):
    """Return (a,b) with a²+b² = −1 mod the odd modulus, without factoring."""
    return gaussian_norm_pair(modulus, -1, rng)


def _doubled_units():
    for index in range(4):
        for sign in (-2, 2):
            coordinates = [0] * 4
            coordinates[index] = sign
            yield tuple(coordinates)
    yield from product((-1, 1), repeat=4)


def primitive_four_squares(m, rng, *, split=None):
    """Draw exactly uniformly from primitive integer quadruples of squared norm m.

    The optional reusable split is any pair (a,b) with a²+b² = −1 mod m.
    Random choices use rng.randrange. For an odd squarefree m every quadruple
    is primitive. No factorization of m is used.
    """
    if m < 1 or m % 2 == 0:
        raise ValueError("the norm must be positive and odd")
    a, b = norm_minus_one(m, rng) if split is None else split
    if (a * a + b * b + 1) % m:
        raise ValueError("the split must have norm −1 modulo m")
    while True:
        x, y = rng.randrange(m), rng.randrange(m)
        if gcd(x, y, m) == 1:
            break
    # The quaternion below maps to [[2x, 0], [2y, 0]] modulo m in the
    # splitting i=[[0,1],[-1,0]], j=[[a,b],[b,-a]]. Its columns span the
    # uniform line represented by (x,y). A left gcd generates its right ideal.
    coordinates = (x, -y, -a * x - b * y, -b * x + a * y)
    lifted = tuple(2 * (n % m) for n in coordinates)
    generator = _left_gcd((2 * m, 0, 0, 0), lifted)
    integral = []
    for unit in _doubled_units():
        associate = _multiply(generator, unit)
        if all(n % 4 == 0 for n in associate):
            integral.append(tuple(n // 4 for n in associate))
    if len(integral) != 8:
        raise ArithmeticError("an odd-norm Hurwitz ideal must have eight integral generators")
    return integral[rng.randrange(len(integral))]
