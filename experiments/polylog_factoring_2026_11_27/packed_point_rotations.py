"""Rational rotations between integral points on the packed ternary sphere."""

from collections import defaultdict
from dataclasses import dataclass
from math import gcd, isqrt

if __package__:
    from .four_square_sampler import _multiply
    from .packed_ternary_decision import _packing_parameters
else:
    from four_square_sampler import _multiply
    from packed_ternary_decision import _packing_parameters


def _odd_part(n):
    return n // (n & -n)


def _norm(v):
    return sum(c * c for c in v)


def _cross(u, v):
    a, b, c = u
    x, y, z = v
    return b * z - c * y, c * x - a * z, a * y - b * x


def _primitive_quaternion(h):
    divisor = gcd(*h)
    if not divisor:
        raise ValueError("the zero quaternion cannot define a rotation")
    if next(c for c in h if c) < 0:
        divisor = -divisor
    return tuple(c // divisor for c in h)


def rotate_integrally(h, w):
    """Return the rotated point if all three coordinates are integral."""
    a, b, c, d = h
    norm = _norm(h)
    if not norm:
        raise ValueError("the zero quaternion cannot define a rotation")
    numerator = _multiply(_multiply(h, (0, *w)), (a, -b, -c, -d))[1:]
    if any(t % norm for t in numerator):
        return None
    return tuple(t // norm for t in numerator)


def _frame_quaternion(columns, denominator):
    # For angle different from pi use (1+tr(R), R32−R23, R13−R31, R21−R12).
    first, second, third = columns
    real = denominator + first[0] + second[1] + third[2]
    if real:
        return _primitive_quaternion(
            (real, second[2] - third[1], third[0] - first[2], first[1] - second[0])
        )
    # At a half-turn every nonzero column of R+I is a vector on its axis.
    for i, column in enumerate(columns):
        axis = tuple(t + (denominator if i == j else 0) for j, t in enumerate(column))
        if any(axis):
            return _primitive_quaternion((0, *axis))
    raise ArithmeticError("invalid rotation columns")


@dataclass(frozen=True)
class PointTransporter:
    source: tuple
    target: tuple
    prime: int
    root: int
    last_coefficient: int
    first_quaternion: tuple
    second_quaternion: tuple

    def at(self, c, b):
        """Give the primitive Euler quaternion and reduced denominator at [c:b]."""
        if gcd(c, b) != 1:
            raise ValueError("expected coprime projective coordinates")
        h = _primitive_quaternion(
            tuple(c * i + b * j for i, j in zip(self.first_quaternion, self.second_quaternion))
        )
        value = self.prime * c * c + 2 * self.root * c * b + self.last_coefficient * b * b
        return h, _odd_part(value)

    def coprime_denominator_example(self):
        """Transport to the supplied point with denominator coprime to the input."""
        c = 0
        if self.last_coefficient % self.prime == 0:
            c = -1 if self.root > 0 else 1
        return self.at(c, 1)


def packed_point_transporter(n, x, u, w):
    """Classify all rational rotations to a supplied noncoordinate packed point.

    The source is (x+QN,1,0), and the target w is congruent to the supplied
    first column u mod Q. Constructing the transporter requires knowing w.
    """
    scale, m = _packing_parameters(n, x)
    if len(u) != 3 or len(w) != 3 or _norm(u) != x * x + 1 or _norm(w) != m:
        raise ValueError("unexpected first-vector or target norm")
    if any((i - j) % scale for i, j in zip(u, w)):
        raise ValueError("the target must lie in the supplied residue class")
    v = tuple((i - j) // scale for i, j in zip(u, w))
    if _norm(v) != n * n or sum(i * j for i, j in zip(u, v)) != -n * x:
        raise ValueError("target does not unpack to the mixed Gram")
    content = gcd(*v)
    r = n // content
    if r <= 1 or content <= 1 or n % content:
        raise ValueError("a noncoordinate target on the distinct-semiprime family is required")
    normal = _cross(u, v)
    columns = (
        tuple(-i // content for i in v),
        tuple((n * i + x * j) // content for i, j in zip(u, v)),
        tuple(i // content for i in normal),
    )
    h0 = _frame_quaternion(columns, r)
    if _norm(h0) not in (r, 2 * r, 4 * r):
        raise ArithmeticError("primitive frame quaternion has unexpected norm")
    source = (x + scale * n, 1, 0)
    product = _multiply(h0, (0, *source))
    index = next(i for i, t in enumerate(h0) if t % r)
    root = -product[index] * pow(h0[index], -1, r) % r
    if 2 * root > r:
        root -= r
    numerator = tuple(a * root + b for a, b in zip(h0, product))
    if any(t % r for t in numerator) or (root * root + m) % r:
        raise ArithmeticError("failed to saturate the quaternion transporter at the odd prime")
    h1 = tuple(t // r for t in numerator)
    return PointTransporter(source, tuple(w), r, root, (root * root + m) // r, h0, h1)


def small_prime_neighbors(point, prime):
    """Enumerate the integral images of one point under all rotations of denominator prime.

    This is an input-only operation for a small odd prime, using its three
    primitive Euler quaternion shells. Work grows polynomially in the prime.
    """
    if prime < 3 or prime % 2 == 0 or any(prime % i == 0 for i in range(3, isqrt(prime) + 1, 2)):
        raise ValueError("an odd prime denominator is required")
    pairs = defaultdict(list)
    for a in range(-isqrt(4 * prime), isqrt(4 * prime) + 1):
        bound = isqrt(4 * prime - a * a)
        for b in range(-bound, bound + 1):
            pairs[a * a + b * b].append((a, b))
    result = set()
    for norm in (prime, 2 * prime, 4 * prime):
        for length, left in pairs.items():
            for i, j in left:
                for k, l in pairs.get(norm - length, ()):
                    if gcd(i, j, k, l) != 1 or next(t for t in (i, j, k, l) if t) < 0:
                        continue
                    target = rotate_integrally((i, j, k, l), point)
                    if target is not None:
                        result.add(target)
    return result
