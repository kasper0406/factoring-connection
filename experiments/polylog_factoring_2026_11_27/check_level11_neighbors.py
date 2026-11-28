"""Check the two-ideal repair against enumerated shells and exact CRT row laws."""

from collections import Counter
from fractions import Fraction
from math import gcd
from random import Random

from sympy import Matrix

if __package__:
    from .check_level11_marginals import eta_coefficients, image_line, norm_points
    from .check_level11_sampler import projective_lines
    from .level11_neighbors import norm_two_translate, repair_line
    from .level11_sampler import _coordinate_inverse, _inverse_splitting, norm
else:
    from check_level11_marginals import eta_coefficients, image_line, norm_points
    from check_level11_sampler import projective_lines
    from level11_neighbors import norm_two_translate, repair_line
    from level11_sampler import _coordinate_inverse, _inverse_splitting, norm


def _forward(n, inverse):
    return (inverse.adjugate() * pow(int(inverse.det()), -1, n)).applyfunc(lambda v: int(v) % n)


def _contains(n, line, matrix):
    u, v = line
    a, b, c, d = matrix
    return (u*c - v*a) % n == 0 and (u*d - v*b) % n == 0


def _check_repairs(n, inverse, lines, exhaustive):
    forward = _forward(n, inverse)
    mapping = {point: tuple(forward * Matrix(point)) for point in norm_points(n)
               if gcd(*point) == 1} if exhaustive else {}
    count = twists = 0
    for line in lines:
        point, twisted = repair_line(n, inverse, line)
        target = norm_two_translate(n, inverse, line) if twisted else line
        assert norm(point) == n and gcd(*point) == 1
        assert _contains(n, target, tuple(forward * Matrix(point)))
        if exhaustive:
            assert point in mapping
            assert twisted == (not any(_contains(n, line, v) for v in mapping.values()))
        count += 1
        twists += twisted
    return count, twists


def _check_rows(p, q, coefficients):
    n = p*q
    forward = _forward(n, _inverse_splitting(n, Random(n)))
    rows = forward.tolist()
    hits = {(image_line(v, rows, p), image_line(v, rows, q)) for v in norm_points(n)}
    action = tuple(forward * Matrix((1, 0, 1, 0)))

    def translate(slope, prime):
        x, y = (1, slope) if slope != prime else (0, 1)
        u = int(action[0]*x + action[1]*y) % prime
        v = int(action[2]*x + action[3]*y) % prime
        return v*pow(u, -1, prime) % prime if u else prime

    for first, second, flip in ((p, q, False), (q, p, True)):
        principal = {image_line(v, rows, first) for v in norm_points(first)}
        for x in range(first + 1):
            gx = translate(x, first)
            counts = Counter()
            for y in range(second + 1):
                gy = translate(y, second)
                original = (y, x) if flip else (x, y)
                shifted = (gy, gx) if flip else (gx, gy)
                counts[original in hits, shifted in hits] += 1
            pi = Fraction(3, 5)
            slope = Fraction(coefficients[second], second + 1)
            a = pi + slope * (int(x in principal) - pi)
            b = pi + slope * (int(gx in principal) - pi)
            expected = {(False, False): 0, (False, True): 1-a,
                        (True, False): 1-b, (True, True): a+b-1}
            assert all(counts[key] == (second + 1) * v for key, v in expected.items())
    return (p + 1)*(q + 1)


def main():
    rng = Random(114407)
    count = twists = 0
    for n in range(1, 44, 2):
        if gcd(n, 22) != 1:
            continue
        checks, repairs = _check_repairs(n, _inverse_splitting(n, rng), projective_lines(n), True)
        count += checks
        twists += repairs

    # Conjugate the supplied splitting to ensure the action is expressed in
    # the caller's coordinates, including when j has a different matrix.
    n = 15
    basis = _forward(n, _inverse_splitting(n, rng))
    change = Matrix([[1, 1], [0, 1]])
    omega, j = (change * basis.col(i).reshape(2, 2) * change.inv() for i in (1, 2))
    alternate = _coordinate_inverse(n, omega, j)
    extra, _ = _check_repairs(n, alternate, projective_lines(n), True)

    large = 0
    for n in (5**8, 131*137, (2**31-1)*(2**19-1), 2**127-1):
        inverse = _inverse_splitting(n, rng)
        lines = [(1, rng.randrange(n)) for _ in range(4)]
        checks, _ = _check_repairs(n, inverse, lines, False)
        large += checks
    coefficients = eta_coefficients(29)
    cells = sum(_check_rows(p, q, coefficients) for p, q in ((3, 5), (5, 13), (7, 19), (19, 29)))
    print(f'{count} small lines ({twists} repaired), {extra} in a conjugated splitting, '
          f'{large} larger inputs, {cells} CRT cells checked in both directions')


if __name__ == '__main__':
    main()
