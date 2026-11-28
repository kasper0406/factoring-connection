"""Check the binary denominator form against integer Euler matrices and sphere images."""

from collections import Counter
from math import gcd, isqrt

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .four_square_sampler import _multiply
    from .packed_point_rotations import (
        packed_point_transporter,
        rotate_integrally,
        small_prime_neighbors,
    )
    from .packed_ternary_decision import _packing_parameters
else:
    from check_mixed_ternary_embeddings import _all_frames
    from four_square_sampler import _multiply
    from packed_point_rotations import (
        packed_point_transporter,
        rotate_integrally,
        small_prime_neighbors,
    )
    from packed_ternary_decision import _packing_parameters


def matrix_denominator(h):
    a, b, c, d = h
    norm = a * a + b * b + c * c + d * d
    numerators = []
    for unit in ((0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)):
        numerators.extend(_multiply(_multiply(h, unit), (a, -b, -c, -d))[1:])
    return norm // gcd(norm, *numerators)


def check_frames():
    parity = Counter()
    exercised_alternative = points = rotations = 0
    for n, x in ((15, 0), (21, 1), (35, 0), (65, 4), (77, 1), (91, 0), (143, 0)):
        scale, m = _packing_parameters(n, x)
        for u, v in _all_frames(n, x):
            if gcd(*v) == n:
                continue
            points += 1
            w = tuple(a - scale * b for a, b in zip(u, v))
            form = packed_point_transporter(n, x, u, w)
            r, beta, c2 = form.prime, form.root, form.last_coefficient
            assert r == n // gcd(*v) and beta * beta + m == r * c2 and 2 * abs(beta) <= r <= c2
            multiplier = sum(t * t for t in form.first_quaternion) // r
            parity[multiplier] += 1
            for c in range(-3, 4):
                for b in range(4):
                    if gcd(c, b) != 1:
                        continue
                    h, den = form.at(c, b)
                    assert matrix_denominator(h) == den
                    assert rotate_integrally(h, form.source) == w
                    assert den == r if not b else 2 * r * den >= m * b * b
                    rotations += 1
            h, den = form.coprime_denominator_example()
            assert rotate_integrally(h, form.source) == w
            assert gcd(den, n) == 1 and 2 * r * den >= m and den * r < m + r * r
            exercised_alternative += int(c2 % r == 0)
    return points, rotations, dict(parity), exercised_alternative


def check_neighbors():
    n, x = 77, 1
    scale, m = _packing_parameters(n, x)
    source = (x + scale * n, 1, 0)
    informative = {}
    for u, v in _all_frames(n, x):
        if gcd(*v) != n:
            informative[tuple(i - scale * j for i, j in zip(u, v))] = gcd(*v)
    counts = []
    for prime in (3, 5, 7, 11, 13, 17, 19, 23):
        neighbors = small_prime_neighbors(source, prime)
        assert all(sum(t * t for t in w) == m for w in neighbors)
        matches = set(informative) & neighbors
        assert matches == (set(informative) if prime == 11 else set())
        assert all(informative[w] == 7 for w in matches)
        # Each root of X² = -M modulo the small prime gives a signed-cube
        # orbit of neighbors. The lower bound M/2 excludes point stabilizers
        # from identifying two such small-denominator rotations here.
        roots = sum((a * a + m) % prime == 0 for a in range(prime))
        assert len(neighbors) == 24 * roots
        counts.append((prime, len(neighbors), len(matches)))
    u, v = (1, 1, 0), (-14, -63, 42)
    w = tuple(i - scale * j for i, j in zip(u, v))
    form = packed_point_transporter(n, x, u, w)
    assert (form.prime, 2 * form.root, form.last_coefficient) == (11, 6, 3107261354)
    h, den = form.coprime_denominator_example()
    assert rotate_integrally(h, form.source) == w
    assert den == 1553630677 == (m + 21) // 22
    assert all(den % d for d in range(2, isqrt(den) + 1))
    return counts, den


if __name__ == "__main__":
    print("Packed points, transports, Euler norm multipliers, adjusted roots:", check_frames())
    print("Prime/neighbor/packed-hit counts and coprime denominator:", check_neighbors())
