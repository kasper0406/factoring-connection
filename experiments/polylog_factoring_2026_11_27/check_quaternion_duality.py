"""Check self-pairing a primitive four-square point with a twisted conjugate."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from random import Random

from sympy import Matrix

if __package__:
    from .level11_sampler import _inverse_splitting, _primitive_pair, norm as level11_norm
    from .four_square_sampler import _left_gcd, _multiply, norm_minus_one, primitive_four_squares
else:
    from level11_sampler import _inverse_splitting, _primitive_pair, norm as level11_norm
    from four_square_sampler import _left_gcd, _multiply, norm_minus_one, primitive_four_squares


def conjugate(x):
    return x[0], -x[1], -x[2], -x[3]


def norm(x):
    return sum(t * t for t in x)


def dual_gcd(n, point, twist):
    """Return the norm of the common right divisor of point and conjugate(point)*twist.

    This scalar formula applies to squarefree odd n and norm(point)=n,
    including point-dependent and singular twists.
    """
    return gcd(n, sum(x * v for x, v in zip(point, twist)))


def right_gcd_norm(point, twist):
    # Conjugation converts a common right divisor into a common left divisor.
    left = conjugate(point)
    right = _multiply(conjugate(twist), point)
    g = _left_gcd(tuple(2 * x for x in left), tuple(2 * x for x in right))
    return norm(g) // 4


def sample_dual(n, twist, rng, *, split=None):
    point = primitive_four_squares(n, rng, split=split)
    return point, dual_gcd(n, point, twist)


def points(n):
    pairs = defaultdict(list)
    limit = isqrt(n)
    for x in range(-limit, limit + 1):
        r = isqrt(n - x * x)
        for y in range(-r, r + 1):
            pairs[x * x + y * y].append((x, y))
    for squared, left in pairs.items():
        for x, y in left:
            for z, w in pairs.get(n - squared, ()):
                yield x, y, z, w


def unit_twist(n, rng):
    while True:
        v = tuple(rng.randrange(n) for _ in range(4))
        if gcd(n, norm(v)) == 1:
            return v


def projective(point, p):
    pivot = next(x for x in point if x % p)
    inverse = pow(pivot, -1, p)
    return tuple(x * inverse % p for x in point)


FOUR_SQUARE_TAIL = {(0, 2, 0, 0): 1, (0, 0, 2, 0): 1, (0, 0, 0, 2): 1}
LEVEL11_TAIL = {(1, 1, 0, 0): 1, (0, 2, 0, 0): 3, (0, 0, 2, 0): 1,
                (0, 0, 1, 1): 1, (0, 0, 0, 2): 3}


def quadric_remainder(n, terms, tail=FOUR_SQUARE_TAIL):
    """Reduce a polynomial modulo n and a monic quadratic x0² + tail."""
    remainder = {monomial: value % n for monomial, value in terms.items() if value % n}
    largest_power = max((m[0] for m in remainder), default=0)
    for power in range(largest_power, 1, -1):
        for m in [m for m in remainder if m[0] == power]:
            value = remainder.pop(m)
            for (a, b, c, d), coefficient in tail.items():
                target = (m[0] - 2 + a, m[1] + b, m[2] + c, m[3] + d)
                remainder[target] = (remainder.get(target, 0) - coefficient * value) % n
    return {m: value for m, value in remainder.items() if value}


def point_twist(point):
    return tuple(scale * x for scale, x in enumerate(point, 1))


def projective_points(p):
    for pivot in range(4):
        for tail in product(range(p), repeat=3 - pivot):
            yield (0,) * pivot + (1,) + tail


def check_homogeneous():
    # Q*x1 + c*x0*x2². Its restriction vanishes at exactly those prime
    # factors of 391 dividing c. Also test monic substitutions above degree 3.
    for tail in (FOUR_SQUARE_TAIL, LEVEL11_TAIL):
        q = {(2, 0, 0, 0): 1, **tail}
        q_times_x1 = {(a, b + 1, c, d): v for (a, b, c, d), v in q.items()}
        for c, expected in ((17, 17), (23, 23), (391, 391), (1, 1)):
            terms = {**q_times_x1, (1, 0, 2, 0): c}
            assert gcd(391, *quadric_remainder(391, terms, tail).values()) == expected
        shifted = {(a + 4, b, c, d): v for (a, b, c, d), v in q_times_x1.items()}
        assert quadric_remainder(391, shifted, tail) == {}

    for r in (5, 7, 13):
        cells = [(x, y, z, w) for x, y, z, w in projective_points(r)
                 if (x*x + x*y + 3*y*y + z*z + z*w + 3*w*w) % r == 0]
        assert len(cells) == (r + 1)**2
        # A plane section and a union of two plane sections.
        assert 0 < sum(x == 0 for x, y, z, w in cells) <= 2 * (r + 1)
        assert 0 < sum(x*y % r == 0 for x, y, z, w in cells) <= 4 * (r + 1)

    # A sample-dependent linear twist, whose norm can be singular even
    # when all four point coordinates are units.
    for point, expected_norm_gcd, expected_right_gcd in (
        ((-2, -19, -5, -1), 1, 23),
        ((-11, -10, -7, -11), 23, 17),
        ((-5, -14, -13, -1), 391, 1),
    ):
        twist = point_twist(point)
        assert norm(point) == 391 and all(gcd(x, 391) == 1 for x in point)
        assert gcd(391, norm(twist)) == expected_norm_gcd
        assert dual_gcd(391, point, twist) == right_gcd_norm(point, twist) == expected_right_gcd


def level11_modular_probe(n, inverse, blocks, evaluate, rng):
    """Test a circuit on independent modular rank-one matrices in norm coordinates."""
    inputs = []
    for _ in range(blocks):
        a, b = _primitive_pair(n, rng)
        c, d = _primitive_pair(n, rng)
        matrix_entries = Matrix([a*c, a*d, b*c, b*d])
        point = tuple(int(z) % n for z in inverse * matrix_entries)
        inputs.append(point)
    return gcd(n, evaluate(*inputs))


def check_joint_probe():
    # A planted identity modulo 17 must also be visible without obtaining any
    # integer points of norm 391. Modulo 23 its restriction is nonzero.
    def planted(x, y):
        dot = x[0]*y[0] + x[1]*y[1]
        return level11_norm(x)*y[0]**2 + 17*dot**2

    rng = Random(20261128)
    inverse = _inverse_splitting(391, rng)
    observed = Counter(level11_modular_probe(391, inverse, 2, planted, rng) for _ in range(128))
    assert set(observed) <= {17, 391} and observed[17] > 0
    for r in (3, 5, 7):
        # Determinant of the first columns of two rank-one matrices.
        # It is nonzero iff both columns are nonzero and the image lines differ.
        lines = [(1, i) for i in range(r)] + [(0, 1)]
        columns = [(a*c % r, b*c % r) for a, b in lines for c, _ in lines]
        nonzero = sum((a*d-b*c) % r != 0 for a, b in columns for c, d in columns)
        assert Fraction(nonzero, len(columns)**2) == Fraction(r, r+1)**3
    return observed[17]


def check_moments(p, q):
    # Each projective rank-one matrix is one kernel/image pair. Each
    # projective norm-unit quaternion acts on the local projective line.
    cells = Counter(projective(x, p) for x in points(p * q))
    size = sum(cells.values())
    collisions = Fraction(sum(c * c for c in cells.values()), size * size)
    line_count = p + 1
    hit_counts = []
    for twist in projective_points(p):
        if norm(twist) % p:
            hit_counts.append(sum(
                count for x, count in cells.items()
                if sum(a * b for a, b in zip(x, twist)) % p == 0
            ))
    assert len(hit_counts) == p * (p * p - 1)
    mean = Fraction(sum(hit_counts), size * len(hit_counts))
    second = Fraction(sum(h * h for h in hit_counts), size * size * len(hit_counts))
    assert mean == Fraction(1, line_count)
    assert second - mean * mean == (collisions - mean * mean) / (line_count - 1)


def run():
    check_homogeneous()
    # Prechosen from n before point selection. All four point coordinates
    # in each fixture are units, including the two proper factor witnesses.
    twist = unit_twist(391, Random('twists:391'))
    assert twist == (13, 9, 57, 122)
    witnesses = (
        ((-19, -5, -2, -1), 1), ((-19, 5, 2, 1), 17),
        ((-14, -13, -1, -5), 23), ((-5, -19, 2, 1), 391),
    )
    for point, expected in witnesses:
        assert norm(point) == 391 and all(gcd(x, 391) == 1 for x in point)
        assert right_gcd_norm(point, twist) == dual_gcd(391, point, twist) == expected
    for p, q in ((3, 5), (5, 7), (7, 11), (11, 13)):
        check_moments(p, q)

    # Fix each twist from n before drawing points. The factors only define
    # the benchmark family and verify the known number of representations.
    primes = (17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)
    means = [Fraction(0) for _ in range(4)]
    moduli = checks = 0
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            if q >= 2 * p:
                break
            n = p * q
            twists = ((1, 0, 0, 0), (1, 2, 3, 4), unit_twist(n, Random(f'twists:{n}')))
            hits = [0] * len(means)
            size = 0
            for point in points(n):
                size += 1
                for j, v in enumerate(twists):
                    if 1 < dual_gcd(n, point, v) < n:
                        hits[j] += 1
                if 1 < dual_gcd(n, point, point_twist(point)) < n:
                    hits[-1] += 1
            assert size == 8 * (p + 1) * (q + 1)
            for j, hit in enumerate(hits):
                means[j] += Fraction(hit, size)
            rng = Random(f'points:{n}')
            split = norm_minus_one(n, rng)
            for v in twists:
                for _ in range(8):
                    point, divisor = sample_dual(n, v, rng, split=split)
                    assert divisor == right_gcd_norm(point, v)
                    checks += 1
            for _ in range(8):
                point = primitive_four_squares(n, rng, split=split)
                v = point_twist(point)
                assert dual_gcd(n, point, v) == right_gcd_norm(point, v)
                checks += 1
            moduli += 1
    print(f'{moduli} moduli, {checks} sampled right-gcd cross-checks')
    for label, rate in zip(('coordinate', 'fixed dense', 'input-seeded', 'point-dependent'), means):
        print(f'{label}: mean exact coverage {float(rate / moduli):.6f}')


if __name__ == '__main__':
    print(f"{check_joint_probe()}/128 modular probes detected the planted two-block identity")
    run()
