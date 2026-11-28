"""Compare every sampler fiber with direct enumeration of four-dimensional frames."""

from itertools import combinations
from math import gcd, isqrt
from random import Random

if __package__:
    from .check_four_square_incidence import quadruples
    from .check_ternary_sphere import primitive_points
    from .four_square_sampler import _multiply
    from .quaternion_gram import _lifts, embed_four_quaternion
else:
    from check_four_square_incidence import quadruples
    from check_ternary_sphere import primitive_points
    from four_square_sampler import _multiply
    from quaternion_gram import _lifts, embed_four_quaternion


def gram_of(u, v):
    return sum(x*x for x in u), sum(x*y for x, y in zip(u, v)), sum(y*y for y in v)


def imaginary(u, v):
    real, i, j, k = u
    return _multiply((real, -i, -j, -k), v)[1:]


def frame_content(u, v):
    return gcd(*(u[i]*v[j] - u[j]*v[i] for i, j in combinations(range(4), 2)))


def check_gram(gram):
    a, b, c = gram
    n = isqrt(a*c - b*b)
    expected = set()
    primitive = 0
    right = list(quadruples(c))
    for u in quadruples(a):
        for v in right:
            if sum(x*y for x, y in zip(u, v)) != b:
                continue
            prim = frame_content(u, v) == 1
            if gcd(*imaginary(u, v)) == 1:
                assert prim
                expected.add((u, v))
            primitive += prim
    assert gcd(a, b, c) != 1 or primitive == len(expected)
    lifted = set()
    for point in primitive_points(n):
        pairs = _lifts(a if a % 2 else c, (b, *point))
        assert len(set(pairs)) == 8
        for u, v in pairs:
            assert imaginary(u, v) == point
            if a % 2 == 0:
                u, v = v, u
            assert gram_of(u, v) == gram
            lifted.add((u, v))
    assert lifted == expected
    return len(lifted), primitive


if __name__ == "__main__":
    cases = [(1, 0, 1), (2, 1, 5), (3, 0, 3), (5, -2, 1), (5, 0, 5),
             (10, 5, 5), (9, 0, 25), (9, 3, 26), (17, 4, 1)]
    frames = [check_gram(gram) for gram in cases]
    assert frames[2][1] > frames[2][0]
    print(f"{len(cases)} exact Gram controls: {sum(x for x,_ in frames)} supported frames, "
          f"{sum(y for _,y in frames)} primitive frames. At 3I:", frames[2])
    rng = Random(11281937)
    calls = 0
    for gram in cases + [(1, 0, (1009*1013)**2), (1009, 0, 1009*1013**2)]:
        for _ in range(8):
            pair = embed_four_quaternion(gram, rng)
            assert gram_of(*pair) == gram
            assert frame_content(*pair) == gcd(*imaginary(*pair)) == 1
            calls += 1
    print(f"{calls} native sampler calls passed, including large square determinants.")
