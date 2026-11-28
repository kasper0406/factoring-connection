"""Published pairings, independent chord crossings, and relations on X_0(n)."""

from math import gcd
from random import Random

if __package__:
    from .growing_level_intersection import level_intersection, unit_chord_reduction
else:
    from growing_level_intersection import level_intersection, unit_chord_reduction


def edge_loop(n, a):
    inverse = pow(a, -1, n)
    return inverse, -1, 1 - inverse*a, a


def product(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e + b*g, a*f + b*h, c*e + d*g, c*f + d*h


def chord(n, a, b):
    # The two rays of the first edge meet the cusp-0 circle at a and -1/a.
    # Crossing the oriented chord is an exit or entry into the open arc.
    first = -pow(a, -1, n) % n
    second = -pow(b, -1, n) % n
    if len({a, first, b, second}) < 4:
        return 0
    arc_width = (a - first) % n
    return int((second-first) % n < arc_width) - int((b-first) % n < arc_width)


def check():
    g, h = [edge_loop(11, a) for a in (8, 9)]
    assert level_intersection(11, g, h) == 1  # Stein, Example 3.9: integral basis

    # The handbook gives the first scalar as -1. Its printed matrix lists
    # the first function argument down each column (Magma iteration order).
    basis = [edge_loop(37, a) for a in (29, 22, 12, 18)]
    expected = [[0, -1, 0, -1], [1, 0, -1, -1], [0, 1, 0, 0], [1, 1, 0, 0]]
    assert [[level_intersection(37, g, h) for h in basis] for g in basis] == expected

    chord_cases = 0
    for n in range(2, 48):
        units = [a for a in range(n) if gcd(a, n) == 1]
        for a in units:
            for b in units:
                assert level_intersection(n, edge_loop(n, a), edge_loop(n, b)) == chord(n, a, b)
                chord_cases += 1

    # Build curves with c=N, c=2N, etc., independently of cusp-0 chords.
    rng = Random(737)
    relation_cases = 0
    reductions = 0
    splits = 0
    for n in (1, 2, 4, 7, 11, 15, 25, 37, 44, 63, 91, 551):
        loops = []
        for factor in range(1, 21):
            c = factor*n
            a = rng.randrange(1, 4*c+1)
            while gcd(a, c) != 1:
                a += 1
            d = pow(a, -1, c)
            loops.append((a, (a*d-1)//c, c, d))
        for g, h, k in zip(loops, loops[1:], loops[2:]):
            value = level_intersection(n, g, h)
            assert value == -level_intersection(n, h, g)
            assert level_intersection(n, product(g, k), h) == value + level_intersection(n, k, h)
            assert level_intersection(n, g, product(k, h)) == value + level_intersection(n, g, k)
            if n in (1, 2, 4, 7, 25):
                assert value == 0  # genus zero
            divisor, indices = unit_chord_reduction(n, g)
            if divisor is not None:
                assert 1 < divisor < n and n % divisor == 0
                splits += 1
            else:
                assert value == sum(level_intersection(n, edge_loop(n, a), h) for a in indices)
                reductions += 1
            relation_cases += 1

    # Stabilizers of orders two and three in PSL_2 at level 37.
    for elliptic in ((6, -1, 37, -6), (11, -3, 37, -10)):
        for h in basis:
            assert level_intersection(37, elliptic, h) == level_intersection(37, h, elliptic) == 0

    # Insert a huge full cusp walk between two nonzero test classes.
    power = (1 << 4096) + 1
    for parabolic in ((1, -power, 0, 1), (1, 0, 37*power, 1)):
        assert level_intersection(37, basis[0], product(parabolic, basis[1])) == -1
        assert level_intersection(37, product(basis[0], parabolic), basis[1]) == -1
    n = ((1 << 521)-1) * ((1 << 127)-1)
    total_large = 0
    for a in (8, 13, 17, 29, 93):
        for b in (23, 31, 41, 67, 109):
            value = level_intersection(n, edge_loop(n, a), edge_loop(n, b))
            assert value == chord(n, a, b)
            total_large += abs(value)
    assert total_large > 0
    for p, q in ((7, 13), (11, 13), (19, 29)):
        n = p*q
        indices = (2, 4, 5)
        g = product(product(edge_loop(n, 2), edge_loop(n, 4)), edge_loop(n, 5))
        row_sums = [0]*p
        for b in range(n):
            if gcd(b, n) == 1:
                value = level_intersection(n, g, edge_loop(n, b))
                assert value == sum(chord(n, a, b) for a in indices)
                row_sums[b % p] += value
        assert sum(row_sums) == 0 and max(map(abs, row_sums)) <= 2*len(indices)
    assert reductions and splits
    return chord_cases, relation_cases, reductions, splits, total_large


if __name__ == "__main__":
    chords, relations, reductions, splits, large_nonzero = check()
    print(f"Matched the level-11 basis and the level-37 intersection matrix, {chords} unit chord pairs, "
          f"{relations} group-law checks ({reductions} chord reductions and {splits} proper divisors), "
          f"and 25 large-modulus chord pairs ({large_nonzero} nonzero).")
