"""Cross-check oriented forms, embeddings and scalar roots, including prime powers."""

from collections import Counter, defaultdict
from itertools import permutations, product
from math import gcd
from pathlib import Path
from random import Random
import sys

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.nonsmooth.forms import reduce_form

if __package__:
    from .check_four_square_incidence import prime_powers
    from .check_ternary_sphere import primitive_points
    from .four_square_sampler import norm_minus_one, primitive_four_squares
    from .ternary_embedding import (point_from_residue, residue_from_point, residue_from_quaternion,
                                    residue_matching_norm_minus_one)
    from .ternary_sphere import hopf_point, perpendicular_form
else:
    from check_four_square_incidence import prime_powers
    from check_ternary_sphere import primitive_points
    from four_square_sampler import norm_minus_one, primitive_four_squares
    from ternary_embedding import (point_from_residue, residue_from_point, residue_from_quaternion,
                                   residue_matching_norm_minus_one)
    from ternary_sphere import hopf_point, perpendicular_form


def line_form(n, z):
    """Independent HNF and norm form on Z*z+n Z[i], with positive orientation."""
    h = hermite_normal_form(Matrix([[n, 0, z[0]], [0, n, z[1]]]))
    assert h.det() == n
    u, v = h[:, 0], h[:, 1]
    return reduce_form(tuple(map(int, (u.dot(u), 2*u.dot(v), v.dot(v)))))


def global_sign_matches(n, t, u):
    return t == u or tuple(-c % n for c in t) == u


def rotation_orbit(point):
    orbit = set()
    for axes in permutations(range(3)):
        parity = sum(axes[i] > axes[j] for i in range(3) for j in range(i + 1, 3))
        for signs in product((-1, 1), repeat=3):
            if signs[0]*signs[1]*signs[2] == (-1)**parity:
                orbit.add(tuple(sign*point[axis] for sign, axis in zip(signs, axes)))
    return orbit


def check_radius(n):
    residues = {(a, b) for a in range(n) for b in range(n) if (a*a + b*b + 1) % n == 0}
    forms = {t: line_form(n, t) for t in residues}
    for t in residues:
        point = point_from_residue(n, t)
        assert gcd(*point) == 1 and sum(c*c for c in point) == n*n
        assert global_sign_matches(n, t, residue_from_point(point))
        assert reduce_form(perpendicular_form(point)) == forms[t]
    counts = Counter()
    examples = {}
    orbits = defaultdict(set)
    for point in primitive_points(n):
        t = residue_from_point(point)
        form = reduce_form(perpendicular_form(point))
        assert t in residues and forms[t] == form
        counts[form] += 1
        examples[form] = point
        orbit = rotation_orbit(point)
        assert len(orbit) == (6 if n == 1 else 24)
        orbits[form].add(min(orbit))
    assert len(set(counts.values())) == 1
    assert counts == {f: 6*m for f, m in Counter(forms.values()).items()}
    if n > 1:
        prime_count = len(prime_powers(n))
        assert all(len(parts) == 2**(prime_count - 1) for parts in orbits.values())
    matched = 0
    if n <= 21:
        for a in range(n):
            for b in range(n):
                if gcd(a*a + b*b, n) != 1:
                    continue
                f = line_form(n, (a, b))
                inverse = reduce_form((f[0], -f[1], f[2]))
                for form, point in examples.items():
                    c = residue_matching_norm_minus_one((a, b), point)
                    assert (c is not None) == (form in (f, inverse))
                    if c is not None:
                        assert (c*c + a*a + b*b) % n == 0
                        matched += 1
    return len(residues), sum(counts.values()), matched


def check_seeded():
    rng = Random(618011)
    checked = 0
    for n in (45, 75, 65, 1009*1013, 3**23*5**17):
        split = norm_minus_one(n, rng)
        checked_emb = point_from_residue(n, split)
        assert global_sign_matches(n, split, residue_from_point(checked_emb))
        for _ in range(4):
            while True:
                row = primitive_four_squares(n, rng, split=split)
                point = hopf_point(row)
                if gcd(*point) == 1:
                    break
            t = residue_from_quaternion(row, n)
            assert global_sign_matches(n, t, residue_from_point(point))
            assert reduce_form(perpendicular_form(point)) == line_form(n, t)
            checked += 1
    return checked


if __name__ == "__main__":
    radii = (1, 3, 5, 7, 9, 13, 15, 21, 25, 27, 35, 45, 49, 65)
    checks = [check_radius(n) for n in radii]
    print("Residues, primitive points, matched scalar roots:", tuple(map(sum, zip(*checks))))
    n, z = 15, (2, 5)
    first = point_from_residue(n, z)
    second = point_from_residue(n, tuple(4*c % n for c in z))
    c, d = [residue_matching_norm_minus_one(z, point) for point in (first, second)]
    assert reduce_form(perpendicular_form(first)) == reduce_form(perpendicular_form(second))
    assert second not in rotation_orbit(first)
    assert 1 < gcd(c - d, n) < n
    print(f"Radius 15, points {first}, {second}, roots {c}, {d}, divisor {gcd(c-d,n)}")
    print("Seeded quaternion checks:", check_seeded())
