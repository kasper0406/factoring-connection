"""Exact controls for the six returning ideal classes and for binary powering."""

from collections import Counter
from math import gcd

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .check_ternary_embedding import rotation_orbit
    from .packed_sphere_class_action import packed_class_power
    from .packed_ternary_decision import _packing_parameters
else:
    from check_mixed_ternary_embeddings import _all_frames
    from check_ternary_embedding import rotation_orbit
    from packed_sphere_class_action import packed_class_power
    from packed_ternary_decision import _packing_parameters

from src.nonsmooth.forms import compose, identity, reduce_form, reduced_forms


def returning_classes(n, x):
    scale, m = _packing_parameters(n, x)
    coordinate, informative = set(), set()
    informative_contents = set()
    for u, v in _all_frames(n, x):
        point = tuple(i - scale * j for i, j in zip(u, v))
        content = gcd(*v)
        if content == n:
            coordinate.add(point)
        else:
            informative.add(point)
            informative_contents.add(content)
    assert (len(coordinate), len(informative)) == (24, 48)

    # The factor is supplied only to build the control classes. Find r
    # independently from one of the directly enumerated noncoordinate frames.
    content, = informative_contents
    r = n // content
    beta = next(b for b in range(1, r) if (b * b + m) % r == 0)
    root2 = m % 2
    t2 = reduce_form((2, 2 * root2, (root2 * root2 + m) // 2))
    principal = identity(-4 * m)
    assert compose(t2, t2) == principal and t2 != principal
    classes = {principal, t2}
    for sign in (-1, 1):
        prime_class = reduce_form((r, 2 * sign * beta, (beta * beta + m) // r))
        translated = compose(prime_class, t2)
        assert (prime_class[0], translated[0]) == (r, 2 * r)
        classes.update((prime_class, translated))
    assert len(classes) == 6

    seen = Counter()
    for form in classes:
        reduced, point = packed_class_power(n, x, form)
        assert reduced == form
        orbit = frozenset(rotation_orbit(point))
        assert len(orbit) == 24
        if gcd(reduced[0], n) == 1:
            assert orbit == coordinate
        else:
            assert gcd(reduced[0], n) == r and orbit <= informative
        seen[orbit] += 1
    assert sorted(seen.values()) == [2, 2, 2]
    assert set().union(*seen) == coordinate | informative
    return classes, coordinate | informative


def check_full_class_group():
    n, x = 15, 0
    _, m = _packing_parameters(n, x)
    returns, targets = returning_classes(n, x)
    all_orbits = Counter()
    hits = set()
    for form in reduced_forms(-4 * m):
        reduced, point = packed_class_power(n, x, form)
        assert reduced == form
        all_orbits[min(rotation_orbit(point))] += 1
        if point in targets:
            hits.add(form)
    assert hits == returns and set(all_orbits.values()) == {2}
    return sum(all_orbits.values()), len(all_orbits), len(hits)


def check_large_power():
    n, x = 77, 1
    scale, m = _packing_parameters(n, x)
    assert scale == 7**4 and m == 34179874885
    exponent = 7**12 * (7**2 - 1) * (7**2 - 7)
    expected = (
        (11, 3, (166934, 36086, 206701), (384, 917, 925)),
        (19, 7, (5981, 5498, 5716006), (120, 684, 688)),
    )
    for ell, beta, expected_form, residues in expected:
        form = (ell, 2 * beta, (beta * beta + m) // ell)
        reduced, point = packed_class_power(n, x, form, exponent)
        assert reduced == expected_form
        assert sum(t * t for t in point) == m
        assert tuple(sorted(min(t % scale, -t % scale) for t in point)) == residues
        assert residues != (0, 1, 1)
    assert gcd(19, n * scale) == 1
    return exponent


if __name__ == "__main__":
    fixtures = ((15, 0), (21, 1), (35, 0), (65, 4), (77, 1), (91, 0), (143, 0))
    for n, x in fixtures:
        returning_classes(n, x)
    assert _packing_parameters(143, 0)[1] % 25 == 0
    print("Six returning classes and all 72 packed points on", len(fixtures), "inputs")
    print("Entire N=15 class group (classes, proper orbits, returning classes):",
          check_full_class_group())
    print("Failed uniform-return exponent for the norm-11 and norm-19 classes:",
          check_large_power())
