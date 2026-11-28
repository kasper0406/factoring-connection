"""Independent checks of ternary plane incidences and coprime binary norm padding."""

from collections import Counter, defaultdict
from itertools import product
from math import gcd, isqrt

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .packed_plane_sections import _cross, _dot, homogeneous_plane_factor, plane_partners
    from .packed_ternary_decision import _packing_parameters
else:
    from check_mixed_ternary_embeddings import _all_frames
    from packed_plane_sections import _cross, _dot, homogeneous_plane_factor, plane_partners
    from packed_ternary_decision import _packing_parameters


def check_planes():
    counts = Counter()
    for n, x in ((15, 0), (15, 4), (33, 2), (65, 4)):
        q, _ = _packing_parameters(n, x)
        grouped = defaultdict(list)
        for u, v in _all_frames(n, x):
            grouped[u].append(v)
        for u, vs in grouped.items():
            assert plane_partners(n, x, u, u, _dot(u, u) + q*n*x) is None
            assert plane_partners(n, x, u, u, _dot(u, u) + q*n*x + 1) == ()
            points = {v: tuple(a - q*b for a, b in zip(u, v)) for v in vs}
            for v, point in points.items():
                if gcd(*v) == n:
                    continue
                normals = {_cross(point, t) for t in product((0, 1), repeat=3)}
                normals.add(_cross(u, v))
                for raw in normals:
                    content = gcd(*raw)
                    if not content:
                        continue
                    w = tuple(a // content for a in raw)
                    assert homogeneous_plane_factor(n, x, u, w) in (gcd(*v), n // gcd(*v))
                    assert gcd(n, _dot(u, w)) in (n, gcd(*v))
                    counts["squared normal" if gcd(n, _dot(u, w)) == n else "dot product"] += 1
                for w in ((1, 2, -3), (0, 4, 1), (1, 0, 0)):
                    for height in (_dot(w, point), _dot(w, point) + q):
                        observed = plane_partners(n, x, u, w, height)
                        expected = {v0 for v0, w0 in points.items() if _dot(w, w0) == height}
                        assert (set(vs) if observed is None else set(observed)) == expected
                        counts["affine"] += 1
    return counts


def norm_points(radius, d):
    points = set()
    for b in range(-radius, radius + 1):
        remaining = radius*radius - d*b*b
        if remaining >= 0 and isqrt(remaining)**2 == remaining:
            points.update((a, b) for a in {isqrt(remaining), -isqrt(remaining)})
    return points


def check_padding():
    total = 0
    for d, radius, pad in product((1, 2), (15, 33), (7, 13, 26)):
        factors = norm_points(radius, d), norm_points(pad, d)
        product_counts = Counter()
        for (a, b), (c, e) in product(*factors):
            z = (a*c - d*b*e, a*e + b*c)
            assert gcd(radius, a, b) == gcd(radius, *z)
            product_counts[z] += 1
        assert set(product_counts) == norm_points(radius*pad, d)
        assert set(product_counts.values()) == {4 if d == 1 else 2}
        total += sum(product_counts.values())
    return total


if __name__ == "__main__":
    print(check_planes(), "binary norm product pairs", check_padding())
