"""Exhaustive enumeration of split polynomials prod (X - a), a in [-B, B], |multiset| = M."""
from __future__ import annotations

from typing import Iterator, Tuple


def split_polys(B: int, M: int) -> Iterator[Tuple[tuple, tuple]]:
    """Yield (roots, coeffs) for every multiset of size M from [-B, B]; roots nondecreasing,
    coeffs = (c_0, ..., c_{M-1}) of the monic polynomial (the leading 1 omitted)."""
    vals = list(range(-B, B + 1))

    def rec(start, depth, poly, roots):
        if depth == M:
            yield tuple(roots), tuple(poly[:-1])
            return
        for idx in range(start, len(vals)):
            a = vals[idx]
            # poly * (X - a)
            new = [0] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new[i + 1] += c
                new[i] -= a * c
            roots.append(a)
            yield from rec(idx, depth + 1, new, roots)
            roots.pop()

    yield from rec(0, 0, [1], [])


def group_by_top(entries, K: int):
    """Group coefficient tuples by their coefficients in degrees >= K."""
    groups: dict = {}
    for roots, co in entries:
        groups.setdefault(co[K:], []).append((roots, co[:K]))
    return groups
