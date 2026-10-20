"""Integer hull of the hyperbola region {xy >= N} (DIAGNOSTIC; exact integer arithmetic).

The integer hull I(H_N) = conv{(x, y) in Z^2 : x, y >= 1, xy >= N}.  Its bounded boundary is
the lower convex chain of the minimal points (x, ceil(N/x)), 1 <= x <= N.

Facts used in notes/analytic_results.md:
  * every divisor d of N gives a point (d, N/d) ON the hyperbola, hence an extreme point of
    the strictly convex region H_N and therefore a vertex of I(H_N) (proved, trivial);
  * #vertices = Theta(N^{1/3} log N) (Balog--Barany, arXiv 2602.06897; lower bound and
    O(N^{1/3} log N) upper bound: Alcantara--Blanco--Criado--Santos, arXiv 2501.19193).
So any algorithm that walks the hull vertex by vertex costs Omega(N^{1/3} log N).  Hiary's
blocks have the same order of count but are NOT hull edges (measured ratio -> ~1.84 per window).

`hull_vertices(N, X1, X2)` brute-forces the lower hull of the minimal points with x in a
padded window and returns the vertices with X1 <= x <= X2.  O(X2 - X1 + padding) time; it
is a measurement tool, not an algorithm (it touches every x).
"""
from __future__ import annotations

from typing import List, Tuple

import gmpy2


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def _lower_hull(points: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Andrew monotone chain, lower hull, points sorted by x (distinct x)."""
    h: List[Tuple[int, int]] = []
    for p in points:
        while len(h) >= 2:
            (x1, y1), (x2, y2) = h[-2], h[-1]
            # cross((p2 - p1), (p - p1)) <= 0  => p2 not strictly convex => pop
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                h.pop()
            else:
                break
        h.append(p)
    return h


def hull_vertices(N: int, X1: int, X2: int, pad: int | None = None) -> List[Tuple[int, int]]:
    """Vertices (x, y) of the integer hull of {xy >= N} with X1 <= x <= X2."""
    if pad is None:
        c = int(gmpy2.iroot(gmpy2.mpz(N), 3)[0]) or 1
        pad = 64 * X2 // c + 100
    lo = max(1, X1 - pad)
    hi = X2 + pad
    pts = [(x, ceil_div(N, x)) for x in range(lo, hi + 1)]
    h = _lower_hull(pts)
    return [(x, y) for (x, y) in h if X1 <= x <= X2]


def hull_vertex_count(N: int, X1: int, X2: int, pad: int | None = None) -> int:
    return len(hull_vertices(N, X1, X2, pad))
