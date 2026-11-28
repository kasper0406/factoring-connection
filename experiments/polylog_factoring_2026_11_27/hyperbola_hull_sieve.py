"""Stream (product, x, y) over a product interval and a fixed positive x interval."""

from dataclasses import dataclass
from heapq import heapify, heappop, heappush
from math import gcd


def _add(a, k, b):
    return a[0]+k*b[0], a[1]+k*b[1]


def _det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def _steps(n, hi, start, direction):
    """Inclusive integer interval of ray parameters inside xy>=n, x<=hi, y>=1."""
    x, y = start
    dx, dy = direction  # Both rays move rightward, downward, or both.
    if x > hi or y < 1:
        return None
    limits = []
    if dx:
        limits.append((hi-x)//dx)
    if dy:
        limits.append((y-1)//(-dy))
    end = min(limits)

    def product(k):
        return (x+k*dx)*(y+k*dy)

    if dx and dy:
        # Adjacent integers nearest the vertex maximize this concave quadratic.
        floor_peak = (x*dy+y*dx)//(-2*dx*dy)
        peak = min(end, max(0, floor_peak))
        if peak < end and product(peak+1) > product(peak):
            peak += 1
    else:
        peak = end if dx else 0
    if product(peak) < n:
        return None
    a, b = 0, peak
    while a < b:
        mid = (a+b)//2
        if product(mid) >= n:
            b = mid
        else:
            a = mid+1
    first = a
    a, b = peak, end
    while a < b:
        mid = (a+b+1)//2
        if product(mid) >= n:
            a = mid
        else:
            b = mid-1
    return first, a


def _next_vertex(n, hi, point):
    """ABCS Algorithm 2 on the hyperbola clipped by the wall x=hi."""
    inside, outside = (1, 0), (0, -1)
    while True:
        _, last = _steps(n, hi, _add(point, 1, inside), outside)
        inside = _add(inside, last, outside)
        steps = _steps(n, hi, _add(point, 1, outside), inside)
        if steps is None:
            _, last = _steps(n, hi, point, inside)
            return _add(point, last, inside)
        first, _ = steps
        outside = _add(outside, first-1, inside)


def initial_hull(n, lo, hi):
    """All extreme vertices of the lower integer hull, including both end columns."""
    if n < 1 or not 1 <= lo <= hi:
        raise ValueError("expected positive n and 1 <= lo <= hi")
    vertices = [(lo, (n+lo-1)//lo)]
    while vertices[-1][0] < hi:
        vertices.append(_next_vertex(n, hi, vertices[-1]))
    return vertices


def _primitive(v):
    g = gcd(*v)
    if not g:
        raise ValueError("expected a nonzero ray")
    return v[0]//g, v[1]//g


def corner_chain(u, v):
    """Extreme points of the punctured lattice cone from ray u to ray v."""
    u, v = _primitive(u), _primitive(v)
    det = _det(u, v)
    if det == 0:
        raise ValueError("expected linearly independent rays")
    if det < 0:
        return list(reversed(corner_chain(v, u)))
    a, b = u
    if abs(a) > 1:
        z = pow(-b, -1, abs(a))
        adjacent = (z, (1+b*z)//a)
    else:
        adjacent = (0, a) if a else (-b, 0)
    # det(u,adjacent)=1. Choose its representative with offset in [0,det).
    quotient, offset = divmod(_det(adjacent, v), det)
    current = _add(adjacent, -quotient, u)
    previous = u
    corners = [u]
    while offset:
        difference = det-offset
        if offset >= difference:
            # Skip a whole run of collinear Hilbert basis points.
            count = offset//difference
            step = _add(current, -1, previous)
            previous, current = _add(current, count-1, step), _add(current, count, step)
            offset -= count*difference
            det = offset+difference
        else:
            corners.append(current)
            multiple = (det+offset-1)//offset
            previous, current = current, _add((-previous[0], -previous[1]), multiple, current)
            det, offset = offset, multiple*offset-det
    corners.append(v)
    return corners


@dataclass(slots=True, eq=False)
class _Node:
    point: tuple
    left: object = None
    right: object = None


def divisor_events(start, stop, lo=1, hi=None, *, stats=None):
    """Yield every (n,x,y) with start<=n=x*y<=stop and lo<=x<=hi in (n,x) order.

    All inputs are positive integers. The default x interval includes both
    orientations of every divisor pair. Empty product levels require no loop.
    """
    hi = stop if hi is None else hi
    if not 1 <= start <= stop or not 1 <= lo <= hi:
        raise ValueError("expected 1 <= start <= stop and 1 <= lo <= hi")
    if stats is None:
        stats = {}
    stats.update(initial_vertices=0, events=0, inserted_vertices=0, max_frontier=0)
    if lo == hi:
        stats.update(initial_vertices=1, max_frontier=1)
        for y in range((start+lo-1)//lo, stop//lo+1):
            stats['events'] += 1
            yield lo*y, lo, y
        return
    frontier = []
    last = None
    for point in initial_hull(start, lo, hi):
        node = _Node(point, last)
        if last is not None:
            last.right = node
        frontier.append((point[0]*point[1], point[0], node))
        last = node
    heapify(frontier)
    stats.update(initial_vertices=len(frontier), max_frontier=len(frontier))
    while frontier[0][0] <= stop:
        product, x, node = heappop(frontier)
        point = node.point
        stats['events'] += 1
        yield product, x, point[1]
        left, right = node.left, node.right
        u = _add(left.point, -1, point) if left is not None else (0, 1)
        v = _add(right.point, -1, point) if right is not None else (0, 1)
        last = left
        for offset in corner_chain(u, v):
            new = _add(point, 1, offset)
            if (left is not None and new == left.point) or (right is not None and new == right.point):
                continue
            fresh = _Node(new, last)
            if last is not None:
                last.right = fresh
            last = fresh
            heappush(frontier, (new[0]*new[1], new[0], fresh))
            stats['inserted_vertices'] += 1
        if last is not None:
            last.right = right
        if right is not None:
            right.left = last
        node.left = node.right = None
        stats['max_frontier'] = max(stats['max_frontier'], len(frontier))
