"""Batch a balanced-semiprime prime-bit query inside its exact public size interval."""

from math import gcd, isqrt
from typing import NamedTuple

if __package__:
    from .affine_prime_window import affine_divisors, window_radius
    from .sqrt_prime_marks import _DecimalRing
else:
    from affine_prime_window import affine_divisors, window_radius
    from sqrt_prime_marks import _DecimalRing
from src.trackC.polyback import product_tree_linear, multipoint_eval


class Cube(NamedTuple):
    base: int
    bits: tuple
    first: int
    stop: int

    @property
    def size(self):
        return self.stop - self.first


def _spread(rank, bits):
    """Scatter consecutive rank bits into distinct, increasing integer bit weights."""
    value = 0
    while rank:
        least = rank & -rank
        value |= bits[least.bit_length()-1]
        rank -= least
    return value


def _make_cube(base, free, lower, upper):
    bits = tuple(1 << j for j in range(1, free.bit_length()) if free & (1 << j))
    def lower_bound(value):
        lo, hi = 0, 1 << len(bits)
        while lo < hi:
            mid = (lo+hi)//2
            if base + _spread(mid, bits) < value:
                lo = mid+1
            else:
                hi = mid
        return lo
    return Cube(base, bits, lower_bound(lower), lower_bound(upper+1))


def _cube_factor(n, cube, stats):
    """Find a nonunit in this rank interval using one shared low-digit polynomial."""
    if cube.size == 0:
        return None
    low_count = (cube.size.bit_length()-1)//2
    stride = 1 << low_count
    first_full = (cube.first+stride-1)//stride
    stop_full = cube.stop//stride
    # At most two partial rows. A singleton cube has one full row of length one.
    for edge in (range(cube.first, min(cube.stop, first_full*stride)),
                 range(max(cube.first, first_full*stride, stop_full*stride), cube.stop)):
        for rank in edge:
            stats["edge_gcds"] += 1
            divisor = gcd(n, cube.base + _spread(rank, cube.bits))
            if divisor > 1:
                return divisor
    if first_full >= stop_full:
        return None

    offsets = [0]
    for bit in cube.bits[:low_count]:
        offsets += [x+bit for x in offsets]
    points = [cube.base + _spread(high, cube.bits[low_count:])
              for high in range(first_full, stop_full)]
    ring = _DecimalRing(n)
    polynomial = product_tree_linear(ring, offsets)[-1][0]
    values = multipoint_eval(ring, polynomial, points)
    stats["polynomial_degree"] = max(stats["polynomial_degree"], stride)
    stats["block_values"] += len(points)
    for point, value in zip(points, values):
        divisor = gcd(n, value)
        if 1 < divisor < n:
            return divisor
        if divisor == n:
            stats["saturated_rows"] += 1
            for offset in offsets:
                stats["refinement_gcds"] += 1
                divisor = gcd(n, point+offset)
                if 1 < divisor < n:
                    return divisor
            raise ArithmeticError("saturated row has no proper factor in the balanced interval")
    return None


def _cost(cubes):
    return sum(isqrt(cube.size-1)+1 for cube in cubes if cube.size)


def _query_bounds(n, mask):
    if n < 15 or n % 2 == 0:
        raise ValueError("expected an odd balanced semiprime")
    allowed = (1 << (n.bit_length()-1))-2
    if mask < 0 or mask & ~allowed:
        raise ValueError("expected dyadic positions from 1 through bitlength(n)-2")
    lower, upper = isqrt((n-1)//2)+1, isqrt(2*n)
    active = (1 << upper.bit_length())-2
    return lower, upper, active


def submask_zero(n, mask, *, stats=None):
    """Decide whole-N divisibility by the dyadic ratios selected by mask.

    Requires n=pq, distinct odd primes p<q<2p. No supplied factor or prime
    test is used. Report the factor in stats whenever a proper gcd was found.
    """
    lower, upper, active = _query_bounds(n, mask)
    selected = mask & active
    accepted = [_make_cube(1, active ^ selected, lower, upper)]
    # Partition the complement by its first selected bit equal to one.
    rejected = []
    available = active
    while selected:
        bit = selected & -selected
        available ^= bit
        rejected.append(_make_cube(1+bit, available, lower, upper))
        selected ^= bit

    search_accepted = _cost(accepted) <= _cost(rejected)
    cubes = accepted if search_accepted else rejected
    work = {} if stats is None else stats
    work.update(avoiding=accepted[0].size, hitting=sum(c.size for c in rejected),
                free=(active & ~mask).bit_count(), side="avoiding" if search_accepted else "hitting",
                polynomial_degree=0, block_values=0, edge_gcds=0, saturated_rows=0,
                refinement_gcds=0, divisor=None)
    for cube in cubes:
        divisor = _cube_factor(n, cube, work)
        if divisor is not None:
            work["divisor"] = divisor
            return bool(divisor & mask and (n//divisor) & mask)
    # If there is no avoiding prime then both hit. If there is no hitting
    # prime then neither hits and the whole-modulus test is false.
    return search_accepted


def submask_window_zero(n, mask, *, stats=None):
    """Search assignments outside the densest run of free bits of certified width."""
    lower, upper, active = _query_bounds(n, mask)
    free = active & ~mask
    work = {} if stats is None else stats
    work.update(free=free.bit_count(), grouped=0, outside=0, intervals=0,
                maximum_dimension=0, divisor=None)
    avoiding_count = _make_cube(1, free, lower, upper).size
    if avoiding_count == 0:
        return True
    if free == active:
        return False

    radius_limit = max(2, window_radius(n))
    # For each possible low weight find the largest admissible inner set.
    # Ties use the shorter coefficient interval to reduce LLL dimensions.
    remaining = free
    inner, stride, width = 0, 1, 0
    while remaining:
        low = remaining & -remaining
        span = (2*radius_limit).bit_length()
        possible = remaining & ((low << span)-1)
        if possible//low > 2*radius_limit:
            possible ^= 1 << (possible.bit_length()-1)
        current_width = possible//low
        if (possible.bit_count(), -current_width) > (inner.bit_count(), -width):
            inner, stride, width = possible, low, current_width
        remaining ^= low

    outer = _make_cube(1, free ^ inner, lower-width*stride, upper)
    inverse = pow(stride, -1, n)
    work.update(grouped=inner.bit_count(), outside=outer.size, width=width, stride=stride)
    for rank in range(outer.first, outer.stop):
        base = outer.base + _spread(rank, outer.bits)
        lo = max(0, -((base-lower)//stride))
        hi = min(width, (upper-base)//stride)
        if lo > hi:
            continue
        midpoint = (lo+hi)//2
        radius = hi-midpoint
        stats_window = {}
        work["intervals"] += 1
        candidates = affine_divisors(n, (base*inverse+midpoint) % n, radius, lower,
                                     stats=stats_window)
        work["maximum_dimension"] = max(work["maximum_dimension"], stats_window["dimension"])
        for divisor in candidates:
            if 1 < divisor < n:
                work["divisor"] = divisor
                return bool(divisor & mask and (n//divisor) & mask)
    return True
