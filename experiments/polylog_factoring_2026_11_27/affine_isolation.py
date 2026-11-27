"""Nonadaptive mixed-root isolation using public weighted disk queries."""

from math import gcd, isqrt
from random import Random


def query_form(n, scale, offset, width):
    radius = width // 2
    shift = radius - offset
    leading = (n + 1)**4
    cross = 2*shift % n
    cross += n*((cross - 1) % (n+1))
    padding = n*(n+1)
    max_radius = n//8
    common_base = leading*(max_radius+1)**2
    cutoff = (2*common_base + padding*max_radius + 2*n - 1) // n
    base = n*cutoff - leading*radius**2 - padding*radius - n
    last = base + (shift*shift-scale*scale-base) % n
    return (leading, cross, last), cutoff


def decode(n, scale, offset, width, raw):
    radius = width // 2
    moment = raw * pow(2, -1, n)
    inverse = pow(scale, -1, n)
    if width % 2 == 0:
        endpoint = (width-offset)*inverse % n
        if (endpoint**2-1) % n == 0:
            moment -= width-radius
    for sign in (1, -1):
        position = (offset+sign*scale) % n
        if position < width:
            moment -= position-radius
    position = (moment+radius) % n
    if position >= width:
        return None
    root = (position-offset)*inverse % n
    if root in (1, n-1) or (root*root-1) % n:
        return None
    divisor = gcd(root-1, n)
    assert 1 < divisor < n
    return divisor


def make_batch(n, rng, repeats):
    widths = [1 << k for k in range(n.bit_length()) if (1 << k) <= n//4]
    batch = []
    for width in widths:
        for _ in range(repeats):
            while True:
                scale = rng.randrange(1, n)
                divisor = gcd(scale, n)
                if divisor == 1:
                    break
                if 1 < divisor < n:
                    return divisor, []
            batch.append((scale, rng.randrange(n), width))
    return None, batch


def find_factor(n, rng, disk_moment, repeats=None):
    if repeats is None:
        repeats = min(16 * (n.bit_length() + 1), 65_536)
    divisor, batch = make_batch(n, rng, repeats)
    if divisor is not None:
        return divisor
    for scale, offset, width in batch:
        form, cutoff = query_form(n, scale, offset, width)
        raw = disk_moment(n, form, cutoff)
        divisor = decode(n, scale, offset, width, raw)
        if divisor is not None:
            return divisor
    return None


def disk_reference(n, form, cutoff):
    a, b, c = form
    height = n*cutoff
    discriminant = 4*a*c-b*b
    max_y = isqrt(4*a*height//discriminant)
    result = 0
    for y in range(-max_y, max_y+1):
        bound = isqrt(4*a*height-discriminant*y*y)
        low = -((bound+b*y)//(2*a))
        high = (bound-b*y)//(2*a)
        for x in range(low, high+1):
            value = a*x*x+b*x*y+c*y*y
            if 0 < value <= height and value % n == 0:
                result += x*y
    return result % n


if __name__ == "__main__":
    rng = Random(261127)
    checked = isolated = 0
    for n in (15, 21, 45, 105, 315, 385, 945, 1155, 1575, 3465):
        roots = [r for r in range(1, n) if (r*r-1) % n == 0]
        found = 0
        cutoffs = set()
        for width in [1 << k for k in range(n.bit_length()) if (1 << k) <= n//4]:
            for _ in range(200):
                scale = rng.randrange(1, n)
                while gcd(scale, n) != 1:
                    scale = rng.randrange(1, n)
                offset = rng.randrange(n)
                form, cutoff = query_form(n, scale, offset, width)
                a, b, c = form
                cutoffs.add(cutoff)
                assert gcd(gcd(a, b), c) == 1 and 3*cutoff**2 < 4*a*c-b*b
                moment = disk_reference(n, form, cutoff)
                radius = width//2
                positions = [(scale*r+offset) % n for r in roots]
                expected = 2*sum(z-radius for z in positions if z <= 2*radius) % n
                assert moment == expected
                divisor = decode(n, scale, offset, width, moment)
                mixed = [(scale*r+offset) % n for r in roots if r not in (1, n-1)]
                if sum(z < width for z in mixed) == 1:
                    assert divisor is not None
                    found += 1
                    isolated += 1
                checked += 1
        assert found and len(cutoffs) == 1
    print(f"{checked} weighted disk queries, {isolated} isolated roots")
