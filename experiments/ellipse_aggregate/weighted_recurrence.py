"""A fixed-cap reduction using five weighted residues of one public disk."""

import json
from argparse import ArgumentParser
from collections import Counter
from math import isqrt

from factor import branch_from_n, residue_disk, split


def squares(k):
    return k * (k + 1) * (2 * k + 1) // 6


def fourths(k):
    return k * (k + 1) * (2 * k + 1) * (3 * k * k + 3 * k - 1) // 30


def raw_moments(br, work):
    """Baseline oracle. The join determines every positive match in the cap."""
    disk = residue_disk(br, br.cap, work)
    j1 = j2 = y1 = 0
    for i, j in enumerate(disk.indices[1:], 1):
        count = disk.counts[i] - disk.counts[i - 1]
        ysq = disk.moments[i] - disk.moments[i - 1]
        j1 += j * count
        j2 += j * j * count
        y1 += j * ysq
    return tuple(x % br.n for x in (disk.counts[-1], j1, j2, disk.moments[-1], y1))


def decode(br, raw):
    n = br.n
    k = isqrt(br.cap)
    principal = (
        k,
        squares(k),
        fourths(k),
        br.ell**2 * squares(k),
        br.ell**2 * fourths(k),
    )
    count, j1, j2, y0, y1 = ((r - t) % n for r, t in zip(raw, principal, strict=True))
    if not 1 <= count <= max(2, k):
        raise ArithmeticError("invalid mixed count")
    if count >= 3:
        unit = pow(squares(count), -1, n)
        a, weight = j1 * unit % n, y0 * unit % n
        route = "single_ray"
    elif count == 1:
        a, weight = j1, y0
        route = "single_point"
    else:
        gap2 = (2 * j2 - j1**2) % n
        gap = isqrt(gap2)
        if gap * gap != gap2 or (j1 - gap) % 2:
            raise ArithmeticError("invalid norm moments")
        a = (j1 - gap) // 2
        if gap:
            weight = ((a + gap) * y0 - y1) * pow(gap, -1, n) % n
            route = "distinct_norms"
        else:
            weight = y0
            route = "tied_minimum"
    if not 1 <= a <= br.cap:
        raise ArithmeticError("invalid decoded first shell")
    return a, weight, count, route


def fixed_cap_factor(n):
    br, d = branch_from_n(n)
    if d is not None:
        return {"N": n, "factors": sorted((d, n // d)), "route": "branch_gcd"}
    work = Counter()
    raw = raw_moments(br, work)
    a, weight, count, route = decode(br, raw)
    d = split(br, a, weight)
    if not 1 < d < n or n % d:
        raise ArithmeticError("failed divisor verification")
    return {
        "N": n,
        "factors": sorted((d, n // d)),
        "A": br.cap,
        "a": a,
        "S": weight,
        "mixed_points": count,
        "route": route,
        "raw": raw,
        "work": dict(work),
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    print(json.dumps(fixed_cap_factor(parser.parse_args().N), indent=2))
