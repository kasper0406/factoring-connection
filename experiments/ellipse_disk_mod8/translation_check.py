"""Check the translation-parity obstruction on full periods and actual ellipse cutoffs."""

import argparse
import json
from collections import Counter
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path


def mixed_weight(u, v, n):
    """Integer completed mixed weight, computed from public residues."""
    if u * v % n:
        return 0
    common = gcd(gcd(u, v), n)
    multiplicity = 1 if common == 1 else (4 if common == n else 2)
    return multiplicity - (u % n == 0) - (v % n == 0)


def toggle(odd, key):
    if key in odd:
        odd.remove(key)
    else:
        odd.add(key)


def odd_cosets(offsets, p, q):
    left, right = set(), set()
    for u, v in offsets:
        toggle(left, (u % p, v % q))
        toggle(right, (u % q, v % p))
    return left, right


def direct_period(p, q):
    """Count every anchor independently with a bitset for each translate."""
    n = p * q
    universe = [(u, v) for u in range(n) for v in range(n)]
    origin = (0, 0)
    base = [(u, v) for u, v in universe if mixed_weight(u, v, n) % 2]
    cache = {}

    def table(offset):
        if offset not in cache:
            x, y = offset
            cache[offset] = sum(1 << (((u + x) % n) * n + (v + y) % n) for u, v in base)
        return cache[offset]

    checks = Counter()

    def check(offsets, family):
        state = 0
        for point in offsets:
            state ^= table(point)
        left, right = odd_cosets(offsets, p, q)
        a, b = len(left), len(right)
        assert state.bit_count() == n * (a + b) - 2 * a * b
        if offsets and len(offsets) < n and state == 0:
            u0, v0 = offsets[0]
            assert any(
                1 < gcd(coord, n) < n
                for u, v in offsets[1:]
                for coord in (u - u0, v - v0)
            )
        checks[family] += 1
        return state

    check([], "empty")
    check([origin], "singleton")
    primitive = [(u, v) for u, v in base if gcd(gcd(u, v), n) == 1]
    phi = (p - 1) * (q - 1)
    assert len(primitive) == 2 * phi
    for shear in range(n):
        surviving = sum(u * (v - shear * u) % n == 0 for u, v in primitive)
        assert surviving == phi * (int(shear % p == 0) + int(shear % q == 0))
        checks["all_shear_parameters"] += 1
    for shift in universe[1:]:
        assert check([origin, shift], "all_nonzero_reflection_centers") != 0
    for step in (1, 2, n - 1):
        for length in range(1, n + 1):
            offsets = [(i, (step * i) % n) for i in range(length)]
            assert check(offsets, "progression").bit_count() == 2 * length * (
                n - length
            )
    # These directions use the factors. They certify sharpness of the obstruction.
    h, k = (p, q), (q, p)
    rectangle = [origin, h, k, ((p + q) % n, (p + q) % n)]
    assert check(rectangle, "factored_parallelogram") == 0
    for u, v in universe:
        weights = [mixed_weight(u - x, v - y, n) for x, y in rectangle]
        assert weights[0] - weights[1] - weights[2] + weights[3] == 0
        assert (4 * sum(weights)) % 8 == 0
    checks["integer_signed_parallelogram"] += 1
    for subset in combinations(rectangle, 3):
        assert check(subset, "proper_triangle") != 0
    # All six pairs, four triples, and the full set of four unit-axis offsets.
    square = [origin, (1, 0), (0, 1), (1, 1)]
    for length in (2, 3, 4):
        for subset in combinations(square, length):
            check(subset, "unit_square_subsets")
    return {"N": n, "checks": dict(checks), "anchors_per_check": n * n}


def audit_disk(row, cutoff, mask):
    """Physical query at anchor zero versus factor-backed coset accounting."""
    n, ell, c, delta = (row[k] for k in ("N", "ell", "c", "delta"))
    p, q = row["factors"]
    assert p * q == n
    bound = 4 * n * cutoff
    odd_left, odd_right = set(), set()
    candidates = direct = 0
    for y in range(1, isqrt(bound // delta) + 1):
        if y & mask != mask:
            continue
        # w == delta*y mod 2 includes the physical lattice condition.
        start = 1 if (delta * y) % 2 else 2
        for w in range(start, isqrt(bound - delta * y * y) + 1, 2):
            candidates += 1
            u, v = ell * w - c * y, ell * w + c * y
            if (w * w + delta * y * y) % (4 * n) == 0 and u % n and v % n:
                direct += 1
            toggle(odd_left, (u % p, v % q))
            toggle(odd_right, (u % q, v % p))
    assert (direct % 2) == (((0, 0) in odd_left) ^ ((0, 0) in odd_right))
    a, b = len(odd_left), len(odd_right)
    no_collision = a == b == candidates
    if 2 * cutoff <= row["a"]:
        assert no_collision
    return candidates, direct, no_collision


def triangle_control(row):
    """Public triangle count versus the factor-backed threshold at every cutoff."""
    n = row["N"]
    q = max(row["factors"])
    lo, hi = isqrt(n) + 1, isqrt(2 * n)
    hits = 0
    for v in range(1, hi + 1):
        hits += sum(u * v % n == 0 for u in range(1, v))
        if v >= lo:
            assert hits == int(v >= q)
    return hi - lo + 1


def dyadic_control(n):
    """Exhibit a physical point in every child cell at one unmasked cutoff."""
    c = isqrt(4 * n)
    if 4 * n - c * c <= 4:
        c -= 1
    delta = 4 * n - c * c
    side = isqrt(4 * n // (delta + 9))
    stride = 1 << (side.bit_length() - 1)
    anchors = set()
    scale = pow(stride, -1, n)
    for x in range(1, stride + 1):
        for y in range(1, stride + 1):
            w = 2 * x + (delta % 2) * y
            assert 0 < w * w + delta * y * y <= 4 * n
            # The standard residues of (x,y) give this child's translation.
            ry = y % stride
            rw = 2 * (x % stride) + (delta % 2) * ry
            anchors.add(((rw - c * ry) * scale % n, (rw + c * ry) * scale % n))
    assert len(anchors) == stride * stride
    return {
        "N": n,
        "delta": delta,
        "B": 1,
        "stride": stride,
        "distinct_nonempty_children": len(anchors),
    }


def run():
    source = json.loads(Path(__file__).with_name("results.json").read_text())
    groups = {}
    collision_pair = []
    for key in ("records", "tie_controls"):
        counts = Counter()
        for row in source[key]:
            counts["triangle_threshold_cutoffs"] += triangle_control(row)
            first = row["a"]
            for cutoff in sorted({first // 2, first - 1, first}):
                for mask in (0, 1, 3):
                    size, direct, injective = audit_disk(row, cutoff, mask)
                    counts["queries"] += 1
                    counts["candidate_coordinates"] += size
                    counts["parity_one_at_origin"] += direct % 2
                    if size:
                        counts["nonempty_stencils"] += 1
                        counts["nonempty_injective_in_both_quotients"] += injective
                        if 2 * cutoff <= first:
                            counts["nonempty_through_half_minimum"] += 1
                    if cutoff == first and mask == 0:
                        assert direct == row["mixed_orbits"]
                        counts["first_shell_singletons"] += direct == 1
                        counts["first_shell_ties"] += direct == 2
                        if key == "records" and row["N"] in (11009, 28829):
                            collision_pair.append(
                                {
                                    "N": row["N"],
                                    "B": cutoff,
                                    "candidates": size,
                                    "mixed": direct,
                                    "injective": injective,
                                }
                            )
        groups[key] = dict(branches=len(source[key]), **counts)
    return {
        "full_period": [direct_period(p, q) for p, q in ((3, 5), (3, 7), (5, 7))],
        "disks": groups,
        "dyadic_descent": [
            dyadic_control(n) for n in (11009, 28829, 438751, 6494011, 38942974393)
        ],
        "historical_pair": collision_pair,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report, end="")
