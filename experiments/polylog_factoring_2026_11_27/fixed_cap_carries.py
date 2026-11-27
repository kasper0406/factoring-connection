"""Exact floor-square formulations of the fixed-cap absolute moment."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd, isqrt
import json


def odd_primes(limit):
    return [
        p
        for p in range(3, limit + 1, 2)
        if all(p % divisor for divisor in range(3, isqrt(p) + 1, 2))
    ]


def public_branch(n):
    c = isqrt(4 * n)
    delta = 4 * n - c * c
    if delta <= 4:
        c -= 1
        delta = 4 * n - c * c
    divisor = gcd(c, n)
    return c, delta, isqrt(delta // 3), divisor


def triangular_root(value):
    """The triangular number at floor(sqrt(value)), with zero below zero."""
    if value < 0:
        return 0
    root = isqrt(value)
    return root * (root + 1) // 2


def cumulative_y(height, delta):
    """Sum w over positive lattice points w^2+delta*y^2 <= height."""
    if height < delta + 1:
        return 0
    return sum(
        triangular_root(height - delta * y * y)
        for y in range(1, isqrt(height // delta) + 1)
    )


def cumulative_w(height, delta):
    """The same cumulative sum after exchanging the two summations."""
    if height < delta + 1:
        return 0
    return sum(
        w * isqrt((height - w * w) // delta)
        for w in range(1, isqrt(height - delta) + 1)
    )


def floor_difference_trace(n):
    """Compute the public fixed-cap trace by exact floor-square differences."""
    _, delta, cap, divisor = public_branch(n)
    if divisor > 1:
        return divisor, None, Counter()
    work = Counter()
    trace = 0
    for level in range(1, cap + 1):
        height = 4 * n * level
        maximum_y = isqrt((height - 1) // delta)
        work["levels"] += 1
        work["floor_square_terms"] += maximum_y
        for y in range(1, maximum_y + 1):
            value = height - delta * y * y
            root = isqrt(value)
            if root > isqrt(value - 1):
                trace += root
                work["boundary_hits"] += 1
    return gcd(trace, n), trace % n, work


def check():
    cumulative_checks = 0
    for delta in range(1, 48):
        for height in range(0, 240):
            first = cumulative_y(height, delta)
            assert first == cumulative_w(height, delta)
            shell = sum(
                w
                for y in range(1, isqrt(height // delta) + 1)
                for w in [isqrt(height - delta * y * y)]
                if w > 0 and w * w + delta * y * y == height
            )
            assert first - cumulative_y(height - 1, delta) == shell
            cumulative_checks += 1

    branches = hits = terms = 0
    for p, q in combinations(odd_primes(151), 2):
        if q > 2 * p:
            continue
        n = p * q
        _, delta, cap, divisor = public_branch(n)
        if divisor > 1:
            continue
        expected = sum(
            w
            for y in range(1, isqrt(4 * n * cap // delta) + 1)
            for w in range(1, isqrt(4 * n * cap - delta * y * y) + 1)
            if (w * w + delta * y * y) % (4 * n) == 0
        )
        _, trace, work = floor_difference_trace(n)
        assert trace == expected % n
        branches += 1
        hits += work["boundary_hits"]
        terms += work["floor_square_terms"]

    n = 38_942_974_393
    _, trace, work = floor_difference_trace(n)
    assert trace == 69_529_082
    return {
        "cumulative_reciprocity_checks": cumulative_checks,
        "small_balanced_branches": branches,
        "small_boundary_hits": hits,
        "small_floor_square_terms": terms,
        "scaling_control": {"N": n, "trace": trace, **work},
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int, nargs="?")
    args = parser.parse_args()
    if args.N is None:
        print(json.dumps(check(), indent=2))
    else:
        divisor, trace, work = floor_difference_trace(args.N)
        print(json.dumps({"N": args.N, "gcd": divisor, "trace": trace, **work}, indent=2))
