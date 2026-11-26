"""Controls for the native centered count, adaptive descent, and local-solubility test."""

import json
from argparse import ArgumentParser
from collections import Counter
from functools import cache
from math import gcd, isqrt
from pathlib import Path

from centered import (
    basis,
    guided_factor,
    native_oracle,
    norm,
    principal,
    radial_count,
    reduce_form,
    restrict,
)


def direct(n, f, cutoff):
    """Independent integer row intervals. Retain exact signed counts including axes."""
    a, b, c = f
    height = n * cutoff
    delta = 4 * a * c - b * b
    ymax = isqrt(4 * a * height // delta)
    count = 0
    for y in range(-ymax, ymax + 1):
        root = isqrt(4 * a * height - delta * y * y)
        lo = -((root + b * y) // (2 * a))
        hi = (root - b * y) // (2 * a)
        for x in range(lo, hi + 1):
            if (x or y) and norm(f, (x, y)) % n == 0:
                count += 1
    return count


def check_branch(row):
    n, ell, c, delta, first = (row[k] for k in ("N", "ell", "c", "delta", "a"))
    native, work = native_oracle()
    checks = Counter()

    @cache
    def checked(n, f, cutoff):
        expected = direct(n, f, cutoff)
        actual = native(n, f, cutoff)
        assert actual == expected % 16, (n, f, cutoff, expected, actual)
        checks["distinct_residues_checked"] += 1
        checks["positive_cutoff_nonprincipal_requests"] += cutoff > 0 and f[0] > 1
        checks["exact_counts_at_least_16"] += expected >= 16
        return actual

    answer = guided_factor(n, ell, c, delta, checked)
    assert answer["factors"] == sorted(row["factors"]) and answer["first"] == first
    guided_nodes = work["nodes"]
    hits = []
    for cutoff in (first - 1, first):
        residue = checked(n, principal(delta), cutoff)
        hits.append((residue - 4 * isqrt(cutoff)) % 16)
    assert hits == [0, 4 * row["mixed_orbits"]]
    return {
        "N": n,
        "ell": ell,
        **answer,
        "mixed_signed_residues": hits,
        "guided_oracle_calls": checked.cache_info().hits
        + checked.cache_info().misses
        - 2,
        "guided_nodes": guided_nodes,
        **checks,
    }


def class_census(n):
    # N-only principal branch. An inert 2 makes every r+omega a unit modulo 2.
    c = isqrt(4 * n)
    if 4 * n - c * c <= 4:
        c -= 1
    delta = 4 * n - c * c
    assert gcd(n, c) == 1 and delta % 8 == 3
    m = 1 << ((isqrt(n) // 2).bit_length() - 1)
    f = principal(delta)
    r0 = ((c - 1) // 2) % 8
    classes = set()
    # Check a concrete higher-power congruence independently by enumerating squares.
    modulus = 2048
    squares = {i * i % modulus for i in range(1, modulus, 2)}
    for r in range(r0, m, 8):
        assert norm(f, (r, 1)) <= n and norm(f, (r - m, 1)) <= n
        assert (n * pow(norm(f, (r, 1)), -1, modulus)) % modulus in squares
        a, b, d = reduce_form(restrict(f, *basis("y", m, r)))
        assert gcd(gcd(a, b), d) == 1
        classes.add((a, abs(b), d))
    assert len(classes) == m // 8
    # Factor-free reference on the single audit shell N; this is outside the evaluator.
    represented, mixed = set(), set()
    for y in range(1, isqrt(4 * n // delta) + 1, 2):
        radicand = 4 * n - delta * y * y
        w = isqrt(radicand)
        if w * w != radicand:
            continue
        for signed_w in (w, -w):
            x = (signed_w - y) // 2
            residue = x * pow(y, -1, m) % m
            if residue % 8 == r0:
                represented.add(residue)
                if (signed_w - c * y) % n and (signed_w + c * y) % n:
                    mixed.add(residue)
    return {
        "N": n,
        "delta": delta,
        "stride": m,
        "inequivalent_even_up_to_reflection": len(classes),
        "represent_N": len(represented),
        "contain_a_mixed_representation_of_N": len(mixed),
        "mixed_residues_in_selected_set": sorted(mixed),
    }


def run():
    source = json.loads(Path(__file__).with_name("results.json").read_text())
    report = {}
    for key in ("records", "tie_controls"):
        rows = [check_branch(row) for row in source[key]]
        report[key] = {
            "branches": len(rows),
            "distinct_residues_checked": sum(
                r["distinct_residues_checked"] for r in rows
            ),
            "positive_cutoff_nonprincipal_requests": sum(
                r["positive_cutoff_nonprincipal_requests"] for r in rows
            ),
            "exact_counts_at_least_16": sum(
                r["exact_counts_at_least_16"] for r in rows
            ),
            "max_guided_oracle_calls": max(r["guided_oracle_calls"] for r in rows),
            "max_guided_levels": max(r["levels"] for r in rows),
            "total_guided_nodes": sum(r["guided_nodes"] for r in rows),
            "singletons": sum(r["mixed_signed_residues"][1] == 4 for r in rows),
            "ties": sum(r["mixed_signed_residues"][1] == 8 for r in rows),
            "examples": rows[:2],
        }
    controls = 0
    for f in ((1, 0, 1), (4, 0, 4), (3, 2, 6), (1, 17, 75)):
        for n in (15, 35):
            for cutoff in (0, 1, 2, 4, n):
                assert (
                    radial_count(n, f, cutoff, Counter()) == direct(n, f, cutoff) % 16
                )
                controls += 1
    report["extra_disk_checks"] = controls
    report["inequivalent_locally_soluble"] = [
        class_census(n) for n in (11009, 28829, 38942974393, 33119086993)
    ]
    return report


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--census", type=int, help="check one inert-2 N-only fixture")
    args = parser.parse_args()
    report = run() if args.census is None else class_census(args.census)
    result = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")
