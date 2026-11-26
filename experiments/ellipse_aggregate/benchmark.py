"""Independent small controls and input-only timings for the two moment evaluators."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import isqrt
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter
import json
import platform

from factor import factor, residue_disk
from row_batches import Branch, disk_query


def odd_primes(limit):
    return [
        p
        for p in range(3, limit + 1, 2)
        if all(p % d for d in range(3, isqrt(p) + 1, 2))
    ]


def trial(n):
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return [d, n // d]
    raise ValueError("expected an odd composite")


def next_prime(n):
    v = n + 1 if n % 2 == 0 else n + 2
    while any(v % d == 0 for d in range(3, isqrt(v) + 1, 2)):
        v += 2
    return v


def timed(f, *args):
    start = perf_counter()
    result = f(*args)
    return result, perf_counter() - start


def compare_prefixes(br):
    disk = residue_disk(br, br.cap, Counter())
    sat = hits = 0
    for B in range(1, br.cap + 1):
        count, moment = disk.prefix(B)
        ref = disk_query(br, B)
        assert 4 * count % 16 == ref["signed_count_mod16"]
        assert moment * br.c * br.c % br.n == ref["positive_second_moment_mod_N"]
    # At the full rectangular width the naive polynomial match vanishes
    # modulo 4N on every y; count how many y have a hit even in that rectangle.
    mod = 4 * br.n
    X = isqrt(mod * br.cap)
    squares = {w * w % mod for w in range(1, X + 1)}
    for y in range(1, isqrt(mod * br.cap // br.delta) + 1):
        target = (-br.delta * y * y) % mod
        if target in squares:
            hits += 1
        else:
            sat += 1
    return br.cap, sat, hits


def product_control(br):
    # Explicitly check the product test on both truth values, on the B=1 rectangle.
    mod = 4 * br.n
    residues = [w * w % mod for w in range(1, isqrt(mod) + 1)]
    tests = absent = 0
    for y in range(1, isqrt(mod // br.delta) + 1):
        t = (-br.delta * y * y) % mod
        value = 1
        for r in residues:
            value = value * (t - r) % mod
        assert value == 0
        tests += 1
        absent += t not in residues
    return tests, absent


def small_controls():
    pairs = [(p, q) for p, q in combinations(odd_primes(241), 2) if q <= 2 * p]
    gcd_paths = adjusted = 0
    for p, q in pairs:
        n = p * q
        runs = [factor(n, m) for m in ("shell", "join")]
        assert all(r["factors"] == [p, q] for r in runs)
        if runs[0]["route"] == "branch_gcd":
            gcd_paths += 1
        else:
            assert runs[0]["a"] == runs[1]["a"] and runs[0]["S"] == runs[1]["S"]
        adjusted += 4 * n - isqrt(4 * n) ** 2 <= 4
    cases = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "ellipse_row_batches"
            / "public_cases.json"
        ).read_text()
    )
    checked = ys_empty = ys_match = products = false_positives = 0
    for params in cases:
        br = Branch(*params)
        ncount, empty, hit = compare_prefixes(br)
        checked += ncount
        ys_empty += empty
        ys_match += hit
        tests, absent = product_control(br)
        products += tests
        false_positives += absent
    return dict(
        exhaustive_balanced_pairs=len(pairs),
        unusable_discriminant_adjustments=adjusted,
        branch_gcd_shortcuts=gcd_paths,
        source_branches=len(cases),
        independent_disk_prefixes=checked,
        cap_rectangle_y_with_match=ys_match,
        cap_rectangle_y_without_match=ys_empty,
        explicit_B1_product_tests=products,
        explicit_B1_product_false_positives=false_positives,
    )


def scaling():
    rng = Random(20261127)
    out = []
    for bits in (10, 12, 14, 16, 18):
        for case in range(3):
            # Fixed recipe; no selection by discriminant, minimum, hit or timing.
            p = next_prime(rng.randrange(2 ** (bits - 1), 3 * 2 ** (bits - 2)))
            q = next_prime(rng.randrange(p + 1, 2 * p - 2))
            n = p * q
            assert p < q <= 2 * p
            trialruns = [timed(trial, n) for _ in range(5)]
            assert all(f == [p, q] for f, t in trialruns)
            row = dict(
                N=n, bits=n.bit_length(), trial_seconds=median(t for _, t in trialruns)
            )
            for m in ("shell", "join"):
                r, t = timed(factor, n, m)
                assert r["factors"] == [p, q]
                row[m] = dict(seconds=t, **r)
            assert row["shell"].get("S") == row["join"].get("S")
            out.append(row)
            print(
                n,
                n.bit_length(),
                row["shell"].get("a"),
                round(row["shell"]["seconds"], 3),
                round(row["join"]["seconds"], 3),
                flush=True,
            )
    return out


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = small_controls()
    print(summary, flush=True)
    rows = scaling()
    args.output.write_text(
        json.dumps(
            dict(python=platform.python_version(), controls=summary, scaling=rows),
            indent=2,
        )
        + "\n"
    )
