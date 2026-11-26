"""Blind class and rectangle comparisons on the same primitive nonsquare band."""

import json
import os
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt, prod
from pathlib import Path
from statistics import median
from time import perf_counter

from benchmark import Branch, independent
from class_prefix_check import forms
from conductor_band import (
    branch_from_n,
    candidates,
    canonical,
    compose,
    evaluate,
    public_factorization,
    reduce_form,
    verify_levels,
)
from rectangle_rays import Rectangles
from sympy import primerange

HERE = Path(__file__).resolve().parent


def physical_reference(br, hits, cutoff):
    # Use the physical coordinates and the original definition of a removable h.
    even = br.delta % 4 == 0
    d, level = (br.delta // 4, br.n) if even else (br.delta, 4 * br.n)
    excluded = prod(p for p in range(2, isqrt(d) + 1) if d % (p * p) == 0)
    total = 0
    for t, y in hits:
        if 4 * t <= cutoff or t > cutoff or isqrt(t) ** 2 == t:
            continue
        x = isqrt(level * t - d * y * y)
        if not any(
            x % h == 0 and t % (h * h) == 0 and gcd(h, excluded) == 1
            for h in range(2, isqrt(t) + 1)
        ):
            total += 1
    return total


def class_vectors(f, cutoff, excluded):
    """Independent literal square-box enumeration of one reduced form."""
    a, b, c = f
    high_x, high_y = isqrt(2 * cutoff // a), isqrt(2 * cutoff // c)
    result = []
    for y in range(high_y + 1):
        for x in range(-high_x if y else 1, high_x + 1):
            t = a * x * x + b * x * y + c * y * y
            if 4 * t <= cutoff or t > cutoff or isqrt(t) ** 2 == t:
                continue
            g = gcd(x, y)
            while (d := gcd(g, excluded)) > 1:
                g //= d
            if g == 1:
                result.append((x, y))
    return result


def bad_product(delta):
    d = delta // 4 if delta % 4 == 0 else delta
    return prod(p for p, e in public_factorization(d, Counter()).items() if e >= 2)


def all_class_checks(gp):
    stats = Counter()
    for delta in range(5, 1001):
        if delta % 4 not in (0, 3):
            continue
        fs = {canonical(f) for f in forms(delta)}
        cap = isqrt(delta // 3)
        excluded = bad_product(delta)
        # Selection depends on the discriminant alone. No N operations run here.
        br = Branch(0, 1, 0, delta)
        for b in {cap, cap // 2}:
            selected, prime = candidates(br, b, gp, Counter())
            want = set()
            for f in fs:
                square = compose(f, f)
                events = class_vectors(square, b, excluded)
                directions = {(x // gcd(x, y), y // gcd(x, y)) for x, y in events}
                assert len(directions) <= (1 if 4 * b * b < delta else 3)
                if len(directions) > 1:
                    stats["multi_direction_labels"] += 1
                if len(events) & 1:
                    want.add(f)
                stats["labels"] += 1
            assert want == set(selected), (delta, b, want ^ set(selected))
            stats["lists"] += 1
            stats["lists_with_descent"] += prime > 1
    return stats


def timing_comparison(br, b, expected, gp):
    """Rotate test order over five repetitions, with fresh evaluator objects."""
    tags = [
        "rectangle",
        "all_levels",
        "inert_first",
        "geometry_levels",
        "class_membership",
    ]
    times = {tag: [] for tag in tags}
    for run in range(len(tags)):
        for tag in tags[run:] + tags[:run]:
            start = perf_counter()
            if tag == "geometry_levels":
                answer = verify_levels(br, b)["parity"]
            elif tag == "class_membership":
                answer = evaluate(br, b, gp)["parity"]
            else:
                rect = Rectangles(br)
                levels = None if tag == "rectangle" else range(b // 4 + 1, b + 1)
                if tag == "inert_first":
                    levels = [
                        t
                        for t in levels
                        if isqrt(t) ** 2 != t and not rect.primitive_rule(t)[1]
                    ]
                answer = rect.count(
                    b // 4 + 1,
                    b,
                    modulus=2,
                    base_nonsquares=True,
                    proposed_levels=levels,
                )
            times[tag].append(perf_counter() - start)
            assert answer == expected
    return {tag: median(values) for tag, values in times.items()}


def compare(br, b, gp, *, other_baselines=False):
    ans = evaluate(br, b, gp)
    start = perf_counter()
    rectangles = Rectangles(br)
    exact = rectangles.count(b // 4 + 1, b, modulus=None, base_nonsquares=True)
    end = perf_counter()
    direct = verify_levels(br, b)
    assert direct["exact_count"] == exact and ans["parity"] == exact % 2
    baselines = {}
    if other_baselines:
        # Same target-root verifier, with and without a cheap inert-prime screen.
        for tag in ("all_levels", "inert_first"):
            start_control = perf_counter()
            control = Rectangles(br)
            proposed = range(b // 4 + 1, b + 1)
            if tag == "inert_first":
                proposed = [
                    t
                    for t in proposed
                    if isqrt(t) ** 2 != t and not control.primitive_rule(t)[1]
                ]
            got = control.count(
                b // 4 + 1,
                b,
                modulus=None,
                base_nonsquares=True,
                proposed_levels=proposed,
            )
            baselines[tag] = {
                "seconds": perf_counter() - start_control,
                "work": control.metrics(),
            }
            assert got == exact
    return {
        **ans,
        "B": b,
        "level_route": direct,
        "same_verifier_controls": baselines,
        "exact_rectangle_count": exact,
        "rectangle_seconds": end - start,
        "rectangle_work": rectangles.metrics(),
    }


def cohort(rows, gp):
    stats = Counter()
    for row in rows:
        br = Branch(row["N"], row["ell"], row["c"], row["delta"])
        hits = independent(br)
        for b in {max(0, row["a"] - 1), row["a"], br.cap}:
            ans = compare(br, b, gp)
            assert ans["exact_rectangle_count"] == physical_reference(br, hits, b)
            if ans["factor_witness"]:
                assert ans["factor_witness"] in row["factors"]
            stats["nonempty_even"] += (
                ans["exact_rectangle_count"] > 0 and ans["parity"] == 0
            )
            stats["queries"] += 1
            stats["positive_bits"] += ans["parity"]
        stats["branches"] += 1
    return stats


def run(gp):
    orig = json.loads((HERE / "results.json").read_text())
    held = json.loads((HERE / "class_selection_results.json").read_text())["holdout"]
    histogram = Counter()
    for tier in held.values():
        for row in tier["rows"]:
            br, shortcut = branch_from_n(row["N"])
            assert shortcut is None
            fac = public_factorization(br.delta, Counter())
            histogram[
                max((p for p, e in fac.items() if p != 2 and e >= 2), default=1)
            ] += 1
    rows = [tier["rows"][0] for tier in held.values()]
    rows.extend([held["18"]["rows"][7], held["22"]["rows"][16]])
    rows.append({"N": 140000556000551, "p": 10000019, "q": 14000029})
    cases = []
    for row in rows:
        br, shortcut = branch_from_n(row["N"])
        assert shortcut is None
        # Only independent labels and benchmark cutoff selection use fixture factors.
        f = reduce_form((row["p"], br.c, row["q"]))
        square = compose(f, f)
        first = square[0]
        measurements = [
            compare(br, b, gp, other_baselines=True) for b in (first - 1, first)
        ]
        for item in measurements:
            assert item["exact_rectangle_count"] == len(
                class_vectors(square, item["B"], bad_product(br.delta))
            )
        fac = public_factorization(br.delta, Counter())
        prime_discriminant = list(fac.items()) == [(br.delta, 1)]
        split_primes = (
            sum(
                pow(-br.delta, (p - 1) // 2, p) == 1
                for p in primerange(max(3, first // 4 + 1), first + 1)
            )
            if prime_discriminant
            else None
        )
        cases.append(
            {
                "N": br.n,
                "delta": br.delta,
                "c": br.c,
                "first": first,
                "split_odd_primes_in_upper_band_if_prime_discriminant": split_primes,
                "cutoffs": measurements,
                "median_seconds_at_first": timing_comparison(
                    br, first, measurements[1]["parity"], gp
                )
                if prime_discriminant or row["N"] == 140000556000551
                else None,
            }
        )
    return {
        "all_classes": all_class_checks(gp),
        "small": cohort(orig["records"] + orig["tie_controls"], gp),
        "heldout_conductor_steps": histogram,
        "cases": cases,
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--gp", default=os.environ.get("GP", "gp"))
    opt = parser.parse_args()
    text = json.dumps(run(opt.gp), indent=2) + "\n"
    if opt.output:
        opt.output.write_text(text)
    else:
        print(text, end="")
