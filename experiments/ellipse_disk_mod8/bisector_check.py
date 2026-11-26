"""Exact controls, two slope groups, coverage samples, and matched single-shell timings."""

import json
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt
from pathlib import Path
from random import Random
from time import perf_counter

from bisector import (
    WindowProducts,
    block_summary,
    branch_from_n,
    discriminant_divisors,
    linear_hits,
    mixed_points,
)
from flint import fmpz
from row_batches import Branch


def all_points(br, j):
    for y in range(1, isqrt(4 * br.n * j // br.delta) + 1):
        value = 4 * br.n * j - br.delta * y * y
        w = isqrt(value)
        if w and w * w == value and (w - br.delta * y) % 2 == 0:
            yield w, y


def reference(br, j=1):
    return sorted(
        (w, y)
        for w, y in all_points(br, j)
        if (br.ell * w - br.c * y) % br.n and (br.ell * w + br.c * y) % br.n
    )


def verify(br):
    prep, search = Counter(), Counter()
    points = mixed_points(br, prep, search)
    full = sorted(all_points(br, 1))
    assert sorted(points + [(br.c, br.ell)]) == full, (br, points)
    signed_x = [((sw - (br.delta % 2) * y) // 2, y) for w, y in full for sw in (w, -w)]
    for r in range(8):
        group = block_summary(br, points, 32, r, 8)
        assert (
            group["signed_mod16"]
            == 2 * sum((x - r * y) % 8 == 0 for x, y in signed_x) % 16
        )
        for child in group["children"]:
            t = child["residue"]
            assert (
                child["signed_mod16"]
                == 2 * sum((x - t * y) % 16 == 0 for x, y in signed_x) % 16
            )
    for w, y in points:
        d = gcd((br.c * y - br.ell * w) // 2, br.n)
        assert 1 < d < br.n and br.n % d == 0
    return len(points)


def timed(n):
    start = perf_counter()
    br, shortcut = branch_from_n(n)
    assert not shortcut
    prep, search = Counter(), Counter()
    points = mixed_points(br, prep, search)
    group = block_summary(br, points, 65536, (br.c - br.delta % 2) // 2, 8)
    kernel_seconds = perf_counter() - start
    start = perf_counter()
    expected = reference(br)
    assert points == expected
    row_seconds = perf_counter() - start
    return {
        "N": n,
        "bits": n.bit_length(),
        "delta": br.delta,
        "mixed_points": points,
        "raw_mod16": 4 * (len(points) + 1) % 16,
        "group": group,
        "discriminant_work": prep,
        "search_work": search,
        "reference_rows": isqrt(4 * n // br.delta),
        "bisector_seconds": kernel_seconds,
        "reference_seconds": row_seconds,
    }


def constructed_positive(k):
    # Fixture generation only. The two prime norms use successive odd real coordinates.
    x = (1 << k) + 1
    p = x * x // 2
    p += 1 - p % 2
    while not (fmpz(p).is_prime() and fmpz(p + x + 1).is_prime()):
        p += 2
    q = p + x + 1
    return p * q


def draw_prime(rng, pscale):
    while True:
        z = rng.randrange(pscale + 1, 2 * pscale, 2)
        if fmpz(z).is_prime():
            return z


def coverage():
    rows = []
    for power in (8, 12, 16, 20):
        rng, seen = Random(8294 + power), set()
        hits = Counter()
        while len(seen) < 64:
            p, q = sorted((draw_prime(rng, 1 << power), draw_prime(rng, 1 << power)))
            if p == q or p * q in seen:
                continue
            n = p * q
            seen.add(n)
            first = []
            for offset in range(3):
                c = isqrt(4 * n) - offset
                delta = 4 * n - c * c
                if delta <= 4 or gcd(n, c * delta) != 1:
                    first.append(None)
                    continue
                br = Branch(n, 1, c, delta)
                value = next(
                    (j for j in range(1, min(8, br.cap) + 1) if reference(br, j)), None
                )
                first.append(value)
                if offset == 0:
                    assert bool(verify(br)) == (value == 1)
            if any(v == 1 for v in first):
                hits["any_of_three_B1"] += 1
            if any(v is not None for v in first):
                hits["any_of_three_through_8_capped"] += 1
            hits["first_branch_B1"] += first[0] == 1
        rows.append(
            {
                "prime_scale_power": power,
                "pairs": len(seen),
                **{
                    k: hits[k]
                    for k in (
                        "first_branch_B1",
                        "any_of_three_B1",
                        "any_of_three_through_8_capped",
                    )
                },
            }
        )
    return rows


def run():
    fixture = json.loads(Path(__file__).with_name("results.json").read_text())
    cohorts = []
    for name in ("records", "tie_controls"):
        counter = Counter()
        for r in fixture[name]:
            br = Branch(r["N"], r["ell"], r["c"], r["delta"])
            count = verify(br)
            assert bool(count) == (r["a"] == 1)
            counter["branches"] += 1
            counter["positive"] += bool(count)
        cohorts.append({"name": name, **counter})
    small, even = Counter(), Counter()
    primes = [i for i in range(3, 242, 2) if fmpz(i).is_prime()]
    for p in primes:
        for q in primes:
            if not p < q <= 2 * p:
                continue
            n = p * q
            br, d = branch_from_n(n)
            assert not d
            small["branches"] += 1
            small["positive"] += bool(verify(br))
            if q > 47:
                continue
            for ell in (2, 4, 6, 8):
                for c in range(2, isqrt(4 * n) + 1, 2):
                    delta, r = divmod(4 * n - c * c, ell * ell)
                    if (
                        r
                        or delta <= 4
                        or delta % 4 not in (0, 3)
                        or gcd(c * ell * delta, n) > 1
                    ):
                        continue
                    even["branches"] += 1
                    even["positive"] += bool(verify(Branch(n, ell, c, delta)))
    windows_checked = 0
    for n in (15, 35, 77):
        work = Counter()
        windows = WindowProducts(n, work)
        for l, s in ((1, 1), (9, 3), (21, 13), (89, 17)):
            for ratio in range(1, n):
                if gcd(ratio, n) != 1:
                    continue
                expected = {
                    (i, j)
                    for i in range(1, l + 1)
                    for j in range(1, s + 1)
                    if gcd(i - ratio * j, n) > 1
                }
                actual = list(linear_hits(n, l, s, ratio, windows, work))
                assert len(actual) == len(expected) and set(actual) == expected
                windows_checked += 1
    for d in range(1, 1025):
        assert discriminant_divisors(d, Counter()) == [
            k for k in range(1, d + 1) if d % k == 0
        ]
    named = [timed(n) for n in (33119086993, 38942974393)]
    assert [row["group"]["mixed_residues"] for row in named] == [[12122], []]
    # These are intentionally positive; the even physical y is outside the odd-y chart.
    series = [timed(constructed_positive(k)) for k in (12, 16, 20, 23)]
    assert all(len(row["mixed_points"]) == 1 for row in series)
    return {
        "repository": cohorts,
        "small_N_only": small,
        "even_ell": even,
        "window_controls": windows_checked,
        "divisor_controls": 1024,
        "named_fixtures": named,
        "constructed_positive_series": series,
        "coverage": coverage(),
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")
