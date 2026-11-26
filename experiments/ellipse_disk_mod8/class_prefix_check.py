"""Ground truth on either side of a mixed minimum and exact public-class catalogues."""

import json
import os
import sys
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / "src"),
    str(ROOT / "experiments/nov26_consolidation"),
    str(ROOT / "experiments/ellipse_row_batches"),
]

from check_identities import forms
from class_prefix import (
    branch_from_n,
    candidate_lifts,
    canonical,
    decision,
    represented_factors,
)
from row_batches import Branch

from genusleak.forms import compose, reduce_form


def catalogues(gp):
    orders = classes = 0
    for delta in range(5, 501):
        if delta % 4 not in (0, 3):
            continue
        cap = isqrt(delta // 3)
        expected = list(forms(delta))
        for low, high in ((1, cap), (max(1, cap // 2), cap)):
            want = {canonical(f) for f in expected if low <= compose(f, f)[0] <= high}
            got = candidate_lifts(delta, low, high, gp, Counter())
            assert set(got) == want, (
                delta,
                low,
                high,
                want.difference(got),
                set(got).difference(want),
            )
        orders += 1
        classes += len(expected)
    return {"orders": orders, "proper_classes": classes}


def check_cohort(rows, gp):
    stats = Counter()
    for row in rows:
        br = Branch(row["N"], row["ell"], row["c"], row["delta"])
        a = row["a"]
        for low, high, positive in (
            (1, a - 1, False),
            (1, a, True),
            (a, a, True),
            (1, br.cap, True),
        ):
            answer = decision(br, low, high, gp)
            if positive:
                assert answer["factor"] in row["factors"]
            else:
                assert answer["factor"] is None
            stats["decisions"] += 1
        stats["branches"] += 1
        stats["nonsquare_minima"] += isqrt(a) ** 2 != a
    return stats


def public_line_checks(br, candidates):
    """Reduce both known norm-N sublattices after dividing their form by N."""
    n = br.n
    t = br.c * pow(br.ell, -1, n) % n
    for a, b, d in candidates:
        for root in (t, -t):
            r = (root - b) * pow(2 * a, -1, n) % n
            last, rem = divmod(a * r * r + b * r + d, n)
            assert rem == 0
            assert reduce_form((n * a, 2 * a * r + b, last)) == (a, b, d)
    return 2 * len(candidates)


def local_frame_checks(br, candidates):
    """Construct unit-determinant frames with prescribed middle coefficient locally."""
    checked = 0
    for a, b, d in candidates:
        for m in (8, 32, 9, 25, 49):
            for x, y, z, t in ((1, 0, 0, 1), (0, 1, -1, 0), (1, 1, 0, 1)):
                val = a * x * x + b * x * y + d * y * y
                if gcd(val, m) == 1:
                    break
            else:
                raise ArithmeticError("primitive form missing a unit")
            middle = 2 * a * x * z + b * (x * t + y * z) + 2 * d * y * t
            g = gcd(2 * val, m)
            assert (br.c - middle) % g == 0
            k = ((br.c - middle) // g) * pow(2 * val // g, -1, m // g)
            z, t = z + k * x, t + k * y
            assert (
                x * t - y * z == 1
                and (2 * a * x * z + b * (x * t + y * z) + 2 * d * y * t - br.c) % m
                == 0
            )
            checked += 1
    return checked


def scaling(gp, n, p, q, *, blind=True):
    br, d = branch_from_n(n)
    assert d is None and p * q == n
    # Control and proposed single-index run only. The blind run receives just N.
    oracle_class = reduce_form((p, br.c, q))
    a = compose(oracle_class, oracle_class)[0]
    results = {}
    intervals = [("given_index", a, a)]
    if blind:
        intervals.append(("blind_cap", 1, br.cap))
    for tag, lo, hi in intervals:
        build, search = Counter(), Counter()
        start = perf_counter()
        candidates = candidate_lifts(br.delta, lo, hi, gp, build)
        middle = perf_counter()
        assert canonical(oracle_class) in candidates
        hits = represented_factors(br, candidates, search)
        end = perf_counter()
        checked = public_line_checks(br, candidates)
        assert set(hits) == {p, q}
        results[tag] = {
            "public_line_checks": checked,
            "local_frame_checks": local_frame_checks(br, candidates),
            "build": build,
            "search": search,
            "build_seconds": middle - start,
            "search_seconds": end - middle,
        }
    return {
        "N": n,
        "delta": br.delta,
        "first": a,
        "cap": br.cap,
        "factors": [p, q],
        **results,
    }


def run(gp):
    base = Path(__file__).parent
    sample = json.loads((base / "results.json").read_text())
    old = json.loads((base.parent / "ellipse_aggregate/results.json").read_text())[
        "scaling"
    ]
    measures = [scaling(gp, row["N"], *row["shell"]["factors"]) for row in old]
    return {
        "catalogues": catalogues(gp),
        "repository": check_cohort(sample["records"], gp),
        "ties": check_cohort(sample["tie_controls"], gp),
        "scaling": measures,
        "larger_controls": [
            scaling(gp, 10000019 * 14000029, 10000019, 14000029),
            scaling(gp, 1000000007 * 1400000023, 1000000007, 1400000023, blind=False),
        ],
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--gp", default=os.environ.get("GP", "gp"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(run(args.gp), indent=2) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")
