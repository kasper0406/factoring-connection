"""Check global projection, rectangle clips, ray cancellation, and factor readout."""

import json
import sys
from argparse import ArgumentParser
from collections import Counter
from math import isqrt
from pathlib import Path
from time import perf_counter

from benchmark import Branch, independent
from rectangle_rays import Rectangles, branch_from_n, recover

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "nov26_consolidation"), str(HERE.parents[1] / "src")]
from check_identities import count as form_count

from genusleak.forms import compose, reduce_form


def mixed_form(br, p, q):
    f = reduce_form((p, br.c, q))
    return compose(f, f)


def held_out(tiers):
    result = []
    for tier in tiers.values():
        for row in tier["rows"][:2]:
            br, shortcut = branch_from_n(row["N"])
            assert shortcut is None
            f = mixed_form(br, row["p"], row["q"])
            counter = Rectangles(br)
            cutoffs = sorted({f[0] - 1, f[0], min(2 * f[0], br.cap), br.cap})
            values = []
            for b in cutoffs:
                exact = isqrt(b) + form_count(f, b)
                assert counter.count(1, b, modulus=2) == exact % 2
                assert counter.count(1, b, modulus=4) == exact % 4
                values.append(exact)
            result.append(
                {"N": br.n, "cutoffs": cutoffs, "reference_positive_counts": values}
            )
    return result


def modes(br, cutoff):
    result = {}
    for label, rays, modulus, screen in (
        ("plain", False, None, False),
        ("rays_exact", True, None, False),
        ("parity_unscreened", True, 2, False),
        ("parity", True, 2, True),
        ("mod4", True, 4, True),
    ):
        start = perf_counter()
        rect = Rectangles(br)
        answer = rect.count(1, cutoff, modulus=modulus, rays=rays, screen=screen)
        result[label] = {
            "answer": answer,
            "seconds": round(perf_counter() - start, 6),
            "work": rect.metrics(),
        }
    exact = result["plain"]["answer"]
    assert (
        result["rays_exact"]["answer"] == exact
        and result["parity_unscreened"]["answer"] == exact % 2
    )
    assert (
        result["mod4"]["answer"] == exact % 4
        and result["parity"]["answer"] == exact % 2
    )
    return result


def controls(old):
    stats = Counter()
    for row in old["records"] + old["tie_controls"]:
        br = Branch(row["N"], row["ell"], row["c"], row["delta"])
        hits = independent(br)  # Direct area enumeration using only the public branch.
        counter = Rectangles(br)
        for hi in range(br.cap + 1):
            for lo, ymin, ymax in (
                (1, 1, None),
                (max(1, hi // 2), 2, 7),
                (max(1, hi), 1, 2 + hi % 5),
            ):
                want = sum(
                    lo <= j <= hi and y >= ymin and (ymax is None or y <= ymax)
                    for j, y in hits
                )
                for rays, modulus in (
                    (False, None),
                    (True, None),
                    (True, 2),
                    (True, 4),
                ):
                    got = counter.count(lo, hi, ymin, ymax, rays=rays, modulus=modulus)
                    assert got == (want if modulus is None else want % modulus), (
                        br,
                        lo,
                        hi,
                        ymin,
                        ymax,
                        got,
                        want,
                    )
                    stats["count_checks"] += 1
        a, d, _ = recover(br)
        assert a == row["a"] and d in row["factors"]
        stats["factor_checks"] += 1
        stats["squareful_D"] += bool(counter.bad)
        stats["tie_cases"] += row["mixed_orbits"] == 2
    return stats


def scaling(old):
    result = []
    for row in old:
        br, shortcut = branch_from_n(row["N"])
        assert shortcut is None
        cutoff = row["first"]
        measured = modes(br, cutoff)
        f = mixed_form(br, *row["factors"])
        assert f[0] == cutoff and measured["plain"]["answer"] == isqrt(
            cutoff
        ) + form_count(f, cutoff)
        result.append({"N": br.n, "delta": br.delta, "cutoff": cutoff, **measured})
    return result


def run():
    old = json.loads((HERE / "results.json").read_text())
    scale = json.loads((HERE / "class_prefix_results.json").read_text())["scaling"]
    previous = json.loads((HERE / "class_selection_results.json").read_text())
    fixture = previous["next_cutoff_fixture"]
    br, shortcut = branch_from_n(fixture["N"])
    assert shortcut is None
    exact = {r["B"]: r["positive_disk"] for r in fixture["disks"]}
    cutoffs = [17, 65, 145, 257, 401, *exact]
    fixed = {}
    for b in cutoffs:
        m = modes(br, b)
        assert m["plain"]["answer"] == exact.get(b, isqrt(b))
        fixed[b] = m
    start = perf_counter()
    a, d, work = recover(br)
    elapsed = perf_counter() - start
    assert a == 582 and 1 < d < br.n and br.n % d == 0
    # Verify the public nonsquare residue separately from timed evaluation.
    witness = (7, 39389)
    assert (witness[1] ** 2 - br.n * witness[0]) % (br.delta // 4) == 0
    return {
        "controls": controls(old),
        "held_out": held_out(previous["holdout"]),
        "scaling": scaling(scale),
        "fixture": {
            "N": br.n,
            "delta": br.delta,
            "cutoffs": fixed,
            "nonsquare_residue": witness,
            "factor": d,
            "first": a,
            "decoder_work": work,
            "decoder_seconds": round(elapsed, 6),
        },
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    opt = parser.parse_args()
    data = json.dumps(run(), indent=2) + "\n"
    if opt.output:
        opt.output.write_text(data)
    else:
        print(data, end="")
