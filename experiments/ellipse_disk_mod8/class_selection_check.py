"""Blind selector versus same-size low minima and the full class list."""

import json
import sys
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt
from pathlib import Path
from random import Random
from time import perf_counter

from class_selection import (
    branch_from_n,
    classes_at,
    denominator_budget,
    propose,
    search_minima,
)
from sympy import isprime

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "src"), str(HERE.parent / "nov26_consolidation")]
from check_identities import forms

from genusleak.forms import reduce_form

SEED = 20261126
SCALES = (10, 14, 18, 22)
PER_SCALE = 24
GAP = 4


def random_prime(rng, scale):
    while True:
        x = rng.randrange(2**scale + 1, 2 ** (scale + 1), 2)
        if isprime(x):
            return x


def fixtures(scale):
    rng = Random(SEED + scale)
    pairs = set()
    while len(pairs) < PER_SCALE:
        p, q = sorted((random_prime(rng, scale), random_prime(rng, scale)))
        if p < q and (p, q) not in pairs:
            pairs.add((p, q))
            yield p, q


def measured_search(br, minima):
    build, search = Counter(), Counter()
    start = perf_counter()
    hits = search_minima(br, minima, build, search)
    return set(hits), {
        **build,
        "values": search["multipoint_values"],
        "seconds": round(perf_counter() - start, 6),
    }


def check(p, q):
    n = p * q
    start = perf_counter()
    br, shortcut = branch_from_n(n)
    if shortcut:
        assert shortcut in (p, q)
        return {"N": n, "shortcut": shortcut}
    work = Counter()
    limit = denominator_budget(n)
    minima, direct = propose(br, limit, GAP, work)
    setup_seconds = perf_counter() - start
    selected, selected_work = measured_search(br, minima)
    small, small_work = measured_search(br, range(1, len(minima) + 1))
    full, full_work = measured_search(br, range(1, br.cap + 1))

    # Truth is used only after the three public searches have finished.
    a = reduce_form((p, -br.c, q))[0]
    truth = {p, q}
    assert full == truth and selected == (truth if a in minima else set())
    assert small == (truth if a <= len(minima) else set()) and direct <= {p}
    assert 0 < a <= br.cap
    return {
        "N": n,
        "p": p,
        "q": q,
        "cap": br.cap,
        "true_minimum": a,
        "limit": limit,
        "largest_proposal": max(minima, default=0),
        "direct": bool(direct),
        "selected": bool(selected),
        "small": bool(small),
        "propose": {**work, "seconds": round(setup_seconds, 6)},
        "selection_search": selected_work,
        "small_search": small_work,
        "full_search": full_work,
    }


def catalogues():
    count = 0
    for delta in range(5, 501):
        if delta % 4 not in (0, 3):
            continue
        cap = isqrt(delta // 3)
        # Check both a sparse index set and the full cap.
        for step in (1, 3):
            keys = range(1, cap + 1, step)
            want = {(a, abs(b), d) for a, b, d in forms(delta) if a in keys}
            assert set(classes_at(delta, keys, Counter())) == want
        count += 1
    return count


def summarize(rows):
    active = [r for r in rows if "shortcut" not in r]
    return {
        "inputs": len(rows),
        "branch_shortcuts": len(rows) - len(active),
        "direct": sum(r["direct"] for r in active),
        "selected": sum(r["selected"] for r in active),
        "combined": sum(r["direct"] or r["selected"] for r in active),
        "equally_many_low_minima": sum(r["small"] for r in active),
        "full": len(active),
        "totals": {
            key: dict(sum((Counter(r[key]) for r in active), Counter()))
            for key in ("propose", "selection_search", "small_search", "full_search")
        },
    }


def nonsquare_cutoffs():
    """Reference values for two adjacent public disks, with shell work charged."""
    n, cutoffs = 140000556000551, (581, 582)
    br, shortcut = branch_from_n(n)
    assert shortcut is None and br.delta % 4 == 0
    d = br.delta // 4
    # Here Q0 = x^2 + d*y^2. Both axes are empty on the punctured disk.
    assert gcd(d, n) == 1 and isqrt(n * max(cutoffs)) < n
    modulus = 8 * 3 * 5 * 7
    squares = {i * i % modulus for i in range(modulus)}
    positive = tests = 0
    results = []
    for j in range(1, max(cutoffs) + 1):
        target = n * j
        last = isqrt((target - 1) // d)
        shell_y = []
        for residue in range(1, modulus + 1):
            if (target - d * residue * residue) % modulus not in squares:
                continue
            for y in range(residue, last + 1, modulus):
                value = target - d * y * y
                x = isqrt(value)
                tests += 1
                if x * x == value:
                    shell_y.append(y)
        positive += len(shell_y)
        if j in cutoffs:
            results.append(
                {
                    "B": j,
                    "positive_disk": positive,
                    "shell_y": sorted(shell_y),
                    "signed_mod16": (4 * positive) % 16,
                    "square_tests_so_far": tests,
                }
            )
    assert [disk["signed_mod16"] for disk in results] == [0, 4]
    (y,) = results[-1]["shell_y"]
    x = isqrt(n * cutoffs[-1] - d * y * y)
    assert 1 < gcd(br.c * y // 2 - x, n) < n
    return {"N": n, "ell": br.ell, "c": br.c, "delta": br.delta, "disks": results}


def run():
    old = json.loads((HERE / "results.json").read_text())["records"]
    controls = [check(*r["factors"]) for r in old]
    holdout = {}
    for scale in SCALES:
        rows = [check(p, q) for p, q in fixtures(scale)]
        holdout[str(scale)] = {"summary": summarize(rows), "rows": rows}
    return {
        "seed": SEED,
        "gap": GAP,
        "next_cutoff_fixture": nonsquare_cutoffs(),
        "catalogue_discriminants": catalogues(),
        "repo_N_only": summarize(controls),
        "collision_N_only": controls[:2],
        "holdout": holdout,
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
