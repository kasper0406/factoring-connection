"""Compare the five raw cap moments and decoded minimum with direct integer shells."""

import json
from collections import Counter
from itertools import combinations
from math import isqrt
from pathlib import Path

from benchmark import odd_primes
from factor import branch_from_n, split
from flint import fmpz
from row_batches import Branch
from weighted_recurrence import decode, fixed_cap_factor, raw_moments
from weighted_recurrence_control import factored_observables
from weighted_recurrence_single import decode_trace, factor_trace, raw_trace


def direct(br):
    raw = [0] * 5
    trace = 0
    first, weight, count = None, 0, 0
    for j in range(1, br.cap + 1):
        t = isqrt(j)
        for y in range(1, isqrt((4 * br.n * j - 1) // br.delta) + 1):
            w2 = 4 * br.n * j - br.delta * y * y
            w = isqrt(w2)
            if w * w != w2:
                continue
            trace += w
            terms = (1, j, j * j, y * y, j * y * y)
            for index, term in enumerate(terms):
                raw[index] += term
            if t * t == j and (w, y) == (br.c * t, br.ell * t):
                continue
            count += 1
            if first is None:
                first = j
            if j == first:
                weight += y * y
    return tuple(x % br.n for x in raw), trace % br.n, (first, weight, count)


def check(br, factor):
    raw, trace, truth = direct(br)
    assert raw == raw_moments(br, Counter())
    oracle_raw, oracle_trace, _ = factored_observables(br, factor, br.n // factor)
    assert oracle_raw == raw
    assert trace == raw_trace(br, Counter()) == oracle_trace
    divisor, tests = decode_trace(br, trace)
    assert divisor in (factor, br.n // factor)
    a, weight, count, route = decode(br, raw)
    assert (a, weight, count) == truth
    d = split(br, a, weight)
    assert 1 < d < br.n and br.n % d == 0
    return route, a, weight, tests


def main():
    small, historical, ties = Counter(), Counter(), Counter()
    trace_tests = Counter()
    for p, q in combinations(odd_primes(241), 2):
        if q <= 2 * p:
            br, d = branch_from_n(p * q)
            if br is None:
                assert d in (p, q)
                small["branch_gcd"] += 1
            else:
                route, _, _, tests = check(br, p)
                small[route] += 1
                trace_tests[tests] += 1
    row_dir = Path(__file__).resolve().parents[1] / "ellipse_row_batches"
    data = json.loads((row_dir / "results.json").read_text())
    collisions = []
    for group, hist in [("records", historical), ("additional_tie_controls", ties)]:
        for record in data[group]:
            br = Branch(*(record[key] for key in ("N", "ell", "c", "delta")))
            route, a, weight, tests = check(br, record["factor"])
            trace_tests[tests] += 1
            assert (a, weight) == (record["a"], record["S"])
            hist[route] += 1
            if group == "records" and len(collisions) < 2:
                collisions.append({"N": br.n, "ell": br.ell, "a": a, "S": weight})
    p = 1329227995784915872903807060280376091
    q = 1661534994731144841129758825350457951
    assert fmpz(p).is_prime() and fmpz(q).is_prime() and p < q <= 2 * p
    br, direct_factor = branch_from_n(p * q)
    assert direct_factor is None
    raw, trace, reductions = factored_observables(br, p, q)
    d, tests = decode_trace(br, trace)
    assert d in (p, q)
    a, weight, count, route = decode(br, raw)
    assert split(br, a, weight) in (p, q)
    control = {
        "bits": br.n.bit_length(),
        "A": br.cap,
        "a": a,
        "count": count,
        "route": route,
        "Gauss_steps": reductions,
        "trace_square_tests": tests,
    }
    print(
        json.dumps(
            {
                "small": small,
                "historical": historical,
                "additional_ties": ties,
                "factor_supplied_control": control,
                "trace_tests": trace_tests,
                "trace_scaling_example": factor_trace(38942974393),
                "collision_controls": collisions,
                "scaling_example": fixed_cap_factor(38942974393),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
