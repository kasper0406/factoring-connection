"""Unramified cubic characters and quartic parity on growing discriminants."""

import json
from argparse import ArgumentParser
from collections import Counter
from itertools import combinations, islice
from math import gcd, isqrt
from pathlib import Path

from class_trace_check import balanced_cases
from flint import fmpz_mod_poly_ctx
from sympy import divisors, integer_nthroot, isprime, primerange

from genusleak.forms import identity, inverse, power, prime_form, reduce_form

# Coefficients s,a,b of x³ + sx² + ax + b; public calibration primes.
FIELDS_4027 = (
    ((0, -44, 113), 13),
    ((0, -8, 15), 17),
    ((1, -7, 12), 19),
    ((0, 10, 1), 29),
)


def discriminant(cubic):
    s, a, b = cubic
    return s * s * a * a - 4 * a**3 - 4 * s**3 * b - 27 * b * b + 18 * s * a * b


def irreducible(cubic):
    s, a, b = cubic
    return b != 0 and all(
        ((x + s) * x + a) * x + b != 0
        for d in divisors(abs(b))
        for x in (int(d), -int(d))
    )


def quartic(cubic):
    """Integer coefficients in ascending order; discard primes dividing 6abD after depression."""
    s, a, b = cubic
    if s:
        a, b = 9 * a - 3 * s * s, 27 * b - 9 * s * a + 2 * s**3
    if b == 1:
        return [-4 * a, -8, 0, 0, 1], 6 * a
    return [-3 * b * b - 4 * a**3, -8 * b * b, -6 * b * b, 0, b * b], 6 * a * b


def local_degree(coeffs, r):
    ring = fmpz_mod_poly_ctx(r)
    f = ring(coeffs)
    x = ring.gen()
    return f.gcd(x.pow_mod(r, f) - x).degree()


def cubic_local(cubic, r):
    s, a, b = cubic
    return local_degree([b, a, s, 1], r)


def canonical(f):
    return min(f, inverse(f))


def reduced_classes(d):
    reduced = set()
    for a in range(1, isqrt(d // 3) + 1):
        for b in range(-a, a + 1):
            c, rem = divmod(b * b + d, 4 * a)
            if rem == 0 and a <= c and gcd(a, b, c) == 1:
                reduced.add(reduce_form((a, b, c)))
    return reduced


def calibrate():
    d = 4027
    unit = identity(-d)
    reduced = reduced_classes(d)
    assert len(reduced) == 9 and all(power(f, 3) == unit for f in reduced)
    labels = {}
    for i, (f, r) in enumerate(FIELDS_4027):
        assert discriminant(f) == -d and irreducible(f)
        profile = tuple(cubic_local(g, r) for g, _ in FIELDS_4027)
        assert profile == tuple(3 if j == i else 0 for j in range(4))
        labels[i] = canonical(prime_form(-d, r))
    assert {unit, *labels.values()} == {canonical(f) for f in reduced}
    return labels


def decode_4027(parity_query, labels):
    """Recover the inverse class pair from at most three calls to a quartic-parity oracle."""

    def bit(i):
        return parity_query(tuple(quartic(FIELDS_4027[i][0])[0]))

    first, second = bit(0), bit(1)
    if not first and not second:
        return identity(-4027), 2
    if not first:
        return labels[0], 2
    if not second:
        return labels[1], 2
    return (labels[3] if bit(2) else labels[2]), 3


def decoder_controls(start, limit, labels):
    histogram, queries = Counter(), Counter()
    profiles = Counter()
    for p, q, c in islice(balanced_cases(4027, start), limit):
        n = p * q
        assert 4 * n == c * c + 4027
        answers = {}
        for cubic, _ in FIELDS_4027:
            h, bad = quartic(cubic)
            assert gcd(bad * 4027, n) == 1
            cp, cq = [cubic_local(cubic, r) for r in (p, q)]
            nr = [local_degree(h, r) for r in (p, q)]
            assert (
                cp == cq
                and cp in (0, 3)
                and all(k in ((0, 4) if cp else (1,)) for k in nr)
            )
            answers[tuple(h)] = (nr[0] * nr[1]) % 2
            profiles[nr[0] * nr[1]] += 1
        actual, calls = decode_4027(answers.__getitem__, labels)
        expected = canonical(reduce_form((p, c, q)))
        assert actual == expected and canonical(power(actual, 2)) == actual
        queries[calls] += 1
        histogram[str(expected)] += 1
    assert sum(queries.values()) == limit and len(histogram) == 5
    return {
        "min_prime": start,
        "inputs": limit,
        "decoded_classes": dict(histogram),
        "adaptive_query_counts": dict(queries),
        "four_control_root_counts": dict(profiles),
    }


def catalogue(max_gap):
    # D fundamental prime; permit any square multiplier and the nearest odd c.
    rows = {}
    for a in range(2, int(integer_nthroot(max_gap // 4, 3)[0]) + 1, 2):
        for b in range(1, isqrt((max_gap - 4 * a**3) // 27) + 1, 2):
            d = 4 * a**3 + 27 * b * b
            cubic = (0, a, b)
            if isprime(d) and irreducible(cubic):
                for ell in range(1, isqrt(max_gap // d) + 1, 2):
                    rows.setdefault(d * ell * ell, (cubic, ell))
    return rows


def coverage(P):
    cat = catalogue(16 * P)
    hits = Counter()
    total = 0
    used = set()
    tested = {}
    for p, q in combinations(primerange(P, 2 * P), 2):
        n = p * q
        c = isqrt(4 * n)
        if c % 2 == 0:
            c -= 1
        gap = 4 * n - c * c
        total += 1
        if gap not in cat:
            continue
        cubic, ell = cat[gap]
        h, bad = quartic(cubic)
        d = -discriminant(cubic)
        assert c * c + ell * ell * d == 4 * n
        if gcd(c * bad * ell, n) > 1:
            hits["direct_factor_or_excluded"] += 1
            continue
        counts = []
        for r in (p, q):
            # Each distinct prime/cubic combination is cheap at these sizes.
            local = (cubic, r)
            if local not in tested:
                tested[local] = cubic_local(cubic, r)
                assert (local_degree(h, r) % 2) == int(tested[local] == 0)
            counts.append(tested[local])
        assert counts[0] == counts[1] and counts[0] in (0, 3)
        hits["bit_1" if counts[0] == 0 else "bit_0"] += 1
        used.add(d)
    return {
        "P": P,
        "balanced_pairs": total,
        "hits": dict(hits),
        "distinct_discriminants_on_hits": len(used),
        "distinct_quartic_prime_checks": len(tested),
    }


def same_character_controls():
    d = 1579
    forms = reduced_classes(d)
    assert len(forms) == 9 and any(power(f, 3) != identity(-d) for f in forms)
    cubic = (0, 4, 7)
    assert irreducible(cubic) and discriminant(cubic) == -d
    h, bad = quartic(cubic)
    rows = []
    for p, q, c in islice(balanced_cases(d), 2):
        n = p * q
        assert gcd(n, d * bad) == 1 and 4 * n == c * c + d
        counts = [local_degree(h, r) for r in (p, q)]
        assert counts == [1, 1] and all(cubic_local(cubic, r) == 0 for r in (p, q))
        square = canonical(power(reduce_form((p, c, q)), 2))
        rows.append(
            {
                "N": n,
                "c": c,
                "control_factors": [p, q],
                "quartic_parity": 1,
                "squared_class": square,
            }
        )
    assert [row["squared_class"][0] for row in rows] == [5, 17]
    return rows


def run():
    labels = calibrate()
    return {
        "cyclic_order_nine_collision": same_character_controls(),
        "rank_two_4027": [
            decoder_controls(start, 96, labels) for start in (3000, 2**21, 2**36)
        ],
        "varying_discriminants": [coverage(p) for p in (256, 1024, 4096)],
    }


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    result = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end="")
