"""Compare the first additive digit at 2 and four successive N-th powers."""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations
from math import gcd
from random import Random

from sympy import isprime, n_order, primerange


def transforms(x):
    return ((x, 1), (-1 - x, 1), (1, x), (-1, 1 + x),
            (-1 - x, x), (-x, 1 + x))


def orbit_mask(n, x, y):
    """Bit i records a proper gcd between y and transform i of x."""
    return sum(1 << i for i, (a, b) in enumerate(transforms(x))
               if gcd(b, n) == 1 and 1 < gcd(b * y - a, n) < n)


def scan(n):
    orbit = [2]
    hits = []
    propagated_hits = []
    three_power = 3
    easy_split = False
    for _ in range(5):
        x = orbit[-1]
        y = pow(x, n, n)
        shifted = pow(x + 1, n, n)
        divisor = gcd(shifted - y - 1, n)
        hits.append(divisor if 1 < divisor < n else None)
        three_power = pow(three_power, n, n)
        propagated = gcd(three_power - y - 1, n)
        propagated_hits.append(propagated if 1 < propagated < n else None)
        # x(x^(N-1)-1) and (x+1)((x+1)^(N-1)-1) reuse both powers.
        tests = (x, x + 1, x * x + x + 1, y - x, shifted - x - 1)
        easy_split |= any(1 < gcd(test, n) < n for test in tests)
        orbit.append(y)
    assert hits[0] == propagated_hits[0]
    return orbit, hits, propagated_hits, easy_split


def count(pairs):
    totals = Counter()
    first_example = None
    for p, q in pairs:
        n = p * q
        orbit, hits, propagated_hits, easy_split = scan(n)
        totals["pairs"] += 1
        for index, divisor in enumerate(hits):
            totals[f"j{index}"] += bool(divisor)
            totals[f"propagated_{index + 1}"] += bool(propagated_hits[index])
        totals["ordinary_controls"] += easy_split
        if not hits[0] and any(propagated_hits[1:]):
            totals["new_propagated"] += 1
            observations = [(orbit[k], r) for k, r in enumerate(propagated_hits, 1) if r]
            if any((y - 2) % r for y, r in observations):
                totals["new_propagated_without_return"] += 1
            if any(all((b * y - a) % r for a, b in transforms(2)) for y, r in observations):
                totals["new_propagated_outside_symmetry"] += 1
        if hits[0] or not any(hits[1:]):
            continue
        totals["new_on_orbit"] += 1
        totals["j1_new"] += bool(hits[1])
        if easy_split:
            continue
        totals["new_after_ordinary"] += 1
        totals["j1_after_ordinary"] += bool(hits[1])
        # Include x_5, the power-map output for the last additive query.
        masks = [orbit_mask(n, x, y) for x, y in combinations(orbit, 2)]
        if any(mask & 1 for mask in masks):
            totals["ordinary_collision"] += 1
        if any(masks):
            totals["any_symmetry_collision"] += 1
            continue
        totals["new_after_all"] += 1
        totals["j1_after_all"] += bool(hits[1])
        totals["new_after_propagated_too"] += not any(propagated_hits)
        if first_example is None:
            j = next(i for i, divisor in enumerate(hits) if divisor)
            first_example = (n, p, q, j, orbit[j], hits[j])
    return dict(totals), first_example


def draw_prime(rng, lo, hi):
    while True:
        value = rng.randrange(lo + 1, hi, 2)
        if isprime(value):
            return value


def separated_pairs(bits):
    rng = Random(20261128 + bits)
    scale = 1 << bits
    for _ in range(1000):
        yield draw_prime(rng, scale, 2 * scale), draw_prime(rng, 3 * scale, 4 * scale)


def order_bins():
    """Classify the factor returned at j=1 using the known diagnostic primes."""
    primes = list(primerange(4096, 8192))
    labels = {}
    rows = {label: Counter() for label in ("1", "2", "3–8", ">8")}
    for r in primes:
        index = (r - 1) // n_order(2, r)
        label = str(index) if index <= 2 else ("3–8" if index <= 8 else ">8")
        labels[r] = label
        rows[label]["all_sides"] += len(primes) - 1
    for p, q in combinations(primes, 2):
        n = p * q
        x = pow(2, n, n)
        divisor = gcd(pow(x + 1, n, n) - pow(x, n, n) - 1, n)
        if divisor == 1 or divisor == n:
            continue
        row = rows[labels[divisor]]
        row["j1"] += 1
        orbit, hits, _, easy_split = scan(n)
        assert hits[1] == divisor
        if (not hits[0] and not easy_split
                and not any(orbit_mask(n, a, b) for a, b in combinations(orbit, 2))):
            row["j1_after_all"] += 1
    return rows


def complete_orbit_certificate():
    n = 31439
    for seed, expected in ((2, 24709), (3, 5996), (5, 8386)):
        x = seed
        orbit = []
        product = 1
        for _ in range(6):
            assert x not in orbit
            orbit.append(x)
            y = pow(x, n, n)
            product = product * (pow(x + 1, n, n) - y - 1) % n
            x = y
        assert x == seed and product == expected and gcd(product, n) == 1
        print("orbit", seed, orbit, "product", product)

    def simple_gcd(x):
        return gcd(x * (x + 1) * (x * x + x + 1), n)

    counts = Counter()
    for x in range(n):
        y = pow(x, n, n)
        divisor = gcd(pow(x + 1, n, n) - y - 1, n)
        assert divisor == simple_gcd(x) == simple_gcd(y)
        counts[divisor] += 1
    assert counts == {1: 30429, 149: 414, 211: 588, n: 8}
    return dict(sorted(counts.items()))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--order-bins", action="store_true",
                      help="only classify j=1 splits at scale 4096 by the order of 2")
    mode.add_argument("--certificate", action="store_true",
                      help="verify full cycles and the all-seed obstruction at N=31439")
    args = parser.parse_args()
    if args.certificate:
        print("all-seed gcds", complete_orbit_certificate())
    elif args.order_bins:
        for label, row in order_bins().items():
            print(label, row["all_sides"], row["j1"], row["j1_after_all"])
    else:
        for scale in (64, 256, 1024, 4096):
            print("dyadic", scale, *count(combinations(primerange(scale, 2 * scale), 2)),
                  flush=True)
        for bits in (20, 32, 48):
            print("separated", bits, *count(separated_pairs(bits)), flush=True)
