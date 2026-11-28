"""Check repeated rough factors in short rational rows, from both directions.

Compare all unit-square targets for distinct primes in [11,22]. The direct
search enumerates short representatives for each target with the Euclidean
shell routine. The determinant search enumerates primitive signed fractions
in the enclosing height range, without using the target shell. Two distinct
fractions may represent the same target modulo N only if N divides their
nonzero determinant. Both passes count unordered pairs of distinct fractions,
even when a pair has several common primes. The factor base is {2,3}, with
the sign available separately, and the listed thresholds are auxiliary.

Run directly or with -m from the repository root. The module uses the same
gmpy2 dependency as the existing shell benchmark.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from math import gcd, prod

if __package__:
    from . import check_rational_height_shell as shell
else:
    import check_rational_height_shell as shell


THRESHOLDS = (3, 7, 13)


@lru_cache(maxsize=None)
def factor_support(value):
    """Prime divisors and odd-valuation prime divisors of a positive integer."""
    primes, odd = set(), set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            primes.add(divisor)
            exponent = 0
            while value % divisor == 0:
                value //= divisor
                exponent += 1
            if exponent % 2:
                odd.add(divisor)
        divisor += 1
    if value > 1:
        primes.add(value)
        odd.add(value)
    return frozenset(primes), frozenset(odd)


def height(fraction):
    u, v = fraction
    return abs(u) * v


def outside_base_mask(odd_primes):
    return sum(1 << p for p in odd_primes if p > 3)


def parity_basis(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            i = vector.bit_length() - 1
            if i not in pivots:
                pivots[i] = vector
                break
            vector ^= pivots[i]
    return pivots


def has_odd_dependency(rows):
    """The extra bit equals one precisely when the row count is odd."""
    odd_sets = (factor_support(height(fraction))[1] for fraction in rows)
    return 0 in parity_basis(1 | outside_base_mask(odd) for odd in odd_sets)


def from_targets(cohort):
    events = {y: set() for y in THRESHOLDS}
    small_rows = {y: set() for y in THRESHOLDS}
    odd_hits = set()
    for n, targets in cohort.items():
        for t in targets:
            candidates = shell.representations(n, t)[0]
            rows = sorted(f for f in candidates if gcd(f[1], n) == 1)
            if has_odd_dependency(rows):
                odd_hits.add((n, t))
            for f in rows:
                max_odd_prime = max(factor_support(height(f))[1], default=1)
                for y in THRESHOLDS:
                    if max_odd_prime <= y:
                        small_rows[y].add((n, t))
            for a, b in combinations(rows, 2):
                common = gcd(height(a), height(b))
                largest = max(factor_support(common)[0], default=1)
                for y in THRESHOLDS:
                    if largest > y:
                        events[y].add((n, t, a, b))
    return events, small_rows, odd_hits


def from_fractions(cohort):
    maximum = max(cohort)
    groups = defaultdict(list)
    for v in range(1, maximum + 1):
        for u in range(-(maximum // v), maximum // v + 1):
            if u and gcd(u, v) == 1:
                for p in factor_support(abs(u) * v)[0]:
                    if p > min(THRESHOLDS):
                        groups[p].append((u, v))

    events = {y: set() for y in THRESHOLDS}
    for p, group in groups.items():
        for a, b in combinations(sorted(group), 2):
            det = a[0] * b[1] - b[0] * a[1]
            assert det
            for n, targets in cohort.items():
                if det % n or max(height(a), height(b)) > n or gcd(a[1] * b[1], n) != 1:
                    continue
                t = a[0] * pow(a[1], -1, n) % n
                assert (b[0] - t * b[1]) % n == 0
                if t in targets:
                    for y in THRESHOLDS:
                        if p > y:
                            events[y].add((n, t, a, b))
    return events


def check_uncovered_batches(rows, y):
    """Check every subset outside the three events for a zero rough parity."""
    batches = {(0, 0)}  # row mask and cumulative outside-base parity
    for index, (_, u, _, lift, (primes, odd)) in enumerate(rows):
        if max(odd, default=1) <= y:
            continue
        forbidden = 0
        for earlier, (_, a, _, other_lift, (other_primes, _)) in enumerate(rows[:index]):
            if lift != other_lift:
                clash = max(primes & other_primes, default=1) > y
            else:
                numerator_odd = factor_support(abs(u))[1] | factor_support(abs(a))[1]
                clash = max(numerator_odd, default=1) <= y
            if clash:
                forbidden |= 1 << earlier
        parity = outside_base_mask(odd)
        new = {(mask | (1 << index), old ^ parity) for mask, old in batches if not mask & forbidden}
        assert all(value for _, value in new)
        batches.update(new)
    return len(batches) - 1


def conic_check(cohort):
    """Separate exact rational lifts from congruences between distinct lifts."""
    counts = {y: defaultdict(int) for y in THRESHOLDS}
    examples = {}
    for n, squares in cohort.items():
        for t in sorted(s for s in squares if (s + 1) % n in squares):
            rows = []
            for shift in (0, 1):
                for u, v in shell.representations(n, t + shift)[0]:
                    if gcd(v, n) == 1:
                        lift = (u - shift * v, v)
                        assert gcd(*lift) == 1 and v > 0
                        rows.append((shift, u, v, lift, factor_support(abs(u) * v)))
            vectors = (outside_base_mask(odd) for _, _, _, _, (_, odd) in rows)
            dependent = len(parity_basis(vectors)) < len(rows)
            single_limit = min(max(odd, default=1) for _, _, _, _, (_, odd) in rows)
            repeated_limit = 1
            pair_limits = []
            for first, second in combinations(rows, 2):
                i, u, v, lift, (primes, odd) = first
                j, a, b, other_lift, (other_primes, other_odd) = second
                determinant = u * b - a * v + (j - i) * v * b
                assert determinant % n == 0 and (determinant == 0) == (lift == other_lift)
                if determinant:
                    repeated_limit = max(repeated_limit, max(primes & other_primes, default=1))
                else:
                    assert i != j and b == v and gcd(u, a) == 1
                    numerator_odd = factor_support(abs(u))[1] | factor_support(abs(a))[1]
                    pair_limits.append(max(numerator_odd, default=1))
                    # Within this two-row batch there is no distinct-lift event.
                    if outside_base_mask(odd) == outside_base_mask(other_odd):
                        for y in THRESHOLDS:
                            if min(max(odd, default=1), max(other_odd, default=1)) > y:
                                assert max(numerator_odd, default=1) <= y
                                counts[y]["pair_only_witnesses"] += 1
                                examples.setdefault(y, (n, t, (i, u, v), (j, a, b)))
            for y in THRESHOLDS:
                stats = counts[y]
                singleton = single_limit <= y
                paired = any(limit <= y for limit in pair_limits)
                repeated = repeated_limit > y
                assert not dependent or singleton or paired or repeated
                stats["uncovered_batches"] += check_uncovered_batches(rows, y)
                stats["targets"] += 1
                stats["dependent"] += dependent
                stats["singleton"] += singleton
                stats["paired"] += paired
                stats["repeated"] += repeated
                stats["paired_only_targets"] += dependent and paired and not singleton and not repeated
    print("Y,conic_targets,dependent,singleton,paired,repeated,paired_only_targets,pair_only_witnesses,uncovered_batches")
    for y, stats in counts.items():
        assert stats["dependent"] and stats["repeated"] and stats["paired"] and stats["uncovered_batches"]
        print(y, *(stats[k] for k in ("targets", "dependent", "singleton", "paired", "repeated",
                                       "paired_only_targets", "pair_only_witnesses", "uncovered_batches")), sep=",")
    assert 3 in examples
    print("A same-lift pair requiring the paired exception within its batch:", examples[3])


def main():
    primes = [p for p in shell.primes_upto(22) if p >= 11]
    cohort = {}
    for p, q in combinations(primes, 2):
        n = p * q
        cohort[n] = {x * x % n for x in range(1, n) if gcd(x, n) == 1}
        assert len(cohort[n]) == (p - 1) * (q - 1) // 4
    events, small_rows, odd_hits = from_targets(cohort)
    assert events == from_fractions(cohort)
    total = sum(len(targets) for targets in cohort.values())
    print("primes", primes, "moduli", sorted(cohort), "targets", total, "odd_hits", len(odd_hits))
    print("Y,pairs,repeat_targets,small_row_targets,union_targets,odd_hits_outside_union")
    for y in THRESHOLDS:
        repeated = {(n, t) for n, t, _, _ in events[y]}
        union = small_rows[y] | repeated
        assert odd_hits <= union
        print(y, len(events[y]), len(repeated), len(small_rows[y]), len(union),
              len(odd_hits - union), sep=",")

    conic_check(cohort)

    # An odd three-row relation can be genuine even when every singleton fails.
    n, t = 247, 23
    rows = [(-5, 32), (-1, 204), (-17, 10)]
    direct = {f for f in shell.representations(n, t)[0] if gcd(f[1], n) == 1}
    assert t in cohort[n] and set(rows) <= direct
    assert has_odd_dependency(rows) and all(not has_odd_dependency([f]) for f in direct)
    assert prod(u for u, _ in rows) * 768 == -prod(v for _, v in rows)


if __name__ == "__main__":
    main()
