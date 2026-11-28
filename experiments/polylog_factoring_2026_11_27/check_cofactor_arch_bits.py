"""Compare Archimedean bit envelopes with all intersected unit squares."""

from fractions import Fraction
from itertools import combinations
from math import isqrt
from random import Random

from sympy import nextprime

if __package__:
    from .cofactor_arch_bits import arch_outcomes, floor_lift
else:
    from cofactor_arch_bits import arch_outcomes, floor_lift


def floor_words(n, lower=0, upper=None):
    """All integer cells on [sqrt(n/2),sqrt(n)] times its reciprocal."""
    upper = n if upper is None else upper
    pairs = set()
    for x in range(max(lower, isqrt(n//2)), min(upper, isqrt(n))+1):
        for y in range(isqrt(n), isqrt(2*n)+1):
            # In z=u² the lower edge for y and the upper edge for x are open.
            open_low = Fraction(n*n, (y+1)**2)
            open_high = (x+1)**2
            lo = max(Fraction(n, 2), lower*lower, x*x, open_low)
            hi = min(n, upper*upper, open_high, Fraction(n*n, y*y))
            if lo < hi or (lo == hi and open_low < lo < open_high):
                pairs.add((x, y))
    return pairs


def verify_envelopes():
    cases = 0
    for n in [*range(1, 150), 225, 242, 675, 9409, 16383]:
        labels = floor_words(n)
        for mask in range(1 << isqrt(2*n).bit_length()):
            expected = frozenset(bool((x & mask) and (y & mask)) for x, y in labels)
            assert arch_outcomes(n, mask) == expected, (n, mask, expected)
            cases += 1
    return cases


def verify_restricted():
    rng = Random(91)
    count = 0
    for n in [*range(1, 260), 2541, 10007]:
        for _ in range(5):
            lower, upper = sorted(rng.randrange(isqrt(n)+4) for _ in range(2))
            labels = floor_words(n, lower, upper)
            for mask in range(1 << isqrt(2*n).bit_length()):
                expected = frozenset(bool((x & mask) and (y & mask)) for x, y in labels)
                assert arch_outcomes(n, mask, lower=lower, upper=upper) == expected
                count += 1
    return count


def greedy_savings(p, q):
    n = p*q
    prefix = selected = 0
    forced_yes = forced_no = ambiguous = 0
    for j in range(isqrt(2*n).bit_length()-1, 0, -1):
        bit = 1 << j
        mask = selected | bit
        choices = arch_outcomes(n, mask, lower=prefix, upper=prefix+2*bit-1)
        answer = bool((p & mask) and (q & mask))
        assert answer in choices and answer == bool(p & bit)
        forced_yes += choices == {True}
        forced_no += choices == {False}
        ambiguous += len(choices) == 2
        if answer:
            prefix |= bit
        else:
            selected = mask
    assert prefix | 1 == p
    return forced_yes, forced_no, ambiguous


def verify_greedy():
    rng = Random(7096)
    rows = []
    for bits in (16, 32, 64):
        results = []
        while len(results) < 60:
            p, q = sorted(int(nextprime(rng.randrange(2**(bits-1), 2**bits-64)))
                          for _ in range(2))
            if p != q:
                results.append(greedy_savings(p, q))
        rows.append((bits, *(sum(column) for column in zip(*results))))
    return rows, greedy_savings(2**127-1, 162*2**120+1)


def verify_lifts():
    rng = Random(71067)
    count = exact = 0
    for k in (1, 2, 3, 4, 8, 12, 24, 40):
        modulus = 1 << k
        for _ in range(60):
            n = rng.randrange(4096*modulus**6, 8192*modulus**6) | 1
            a = rng.randrange(1, modulus, 2)
            b = n*pow(a, -1, modulus) % modulus
            x, y = floor_lift(n, modulus, a, b)
            assert x % modulus == a and y % modulus == b
            assert 0 < x <= y <= 2*x and n//x == y
            assert 0 <= n-x*y < x and (n-x*y) % modulus == 0
            exact += n == x*y
            count += 1
    assert floor_lift(77, 16, 7, 11) == (7, 11)
    assert floor_lift(77, 16, 1, 13) is None
    return count, exact


def check_large():
    p, q = 2**127-1, 162*2**120+1
    n = p*q
    # Active positions start at j=1. Positions above the largest possible
    # factor are structurally zero and are excluded from these counts.
    bits = [1 << j for j in range(1, isqrt(2*n).bit_length())]
    summary = []
    for length in (1, 2):
        forced_yes = forced_no = mixed = max_cuts = 0
        for active in combinations(bits, length):
            mask = sum(active)
            stats = {}
            outcome = arch_outcomes(n, mask, stats=stats)
            assert bool((p & mask) and (q & mask)) in outcome
            max_cuts = max(max_cuts, stats["cut_points"])
            forced_yes += outcome == {True}
            forced_no += outcome == {False}
            mixed += len(outcome) == 2
        summary.append((length, forced_yes, forced_no, mixed, max_cuts))
    return summary


if __name__ == "__main__":
    print("Small exact envelopes checked:", verify_envelopes())
    print("Restricted-interval exact envelopes checked:", verify_restricted())
    print("Greedy prefix tests (prime bits, forced yes, forced no, either) and 255-bit fixture:", verify_greedy())
    print("Certified rounded lifts and exact products among them:", verify_lifts())
    print("255-bit sample (query size, forced yes, forced no, either, maximum cut points):", check_large())
