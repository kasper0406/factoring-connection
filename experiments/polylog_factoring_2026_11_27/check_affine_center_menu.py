"""Count the balanced semiprimes a prelisted menu of affine windows can split.

The center is A(N)/B(N) in Z/N: a proper denominator gcd is a split and
a denominator divisible by N is unusable. On a unit denominator the constant
terms determine the center. The candidate-prime budget uses only those terms.
"""

from collections import Counter
from itertools import combinations
from math import comb, gcd

from sympy import primerange

if __package__:
    from .affine_prime_window import prime_window_divisor, window_radius
else:
    from affine_prime_window import prime_window_divisor, window_radius


def menu(scale):
    large = 2**(scale.bit_length()**2) + 17
    return (
        ((2*scale+1,), (2,)),
        ((1,), ((scale+3)*(scale+7),)),
        ((-1,), (1,)),
        ((-2,), (1,)),
        ((large+2*scale-1,), (large,)),
        ((2*scale+1, -5, large), (2, 0, -7)),
        ((1, 3), (0, 7)),  # denominator is always divisible by N
    )


def estimate_menu(primes, scale, radius, formulas):
    """Exact prime cover and a cardinality bound from coefficient sizes alone."""
    covered = set()
    budget = 0
    for numerator, denominator in formulas:
        a, b = numerator[0], denominator[0]
        terms = [b] + [a+b*y for y in range(-radius, radius+1)]
        for term in filter(None, terms):
            magnitude = abs(term)
            local = {r for r in primes if magnitude % r == 0}
            # k distinct prime divisors in [R,2R] require R^k <= |term|.
            count, power = 0, 1
            while power*scale <= magnitude:
                count += 1
                power *= scale
            assert len(local) <= count
            budget += count
            covered.update(local)
    return covered, budget


def evaluate(coefficients, n):
    result = 0
    for coefficient in reversed(coefficients):
        result = result*n + coefficient
    return result


def observe(n, numerator, denominator):
    divisor = gcd(denominator, n)
    if divisor != 1:
        return ('denominator' if divisor < n else 'undefined'), divisor
    x = numerator*pow(denominator, -1, n) % n
    return 'window', prime_window_divisor(n, x)


def check_band(scale):
    primes = list(primerange(scale, 2*scale+1))
    formulas = menu(scale)
    upper = window_radius(4*scale**2)
    assert 2*upper < scale
    covered, budget = estimate_menu(primes, scale, upper, formulas)
    hits, counts = set(), Counter()
    one_covered, one_budget = estimate_menu(primes, scale, upper, formulas[:1])
    assert 0 < len(one_covered) <= one_budget < len(primes)//2
    for p, q in combinations(primes, 2):
        n = p*q
        radius = window_radius(n)
        outcomes = []
        for numerator, denominator in formulas:
            a, b = numerator[0], denominator[0]
            mode, value = observe(n, evaluate(numerator, n), evaluate(denominator, n))
            assert (mode, value) == observe(n, a, b)
            if mode == 'window':
                offsets = range(-radius, radius+1)
                found = [r for r in (p, q) if any((a+b*y) % r == 0 for y in offsets)]
                assert value == (p if p in found else 1) * (q if q in found else 1)
                if b and a % b == 0 and abs(a//b) <= upper:
                    # An enclosing-window integer center either hits both or neither.
                    assert value == (n if abs(a//b) <= radius else 1)
                    counts['zero term, both' if value == n else 'zero term, neither'] += 1
            proper = 1 < value < n
            if proper:
                assert value in covered and value in (p, q)
                hits.add(n)
                counts[mode+' splits'] += 1
            outcomes.append(proper)
            counts['formula/pair checks'] += 1
        if outcomes[0]:
            assert p in one_covered or q in one_covered
            counts['first-center splits'] += 1
        # A modulus-dependent choice from the menu is bounded by the same union.
        index = (n ^ (n >> 3)) % len(formulas)
        if outcomes[index]:
            assert n in hits
            counts['selected splits'] += 1
    size = len(primes)
    touched = comb(size, 2)-comb(size-len(covered), 2)
    assert len(hits) <= touched <= len(covered)*(size-1) <= budget*(size-1)
    assert counts['first-center splits'] <= one_budget*(size-1) < comb(size, 2)
    counts['first-center prime budget'] = one_budget
    assert 0 < len(hits) < comb(size, 2) and counts['window splits'] and counts['denominator splits']
    return size, len(covered), budget, len(hits), touched, counts


if __name__ == '__main__':
    for scale in (64, 192):
        size, touched_primes, budget, split_pairs, touched_pairs, counts = check_band(scale)
        print(f'R={scale}: {size} primes, {split_pairs}/{comb(size, 2)} split pairs; '
              f'{touched_primes} candidate primes, {touched_pairs} touched pairs, prime budget {budget}')
        print(dict(counts))
