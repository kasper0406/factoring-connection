"""Prime support of nonempty binary prefixes and remainders, graded on known factors."""

from collections import Counter
from functools import lru_cache
from math import gcd, isqrt, prod

from sympy import divisors, factorint, primerange

if __package__:
    from .check_root_scale_order_hints import prime_pairs
else:
    from check_root_scale_order_hints import prime_pairs


def binary_halves(n):
    """(H, floor(N/H), N mod H) for every split into two nonempty bit strings."""
    return [(1 << i, n >> i, n % (1 << i)) for i in range(1, n.bit_length())]


def partner_bounds(p, radius, h, m):
    """Upper bounds for integer partners giving positive multiples of m."""
    assert gcd(p, h*m) == 1
    first = p*radius//h
    last = p*(2*radius-1)//h
    periods, leftover = divmod(radius, h*m)
    quotient = min(radius, periods*h+min(leftover, h), (last//m)*((h+p-1)//p))
    remainder = min(radius, ((h-1)//m)*((radius+h-1)//h),
                    radius//m + last-first+1)
    return quotient, remainder


def check_conditional_counts():
    cases = integer_slopes = 0
    for radius in (3, 16, 32, 64, 128):
        primes = tuple(map(int, primerange(radius, 2*radius)))
        for p in primes:
            moduli = [m for m in divisors(p*p-1) if m <= 2*radius+1]
            for k in range(1, (4*radius**2).bit_length()):
                h = 1 << k
                digits = [divmod(p*q, h) for q in range(radius, 2*radius)]
                for m in moduli:
                    high = sum(a > 0 and a % m == 0 for a, b in digits)
                    low = sum(b > 0 and b % m == 0 for a, b in digits)
                    hi_bound, lo_bound = partner_bounds(p, radius, h, m)
                    assert high <= hi_bound and low <= lo_bound
                    assert m * max(hi_bound, lo_bound)**2 <= (8*radius)**2
                    if 2*h >= radius:
                        assert m*hi_bound <= 8*radius
                    if h <= 2*radius:
                        assert m*lo_bound <= 8*radius
                    if m >= 2:
                        high_slopes = []
                        low_slopes = []
                        for q in primes:
                            a, b = digits[q-radius]
                            if a > 0 and a % m == 0:
                                k = a//m
                                assert gcd(k, q) == 1
                                high_slopes.append((k, q))
                            if b > 0 and b % m == 0:
                                assert m % 2 and (a*h-p*q) % m == 0
                                if a % q == 0:
                                    assert m == q
                                    integer_slopes += 1
                                low_slopes.append((a, q))
                        for i, (a, q) in enumerate(low_slopes):
                            for a1, q1 in low_slopes[:i]:
                                determinant = a*q1-a1*q
                                assert determinant != 0 and determinant % m == 0
                        assert len(high_slopes)*m <= m+4*radius
                        assert len(low_slopes)*m <= m+4*radius
                    cases += 1
    assert integer_slopes > 0
    return cases, integer_slopes


def check_composite_partners():
    radius, p, h, b = 1 << 20, 1_863_671, 4096, 455
    m = 2*3*19*227
    assert p*b == h*(b*b-1)+1 and (b*b-1) % m == 0 and (p*p-1) % m == 0
    partners = range((radius+b-1)//b, min(h, (2*radius-1)//b+1))
    for j in partners:
        q = b*j
        assert radius <= q < 2*radius and (p*q//h) % m == 0
    assert len(partners)*m > m+4*radius
    return len(partners), m


@lru_cache(maxsize=2048)
def order_powers(p):
    # Factoring the two known prime neighbors is only for external grading.
    factors = Counter(factorint(p-1))
    factors.update(factorint(p+1))
    return {int(base): int(base**exponent) for base, exponent in factors.items()}


def supporting_bases(digit, powers):
    assert digit > 0
    return {base for base in powers if digit % base == 0}


def agreement(powers, bases):
    return prod(powers[base] for base in bases)


def check_pool_witnesses():
    cases = qualifying = small = large = 0
    for radius in (16, 32, 64, 128):
        primes = tuple(map(int, primerange(radius, 2*radius)))
        for p in primes:
            powers = order_powers(p)
            for q in primes:
                if q == p:
                    continue
                n = p*q
                digits = [digit for _, high, low in binary_halves(n) for digit in (high, low)]
                bases = {base for digit in digits for base in supporting_bases(digit, powers)}
                supported = agreement(powers, bases)
                # Independent reference for radical-only coverage. Raising
                # the product retains the full valuation of every shared prime.
                bound = (p*p-1).bit_length()
                assert supported == gcd(pow(prod(digits), bound, p*p-1), p*p-1)
                cases += 1
                if supported**4 < n:
                    continue
                qualifying += 1
                if any(base > 7 for base in bases):
                    large += 1
                    continue
                d = 1
                previous = 1
                for base in sorted(bases):
                    power = powers[base]
                    while d*d < radius and power > 1:
                        previous = d
                        d *= base
                        power //= base
                assert previous**2 < radius <= d*d
                assert d < 7*(isqrt(radius)+1) and (p*p-1) % d == 0
                small += 1
    assert qualifying == small+large and min(small, large) > 0
    return cases, qualifying, small, large


def benchmark(count=128):
    for bits in (16, 32, 48, 64):
        totals = [0, 0, 0, 0]
        for p, q in prime_pairs(bits, count, seed_offset=4_000_037):
            n = p*q
            splits = binary_halves(n)
            middle = 1 << (n.bit_length()//2)
            best = [1, 1, 1, 1]
            for r in (p, q):
                powers = order_powers(r)
                pool = set()
                for h, high, low in splits:
                    for digit in (high, low):
                        bases = supporting_bases(digit, powers)
                        mass = agreement(powers, bases)
                        if h == middle:
                            best[0] = max(best[0], mass)
                        best[1] = max(best[1], mass)
                        pool.update(bases)
                largest_two = sorted((powers[ell] for ell in pool), reverse=True)[:2]
                best[2] = max(best[2], prod(largest_two))
                best[3] = max(best[3], agreement(powers, pool))
            for i, mass in enumerate(best):
                totals[i] += mass**4 >= n
        yield bits, count, totals


if __name__ == '__main__':
    cases, exceptions = check_conditional_counts()
    print(f'{cases} conditional counts checked for both halves, with all integer partners and '
          f'again with prime partners. Found {exceptions} integral remainder slopes.')
    composite, m = check_composite_partners()
    print(f'{composite} composite partners give divisible prefixes in the comparison at m={m}.')
    rows, qualifying, small, large = check_pool_witnesses()
    print(f'{rows} ordered prime pairs checked. {qualifying} qualify. '
          f'{small} use the smooth witness and {large} have a supported prime above seven.')
    print('prime bits  moduli  middle half  any one half  best two pooled bases  full pool')
    for bits, count, totals in benchmark():
        print(f'{bits:10} {count:7} {totals[0]:12} {totals[1]:13} {totals[2]:22} {totals[3]:10}',
              flush=True)
