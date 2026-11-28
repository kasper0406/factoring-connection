"""Raw agreement of hidden factors with prime-local functions of the public product."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, isqrt, prod
from random import Random

from sympy import divisors, primerange

if __package__:
    from .check_root_scale_order_hints import prime_pairs
else:
    from check_root_scale_order_hints import prime_pairs


def public_anchors(n, primes, multipliers=None):
    """One local residue per prime. Random nonzero multipliers randomize each anchor."""
    if multipliers is None:
        return tuple(pow(n, 3, ell) for ell in primes)
    return tuple(c * pow(n, 2, ell) % ell for ell, c in zip(primes, multipliers, strict=True))


def agreement(prime, primes, anchors):
    """The factor is used to grade the received word after public construction."""
    return prod(ell for ell, a in zip(primes, anchors, strict=True) if prime % ell == a)


def exact_hit_probability(primes, hits):
    """Exact probability of precisely this many agreements with random unit anchors."""
    if not 0 <= hits <= len(primes):
        return Fraction(0)
    # In the numerator of the probability generating function a miss at ell
    # has weight ell-2 and a hit has weight one.
    counts = [1] + [0] * hits
    for ell in primes:
        counts = [count * (ell-2) + (counts[i-1] if i else 0)
                  for i, count in enumerate(counts)]
    return Fraction(counts[hits], prod(ell-1 for ell in primes))


def check_crt_words():
    """Exhaust the good unit pairs mod d and the finite dyadic counting bound."""
    checked = 0
    for primes in ((2, 3, 5), (5, 7, 11), (3, 7, 13)):
        for size in range(1, len(primes)+1):
            for subset in combinations(primes, size):
                modulus = prod(subset)
                units = [u for u in range(1, modulus) if gcd(u, modulus) == 1]
                # Fix each possible product residue. Exactly one first factor
                # equals the table word, then the partner is uniquely fixed.
                for u in units:
                    word = public_anchors(u, subset)
                    matching_pairs = [(x, u * pow(x, -1, modulus) % modulus) for x in units
                                      if agreement(x, subset, word) == modulus]
                    assert len(matching_pairs) == 1
                    x, y = matching_pairs[0]
                    assert x*y % modulus == u
                    checked += 1
    radius = 128
    primes = (3, 5, 7, 11)
    factors = tuple(map(int, primerange(radius, 2*radius)))
    # A product which crosses this threshold has a first crossing subproduct.
    # The largest prime in the subproduct is its last coordinate.
    threshold = isqrt(radius) + 1
    witnesses = [subset for size in range(1, len(primes)+1)
                 for subset in combinations(primes, size)
                 if prod(subset[:-1]) < threshold <= prod(subset)]
    total_hits = 0
    witness_count = 0
    integer_upper_bound = 0
    for subset in witnesses:
        modulus = prod(subset)
        count = sum(agreement(p, subset, public_anchors(p*q, subset)) == modulus
                    for p in factors for q in factors if p != q)
        upper_bound = prod(ell-1 for ell in subset) * ((radius+modulus-1)//modulus)**2
        assert count <= upper_bound
        witness_count += count
        integer_upper_bound += upper_bound
    for p in factors:
        for q in factors:
            if p != q:
                total_hits += agreement(p, primes, public_anchors(p*q, primes)) >= threshold
    assert 0 < total_hits <= witness_count <= integer_upper_bound
    return checked, total_hits, witness_count, integer_upper_bound


def check_exclusive_law():
    """Compare the disjoint-event formula with all public randomized tables."""
    primes = (5, 7, 11, 13)
    p, q = 101, 113
    counts = [0] * (len(primes)+1)
    threshold = 35
    tail_hits = 0
    for multipliers in product(*(range(1, ell) for ell in primes)):
        word = public_anchors(p*q, primes, multipliers)
        matching = [ell for ell, a in zip(primes, word, strict=True) if p % ell == a]
        counts[len(matching)] += 1
        tail_hits += prod(matching) >= threshold
    total = sum(counts)
    probabilities = [exact_hit_probability(primes, k) for k in range(len(counts))]
    assert probabilities == [Fraction(count, total) for count in counts]
    # The two smallest primes have product 35, so the product threshold
    # here is equivalent to at least two hits.
    assert sum(probabilities[2:]) == Fraction(tail_hits, total)
    for k, probability in enumerate(probabilities):
        # Each k-subset event is exclusive. Other primes must all miss.
        simple_bound = (comb(len(primes), k) * probabilities[0]
                        / (max(primes)-2)**k)
        assert probability >= simple_bound
    return total, probabilities[2], Fraction(tail_hits, total)


def check_square_remainders():
    identities = residue_classes = 0
    for radius in (32, 64, 128):
        primes = tuple(map(int, primerange(radius, 2*radius)))
        for p in primes:
            for q in primes:
                if p == q:
                    continue
                n = p*q
                for shift in (-3, -1, 0, 1, 4):
                    anchor = isqrt(n) + shift
                    if gcd(n, anchor) == 1:
                        assert gcd(n-anchor**2, p-anchor) == gcd(q-p, p-anchor)
                        identities += 1
        for v in range(1-radius, radius):
            if not v:
                continue
            band = 1 << (abs(v).bit_length()-1)
            points = range(max(radius, radius-v), min(2*radius, 2*radius-v))
            rounded = [r-isqrt(r*(r+v)) for r in points]
            for d in divisors(abs(v)):
                classes = Counter(y % d for y in rounded)
                # On this box the derivative magnitude is at least
                # T²/(32R²) and at most T²/(2R²). Rounding gives at most
                # (2+T²/(2Rd)) (1+32R²/T²) preimages per class.
                for count in classes.values():
                    assert count*d*band**2 <= 128*(radius*band**2+radius**2*d)
                residue_classes += len(classes)
    return identities, residue_classes


def benchmark(count=256):
    for bits in (16, 24, 32, 40):
        # Use B=(2b)^2, of order (log N)^2 on this dyadic cohort.
        cutoff = (2*bits)**2
        primes = tuple(map(int, primerange((cutoff+1)//2, cutoff+1)))
        assert primes[-1] < 1 << (bits-1)
        rng = Random(2026112901 + bits)
        fixed_hits = random_hits = 0
        for p, q in prime_pairs(bits, count, seed_offset=2_000_009):
            n = p*q
            fixed = public_anchors(n, primes)
            randomized = public_anchors(n, primes, [rng.randrange(1, ell) for ell in primes])
            fixed_hits += agreement(p, primes, fixed)**4 >= n
            random_hits += agreement(p, primes, randomized)**4 >= n
        yield bits, len(primes), count, fixed_hits, random_hits


if __name__ == '__main__':
    rows, actual, union, bound = check_crt_words()
    print(f'{rows} unit product residues checked; dyadic ordered prime pairs: '
          f'{actual} exceed the threshold, {union} witness incidences, integer upper bound {bound}.')
    tables, exact_two, tail = check_exclusive_law()
    print(f'{tables} public randomized tables: exactly two hit with probability {exact_two}; '
          f'product at least 35 with probability {tail}.')
    identities, classes = check_square_remainders()
    print(f'{identities} square-remainder gcd identities and {classes} rounded residue classes checked.')
    print('prime bits  band primes  moduli  cubic-table hits  randomized-table hits')
    for bits, size, count, fixed, randomized in benchmark():
        print(f'{bits:10} {size:12} {count:7} {fixed:16} {randomized:22}')
