"""Unplanted test of floor(N^alpha) as an exponent sharing a large part of r±1."""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt
from random import Random

from sympy import divisors, integer_nthroot, isprime, primerange

if __package__:
    from .weighted_order_hint import minimum_hint_bound, weighted_order_divisors
else:
    from weighted_order_hint import minimum_hint_bound, weighted_order_divisors


SCALES = ((Fraction(1, 2), Fraction(9, 25)),
          (Fraction(2, 3), Fraction(5, 12)),
          (Fraction(3, 4), Fraction(7, 16)))
SEED = 2026112804
HOLDOUT_OFFSET = 1_000_003


def public_exponent(n, alpha):
    """The constructor takes only N and a fixed public exponent scale."""
    if alpha == Fraction(1, 2):
        return isqrt(n)
    return int(integer_nthroot(n**alpha.numerator, alpha.denominator)[0])


def _prime(rng, radius):
    # Uniform odd candidates followed by rejection. All random draws here
    # are below 2^40, in SymPy's deterministic isprime range.
    while True:
        candidate = rng.randrange(radius+1, 2*radius, 2)
        if isprime(candidate):
            return candidate


def prime_pairs(bits, count, *, seed_offset=0, exclude=()):
    radius = 1 << (bits-1)
    rng = Random(SEED + bits + seed_offset)
    seen = set(exclude)
    target = len(seen) + count
    while len(seen) < target:
        p, q = sorted((_prime(rng, radius), _prime(rng, radius)))
        if p != q and (p, q) not in seen:
            seen.add((p, q))
            yield p, q


def grade(n, e, prime, offset, gamma):
    common = gcd(e, prime+offset)
    return common**gamma.denominator >= n**gamma.numerator


def benchmark(count=256):
    for bits in (16, 24, 32, 40):
        pairs = list(prime_pairs(bits, count))
        for alpha, gamma in SCALES:
            minus = plus = either = 0
            for p, q in pairs:
                n = p*q
                e = public_exponent(n, alpha)
                left = any(grade(n, e, r, -1, gamma) for r in (p, q))
                right = any(grade(n, e, r, 1, gamma) for r in (p, q))
                minus += left
                plus += right
                either += left or right
            yield bits, alpha, gamma, count, minus, plus, either


def check_dyadic_counts():
    rows = hits = maximum_multiplicity = 0
    for radius in (16, 32, 64, 128, 256):
        primes = list(map(int, primerange(radius, 2*radius)))
        for r in primes:
            for alpha, gamma in SCALES:
                # Count every integer q, including nonprimes, for the
                # multiplicity bound. Grading below uses distinct prime pairs.
                values = [public_exponent(r*q, alpha) for q in range(radius, 2*radius)]
                multiplicity = max(Counter(values).values())
                assert multiplicity <= 3
                maximum_multiplicity = max(maximum_multiplicity, multiplicity)
                low, high = values[0], values[-1]
                power = (radius*radius)**gamma.numerator
                floor, exact = integer_nthroot(power, gamma.denominator)
                minimum = int(floor) + int(not exact)
                for offset in (-1, 1):
                    large_divisors = [int(d) for d in divisors(r+offset) if d >= minimum]
                    marks = sum(grade(r*q, values[q-radius], r, offset, gamma) for q in primes if q != r)
                    multiple_count = sum(high//d - (low-1)//d for d in large_divisors)
                    width_bound = sum(1 + (high-low)//d for d in large_divisors)
                    assert marks <= multiplicity*multiple_count <= 3*width_bound
                    rows += 1
                    hits += marks
    return rows, hits, maximum_multiplicity


def public_sqrt_replay(n):
    """Return (factor, G, lattice calls) using public N only. Stop at the first verified split."""
    exponent = isqrt(n)
    common = gcd(n, exponent)
    if 1 < common < n:
        return common, None, 0
    low, high = isqrt((n-1)//2)+1, isqrt(2*n)
    radius = (high-low)//(2*n.bit_length()) + 1
    bound = minimum_hint_bound(n, exponent, 4, center=low+radius, radius=radius)
    if bound is None:
        return None, None, 0
    calls = 0
    for center in range(low+radius, high+radius+1, 2*radius+1):
        for sign in (1, -1):
            factors = weighted_order_divisors(n, exponent, bound, sign=sign, center=center,
                                             radius=radius, max_dimension=4)
            calls += 1
            if factors:
                return factors[0], bound, calls
    return None, bound, calls


def replay_cohort(bits, count=256, *, holdout=False):
    splits = []
    calls = 0
    original = set(prime_pairs(bits, 256)) if holdout else set()
    seed_offset = HOLDOUT_OFFSET if holdout else 0
    for p, q in prime_pairs(bits, count, seed_offset=seed_offset, exclude=original):
        assert (p, q) not in original
        n = p*q
        factor, bound, used = public_sqrt_replay(n)
        calls += used
        # Only now consult the synthetic factors, to grade the actual result
        # against complete coverage of the two signs in the public boxes.
        e = isqrt(n)
        expected = gcd(e, n) > 1 or (bound is not None and any(
            gcd(e, prime+offset) >= bound for prime in (p, q) for offset in (-1, 1)))
        assert (factor is not None) == expected
        if factor is not None:
            assert factor in (p, q)
            splits.append((n, factor, n//factor, bound))
    return splits, calls


def report_replay(bits, *, holdout=False):
    splits, calls = replay_cohort(bits, holdout=holdout)
    cohort = 'held-out' if holdout else 'original'
    print(f'{bits}-bit {cohort} primes: {len(splits)}/256 public splits, {calls} calls; '
          f'(N, factor, cofactor, G) = {splits}', flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--replay-bits', type=int, choices=(16, 24, 32, 40),
                        help='run only this cohort of the bounded four-dimensional replay')
    parser.add_argument('--holdout', action='store_true',
                        help='run the fixed replay on disjoint inputs with seed offset +1000003')
    args = parser.parse_args()
    if args.replay_bits is not None:
        report_replay(args.replay_bits, holdout=args.holdout)
        raise SystemExit(0)
    if args.holdout:
        for bits in (16, 24, 32, 40):
            report_replay(bits, holdout=True)
        raise SystemExit(0)
    for alpha, gamma in SCALES:
        assert 4*gamma**2 > alpha and alpha-Fraction(1, 2)-gamma < 0
    rows, hits, max_multiplicity = check_dyadic_counts()
    print(f'Exhausted {rows} fixed-prime, sign and scale rows, with {hits} qualifying partners counted '
          f'and maximum integer preimage size {max_multiplicity}.')
    print('prime bits  alpha  gamma  pairs  r-1 hits  r+1 hits  either')
    for bits, alpha, gamma, count, minus, plus, either in benchmark():
        print(f'{bits:10}  {str(alpha):>5}  {str(gamma):>5}  {count:5}  {minus:8}  {plus:8}  {either:6}')
    print('The bit-length window cover for the following replay was chosen after inspecting the two '
          'square-root criterion hits. All 1,024 original inputs are replayed.', flush=True)
    for bits in (16, 24, 32, 40):
        report_replay(bits)
