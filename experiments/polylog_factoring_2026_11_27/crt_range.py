"""Evaluate the retained reciprocal-range slopes with a supplied factorization.

For t distinct prime factors and centered slope cutoff B, prepare and evaluate
O((B + 1) 2**ceil(t/2)) half configurations, up to polynomial factors in
the input bit length, numerical precision and number of denominator intervals.
Local square roots use expected polynomial time. Factoring is not included.
The CLI checks the supplied prime bases with SymPy's isprime.
For large prime bases it is a probable-prime check. The complexity and the
local-root formulas assume the supplied bases are prime.
"""

from argparse import ArgumentParser
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass
from math import gcd, prod
from random import Random
from secrets import randbelow
from time import perf_counter

from mpmath import mp
from sympy import isprime

from reciprocal_range import digit_intersection, quotient_interval, scan_spectrum, value_at


def prime_square_root(target, prime, work):
    """Return one square root modulo a supplied odd prime, or None."""
    target %= prime
    work['local_prime_tests'] += 1
    if target == 1:
        return 1
    if pow(target, (prime - 1) // 2, prime) != 1:
        return None
    if prime % 4 == 3:
        return pow(target, (prime + 1) // 4, prime)

    odd_part = prime - 1
    order_power = 0
    while odd_part % 2 == 0:
        odd_part //= 2
        order_power += 1

    while True:
        candidate = randbelow(prime - 1) + 1
        work['local_nonresidue_probes'] += 1
        if pow(candidate, (prime - 1) // 2, prime) == prime - 1:
            break

    generator = pow(candidate, odd_part, prime)
    root = pow(target, (odd_part + 1) // 2, prime)
    residue = pow(target, odd_part, prime)
    while residue != 1:
        power = residue * residue % prime
        index = 1
        while power != 1:
            power = power * power % prime
            index += 1
        adjustment = pow(generator, 1 << (order_power - index - 1), prime)
        root = root * adjustment % prime
        generator = adjustment * adjustment % prime
        residue = residue * generator % prime
        order_power = index
    return root


def inverse_square_root(target, prime, modulus, work):
    """Lift a unit square root to modulus p^(2e), then invert it."""
    root = prime_square_root(target, prime, work)
    if root is None:
        return None
    precision = prime
    while precision < modulus:
        precision = min(precision * precision, modulus)
        root = (root - (root * root - target) * pow(2 * root, -1, precision)) % precision
        work['hensel_steps'] += 1
    return pow(root, -1, modulus)


def half_configurations(projections, n, scale):
    modulus = n * n
    roots = [0]
    for projection in projections:
        roots = [(root + sign * projection) % modulus
                 for root in roots for sign in (-1, 1)]
    return [(scale * root % n, root) for root in roots]


def lift_segments(n, start, width, ranges):
    """Split the cyclic lower window into regions of constant upper lifts."""
    start %= n
    if start + width <= n:
        window = [(start, start + width)]
    else:
        window = [(0, start + width - n), (start, n)]
    segments = []
    for low, high in ranges:
        if not 0 <= high - low <= n * n:
            raise ValueError('Each denominator interval must have length in [0, N²]')
        for left, right in window:
            cuts = [left] + sorted({bound % n for bound in (low, high)
                                    if left < bound % n < right}) + [right]
            for first_digit, stop_digit in zip(cuts, cuts[1:]):
                first, count = digit_intersection(n, first_digit, low, high)
                if count:
                    segments.append((first_digit, stop_digit, first, count))
    return segments


def phase(value, modulus):
    centered = (value + modulus // 2) % modulus - modulus // 2
    return mp.expjpi(2 * mp.mpf(centered) / modulus)


def lift_weight(n, slope, first, count):
    if slope == 0:
        return mp.mpf(count) / n
    tail = (slope * count + n // 2) % n - n // 2
    if tail == 0:
        return mp.mpf(0)
    angle = mp.mpf(slope) / n
    scaled_sine = mp.sinpi(angle) / angle
    magnitude = mp.sinpi(mp.mpf(tail) / n) / (slope * scaled_sine)
    return magnitude * phase(tail - slope + 2 * slope * first, 2 * n)


@dataclass
class SlopeTable:
    slope: int
    left: list
    right: list
    right_digits: list


class CrtRange:
    """Cache the half-root ordering across upper digits and interval queries."""

    def __init__(self, factors, scale=1, cutoff=8, digits=25):
        if (not factors or len({p for p, _ in factors}) != len(factors)
                or any(p < 3 or p % 2 == 0 or e < 1 or not isprime(p)
                       for p, e in factors)):
            raise ValueError('Supply distinct odd prime bases with positive exponents')
        self.n = prod(p**e for p, e in factors)
        if gcd(scale, self.n) != 1 or cutoff < 0 or digits < 1:
            raise ValueError('Expected a unit scale, nonnegative cutoff and positive precision')
        self.scale = scale % self.n
        self.modulus = self.n * self.n
        cutoff = min(cutoff, self.n // 2)
        self.digits = digits + len(factors) + len(str(2 * cutoff + 1)) + 10
        self.tables = []
        self.setup_work = Counter()
        local = []
        for prime, exponent in factors:
            mod = prime**(2 * exponent)
            complement = self.modulus // mod
            local.append((prime, mod, complement * pow(complement, -1, mod)))

        for slope in range(-cutoff, cutoff + 1):
            self.setup_work['slope_candidates'] += 1
            target = 1 - slope
            if gcd(target, self.n) != 1:
                self.setup_work['nonunit_slopes'] += 1
                continue
            projections = []
            for prime, mod, idempotent in local:
                root = inverse_square_root(target, prime, mod, self.setup_work)
                if root is None:
                    break
                projections.append(root * idempotent % self.modulus)
            if len(projections) != len(local):
                self.setup_work['nonsquare_slopes'] += 1
                continue
            middle = (len(local) + 1) // 2
            left = half_configurations(projections[:middle], self.n, self.scale)
            right = sorted(half_configurations(projections[middle:], self.n, self.scale))
            self.setup_work['left_configurations'] += len(left)
            self.setup_work['right_configurations'] += len(right)
            self.tables.append(SlopeTable(slope, left, right, [a for a, _ in right]))

    def value(self, digit, start, width, ranges):
        """Sum normalized retained terms, adding the supplied half-open intervals."""
        n = self.n
        if not 0 < width <= n:
            raise ValueError('Window width must be in [1,N]')
        segments = lift_segments(n, start, width, ranges)
        upper_scale = self.scale + n * (digit % n)
        work = Counter(lift_segments=len(segments))
        with mp.workdps(self.digits + len(str(max(1, len(segments))))):
            total = mp.mpc(0)
            for table in self.tables:
                total += self._slope_value(table, segments, upper_scale, work)
        return total, work

    def _slope_value(self, table, segments, upper_scale, work):
        n, modulus, b = self.n, self.modulus, table.slope
        weights = [(low, high, lift_weight(n, b, first, count))
                   for low, high, first, count in segments]
        weights = [(low, high, weight) for low, high, weight in weights if weight]
        if not weights:
            return mp.mpc(0)

        # For r=a+Nh and (1−b)r²=C² mod N², f(r)=(2−b)r and
        # f(r)−f(a)=Nbh. Thus f(a)=(2−2b)r+ba mod N².
        multiplier = (2 - 2 * b) * upper_scale
        prefix = [mp.mpc(0)]
        for digit, root in table.right:
            prefix.append(prefix[-1] + phase(multiplier * root + b * digit, modulus))
        work['right_phases'] += len(table.right)
        carry_phase = phase(-b, n)

        total = mp.mpc(0)
        for digit, root in table.left:
            subtotal = mp.mpc(0)
            for low, high, lift in weights:
                for shift, factor in ((-digit, 1), (n - digit, carry_phase)):
                    first = bisect_left(table.right_digits, low + shift)
                    stop = bisect_left(table.right_digits, high + shift)
                    work['range_searches'] += 1
                    if first < stop:
                        subtotal += (prefix[stop] - prefix[first]) * factor * lift
                        work['nonempty_ranges'] += 1
            if subtotal:
                total += phase(multiplier * root + b * digit, modulus) * subtotal
                work['left_phases'] += 1
        return total


def check():
    rng = Random(261127)
    root_checks = 0
    for prime in (3, 5, 7, 13, 17, 41, 97):
        modulus = prime**6
        squares = {x*x % prime for x in range(1, prime)}
        for target in range(-prime, prime):
            if target % prime == 0:
                continue
            root = inverse_square_root(target, prime, modulus, Counter())
            if target % prime in squares:
                assert target * root**2 % modulus == 1
            else:
                assert root is None
            root_checks += 1
    fixtures = [((3, 1), (5, 1)), ((3, 2), (5, 1)), ((7, 1), (11, 1)),
                ((3, 1), (5, 1), (7, 1)), ((3, 2), (5, 1), (7, 1)),
                ((3, 2), (5, 1), (7, 1), (11, 1)), ((101, 1), (109, 1))]
    checks = 0
    for factors in fixtures:
        n = prod(p**e for p, e in factors)
        modulus = n * n
        left, right = quotient_interval(modulus, 2)
        scale = rng.choice([a for a in range(1, n) if gcd(a, n) == 1])
        queries = [(n - 4, min(29, n), [(left, right), (1 - right, 1 - left)]),
                   (0, n, [(-modulus // 2 + 1, modulus // 2 + 1)]),
                   (n - 1, 2, [(-n - 2, n + 7)]),
                   (rng.randrange(n), min(83, n), [(-n * 4 + 2, n * 2 + 5)]),
                   (n - 1, 2, [(0, 0)])]
        for cutoff in (0, 4, 12):
            evaluator = CrtRange(factors, scale, cutoff)
            for start, width, ranges in queries:
                bands = []
                for low, high in ranges:
                    spectrum, _ = scan_spectrum(n, scale, start, width, low, high)
                    bands.append({k: (b, v) for k, (b, v) in spectrum.items()
                                  if min(b, n - b) <= cutoff})
                for digit in (0, 1, n // 3, n - 1):
                    expected = sum(value_at(n, digit, band) for band in bands)
                    observed, _ = evaluator.value(digit, start, width, ranges)
                    assert abs(complex(observed) - expected) < 2e-12, (
                        n, cutoff, start, width, ranges, digit, observed, expected,
                    )
                    checks += 1
    print(f'{root_checks} local root checks and {checks} retained scan-spectrum comparisons passed')


def parse_factor(value):
    prime, _, exponent = value.partition(':')
    return int(prime), int(exponent) if exponent else 1


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('factors', nargs='*', type=parse_factor, metavar='P[:E]')
    parser.add_argument('--scale', type=int, default=1)
    parser.add_argument('--cutoff', type=int, default=8)
    parser.add_argument('--precision', type=int, default=25)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--width', type=int, help='Lower window width, default N')
    parser.add_argument('--range', dest='ranges', nargs=2, type=int, action='append',
                        metavar=('LOW', 'HIGH'), help='Half-open denominator interval, repeatable')
    parser.add_argument('--quotient', type=int, default=2, help='Positive bin used without --range')
    parser.add_argument('--digit', type=int, action='append', help='Upper query digit, repeatable')
    args = parser.parse_args()
    if not args.factors:
        check()
        return
    started = perf_counter()
    try:
        evaluator = CrtRange(args.factors, args.scale, args.cutoff, args.precision)
        if args.ranges is None and args.quotient < 2:
            raise ValueError('Expected a nearest quotient of at least 2')
        ranges = args.ranges or [quotient_interval(evaluator.modulus, args.quotient)]
        width = evaluator.n if args.width is None else args.width
        print(f'N={evaluator.n}, factorization supplied, setup {perf_counter() - started:.3f}s '
              f'{dict(evaluator.setup_work)}')
        for digit in args.digit or [0]:
            started = perf_counter()
            result, work = evaluator.value(digit, args.start, width, ranges)
            print(f'j={digit}, value={mp.nstr(result, args.precision)}, '
                  f'query {perf_counter() - started:.3f}s {dict(work)}')
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
