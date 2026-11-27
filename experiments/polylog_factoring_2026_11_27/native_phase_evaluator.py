"""Evaluate complete fixed cubic and quartic phases without a factorization step."""

from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt
from time import perf_counter

from mpmath import mp
from sympy import integer_nthroot

from stationary_dual import proven_prime


def polynomial(value, family, modulus):
    square = value * value % modulus
    if family == 'cubic':
        return value * (3 - square) % modulus
    return square * (square - 2) % modulus


def public_prime_power(modulus):
    """Recognize a prime power by exact integer roots and a primality proof."""
    for exponent in range(modulus.bit_length() - 1, 1, -1):
        base, exact = integer_nthroot(modulus, exponent)
        if exact:
            return (int(base), exponent) if proven_prime(base) else None
    return None


def character(value, prime):
    return 1 if pow(value % prime, (prime - 1) // 2, prime) == 1 else -1


@dataclass(frozen=True)
class PhasePlan:
    modulus: int
    family: str
    method: str
    setup_checks: int
    weights: tuple[tuple[int, int], ...] | None
    gauss_prime: int | None = None

    def evaluate(self, scale, bits=45):
        """Return the normalized complex sum and the number of phase evaluations."""
        if gcd(scale, self.modulus) != 1:
            raise ValueError('Expected a unit scale')
        if bits < 1:
            raise ValueError('Expected positive precision')
        modulus = self.modulus
        with mp.workprec(modulus.bit_length() + bits + 40):
            if self.weights is None:
                terms = (polynomial(value, self.family, modulus) for value in range(1, (modulus + 1) // 2))
                if self.family == 'cubic':
                    value = 1 + 2 * mp.fsum(mp.cospi(mp.mpf(2 * (scale * exponent % modulus)) / modulus)
                                             for exponent in terms)
                else:
                    value = 1 + 2 * mp.fsum(mp.expjpi(mp.mpf(2 * (scale * exponent % modulus)) / modulus)
                                             for exponent in terms)
                return value / mp.sqrt(modulus), (modulus + 1) // 2

            value = mp.fsum(weight * mp.expjpi(mp.mpf(2 * (scale * exponent % modulus)) / modulus)
                            for exponent, weight in self.weights)
            if self.gauss_prime is not None:
                prime = self.gauss_prime
                value *= character(scale, prime) * (1 if prime % 4 == 1 else 1j)
            return value, len(self.weights)


def prepare(modulus, family):
    """Cache a square-modulus scan or use the public prime-power formula."""
    if modulus < 3 or modulus % 2 == 0 or family not in ('cubic', 'quartic'):
        raise ValueError('Expected an odd modulus and the cubic or quartic phase')
    if family == 'cubic' and modulus % 3 == 0:
        return PhasePlan(modulus, family, 'literal', 0, None)

    power = public_prime_power(modulus)
    if power is not None:
        prime, exponent = power
        roots = (1, -1) if family == 'cubic' else (0, 1, -1)
        weights = Counter()
        for root in roots:
            second_half = -3 * root if family == 'cubic' else 6 * root * root - 2
            weight = character(second_half, prime) if exponent % 2 else 1
            weights[polynomial(root, family, modulus)] += weight
        return PhasePlan(modulus, family, 'public-prime-power', 0, tuple(weights.items()),
                         prime if exponent % 2 else None)

    base = isqrt(modulus)
    if base * base != modulus:
        return PhasePlan(modulus, family, 'literal', 0, None)

    weights = Counter()
    for lower in range(base):
        square = lower * lower % base
        critical = (square == 1 if family == 'cubic' else lower * (square - 1) % base == 0)
        if critical:
            weights[polynomial(lower, family, modulus)] += 1
    return PhasePlan(modulus, family, 'derivative-scan', base, tuple(weights.items()))


def check():
    checks = 0
    methods = Counter()
    for modulus in (15, 35, 45, 125, 343, 625, 1_225, 2_025, 3_125, 30_625):
        for family in ('cubic', 'quartic'):
            prepared = prepare(modulus, family)
            literal = PhasePlan(modulus, family, 'literal', 0, None)
            for scale in (1, 2, modulus - 1):
                actual, _ = prepared.evaluate(scale)
                control, _ = literal.evaluate(scale)
                assert abs(actual - control) < mp.mpf('1e-12')
                methods[prepared.method] += 1
                checks += 1
    print(f'{checks} comparisons against literal sums, {dict(methods)}')

    started = perf_counter()
    prepared = prepare(11_009**2, 'quartic')
    elapsed = perf_counter() - started
    value, phases = prepared.evaluate(17)
    assert abs(abs(value) - mp.mpf('2.709761870264364')) < mp.mpf('1e-12')
    print(f'N=11009 square case: {prepared.setup_checks} derivative checks, '
          f'{phases} cached phases, setup seconds {elapsed:.6g}')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('modulus', nargs='?', type=int, help='Summation modulus, use N² for the square target')
    parser.add_argument('--phase', choices=('cubic', 'quartic'), default='quartic')
    parser.add_argument('--scale', type=int, action='append')
    args = parser.parse_args()
    if args.modulus is None:
        check()
    else:
        plan = prepare(args.modulus, args.phase)
        print(f'{plan.method}: {plan.setup_checks} lower-digit checks')
        for scale in args.scale or [1]:
            value, phases = plan.evaluate(scale)
            print(f'scale={scale}, phases={phases}, normalized magnitude={mp.nstr(abs(value), 16)}')
