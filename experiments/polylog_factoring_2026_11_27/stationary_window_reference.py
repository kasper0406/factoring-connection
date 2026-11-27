"""Preprocess stationary cubic and quartic phases from certified quadratic windows."""

import argparse
from collections import Counter
from dataclasses import dataclass
from math import cos, fsum, gcd, isqrt, prod, sin, tau
from time import perf_counter

from lattice_window import block_square_roots, multiplicity
from stationary_cubic import direct_cubic, direct_quartic, lift_roots


@dataclass(frozen=True)
class StationaryReference:
    n: int
    cubic_exponents: tuple[int, ...]
    quartic_exponents: tuple[tuple[int, int], ...]
    blocks: int
    lattices: int
    direct_points: int
    incidence_tests: int

    @classmethod
    def prepare(cls, n):
        if n < 3 or n % 2 == 0:
            raise ValueError('use an odd modulus greater than one')
        window_limit = max(1, isqrt(n) // 8)
        lattice_multiplicity = multiplicity(n, window_limit) if window_limit > 1 else 0
        roots = []
        blocks = lattices = direct_points = 0
        for start in range(0, n, window_limit + 1):
            blocks += 1
            limit = min(window_limit, n - 1 - start)
            if limit <= 1:
                direct_points += limit + 1
                roots.extend(value for value in range(start, start + limit + 1)
                             if (value * value - 1) % n == 0)
            else:
                lattices += 1
                roots.extend(block_square_roots(n, start, limit, lattice_multiplicity, 1))

        modulus = n * n
        half = pow(2, -1, modulus)
        lifts = []
        for root in roots:
            upper = -((root * root - 1) // n) * pow(2 * root, -1, n) % n
            lift = root + n * upper
            assert lift * lift % modulus == 1
            lifts.append(lift)

        idempotents = sorted((lift + 1) * half % modulus for lift in lifts)
        low_idempotents = [value % n for value in idempotents]
        quartic_exponents = []
        for value in idempotents:
            assert (value * value - value) % modulus == 0
            # An idempotent contained in this support independently selects the
            # sign of each nonzero local critical point.
            weight = sum((value * lower - lower) % n == 0 for lower in low_idempotents)
            quartic_exponents.append((-value % modulus, weight))

        return cls(n, tuple(sorted(2 * lift % modulus for lift in lifts)),
                   tuple(quartic_exponents), blocks, lattices, direct_points,
                   len(idempotents) ** 2)

    def evaluate(self, scale, phase):
        if gcd(scale, self.n) != 1:
            raise ValueError('the scale must be coprime to the modulus')
        modulus = self.n * self.n
        if phase == 'cubic':
            if self.n % 3 == 0:
                raise ValueError('the cubic stationary formula requires a modulus coprime to 3')
            return fsum(cos(tau * ((scale * residue) % modulus) / modulus)
                        for residue in self.cubic_exponents)
        if phase == 'quartic':
            angles = [(tau * ((scale * residue) % modulus) / modulus, weight)
                      for residue, weight in self.quartic_exponents]
            return complex(fsum(weight * cos(angle) for angle, weight in angles),
                           fsum(weight * sin(angle) for angle, weight in angles))
        raise ValueError('use cubic or quartic')


def report(reference, elapsed):
    phases = len(reference.quartic_exponents)
    critical_roots = sum(weight for _, weight in reference.quartic_exponents)
    direct_phase_evaluations = (reference.n * reference.n + 1) // 2
    print(f'N={reference.n} blocks={reference.blocks} lattices={reference.lattices} '
          f'direct_points={reference.direct_points} incidence_tests={reference.incidence_tests} '
          f'weighted_phases={phases} quartic_critical_roots={critical_roots} '
          f'direct_phase_evaluations={direct_phase_evaluations} preprocess_seconds={elapsed:.4f}')


def check(components):
    n = prod(components)
    started = perf_counter()
    reference = StationaryReference.prepare(n)
    elapsed = perf_counter() - started
    modulus = n * n
    direct_square_roots = [value for value in range(n) if (value * value - 1) % n == 0]
    lifted_lower_roots = sorted((residue * pow(2, -1, n)) % n for residue in reference.cubic_exponents)
    assert lifted_lower_roots == direct_square_roots
    assert len(direct_square_roots) == 2 ** len(components)

    critical_points = lift_roots(n, lambda value: value ** 3 - value,
                                 lambda value: 3 * value * value - 1)
    expected_weights = Counter((-value * value) % modulus for value in critical_points)
    assert dict(reference.quartic_exponents) == expected_weights
    assert sum(expected_weights.values()) == 3 ** len(components)

    for scale in (1, 1 + n, modulus - 1):
        quartic = reference.evaluate(scale, 'quartic')
        assert abs(quartic - direct_quartic(n, scale)) < 5e-9
        if n % 3:
            cubic = reference.evaluate(scale, 'cubic')
            assert abs(cubic - direct_cubic(n, scale)) < 5e-9

    report(reference, elapsed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modulus', type=int, help='prepare a reference from this unfactored input')
    parser.add_argument('--scale', type=int, action='append', help='repeat to reuse the cached phases')
    parser.add_argument('--phase', choices=('cubic', 'quartic'), default='quartic')
    args = parser.parse_args()
    if args.modulus is None:
        for components in ((5,), (5, 7), (9, 5), (25, 7), (9, 5, 7), (5, 7, 11), (49, 11), (27, 25)):
            check(components)
        return

    started = perf_counter()
    try:
        reference = StationaryReference.prepare(args.modulus)
        report(reference, perf_counter() - started)
        for scale in args.scale or [1]:
            value = reference.evaluate(scale, args.phase)
            print(f'scale={scale} normalized_sum={value} magnitude={abs(value)}')
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
