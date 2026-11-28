"""Compare the four-square incidence samplers with exact finite counts."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import product
from math import gcd, isqrt
from random import Random
from unittest.mock import patch

if __package__:
    from . import four_square_incidence as incidence
    from . import four_square_sampler as sphere
else:
    import four_square_incidence as incidence
    import four_square_sampler as sphere


def quadruples(norm):
    pairs = defaultdict(list)
    extent = isqrt(norm)
    for a in range(-extent, extent + 1):
        bound = isqrt(norm - a * a)
        for b in range(-bound, bound + 1):
            pairs[a * a + b * b].append((a, b))
    for total, first in list(pairs.items()):
        for left, right in product(first, pairs.get(norm - total, ())):
            yield left + right


def prime_powers(n):
    powers = []
    for prime in range(3, isqrt(n) + 1, 2):
        exponent = 0
        while n % prime == 0:
            n //= prime
            exponent += 1
        if exponent:
            powers.append((prime, exponent))
    if n > 1:
        powers.append((n, 1))
    return powers


def normalized(masses):
    total = sum(masses.values())
    return {key: Fraction(value, total) for key, value in masses.items()}


@cache
def gcd_histogram(n, h):
    return Counter(gcd(n, *values) for values in product(range(n), repeat=h))


@cache
def projected_acceptance(n):
    primitive = [(x, y) for x, y in product(range(n), repeat=2) if gcd(x, y, n) == 1]
    return Fraction(sum(gcd(x, n) == 1 for x, _ in primitive), len(primitive))


class Choices:
    """Replay selected randrange results and verify their requested bounds."""

    def __init__(self, choices):
        self.choices = iter(choices)

    def randrange(self, stop):
        expected, value = next(self.choices)
        assert stop == expected and 0 <= value < stop
        return value

    def finished(self):
        assert next(self.choices, None) is None


def check_factored_content(n):
    factors = prime_powers(n)
    bounds = [exponent * (prime + 1) + prime for prime, exponent in factors]
    frequencies = Counter()
    with patch.object(incidence, "_point_and_scalars", lambda n, g, count, rng: g):
        for tickets in product(*(range(bound) for bound in bounds)):
            rng = Choices(zip(bounds, tickets))
            frequencies[incidence.sample_two_from_factors(factors, rng)] += 1
            rng.finished()
    return normalized(frequencies)


def check_public_content(n, count):
    h = abs(2 - count)
    masses = {}
    for v, frequency in gcd_histogram(n, h).items():
        g = v if count < 2 else n // v
        m = n // g
        # Enumerate all residues for each rejection, including repeated prime powers.
        first_modulus = incidence._other_prime_powers(n, n // v)
        second_modulus = incidence._other_prime_powers(n, m)
        correction = Fraction(gcd_histogram(first_modulus, h)[1], first_modulus**h)
        masses[g] = Fraction(frequency, n**h) * correction * projected_acceptance(second_modulus)
        # Replay an accepted path through the public function for every possible g.
        initial = [(n, v % n)]
        if v != 1:
            initial.extend([(n, 0)] * (h - 1))
        corrections = [(first_modulus, 1 % first_modulus)]
        pair = [(second_modulus, 1 % second_modulus), (second_modulus, 0)]
        rng = Choices(initial + corrections + pair)
        assert incidence.sample_content(n, count, rng) == g
        rng.finished()
    return normalized(masses)


def check_primitive_fibres():
    input_count = 0
    for norm in (3, 5, 9, 15):
        points = {point for point in quadruples(norm) if gcd(*point) == 1}
        residues = list(product(range(norm), repeat=2))
        primitive = [(x, y) for x, y in residues if gcd(x, y, norm) == 1]
        multiplicity = sum(gcd(unit, norm) == 1 for unit in range(norm))
        for a, b in residues:
            if (a * a + b * b + 1) % norm:
                continue
            fibres = Counter()
            for x, y in primitive:
                for index in range(8):
                    rng = Choices([(norm, x), (norm, y), (8, index)])
                    fibres[sphere.primitive_four_squares(norm, rng, split=(a, b))] += 1
                    rng.finished()
            assert set(fibres) == points and set(fibres.values()) == {multiplicity}
            input_count += sum(fibres.values())
    return input_count


def check_conditioned(n, fixed, histograms):
    h = gcd(n, fixed)
    target = {g: g * size for g, size in histograms[n].items() if (fixed * g) % n == 0}
    lifted = {n // h * g: (n // h) * g * size for g, size in histograms[h].items()}
    assert lifted == target
    for g in histograms[h]:
        residues = [t for t in range(h) if (t * g) % h == 0]
        lifts = {t + h * i for t in residues for i in range(n // h)}
        assert lifts == {t for t in range(n) if t * (n // h * g) % n == 0}


def check_draw(n, point, scalars, count):
    assert len(point) == 4 and sum(x * x for x in point) == n * n
    assert len(scalars) == count
    for scalar in scalars:
        assert 0 <= scalar < n and all((scalar * x) % n == 0 for x in point)


def main():
    rng = Random(20261128)
    histograms = {}
    laws = strata = conditional_laws = draws = 0
    for n in range(1, 46, 2):
        histogram = Counter(gcd(*point) for point in quadruples(n * n))
        histograms[n] = histogram
        for count in range(5):
            expected = normalized({g: g**count * size for g, size in histogram.items()})
            actual = check_factored_content(n) if count == 2 else check_public_content(n, count)
            assert actual == expected, (n, count, actual, expected)
            laws += 1
            strata += len(actual)
        for fixed in range(n):
            check_conditioned(n, fixed, histograms)
            conditional_laws += 1

    for n in [*histograms, 169, 945, 1022117]:
        for count in (0, 1, 2, 3, 4, 7):
            for _ in range(4):
                if count == 2:
                    point, scalars = incidence.sample_two_from_factors(prime_powers(n), rng)
                else:
                    point, scalars = incidence.sample_incidence(n, count, rng)
                check_draw(n, point, scalars, count)
                draws += 1
        for fixed in (0, 1, -3, n // 3, 2 * n):
            point, t = incidence.sample_with_fixed_scalar(n, fixed, rng)
            check_draw(n, point, (fixed % n, t), 2)
            draws += 1

    fibres = check_primitive_fibres()
    point_count = sum(sum(histogram.values()) for histogram in histograms.values())
    print(f"{point_count} sphere points, "
          f"{laws} exact laws over {strata} strata, {conditional_laws} fixed-scalar laws, "
          f"{fibres} primitive generator inputs, {draws} seeded draws")


if __name__ == "__main__":
    main()
