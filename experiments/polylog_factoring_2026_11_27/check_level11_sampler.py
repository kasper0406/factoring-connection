"""Exhaustive small-shell checks and seeded trials for the level-11 sampler."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from random import Random

if __package__:
    from . import level11_sampler as sampler
else:
    import level11_sampler as sampler


def direct_shell(n):
    pairs = defaultdict(list)
    limit = isqrt(2*n)
    for x, y in product(range(-limit, limit+1), repeat=2):
        value = x*x + x*y + 3*y*y
        if value <= n:
            pairs[value].append((x, y))
    return {a+b for k, first in pairs.items() for a in first for b in pairs.get(n-k, ())}


def projective_lines(n):
    pairs = {(x, y) for x, y in product(range(n), repeat=2) if gcd(x, y, n) == 1}
    units = [x for x in range(n) if gcd(x, n) == 1]
    while pairs:
        x, y = min(pairs)
        orbit = {(a*x % n, a*y % n) for a in units}
        assert len(orbit) == len(units) and orbit <= pairs
        pairs.difference_update(orbit)
        yield x, y


def check_ideals(n, rng):
    shell = direct_shell(n)
    inverse = sampler._inverse_splitting(n, rng)
    lines = list(projective_lines(n))
    points = []
    for line in lines:
        generators = sampler._line_points(n, inverse, line)
        assert len(generators) in (0, 4)
        points.extend(generators)
    assert len(set(points)) == len(points)
    assert set(points) == {p for p in shell if gcd(*p) == 1}
    assert 3*len(points) >= 4*len(lines)
    return shell, len(points)//4, len(lines)


class Replay:
    """Read one predetermined draw from each requested range."""

    def __init__(self, values):
        self.values = iter(values)

    def randrange(self, *args):
        value = next(self.values)
        assert value in range(*args)
        return value


def proposal_masses(n):
    highest = isqrt(n).bit_length() - 1
    ticket_stop = 1 << (highest+1)
    masses = Counter()
    for ticket in range(1, ticket_stop):
        start = 1 << (highest - (ticket.bit_length()-1))
        for g in range(start, 2*start):
            for toss in range(g*g):
                result = sampler._propose_content(n, Replay((ticket, g, toss)))
                if result is not None:
                    masses[result] += Fraction(1, (ticket_stop-1)*start*g*g)
    return masses


def check_full_law(n, fixtures):
    shell, _, _ = fixtures[n]
    masses = {}
    for g, proposal in proposal_masses(n).items():
        m = n//(g*g)
        _, principal, lines = fixtures[m]
        other = sampler._other_prime_powers(n, m)
        eligible = [(u, v) for u, v in product(range(other), repeat=2) if gcd(u, v, other) == 1]
        correction = Fraction(sum(gcd(u, other) == 1 for u, _ in eligible), len(eligible))
        masses[g] = proposal*correction*Fraction(principal, lines)
    content_counts = Counter(gcd(*point) for point in shell)
    assert set(masses) == set(content_counts)
    assert len({mass/content_counts[g] for g, mass in masses.items()}) == 1
    assert sum(masses.values()) >= Fraction(1, 6)
    return len(content_counts)


def check_ramified(n, rng):
    m, exponent = sampler._remove_eleven(n)
    base = direct_shell(m)
    restored = {sampler._restore_eleven(p, exponent) for p in base}
    shell = direct_shell(n)
    assert len(restored) == len(base) and restored == shell
    primitive = {p for p in shell if gcd(*p) == 1}
    if exponent <= 1:
        assert sampler.sample_primitive(n, rng) in primitive
    else:
        assert not primitive
        try:
            sampler.sample_primitive(n, rng)
        except ValueError as error:
            assert "empty" in str(error)
        else:
            raise AssertionError(f"primitive draws unexpectedly accepted at {n}")
    assert sampler.sample(n, rng) in shell
    return len(shell)


def run():
    rng = Random(11032026)
    norms = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 20, 21, 24, 25,
             27, 32, 35, 40, 45, 48, 49, 60, 63, 65, 75, 81, 84, 91, 105)
    fixtures = {n: check_ideals(n, rng) for n in norms}
    strata = sum(check_full_law(n, fixtures) for n in norms)
    print(f"{len(norms)} norms, {sum(x[2] for x in fixtures.values())} line ideals, "
          f"{sum(len(x[0]) for x in fixtures.values())} shell points, {strata} contents")
    for n in (8, 9, 45, 81, 91):
        shell, principal, lines = fixtures[n]
        print(f"n={n}: {principal}/{lines} principal ideals, "
              f"content counts {dict(sorted(Counter(gcd(*p) for p in shell).items()))}")

    # This includes repeated primes, a three-prime input and two large inputs.
    # Only the integer n is passed to either sampler.
    draws = 0
    for n in (1, 9, 45, 81, 105, 10403, 100160063, 1000036000099, 2**64,
              11 * 2**113 * 1009, (2**127-1)*(2**61-1), (2**127-1)**2):
        for sample in (sampler.sample_primitive, sampler.sample):
            result = sample(n, rng)
            assert sampler.norm(result) == n
            if sample is sampler.sample_primitive:
                assert gcd(*result) == 1
            if n in fixtures:
                assert result in fixtures[n][0]
            draws += 1
    ramified = (11, 22, 33, 44, 55, 99, 121, 176, 242, 363, 1331, 2662)
    restored = sum(check_ramified(n, rng) for n in ramified)
    print(f"{len(ramified)} ramified norms, {restored} shell points checked by bijection")
    observed = Counter(gcd(*sampler.sample(81, rng)) for _ in range(256))
    assert set(observed) == {1, 3, 9}
    for n in (0, -1):
        for sample in (sampler.sample, sampler.sample_primitive):
            try:
                sample(n, rng)
            except ValueError:
                pass
            else:
                raise AssertionError(f"invalid norm {n} accepted")
    print(f"{draws} public draws through 254 bits, plus 256 draws at n=81 "
          f"with contents {dict(sorted(observed.items()))}")


if __name__ == "__main__":
    run()
