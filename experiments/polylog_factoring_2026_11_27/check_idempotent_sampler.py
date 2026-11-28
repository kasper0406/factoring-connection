"""Verify sparse marks, rank weights and conjugacy for the even idempotent sampler."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import gcd, prod
from random import Random
from unittest.mock import patch

from sympy import eye

if __package__:
    from . import idempotent_sampler as sampler
else:
    import idempotent_sampler as sampler


class Unresolved(Exception):
    """A Bernoulli decision not present in the current symbolic tape."""


def check_sparse(limit):
    def replay(numerator, denominator, tape):
        try:
            return next(tape)
        except StopIteration:
            raise Unresolved(Fraction(numerator, denominator)) from None

    outcomes = Counter()
    def visit(tape, mass):
        try:
            marks = tuple(sampler.sparse_integers(limit, iter(tape)))
        except Unresolved as branch:
            probability, = branch.args
            if probability:
                visit(tape+[True], mass*probability)
            if probability < 1:
                visit(tape+[False], mass*(1-probability))
        else:
            outcomes[marks] += mass

    with patch.object(sampler, "bernoulli", replay):
        visit([], Fraction(1))
    assert len(outcomes) == 2**(limit-1) and sum(outcomes.values()) == 1
    for marks, mass in outcomes.items():
        assert mass == prod(Fraction(1 if i in marks else i*i-1, i*i)
                            for i in range(2, limit+1))
    return len(outcomes)


def local_probabilities(p, exponent, dimension):
    lengths, denominator = sampler._rank_intervals(p, exponent, dimension)
    coverage = 1-(1-Fraction(1, p*p))**2
    probabilities = [coverage*Fraction(length, denominator) for length in lengths]
    probabilities[dimension//2] = 1-sum(probabilities)
    return probabilities


def factors(n):
    powers = []
    for p in range(2, n+1):
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            powers.append((p, exponent))
    return powers


def rank2(matrix, p):
    a, b, c, d = matrix
    if (a*d-b*c) % p:
        return 2
    return int(any(t % p for t in matrix))


def check_2x2(n):
    powers = factors(n)
    classes, invertible = defaultdict(set), []
    for matrix in product(range(n), repeat=4):
        a, b, c, d = matrix
        det = a*d-b*c
        if gcd(det, n) == 1:
            invertible.append(matrix)
        errors = a*a+b*c-a, b*(a+d-1), c*(a+d-1), b*c+d*d-d
        if all(t % n == 0 for t in errors):
            ranks = tuple(rank2(matrix, p) for p, _ in powers)
            classes[ranks].add(matrix)

    total = sum(map(len, classes.values()))
    assert len(classes) == 3**len(powers)
    for ranks, matrices in classes.items():
        mass = prod(local_probabilities(p, e, 2)[k] for k, (p, e) in zip(ranks, powers))
        assert mass == Fraction(len(matrices), total)
        if n <= 6:
            # A diagonal representative of this class and every actual conjugate.
            d0, d1 = next((a, d) for a, b, c, d in matrices if b == c == 0)
            orbit = Counter()
            for a, b, c, d in invertible:
                inverse_det = pow(a*d-b*c, -1, n)
                entries = a*d*d0-b*c*d1, a*b*(d1-d0), c*d*(d0-d1), a*d*d1-b*c*d0
                orbit[tuple(t*inverse_det % n for t in entries)] += 1
            assert set(orbit) == matrices and len(set(orbit.values())) == 1
    return total, len(classes)


def check_odd_dimension():
    # At the odd prime 3, exhaust every 3x3 matrix and count parity signs directly.
    signs = Counter()
    for entries in product(range(3), repeat=9):
        if any((sum(entries[3*i+k]*entries[3*k+j] for k in range(3))-entries[3*i+j]) % 3
               for i, j in product(range(3), repeat=2)):
            continue
        a, b, c, d, e, f, g, h, k = (int(i == j)-2*entries[3*i+j]
                                      for i, j in product(range(3), repeat=2))
        signs[(a*(e*k-f*h)-b*(d*k-f*g)+c*(d*h-e*g)) % 3] += 1
    assert set(signs) == {1, 2} and signs[1] == signs[2]
    return dict(signs)


class FixedEntries:
    def __init__(self, entries):
        self.entries = iter(entries)

    def randrange(self, modulus):
        v = next(self.entries)
        assert 0 <= v < modulus
        return v


def run():
    print(f"{check_sparse(8)} sparse masks with exact product probabilities")
    with patch.object(sampler, "sparse_integers", side_effect=[(2, 3, 4, 6, 7), (2, 5, 9)]):
        assert list(sampler.candidate_prime_powers(2**4*3**3*7**2, 2, None)) == [
            (2, 4), (3, 3), (7, 2)
        ]
    totals = [check_2x2(n) for n in (1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15)]
    print(f"{sum(t for t, _ in totals)} idempotents, {sum(s for _, s in totals)} rank profiles "
          f"at {len(totals)} exhaustive moduli")
    basis, inverse = sampler.random_basis(6, 2, FixedEntries([2, 3, 3, 2]))
    assert (basis*inverse).applyfunc(lambda x: x % 6) == eye(2)
    mixed = sampler.sample_idempotent(35, 2, Random(8))
    root = int((eye(2)-2*mixed).det()) % 35
    assert root == 6 and gcd(root-1, 35) == 5
    print(f"3x3 idempotents modulo three by determinant sign: {check_odd_dimension()}")

    rng, profiles = Random(271829), Counter()
    for _ in range(512):
        matrix = sampler.sample_idempotent(6, 2, rng)
        assert all(x % 6 == 0 for x in matrix*matrix-matrix)
        profiles[(rank2(matrix, 2), rank2(matrix, 3))] += 1
    assert len(profiles) == 9
    print(f"512 input-only samples at n=6 reached all nine profiles; "
          f"{sum(v for (r, s), v in profiles.items() if r != s)} had differing ranks")
    for n, d in ((1, 4), (2**12, 4), (3**7*13**2, 8), (77, 16),
                 ((2**127-1)*(2**61-1), 12), (2**256*3**101, 10)):
        matrix = sampler.sample_idempotent(n, d, rng)
        assert matrix.shape == (d, d) and all(x % n == 0 for x in matrix*matrix-matrix)
    for n, d in ((0, 2), (7, 3), (7, 0)):
        try:
            sampler.sample_idempotent(n, d, rng)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid modulus or dimension accepted")
    print("Six further public draws passed, up to dimension 16 and 417-bit modulus.")


if __name__ == "__main__":
    run()
