"""Certify the sampler distribution and contrast against literal modular inversion."""

from collections import defaultdict
from fractions import Fraction
from math import gcd
from random import Random

from clipped_gap_sampler import ClippedGapSampler, GapQuery


class NextChoice(Exception):
    pass


class PrefixRng:
    """Stop at the first unsupplied randrange, allowing exhaustive branch weights."""
    def __init__(self, prefix):
        self.choices = iter(prefix)

    def randrange(self, stop):
        try:
            answer = next(self.choices)
        except StopIteration:
            raise NextChoice(stop) from None
        assert 0 <= answer < stop
        return answer


def exact_sampler_law(sampler):
    result = defaultdict(Fraction)
    pending = [((), Fraction(1))]
    while pending:
        prefix, mass = pending.pop()
        try:
            a = sampler.draw(PrefixRng(prefix))
        except NextChoice as split:
            stop, = split.args
            pending.extend((prefix+(i,), mass/stop) for i in range(stop))
        else:
            result[a] += mass
    return result


def proposal_law(sampler):
    n, c = sampler.n, sampler.c0
    return {a: Fraction(gcd(a*a-c*c, n), n)/sampler.normalizer for a in range(n)}


def literal_contrast(sampler, query, omit_mixed):
    n, c = sampler.n, sampler.c0
    square = n*n
    result = 0
    for a in range(n):
        if (a+query.offset) % n >= query.width or a in (c, n-c) or gcd(a, n) != 1:
            continue
        if omit_mixed and (a*a-c*c) % n == 0:
            continue
        sign = 1 if (query.sign_multiplier*a+query.sign_offset) % n <= (n-1)//2 else -1
        for x in range(a, square, n):
            inverse = pow(x, -1, square)
            at_j = ((x + (c+n*query.j)**2*inverse) % square)//n
            at_zero = ((x + c*c*inverse) % square)//n
            result += sign*((at_j <= (n-1)//2) - (at_zero <= (n-1)//2))
    return Fraction(result, n)


def check_contrast(sampler, law, query, nonstationary):
    estimated_mean = Fraction()
    for a, weight in law.items():
        value = sampler.sample_value(a, query, nonstationary_only=nonstationary)
        assert value in (-sampler.normalizer, 0, sampler.normalizer)
        estimated_mean += weight*value
    assert estimated_mean == literal_contrast(sampler, query, nonstationary)
    return estimated_mean


def run():
    all_queries = 0
    for factors, c in [({3:3}, 2), ({5:2}, 19), ({3:1, 5:1}, 1),
                       ({3:2, 5:1}, 2), ({3:2, 5:1, 7:1}, 4)]:
        sampler = ClippedGapSampler(factors, c)
        n = sampler.n
        law = proposal_law(sampler)
        assert sum(law.values()) == 1
        if n <= 45:
            assert exact_sampler_law(sampler) == law
        for q in [GapQuery(0, 0, n, 2, 4), GapQuery(3, -9, 0, 7, 2),
                  GapQuery(2, n-4, n, 5, -19), GapQuery(n//2, 14, n//2, 3, 7)]:
            for residual in (False, True):
                check_contrast(sampler, law, q, residual)
                all_queries += 1
        if n == 45:
            for a, wanted, residual in [(7, 1, 0), (1, Fraction(-1, 15), Fraction(-1, 15))]:
                query = GapQuery(2, -a, 1, 0, 0)
                assert check_contrast(sampler, law, query, False) == wanted
                assert check_contrast(sampler, law, query, True) == residual
                all_queries += 2
    # Verify aggregation at a larger modulus against the individual seeded draws.
    sampler = ClippedGapSampler({3:90, 5:70}, 19)
    query = GapQuery(23, -7, sampler.n, 37, 91)
    manual_rng = Random(20261128)
    total = sum(sampler.sample_value(sampler.draw(manual_rng), query) for _ in range(128))
    assert sampler.estimate(query, 128, rng=Random(20261128)) == total/128
    print(f'Four exhaustive sampler laws, {all_queries} exact contrasts, '
          f'128 seeded samples on a {sampler.n.bit_length()}-bit modulus.')


if __name__ == '__main__':
    run()
