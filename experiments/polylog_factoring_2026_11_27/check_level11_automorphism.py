"""Check proper gcds from actual target units against the local fixed-line counts."""

from collections import Counter
from itertools import combinations, product
from math import gcd
from random import Random

from sympy import Matrix, kronecker_symbol

if __package__:
    from . import level11_automorphism as probe
    from . import level11_sampler as sampler
    from .check_level11_marginals import eta_coefficients
else:
    import level11_automorphism as probe
    import level11_sampler as sampler
    from check_level11_marginals import eta_coefficients


def crt_lines(p, q):
    """Use the fixture primes only to enumerate all lines for the finite population."""
    p_lines = [(1, i) for i in range(p)] + [(0, 1)]
    q_lines = [(1, j) for j in range(q)] + [(0, 1)]
    inverse = pow(p, -1, q)
    for p_line, q_line in product(p_lines, q_lines):
        yield tuple(a + p*((b-a)*inverse % q) for a, b in zip(p_line, q_line))


def expected_at(r, s, coefficients):
    principal = (3*(s+1) + 2*coefficients[s])//5
    return [edges*(1+int(kronecker_symbol(discriminant, r))) for edges, discriminant in
            ((principal, -4), (s+1-principal, -3))]


def check_pair(p, q, coefficients):
    n = p*q
    inverse = sampler._inverse_splitting(n, Random(n))
    counts = Counter()
    checked_types = set()
    for line in crt_lines(p, q):
        units = probe._unit_numerators(n, inverse, line)
        kind = len(units)//2 - 2  # the four-unit class is 0, the six-unit class is 1
        divisor = probe._integrality_divisor(n, units)
        assert {gcd(n, *u) for u in units if any(u[1:])} == {divisor}
        assert len(set(units)) == len(units) and all(sampler.norm(u) == n*n for u in units)
        counts[kind, divisor] += 1
        if kind not in checked_types:
            checked_types.add(kind)
            basis = sampler._line_basis(n, inverse, line)
            to_coordinates = basis.inv()/n
            for u in units:
                for j in range(4):
                    image = Matrix(probe._product(u, tuple(basis.col(j))))
                    assert all(x.is_Integer for x in to_coordinates*image)
            generators = sampler._line_points(n, inverse, line)
            assert bool(generators) == (kind == 0)
            if generators:
                a = generators[0]
                aj = probe._product(a, (0, 0, 1, 0))
                assert probe._product(aj, probe._conjugate(a)) in units
    total = sum(counts.values())
    assert total == (p+1)*(q+1)
    for j in (0, 1):
        assert counts[j, p] + counts[j, n] == expected_at(p, q, coefficients)[j]
        assert counts[j, q] + counts[j, n] == expected_at(q, p, coefficients)[j]
    both = (1 + kronecker_symbol(-4, p))*(1 + kronecker_symbol(-4, q))
    assert counts[0, n] == both and counts[1, n] == 0
    splits = sum(counts[j, d] for j in (0, 1) for d in (p, q))
    return splits, total


def run():
    primes = (2, 3, 5, 7, 13, 17)
    cohorts = [*combinations(primes, 2), (23, 47)]
    coefficients = eta_coefficients(47)
    pairs = cells = splits = 0
    for p, q in cohorts:
        hits, total = check_pair(p, q, coefficients)
        pairs += 1
        cells += total
        splits += hits
        if (p, q) in ((2, 3), (3, 5), (5, 13), (7, 17), (23, 47)):
            print(f"n={p*q}: {hits}/{total} lines give a proper divisor")
    print(f"{pairs} semiprimes, {cells} line ideals, {splits} exact proper gcds")

    rng = Random(11032027)
    tally = Counter()
    for n in (15, 35, 65, 91, 105, 125, 100160063, 1000036000099):
        for _ in range(12):
            divisor = probe.automorphism_probe(n, rng)
            if divisor is not None:
                assert 1 < divisor < n and n % divisor == 0
                tally[n] += 1
    assert sum(tally.values()) and probe.automorphism_probe(22, rng) == 11
    assert probe.automorphism_probe(11, rng) is None
    assert probe.automorphism_probe(91, Random(2)) == 7
    print(f"96 public trials: {dict(sorted(tally.items()))} verified splits")


if __name__ == "__main__":
    run()
