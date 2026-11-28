"""Exact local counts, CRT checks and public sampling of ordered decompositions."""

from collections import Counter
from itertools import product
from math import gcd, prod
from random import Random

from sympy import eye, zeros

if __package__:
    from .idempotent_sampler import gl_orders
    from .projector_decompositions import (
        _factor_hit, conjugate_profiles, draw_rank_profile, factor_from_decomposition,
        rank_profile_weights, sample_decomposition,
    )
else:
    from idempotent_sampler import gl_orders
    from projector_decompositions import (
        _factor_hit, conjugate_profiles, draw_rank_profile, factor_from_decomposition,
        rank_profile_weights, sample_decomposition,
    )


def rank_mod(point, p, d):
    images = {tuple(sum(point[i*d+j]*v[j] for j in range(d)) % p for i in range(d))
              for v in product(range(p), repeat=d)}
    return [p**j for j in range(d+1)].index(len(images))


def local_counts(p, e, d, s):
    modulus = p**e
    identity = tuple(int(i == j) for i in range(d) for j in range(d))
    zero = (0,) * (d*d)

    def multiply(a, b):
        return tuple(sum(a[i*d+k]*b[k*d+j] for k in range(d)) % modulus
                     for i in range(d) for j in range(d))

    ranks = {a: rank_mod(a, p, d) for a in product(range(modulus), repeat=d*d)
             if multiply(a, a) == a}
    counts = Counter()
    if s == 2:
        for a, rank in ranks.items():
            complement = tuple((u-v) % modulus for u, v in zip(identity, a))
            counts[rank, ranks[complement]] += 1
    else:
        assert s == 3
        for a, left_rank in ranks.items():
            for b, right_rank in ranks.items():
                if multiply(a, b) == zero and multiply(b, a) == zero:
                    c = tuple((u-v-w) % modulus for u, v, w in zip(identity, a, b))
                    counts[left_rank, right_rank, ranks[c]] += 1

    orders = gl_orders(p, e, d)
    weights, dp = rank_profile_weights(p, e, d, s)
    for ranks, count in counts.items():
        assert count * prod(orders[r] for r in ranks) == orders[d]
        assert count * orders[d]**(s-1) == prod(weights[r] for r in ranks)
    assert sum(counts.values()) * orders[d]**(s-1) == dp[s][d]
    return counts


def check_projectors(projectors, n):
    identity = eye(projectors[0].rows)
    assert (sum(projectors, zeros(identity.rows)) - identity).applyfunc(lambda x: x % n).is_zero_matrix
    for i, a in enumerate(projectors):
        for j, b in enumerate(projectors):
            expected = a if i == j else zeros(identity.rows)
            assert (a*b - expected).applyfunc(lambda x: x % n).is_zero_matrix


def run():
    checks = ((2,1,2,2), (3,1,2,2), (2,2,2,2), (2,1,3,3), (3,1,2,3))
    tables = {case: local_counts(*case) for case in checks}
    local_total = sum(sum(c.values()) for c in tables.values())

    # Exhaust all residue pairs at n=12. Extraction groups marks correctly
    # for repeated primes, and the two local Bernoulli marks are independent.
    marks = Counter(tuple(_factor_hit(gcd(12, x, y))) for x, y in product(range(12), repeat=2))
    assert marks == {(): 96, (2,): 32, (3,): 12, (2,3): 4}

    rng = Random(3451)
    profiles = 0
    for p in (2, 3, 13, 37):
        for e in (1, 2):
            for d, s in ((2, 2), (3, 3), (6, 3), (8, 4), (12, 6)):
                weights, dp = rank_profile_weights(p, e, d, s)
                common = (d//s,) * s
                total = dp[s][d]
                central = weights[d//s]**s
                passes = s*s
                # Exact domination needed for the sparse local correction.
                assert central * (p*p)**passes >= total * (p*p-1)**passes
                profile = draw_rank_profile(weights, dp, rng, exclude=common)
                assert sum(profile) == d and profile != common
                profiles += 1

    # Component profiles are imposed in a single common basis. Local zeros
    # and ranks differing by an entire characteristic must remain visible.
    mixed = conjugate_profiles(6, (0, 6), {3: (6, 0)}, rng)
    check_projectors(mixed, 6)
    assert all(int(a.trace()) % 6 == 0 for a in mixed)
    assert factor_from_decomposition(6, mixed) == 2
    equal = conjugate_profiles(35, (1, 2, 1), {}, rng)
    assert factor_from_decomposition(35, equal) is None
    override = conjugate_profiles(60, (1, 1, 1), {4: (0, 2, 1), 5: (3, 0, 0)}, rng)
    check_projectors(override, 60)
    assert [tuple(rank_mod(tuple(map(int, a)), p, 3) for a in override) for p in (2,3,5)] == [
        (0, 2, 1), (1, 1, 1), (3, 0, 0)]
    divisor = factor_from_decomposition(60, override)
    assert 1 < divisor < 60 and 60 % divisor == 0

    # Force one sparse hit to run the thinning and conditional sampler on
    # a large prime. The remaining prime keeps its central rank profile.
    class ForcedHit(Random):
        def __init__(self):
            super().__init__(3107)
            self.tickets = iter([23, 23] + [1, 1]*3 + [0, 0, 0])

        def randrange(self, stop):
            ticket = next(self.tickets, None)
            if ticket is None:
                return super().randrange(stop)
            assert 0 <= ticket < stop
            return ticket

    forced = sample_decomposition(23*53, 2, 2, ForcedHit())
    check_projectors(forced, 23*53)
    assert factor_from_decomposition(23*53, forced) == 23

    # When no constant profile exists, permutation symmetry bounds the
    # probability of a profile match by 1/s, including over different primes.
    odd_shape = local_counts(2, 1, 2, 3)
    other_shape = tables[3, 1, 2, 3]
    matching = sum(count * other_shape[ranks] for ranks, count in odd_shape.items())
    assert 3 * matching <= sum(odd_shape.values()) * sum(other_shape.values())

    samples = 0
    for n, d, s in ((1, 4, 2), (2, 2, 2), (12, 6, 3), (72, 4, 4),
                    (19*37, 6, 3), (1000003*1000033, 4, 2)):
        for _ in range(4):
            result = sample_decomposition(n, d, s, rng)
            check_projectors(result, n)
            factor = factor_from_decomposition(n, result)
            assert factor is None or (1 < factor < n and n % factor == 0)
            samples += 1
    print(f"{local_total} exhaustively counted local decompositions in {len(tables)} cases")
    print(f"{profiles} exact sparse-envelope and conditional-profile checks")
    print(f"{samples} public samples, repeated-prime CRT and characteristic-polynomial decoding passed")


if __name__ == "__main__":
    run()
