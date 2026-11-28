"""Check conditional rank sampling and the exact uniform mass of the free summands."""

from fractions import Fraction
from math import prod
from random import Random

from sympy import eye, zeros

if __package__:
    from . import free_projector_sampler as free
    from .idempotent_sampler import gl_orders
    from .projector_decompositions import factor_from_decomposition
else:
    import free_projector_sampler as free
    from idempotent_sampler import gl_orders
    from projector_decompositions import factor_from_decomposition


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for rank in range(total + 1):
            for tail in compositions(total - rank, parts - 1):
                yield (rank, *tail)


def orbit_sizes(orders, parts):
    return {ranks: orders[-1] // prod(orders[r] for r in ranks)
            for ranks in compositions(len(orders) - 1, parts)}


class FixedTickets:
    def __init__(self, draws):
        self.draws = iter(draws)

    def randrange(self, total):
        expected, ticket = next(self.draws)
        assert expected == total
        return ticket


def check_profile_corrections(n, dimension, parts):
    orders = free._global_gl_orders(n, dimension)
    weights, dp = free.profile_weights_from_orders(orders, parts)
    sizes = orbit_sizes(orders, parts)
    low = dimension // parts
    exceptional = {ranks: size for ranks, size in sizes.items() if max(ranks) - min(ranks) > 1}
    unbalanced = Fraction(sum(exceptional.values()), sum(sizes.values()))
    gate = min(Fraction(1), Fraction(parts * (parts-1), n*n))
    assert unbalanced <= gate
    assert unbalanced == 1 - Fraction(free._balanced_mass(weights, low, parts, dimension), dp[-1][-1])

    # For each exceptional profile derive the complete ordered tree of ticket
    # intervals by enumeration. Check both endpoints of every chosen interval.
    for target in exceptional:
        for last_ticket in (False, True):
            draws = []
            for j, current in enumerate(target):
                by_next = [0] * (dimension+1)
                for ranks in exceptional:
                    if ranks[:j] == target[:j]:
                        by_next[ranks[j]] += prod(weights[r] for r in ranks[j:])
                first = sum(by_next[:current])
                ticket = first + (by_next[current] - 1 if last_ticket else 0)
                draws.append((sum(by_next), ticket))
            assert free._exceptional_profile(weights, dp, low, FixedTickets(draws)) == target
    return len(exceptional)


def free_probability(p, q, dimension, parts):
    local_p = orbit_sizes(gl_orders(p, 1, dimension), parts)
    local_q = orbit_sizes(gl_orders(q, 1, dimension), parts)
    sizes = orbit_sizes(free._global_gl_orders(p*q, dimension), parts)
    assert all(sizes[r] == local_p[r] * local_q[r] for r in sizes)
    return Fraction(sum(sizes.values()), sum(local_p.values()) * sum(local_q.values()))


def run():
    checks = sum(check_profile_corrections(n, d, s)
                 for n in (2, 6, 9, 12, 35, 551)
                 for d, s in ((1, 3), (2, 2), (3, 2), (4, 3), (6, 3), (3, 5)))
    rng = Random(23)
    draws = 0
    for n in (6, 9, 12, 35, 551, (1 << 127) - 1):
        for dimension, parts in ((2, 2), (3, 2), (3, 5), (4, 3)):
            for _ in range(8):
                projectors = free.sample_free_decomposition(n, dimension, parts, rng)
                assert (sum(projectors, zeros(dimension)) - eye(dimension)).applyfunc(lambda x: x % n).is_zero_matrix
                for i, a in enumerate(projectors):
                    for j, b in enumerate(projectors):
                        target = a if i == j else zeros(dimension)
                        assert (a*b - target).applyfunc(lambda x: x % n).is_zero_matrix
                assert factor_from_decomposition(n, projectors) is None
                draws += 1
    print(f"{checks} exceptional profiles checked at both ends of the ticket ranges; {draws} public free draws")
    for dimension, parts in ((2, 2), (3, 2), (4, 3), (6, 3)):
        probability = free_probability(19, 29, dimension, parts)
        print(f"N=551, d={dimension}, s={parts}: free mass {probability} ({float(probability):.8f})")


if __name__ == "__main__":
    run()
