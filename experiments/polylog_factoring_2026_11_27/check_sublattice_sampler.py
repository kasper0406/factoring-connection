"""Independent enumeration of intermediate lattices, local counts and index readouts."""

from collections import Counter
from itertools import combinations_with_replacement, product
from math import gcd, prod
from random import Random

from sympy import ZZ, diag, eye
from sympy.matrices.normalforms import smith_normal_form

if __package__:
    from .sublattice_sampler import (
        _central_mass, draw_type, factor_from_index, factor_from_lattice,
        lattice_from_types, sample_sublattice, type_counts,
    )
else:
    from sublattice_sampler import (
        _central_mass, draw_type, factor_from_index, factor_from_lattice,
        lattice_from_types, sample_sublattice, type_counts,
    )


def hnfs(n, d):
    """Enumerate column-HNF boxes independently, testing nI for membership."""
    divisors = [z for z in range(1, n+1) if n % z == 0]
    for diagonal in product(divisors, repeat=d):
        positions = [(i,j) for i in range(d) for j in range(i+1,d)]
        for entries in product(*(range(diagonal[i]) for i,j in positions)):
            h = diag(*diagonal)
            for (i,j), value in zip(positions, entries):
                h[i,j] = value
            if all(x.is_Integer for x in n*h.inv()):
                yield h


def profile_of(h, p, e):
    smith = smith_normal_form(h, domain=ZZ)
    lengths = []
    for i in range(h.rows):
        z = abs(int(smith[i,i]))
        valuation = 0
        while z % p == 0:
            z //= p
            valuation += 1
        lengths.append(e - valuation)
    assert all(0 <= x <= e for x in lengths)
    return tuple(sum(x >= j for x in lengths) for j in range(1,e+1))


def weights_of(p, e, d):
    transitions, dp = type_counts(p, e, d)
    weights = {}
    for increasing in combinations_with_replacement(range(d+1), e):
        path = (0,) + increasing
        weights[increasing[::-1]] = prod(transitions[s][r] for s,r in zip(path, path[1:]))
    assert sum(weights.values()) == dp[e][0]
    return weights


class PrefixRandom(Random):
    def __init__(self, tickets):
        super().__init__(220)
        self.tickets = iter(tickets)

    def randrange(self, stop):
        ticket = next(self.tickets, None)
        if ticket is None:
            return super().randrange(stop)
        assert 0 <= ticket < stop
        return ticket


def run():
    total_hnfs = 0
    for p,e,d in ((2,1,2), (2,2,2), (2,3,2), (3,2,2), (2,2,3), (2,1,4), (3,1,3)):
        observed = Counter(profile_of(h, p, e) for h in hnfs(p**e, d))
        assert observed == weights_of(p, e, d)
        total_hnfs += sum(observed.values())

    # Force every GL_2(Z/4) element once for every possible subgroup type.
    # All HNFs of each type must occur with exactly the same multiplicity.
    orbit_visits = 0
    for ranks, expected in weights_of(2, 2, 2).items():
        images = Counter()
        for a,b,c,d in product(range(4), repeat=4):
            if gcd(a*d-b*c, 4) == 1:
                h = lattice_from_types(4, 2, 1, {2:ranks}, PrefixRandom((a,b,c,d)))
                images[tuple(h)] += 1
        assert len(images) == expected and len(set(images.values())) == 1
        orbit_visits += sum(images.values())

    rng = Random(109)
    envelopes = 0
    for p in (2, 3, 5, 7, 11, 13, 101):
        for e in (1, 4, 32):
            for d in (2, 4, 10):
                transitions, dp = type_counts(p, e, d)
                total, balanced = dp[e][0], _central_mass(transitions, e, d//2)
                assert balanced * p**(d+4) >= total * (p-1)**(d+4)
                if p >= 11:
                    assert (total-balanced)*(p-1) <= 4*total
                    assert balanced*p**10 >= total*(p-1)**10
                chosen = draw_type(transitions, dp, rng, exclude_central=True)
                assert chosen != (d//2,)*e and list(chosen) == sorted(chosen, reverse=True)
                envelopes += 1

    # All subgroups over a composite give the independent product of local laws.
    for n, factors in ((6, ((2,1),(3,1))), (12, ((2,2),(3,1)))):
        tables = [weights_of(p,e,2) for p,e in factors]
        actual = Counter(tuple(profile_of(h,p,e) for p,e in factors) for h in hnfs(n,2))
        expected = {rs: prod(t[r] for t,r in zip(tables,rs)) for rs in product(*tables)}
        assert actual == expected
        total_hnfs += sum(actual.values())

    # An actual HNF at mixed prime powers whose valuation ratios have the
    # same integer part: v_3(index)/2=3/2 and v_5(index)/3=5/3.
    n = 3**2 * 5**3
    same_floor = lattice_from_types(n, 3, 1, {3:(2,1), 5:(2,1,1)}, rng)
    assert abs(same_floor.det()) == 3**3 * 5**5
    assert factor_from_lattice(n, same_floor) == 75
    assert factor_from_index(n, n**2) is None
    assert factor_from_index(49, 1) == 7

    # On odd-dimensional squarefree inputs, verify success by summing all
    # local masses. Duality pairs the index exponent T with d-T.
    for d in (1, 3, 5):
        for p,q in ((2,3), (3,5), (7,11)):
            a, b = weights_of(p,1,d), weights_of(q,1,d)
            successes = sum(u*v for (r,),u in a.items() for (s,),v in b.items()
                            if factor_from_index(p*q, p**(d-r)*q**(d-s)))
            assert 2*successes >= sum(a.values()) * sum(b.values())

    # A 2D column HNF [[a,b],[0,c]] contains nZ² exactly when a|n,
    # c|n, and ac|bn. Count the index readouts directly on this fixture.
    n = 19 * 29
    divisors = (1, 19, 29, n)
    results = Counter(factor_from_index(n, a*c) for a in divisors for c in divisors
                      for b in range(a) if (b*n) % (a*c) == 0)
    assert sum(results.values()) == 704 and results[None] == 602

    # Two forced sparse hits select the zero subgroup or a cyclic subgroup
    # of order 11 at 11². The latter cannot be a direct summand of (Z/121)^2.
    for last, profile, divisor in ((0,(0,0),121), (1,(1,0),11)):
        forced = sample_sublattice(11**2*13, 2, PrefixRandom([11] + [1]*9 + [0,0,last]))
        assert profile_of(forced,11,2) == profile and profile_of(forced,13,1) == (1,)
        assert factor_from_lattice(11**2*13, forced) == divisor
    assert sample_sublattice(1,4,rng) == eye(4)
    for n,d in ((4,2), (72,4), (11**7*13**2,6), (1000003*1000033,4),
                ((2**127-1)*(2**127+29),6)):
        for _ in range(4):
            lattice = sample_sublattice(n, d, rng)
            assert all(x.is_Integer for x in n*lattice.inv())
            assert 0 < lattice.det() <= n**d
    print(f"{total_hnfs} exhaustive HNFs, {orbit_visits} equal-orbit GL visits")
    print(f"{envelopes} local type/envelope checks, including p=2 and e=32")
    print("N=551: 102 of the 704 intermediate plane lattices split by index")
    print("9 odd-dimensional factor laws, mixed prime-power index and sparse branch passed")


if __name__ == "__main__":
    run()
