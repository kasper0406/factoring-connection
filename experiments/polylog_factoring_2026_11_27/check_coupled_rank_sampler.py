"""Exact laws and input-only tests for the coupled odd-dimensional samplers."""

from collections import Counter
from fractions import Fraction
from math import gcd, isqrt
from random import Random
from unittest.mock import patch

from sympy import GF, ZZ, eye
from sympy.matrices.normalforms import smith_normal_form

if __package__:
    from . import coupled_rank_sampler as sampler
    from .sublattice_sampler import factor_from_lattice
else:
    import coupled_rank_sampler as sampler
    from sublattice_sampler import factor_from_lattice


class FixedDraws:
    def __init__(self, *draws):
        self.draws = iter(draws)

    def randrange(self, *args):
        value = next(self.draws)
        assert value in range(*args)
        return value


def rank_counts(p, d, projector):
    # Count subspaces by the Gaussian binomial recurrence; each has
    # p^(r(d-r)) complements and hence that many projectors onto it.
    row = [1]
    for m in range(1, d+1):
        row = [1] + [row[r-1] + p**r*row[r] for r in range(1,m)] + [1]
    return [v*p**(r*(d-r)) if projector else v for r,v in enumerate(row)]


def radius_law(p, q, d, projector):
    n = p*q
    lo, hi = isqrt((n-1)//2)+1, isqrt(2*n)
    count = hi-lo+1
    tables = [rank_counts(r,d,projector) for r in (p,q)]
    probabilities = [[Fraction(v,sum(row)) for v in row] for row in tables]
    law = Counter()
    if d == 1:
        assert sampler._sample_radii(n,d,projector,FixedDraws()) == {}
        return {(0,0): Fraction(1)}, probabilities
    tickets = list(sampler._gated_tickets(p,q,d,count,projector))
    total = sum(weight for _,_,weight in tickets)
    law[0,0] = Fraction(count-2, count)
    offset = 0
    for i,j,weight in tickets:
        if weight:
            for candidate in (p,q):
                expected = {r:t for r,t in ((p,i),(q,j)) if t}
                assert sampler._sample_radii(n,d,projector,FixedDraws(candidate,offset)) == expected
                assert sampler._sample_radii(n,d,projector,FixedDraws(candidate,offset+weight-1)) == expected
            law[i,j] += Fraction(2*weight, count*total)
        offset += weight
    for i in range(d//2+1):
        for j in range(d//2+1):
            assert law[i,j] == 4*probabilities[0][d//2-i]*probabilities[1][d//2-j]
    return law, probabilities


def exact_check(p, q, d, projector):
    n = p*q
    radius = sampler.window_radius(n)
    width = 2*radius+1
    masks, representatives = Counter(), {}
    for x in range(n):
        divisor = sampler.prime_window_divisor(n, x, radius=radius)
        expected = 1
        for r in (p,q):
            if min(x % r, -x % r) <= radius:
                expected *= r
        assert divisor == expected
        masks[divisor] += 1
        representatives[divisor] = x
    assert masks == {1:(p-width)*(q-width), p:width*(q-width),
                     q:(p-width)*width, n:width*width}

    radii, marginals = radius_law(p,q,d,projector)
    joint = Counter()
    for (i,j), mass in radii.items():
        corrections = {r:t for r,t in ((p,i),(q,j)) if t}
        with patch.object(sampler, "_sample_radii", side_effect=lambda *args: corrections.copy()):
            for divisor, count in masks.items():
                for bit in (0,1):
                    rng = FixedDraws(representatives[divisor], bit)
                    default, ranks = sampler._ranks(n,d,projector,rng)
                    a,b = ranks.get(p,default), ranks.get(q,default)
                    assert (a%2, b%2) == (bit ^ (divisor%p == 0), bit ^ (divisor%q == 0))
                    joint[a,b] += mass*Fraction(count,2*n)
    for r in range(d+1):
        assert sum(joint[r,j] for j in range(d+1)) == marginals[0][r]
        assert sum(joint[j,r] for j in range(d+1)) == marginals[1][r]
    # Within every rank pair both bases are uniform and independent, so
    # the rank-table TV equals TV for full matrices or subspace pairs.
    tv = sum(abs(joint[i,j] - px*qy)
             for i,px in enumerate(marginals[0]) for j,qy in enumerate(marginals[1]))/2
    split = sum(mass for (a,b), mass in joint.items() if a%2 != b%2)
    any_rank_split = sum(mass for (a,b), mass in joint.items() if a != b)
    assert tv == Fraction(abs((p-2*width)*(q-2*width)), 2*n)
    assert split == Fraction(width*(p+q-2*width), n)
    same_radius = sum(mass for (i,j),mass in radii.items() if i == j)
    assert any_rank_split == 1 - (1-split)*same_radius
    return tv, split, any_rank_split


def public_checks():
    calls = 0
    mixed = 0
    recovered = 0
    for p,q,d in ((2,3,3),(3,5,1),(3,5,5),(67,71,3),(1009,1013,3)):
        n = p*q
        for seed in range(8):
            for projector in (False,True):
                default, ranks = sampler._ranks(n,d,projector,Random(seed))
                expected = [ranks.get(r,default) for r in (p,q)]
                if projector:
                    mat = sampler.coupled_projector(n,d,Random(seed))
                    assert all(int(v)%n == 0 for v in mat*mat-mat)
                    smith = smith_normal_form(mat, domain=ZZ)
                    found = any(1 < gcd(n,int(z)) < n for z in smith.diagonal())
                    if n % 2:
                        parity_gcd = gcd(n, int((eye(d)-2*mat).det())-1)
                        assert (1 < parity_gcd < n) == (expected[0]%2 != expected[1]%2)
                else:
                    mat = sampler.coupled_sublattice(n,d,Random(seed))
                    assert mat.det() == p**(d-expected[0])*q**(d-expected[1])
                    found = factor_from_lattice(n,mat) is not None
                assert [mat.to_DM().convert_to(GF(r)).rank() for r in (p,q)] == expected
                assert found == (expected[0] != expected[1])
                calls += 1
                mixed += expected[0]%2 != expected[1]%2
                recovered += found
    print(f"{calls} seeded public draws, {mixed} parity splits, {recovered} rank splits")


if __name__ == "__main__":
    for p,q,d in ((2,3,3),(3,5,1),(3,5,5),(7,11,3),(67,71,3)):
        all_splits = []
        for projector in (False,True):
            tv, split, rank_split = exact_check(p,q,d,projector)
            all_splits.append(str(rank_split))
        print(f"N={p*q}, d={d}, width={2*sampler.window_radius(p*q)+1}, "
              f"full TV={tv}, opposite parity={split}, any rank split {' / '.join(all_splits)}")
    public_checks()
