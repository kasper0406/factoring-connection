"""Exact HNF and type-law controls for general square-root-weighted lattices."""

from collections import Counter
from itertools import product
from math import isqrt
from random import Random
from unittest.mock import patch

from sympy import eye

if __package__:
    from . import weighted_sublattice_general as sampler
    from .check_sublattice_sampler import hnfs, profile_of, weights_of
else:
    import weighted_sublattice_general as sampler
    from check_sublattice_sampler import hnfs, profile_of, weights_of


def add(items):
    a = b = 0
    for x,y in items:
        a,b = a+x,b+y
    return a,b


def multiply(left, right, p):
    a,b = left
    x,y = right
    return a*x+p*b*y, a*y+b*x


def radius_mass(profile, count, p):
    order_power = sum(profile)
    coefficient = count*p**(order_power//2)
    return (0,coefficient) if order_power%2 else (coefficient,0)


def force_choice(masses, p, target):
    """Use a dyadic interval contained in the specified algebraic bucket."""
    total = add(masses)
    lower = add(masses[:target])
    upper = add(masses[:target+1])
    assert sampler._quadratic_sign(*masses[target],p) > 0

    def compare(bound, ticket, scale):
        return sampler._quadratic_sign(bound[0]*scale-total[0]*ticket,
                                       bound[1]*scale-total[1]*ticket,p)

    # Approximate sqrt(p) with integer square roots to find a candidate
    # ticket. Exact endpoint checks certify it before the sampler sees it.
    bits = 8
    while True:
        scale = 1 << bits
        root = isqrt(p << (2*bits))
        denominator = total[0]*scale + total[1]*root
        numerator = (lower[0]+upper[0])*scale + (lower[1]+upper[1])*root
        if denominator:
            ticket = scale*numerator//(2*denominator)
            if (0 <= ticket < scale and compare(lower,ticket,scale) <= 0 and
                    compare(upper,ticket+1,scale) >= 0):
                break
        bits *= 2
    stream = iter(ticket.to_bytes(bits//8,"big"))

    class Bytes:
        def randrange(self, stop):
            assert stop == 256
            return next(stream)

    result = sampler._choose(masses,p,Bytes())
    assert result == target
    return result


def local(p,e,d,enumerate_hnfs):
    center = (d+1)//2
    counts = (Counter(profile_of(h,p,e) for h in hnfs(p**e,d))
              if enumerate_hnfs else weights_of(p,e,d))
    masses = {r:radius_mass(r,c,p) for r,c in counts.items()}
    central = (center,)*e
    transitions, dp, tails = sampler._tables(p,e,d)
    total = add(masses.values())
    assert total == dp[e][0]
    assert masses[central] == sampler._step(tails[e-1],transitions[0][center],center%2,p)
    noncentral = {r[::-1]:v for r,v in masses.items() if r != central}
    exceptional = add(noncentral.values())
    h = sampler.sqrt_mark_probability(p)
    assert sampler._quadratic_sign(total[0]*h.numerator-exceptional[0]*h.denominator,
                                   total[1]*h.numerator-exceptional[1]*h.denominator,p) >= 0

    # Compare every conditional DP branch with the exact prefix sum
    # of full HNF weights, excluding the single central profile.
    prefixes = {r[:length] for r in noncentral for length in range(e)}
    for prefix in prefixes:
        s = prefix[-1] if prefix else 0
        matched = all(r == center for r in prefix)
        choices = sampler._conditional_masses(p,d,transitions,dp,tails,s,e-len(prefix),matched)
        direct = [add(v for r,v in noncentral.items() if r[:len(prefix)+1] == prefix+(t,))
                  for t in range(s,d+1)]
        for a,b in zip(choices,direct):
            assert multiply(a,add(direct),p) == multiply(b,add(choices),p)

    # Exercise every positive-width branch through the real dyadic reader.
    # This also verifies the public gate thinning and the reversed column order.
    real_choose = sampler._choose
    for goal in masses:
        sequence = [int(goal != central)]
        if goal != central:
            path = (0,) + goal[::-1]
            sequence.extend(r-s for s,r in zip(path,path[1:]))
        choices = iter(sequence)
        gate = True

        def choose(weights, prime, _rng):
            nonlocal gate
            if gate:
                # h times the marked conditional exception probability
                # must be exactly the original noncentral mass / total.
                left = multiply(weights[1],total,p)
                right = multiply(add(weights),exceptional,p)
                assert tuple(x*h.numerator for x in left) == tuple(x*h.denominator for x in right)
                gate = False
            with patch.object(sampler,"_choose",real_choose):
                return force_choice(weights,prime,next(choices))

        with patch.object(sampler,"_choose",choose):
            output = sampler._marked_profile(p,e,d,None)
        assert output == (None if goal == central else goal)
    return len(counts), sum(counts.values()) if enumerate_hnfs else 0


def crt_check():
    # Components of sqrt(|H|) in the basis (1, sqrt(2), sqrt(3), sqrt(6)).
    joint = [0]*4
    by_type = Counter()
    for h in hnfs(12,2):
        size = 12**2//int(h.det())
        for i,squarefree in enumerate((1,2,3,6)):
            root = isqrt(size//squarefree)
            if squarefree*root**2 == size:
                joint[i] += root
                by_type[profile_of(h,2,2),profile_of(h,3,1),i] += root
                break
        else:
            raise AssertionError("unexpected square-free part")
    a,b = sampler._tables(2,2,2)[1][2][0]
    c,d = sampler._tables(3,1,2)[1][1][0]
    assert joint == [a*c,b*c,a*d,b*d]
    for left,lcount in weights_of(2,2,2).items():
        for right,rcount in weights_of(3,1,2).items():
            x,y = radius_mass(left,lcount,2)
            v,w = radius_mass(right,rcount,3)
            assert [by_type[left,right,i] for i in range(4)] == [x*v,y*v,x*w,y*w]
    return tuple(joint)


if __name__ == "__main__":
    paths = boxes = 0
    for p,e,d in ((2,1,2),(2,2,2),(3,2,1),(2,1,3),(3,1,3),(5,1,2),
                  (2,3,3),(67,2,4),(257,2,3),(1009,3,3)):
        full = p**e <= 9 and p**e != 8
        checked, count = local(p,e,d,full)
        paths += checked
        boxes += count
    print(f"{paths} full algebraic type paths, {boxes} independently enumerated local HNFs")
    print("N=12, d=2 direct weights by square-free part:", crt_check())
    # Stress the gate bound at higher exponents and both parity choices.
    for p,e,d in product((2,7,61,67,127,1009), (1,8,27), (1,2,5,8)):
        weights,dp,tails = sampler._tables(p,e,d)
        m = (d+1)//2
        a,b = dp[e][0]
        x,y = sampler._step(tails[e-1],weights[0][m],m%2,p)
        h = sampler.sqrt_mark_probability(p)
        assert sampler._quadratic_sign(a*h.numerator-(a-x)*h.denominator,
                                       b*h.numerator-(b-y)*h.denominator,p) >= 0
    rng = Random(473)
    assert sampler.sample_weighted_sublattice(1,3,rng) == eye(3)
    for n,d in ((4,2),(9,3),(12,2),(72,3),(67**2,4),(101*127,3)):
        for _ in range(4):
            lattice = sampler.sample_weighted_sublattice(n,d,rng)
            assert lattice == lattice.upper_triangular()
            assert all(x.is_Integer for x in n*lattice.inv())
    print("72 envelope checks, modulus one, and 24 native HNF draws passed")
