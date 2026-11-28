"""Exact algebraic-law checks and input-only trials for the square-root-size sampler."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import isqrt, prod
from random import Random

if __package__:
    from . import sqrt_weighted_sublattice as sampler
    from .check_determinant_lattices import rank_mod
else:
    import sqrt_weighted_sublattice as sampler
    from check_determinant_lattices import rank_mod


def bracket_sign(vector, p, q):
    """Independent rigorous test using intervals for square roots."""
    if not any(vector):
        return 0
    bits = 4
    while True:
        scale = 1 << bits
        low = high = vector[0]*scale
        for coefficient, radicand in zip(vector[1:], (p, q, p*q)):
            floor = isqrt(radicand*scale*scale)
            low += min(coefficient*floor, coefficient*(floor+1))
            high += max(coefficient*floor, coefficient*(floor+1))
        if low > 0 or high < 0:
            return 1 if low > 0 else -1
        bits *= 2


def explicit_subspaces(p, d):
    zero = (0,)*d
    initial = frozenset([zero])
    spaces, pending = {initial: 0}, [initial]
    vectors = list(product(range(p), repeat=d))
    while pending:
        space = pending.pop()
        for v in vectors:
            if v in space:
                continue
            enlarged = frozenset(tuple((x+i*y) % p for x, y in zip(a, v))
                                 for a in space for i in range(p))
            if enlarged not in spaces:
                spaces[enlarged] = spaces[space]+1
                pending.append(enlarged)
    return spaces


def gaussian(p, d, r):
    return int(prod(Fraction(p**(d-i)-1, p**(r-i)-1) for i in range(r)))


def check_rank_law(p, q, d, hit):
    buckets = sampler._rank_buckets(p, q, d, hit)
    middle = (d+1)//2
    counts, normalizer = {}, [0]*4
    for r, s in product(range(d+1), repeat=2):
        weight = [0]*4
        i = r%2 + 2*(s%2)
        weight[i] = p**(r//2)*q**(s//2)*gaussian(p,d,r)*gaussian(q,d,s)
        counts[r,s] = weight
        normalizer[i] += weight[i]
    assert buckets[-1][1] == tuple(hit.numerator*x for x in normalizer)
    previous = (0,)*4
    per_pair = {}
    for ranks, cumulative in buckets:
        delta = [x-y for x, y in zip(cumulative, previous)]
        assert bracket_sign(delta, p, q) > 0
        if ranks == (middle, middle):
            delta = [x+(hit.denominator-hit.numerator)*z for x, z in zip(delta, normalizer)]
        assert delta == [hit.denominator*x for x in counts[ranks]]
        multiplicity = gaussian(p,d,ranks[0])*gaussian(q,d,ranks[1])
        per_pair[ranks] = tuple(Fraction(x, hit.denominator*multiplicity) for x in delta)
        previous = cumulative
    return per_pair


def replay_buckets(p, q, hit):
    """Drive the actual bit-refinement routine to every rank pair in dimension three."""
    buckets = sampler._rank_buckets(p, q, 3, hit)
    scale = 1 << 256
    roots = (scale, *(isqrt(r*scale*scale) for r in (p, q, p*q)))
    total = buckets[-1][1]
    approximate_total = sum(a*b for a, b in zip(total, roots))

    class Bytes:
        def __init__(self, ticket):
            self.values = iter(ticket.to_bytes(32, 'big'))

        def randrange(self, size):
            assert size == 256
            return next(self.values)

    previous = (0,)*4
    for ranks, current in buckets:
        midpoint = sum((a+b)*root for a, b, root in zip(previous, current, roots))
        ticket = scale*midpoint//(2*approximate_total)
        # Verify the chosen dyadic cell lies strictly inside with independent brackets.
        assert bracket_sign([ticket*z-scale*x for z, x in zip(total, previous)], p, q) > 0
        assert bracket_sign([scale*x-(ticket+1)*z for x, z in zip(current, total)], p, q) > 0
        assert sampler._sample_ranks(p, q, 3, hit, Bytes(ticket)) == ranks
        previous = current
    return len(buckets)


class WindowCenters(Random):
    """Replay the chosen public window centers and use a seeded RNG afterwards."""
    def __init__(self, modulus, centers):
        super().__init__(571)
        self.modulus = modulus
        self.centers = iter(centers)

    def randrange(self, *args):
        if args == (self.modulus,):
            center = next(self.centers, None)
            if center is not None:
                return center
        return super().randrange(*args)


def draw(n, d, rng, p, q):
    stats = {}
    lattice = sampler.sample_sqrt_weighted(n, d, rng, stats=stats)
    ranks = tuple(rank_mod(lattice.tolist(), r) for r in (p, q))
    assert int(lattice.det()) == p**(d-ranks[0])*q**(d-ranks[1])
    assert stats['windows'] <= sampler.WINDOWS
    assert stats['gcd_offsets'] <= sampler.WINDOWS*max(8, n.bit_length()+2)
    return ranks, stats


def run():
    rng = Random(960417)
    for p, q in ((2, 3), (257, 263), (1009, 1321)):
        for _ in range(160):
            vector = [rng.randrange(-1000, 1001) for _ in range(4)]
            assert sampler._field_sign(vector, p, q) == bracket_sign(vector, p, q)
    assert sampler._field_sign([0]*4, 2, 3) == 0
    space2, space3 = explicit_subspaces(2, 3), explicit_subspaces(3, 3)
    law = check_rank_law(2, 3, 3, Fraction(1))
    for a, b in product(space2, space3):
        r, s = space2[a], space3[b]
        vector = [0]*4
        vector[r%2 + 2*(s%2)] = 2**(r//2)*3**(s//2)
        assert law[r,s] == tuple(vector) and (vector[r%2+2*(s%2)]**2
                                              * 2**(r%2)*3**(s%2)) == len(a)*len(b)
    print(f"{len(space2)*len(space3)} subgroup pairs matched exact square-root weights; "
          "481 exact algebraic signs checked")

    extra = 0
    for p, q in ((257, 263), (263, 509), (1009, 1321), (65537, 65539)):
        for d in (1, 3, 9):
            hit = sampler._hit_probability(p, q)
            extra += len(check_rank_law(p, q, d, hit))
    print(f"{extra} further conditional rank buckets matched unconditioned masses")
    assert replay_buckets(2, 3, Fraction(1)) == 16
    assert replay_buckets(257, 263, sampler._hit_probability(257, 263)) == 16

    p, q = 257, 263
    n = p*q
    proper = sum(1 < sampler.prime_window_divisor(n, x) < n for x in range(n))
    assert sampler._hit_probability(p, q) == 1-(1-Fraction(proper,n))**sampler.WINDOWS
    print(f"{proper}/{n} exhaustive centers gave a proper mark")
    observed, window_total = Counter(), 0
    for _ in range(96):
        ranks, stats = draw(n, 3, rng, p, q)
        observed[ranks] += 1
        window_total += stats['windows']
    assert len(observed) > 1
    print(f"96 input-only draws at {n}: {dict(sorted(observed.items()))} using {window_total} windows")
    witness = {}
    lattice = sampler.sample_sqrt_weighted(n, 3, Random(10), stats=witness)
    assert lattice.tolist() == [[67591, 38661, 30721], [0, 263, 116], [0, 0, 1]]
    assert witness['windows'] == 32

    ranks, stats = draw(15, 5, rng, 3, 5)
    assert stats['factor_found'] and stats['windows'] == 0
    p, q = 1009, 1321
    n = p*q
    _, hit_work = draw(n, 5, WindowCenters(n, [p+sampler.window_radius(n)]), p, q)
    assert hit_work['windows'] == hit_work['lattice_windows'] == 1 and hit_work['factor_found']
    ranks, miss_work = draw(n, 5, WindowCenters(n, [0]*sampler.WINDOWS), p, q)
    assert ranks == (3, 3) and miss_work['windows'] == sampler.WINDOWS and not miss_work['factor_found']
    print(f"Certified lattice hit used dimension {hit_work['max_lattice_dimension']} and "
          f"{hit_work['gcd_offsets']} gcd offsets; missed draw used exactly {miss_work['windows']} windows")
    p, q = 16777213, 20000003
    n = p*q
    _, large = draw(n, 3, WindowCenters(n, [p+sampler.window_radius(n)]), p, q)
    assert large['windows'] == large['lattice_windows'] == 1 and large['factor_found']
    print(f"{n.bit_length()}-bit planted center: radius {sampler.window_radius(n)}, "
          f"lattice dimension {large['max_lattice_dimension']}, {large['gcd_offsets']} gcd offset")
    for n, d in ((3, 1), (15, 0), (15, 2), (9, 1), (21, 1)):
        try:
            sampler.sample_sqrt_weighted(n, d, rng)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid input accepted")


if __name__ == "__main__":
    run()
