"""Compare the sparse-noise decoder with literal interval counts and exact binomial valuations."""

from fractions import Fraction
from itertools import product
from math import isqrt
from random import Random

if __package__:
    from .check_fixed_binomial_streams import divisible
    from .fixed_binomial_streams import fixed_menu
    from .noisy_fixed_streams import _seed_candidates, factor_by_noisy_fixed_streams
else:
    from check_fixed_binomial_streams import divisible
    from fixed_binomial_streams import fixed_menu
    from noisy_fixed_streams import _seed_candidates, factor_by_noisy_fixed_streams


def exact_reports(p, q, h):
    bottom = (4*h-1).bit_length()-1
    pairs = {}
    for level, _, t in fixed_menu(p*q,h):
        if level >= bottom:
            pairs.setdefault(level,[]).append(divisible(t,p) and divisible(t,q))
    return {j:tuple(pair) for j,pair in pairs.items()}


def brute_seed(n, reports, h, seed, t, depth, bottom, epsilon):
    """Recompute each complete prefix from scratch, including every sufficiently long interval."""
    seed_depth = min(t,depth)
    for level in range(seed_depth,depth+1):
        mod = 1 << level
        choices = []
        for x in range(seed,mod,1 << seed_depth):
            y = n*pow(x,-1,mod) % mod
            mismatches = []
            for j in range(bottom,level):
                true_pair = (bool((x & y) & (1 << j)), bool(((x*h) & (y*h)) & (1 << j)))
                mismatches.append(true_pair != reports[j])
            admissible = all(sum(mismatches[start:stop]) <= (stop-start)*epsilon
                             for start in range(len(mismatches))
                             for stop in range(start+t,len(mismatches)+1))
            if admissible:
                choices.append(x)
        if len(choices) > 1 << t:
            return []
    return choices


def check_paths():
    n = 157*251
    depth = (isqrt((n-1)//2)+1).bit_length()-1
    counts = 0
    for h in (3,5):
        bottom = (4*h-1).bit_length()-1
        for raw in product(product((False,True), repeat=2), repeat=depth-bottom):
            reports = dict(zip(range(bottom, depth),raw))
            for epsilon in (Fraction(0),Fraction(1,3),Fraction(1)):
                t = bottom
                for seed in range(1,1 << t,2):
                    stats = dict(nodes_tested=0, max_width_seen=1, seeds_abandoned=0)
                    got = _seed_candidates(n,reports,h,seed,t,depth,bottom,epsilon,stats)
                    assert sorted(got) == brute_seed(n,reports,h,seed,t,depth,bottom,epsilon)
                    counts += 1
    # Longer randomized prefixes test overlapping intervals with different starting levels.
    rng = Random(830)
    for _ in range(180):
        h = rng.choice((3,5))
        bottom = (4*h-1).bit_length()-1
        t = rng.randrange(bottom,bottom+3)
        depth = rng.randrange(t,t+7)
        seed = rng.randrange(1,1 << t,2)
        epsilon = Fraction(rng.randrange(4),rng.randrange(4,9))
        reports = {j:(bool(rng.randrange(2)),bool(rng.randrange(2))) for j in range(bottom,depth)}
        stats = dict(nodes_tested=0,max_width_seen=1,seeds_abandoned=0)
        got = _seed_candidates(n,reports,h,seed,t,depth,bottom,epsilon,stats)
        assert sorted(got) == brute_seed(n,reports,h,seed,t,depth,bottom,epsilon)
        counts += 1
    print(f"{counts} seed lists agree with direct interval tests and independent width caps")


def check_stored_and_cap():
    p,q = 2999,3457
    noisy = dict(enumerate(((0,0),(0,0),(1,1),(0,0),(1,0),(1,0),(0,0),(0,0)),3))
    actual = exact_reports(p,q,3)
    assert [j for j in noisy if noisy[j] != actual[j]] == [5]
    assert noisy[5] == tuple(not s for s in actual[5])
    stats = {}
    assert factor_by_noisy_fixed_streams(p*q,noisy,seed_bits=3,error_fraction=Fraction(1,3),
                                          stats=stats) == p
    assert stats["seeds_abandoned"] == 2 and stats["seeds_tested"] == 4
    assert stats["max_width_seen"] == 9
    print(f"stored two-bit flip at one position: {stats['seeds_abandoned']} seeds abandoned, factor {p}")

    # Two nearby errors straddle the seed boundary, so skipping the seed reports would accept p.
    p,q = 1000003,1000033
    reports = exact_reports(p,q,3)
    for j in (4,6):
        reports[j] = tuple(not x for x in reports[j])
    try:
        factor_by_noisy_fixed_streams(p*q,reports,seed_bits=5,error_fraction=Fraction(1,5))
    except ValueError as error:
        assert "no verified balanced divisor" in str(error)
    else:
        raise AssertionError("a seed-crossing burst should exclude the unique balanced prime")

    # With all errors allowed every seed exceeds its own cap on this input.
    try:
        factor_by_noisy_fixed_streams(p*q,reports,seed_bits=3,error_fraction=1,stats=stats)
    except ValueError as error:
        assert "no verified balanced divisor" in str(error)
    else:
        raise AssertionError("all seeds should be overfull")
    assert stats["seeds_abandoned"] == stats["seeds_tested"] == 4


def check_sparse():
    for p,q in ((1009,1559),(1000003,1000033),(850000013,1345000021)):
        for h in (3,5):
            reports = exact_reports(p,q,h)
            t = (4*h-1).bit_length()
            assert factor_by_noisy_fixed_streams(p*q,reports,multiplier=h,seed_bits=t) == p
            # The sparse fixture deliberately alternates single- and double-bit corruptions.
            flips = sorted(reports)[1::11]
            for number,j in enumerate(flips):
                first, second = reports[j]
                reports[j] = (not first, second if number % 2 else not second)
            stats = {}
            assert factor_by_noisy_fixed_streams(p*q,reports,multiplier=h,seed_bits=t,
                                                  error_fraction=Fraction(1,4),stats=stats) == p
            print(f"{(p*q).bit_length()}-bit input, h={h}, corrupted levels {flips}, "
                  f"{stats['nodes_tested']} child tests, verified factor {p}")
    # Empty report menu and depth shorter than the requested seed are valid on small inputs.
    assert factor_by_noisy_fixed_streams(15,{},seed_bits=12) == 3


def check_default_noise():
    # This planted prime falls in the first 100-bit seed, making the large-seed case cheap.
    p, q = 177*2**119+1, 2**127-1
    assert pow(5,(p-1)//2,p) == p-1  # Proth certificate for p
    value = 4
    for _ in range(125):
        value = (value*value-2) % q
    assert value == 0  # Lucas–Lehmer certificate for q
    for h in (3,5):
        reports = exact_reports(p,q,h)
        reports[50] = tuple(not bit for bit in reports[50])
        stats = {}
        assert factor_by_noisy_fixed_streams(p*q,reports,multiplier=h,seed_bits=100,stats=stats) == p
        assert stats["seeds_tested"] == stats["divisions"] == 1
        assert stats["nodes_tested"] == (98 if h==3 else 142)
        print(f"254-bit planted first seed, h={h}, one joint error at level 50, "
              f"default 1% interval bound, {stats['nodes_tested']} child tests")


if __name__ == "__main__":
    check_paths()
    check_stored_and_cap()
    check_sparse()
    check_default_noise()
