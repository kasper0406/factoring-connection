"""Decode two fixed binomial reports with sparse joint errors in every long interval."""

from fractions import Fraction
from math import isqrt
from typing import NamedTuple

if __package__:
    from .fixed_binomial_streams import fixed_menu
else:
    from fixed_binomial_streams import fixed_menu


class _State(NamedTuple):
    x: int
    y: int
    tail: int
    debt: int
    minimum: int


def _observed_pair(x, y, bit, h):
    return (bool(x & y & bit), bool((h*x) & (h*y) & bit))


def _seed_candidates(n, reports, h, seed, t, depth, bottom, epsilon, stats):
    """All surviving final residues for this seed; drop a seed if any width exceeds 2^t."""
    seed_depth = min(t, depth)
    modulus = 1 << seed_depth
    partner = n*pow(seed, -1, modulus) % modulus
    a, d = epsilon.numerator, epsilon.denominator
    tail = 0
    for j in range(bottom, seed_depth):
        mistake = _observed_pair(seed,partner,1 << j,h) != reports[j]
        tail = (tail << 1) | mistake
    # F(l)=d*errors_before_l-a*l, counting only observed levels from bottom.
    length = max(0, seed_depth-bottom)
    states = [_State(seed, partner, tail, d*tail.bit_count()-a*length, 0)]
    cap = 1 << seed_depth
    tail_mask = cap-1
    for j in range(seed_depth, depth):
        bit = 1 << j
        length = j-bottom+1
        following = []
        for state in states:
            # Lift the partner for an unchanged first coordinate; flipping x flips y as well.
            partner_high = (n-state.x*state.y) & bit
            lifted_y = state.y+partner_high
            for x,y in ((state.x,lifted_y),(state.x+bit,lifted_y ^ bit)):
                stats["nodes_tested"] += 1
                mistake = _observed_pair(x,y,bit,h) != reports[j]
                tail = ((state.tail << 1) | mistake) & tail_mask
                debt = state.debt+d*mistake-a
                minimum = state.minimum
                if length >= t:
                    # The newly eligible start is l-t. The last t bits recover F(l-t).
                    eligible = debt-d*tail.bit_count()+a*t
                    minimum = min(minimum, eligible)
                    if debt > minimum:
                        continue
                following.append(_State(x,y,tail,debt,minimum))
                if len(following) > cap:
                    stats["seeds_abandoned"] += 1
                    stats["max_width_seen"] = max(stats["max_width_seen"], len(following))
                    return []
        stats["max_width_seen"] = max(stats["max_width_seen"],len(following))
        states = following
        if not states:
            break
    return [state.x for state in states]


def factor_by_noisy_fixed_streams(n, readings, *, multiplier=3, seed_bits,
                                  error_fraction=Fraction(1,100), stats=None):
    """On n=pq, p<q<2p, verify a factor consistent with a noisy pair at each level.

    Count one error at a level where either or both reports differ. Every interval
    of at least seed_bits observed levels must have at most floor(epsilon*length)
    errors. An overfull seed is abandoned independently. Raises ValueError if no
    retained branch gives a divisor. Reads only the fixed supplied reports.
    """
    if multiplier not in (3,5):
        raise ValueError("expected multiplier 3 or 5")
    bottom = (4*multiplier-1).bit_length()-1
    if not isinstance(seed_bits,int) or seed_bits < bottom:
        raise ValueError("seed_bits must reach the first two-report level")
    epsilon = Fraction(error_fraction)
    if not 0 <= epsilon <= 1:
        raise ValueError("error_fraction must be between zero and one")
    levels = [j for j,h,_ in fixed_menu(n,multiplier) if h == multiplier]
    lower, upper = isqrt((n-1)//2)+1, isqrt(n)
    depth = lower.bit_length()-1
    reports = {}
    for j in levels:
        if j not in readings:
            raise ValueError(f"missing report pair at level {j}")
        pair = tuple(readings[j])
        if len(pair) != 2 or any(value not in (False,True) for value in pair):
            raise ValueError(f"expected two Boolean reports at level {j}")
        reports[j] = pair
    work = {} if stats is None else stats
    work.update(depth=depth, bottom=bottom, report_pairs=len(reports), seeds_tested=0,
                seeds_abandoned=0, max_width_seen=1, nodes_tested=0, divisions=0)
    stride = 1 << depth
    for seed in range(1, 1 << min(seed_bits,depth), 2):
        work["seeds_tested"] += 1
        residues = _seed_candidates(n,reports,multiplier,seed,seed_bits,depth,bottom,epsilon,work)
        for residue in residues:
            candidate = lower + (residue-lower) % stride
            if candidate <= upper:
                work["divisions"] += 1
                if n % candidate == 0:
                    return candidate
    raise ValueError("no verified balanced divisor survived the interval and width bounds")
