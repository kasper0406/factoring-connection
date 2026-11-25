"""First-hit searches. Input is the public quadratic form, prime list and tie coins."""
import math
from collections import Counter
from itertools import combinations

from common import HIGH, K, LOW, exact_factors, exclusion, hatch, matched_divisor

STRIDE = K + 1


def verify(form, u, v, primes, product, excluded, counts, divisor=None):
    counts["norm_tests"] += 1
    value = hatch.norm(form, u, v)
    if not hatch.smooth(value, product):
        return None
    counts["full_factorizations"] += 1
    factors = exact_factors(value, primes)
    if divisor is None:
        counts["divisor_searches"] += 1
        divisor = matched_divisor(factors, excluded)
        if divisor is None:
            return None
    assert LOW < divisor <= HIGH and value % divisor == 0 and math.gcd(divisor, excluded) == 1
    return [u, v, divisor, factors]


def prefix(values):
    result = [0] * STRIDE**2
    for v in range(K):
        row_sum = 0
        for u in range(K):
            row_sum += values[v*K + u]
            i = (v+1)*STRIDE + u+1
            result[i] = result[i-STRIDE] + row_sum
    return result


def box_sum(table, box):
    x, y, X, Y = box
    return table[Y*STRIDE+X] - table[y*STRIDE+X] - table[Y*STRIDE+x] + table[y*STRIDE+x]


def pair_tables(form, primes, explicit, counts):
    """Two independent accumulators for the exact second binomial moment."""
    primitive = bytearray(math.gcd(u, v) == 1 for v in range(1, K+1) for u in range(1, K+1))
    counts["setup_gcds"] = K*K
    weights = [0] * (K*K)
    degree = [0] * (K*K) if not explicit else None
    masks, good = [], []
    excluded = exclusion(form)
    rows = progressions = incidences = rejected = skipped_rows = 0
    for p in primes:
        if excluded % p == 0:
            continue
        counts["root_calls"] += 1
        roots = hatch.line_roots(form, p)
        if not roots:
            continue
        assert len(roots) == 2 and roots[0] != roots[1]
        good.append(p)
        mask = 0
        for v in range(1, K+1):
            rows += 1
            if v % p == 0:
                skipped_rows += 1
                continue
            progressions += 2
            for r in roots:
                for u in range((r*v - 1) % p + 1, K+1, p):
                    idx = (v-1)*K + (u-1)
                    incidences += 1
                    if not primitive[idx]:
                        rejected += 1
                        continue
                    if explicit:
                        mask |= 1 << idx
                    else:
                        weights[idx] += degree[idx]
                        degree[idx] += 1
        if explicit:
            masks.append(mask)

    if explicit:
        occupied = 0
        for first, second in combinations(masks, 2):
            intersection = first & second
            if intersection:
                occupied += 1
            while intersection:
                bit = intersection & -intersection
                weights[bit.bit_length()-1] += 1
                intersection ^= bit
        counts["pair_mask_intersections"] = math.comb(len(masks), 2)
        counts["nonempty_pair_masks"] = occupied
    # Conditional on primitivity, each distinct prime has two allowed slopes out
    # of p+1. Scale the prediction once for exact integer comparisons everywhere.
    denom = math.lcm(*(p+1 for p in good))
    nums = [2 * (denom // (p+1)) for p in good]
    numerator = (sum(nums)**2 - sum(t*t for t in nums)) // 2
    counts.update(good_primes=len(good), prime_rows=rows, v_multiple_skips=skipped_rows,
                  root_progressions=progressions, incidence_gcd_mask_reads=incidences,
                  rejected_nonprimitive_incidence=rejected, prefix_cells=2*K*K)
    wt = prefix(weights)
    counts["primitive_pair_incidences"] = wt[-1]
    return wt, prefix(primitive), numerator, denom**2, primitive


def children(box):
    x, y, X, Y = box
    if X-x >= Y-y:
        mid = (x+X)//2
        return (x, y, mid, Y), (mid, y, X, Y)
    mid = (y+Y)//2
    return (x, y, X, mid), (x, mid, X, Y)


def spatial(form, primes, tie_coins, mode, audit=None):
    counts = Counter()
    product = hatch.primorial(primes)
    excluded = exclusion(form)
    tables = None if mode == "blind" else pair_tables(form, primes, mode == "pairs", counts)
    if audit is not None and tables is not None:
        audit[mode] = tables
    pending = [(1, (0, 0, K, K))]
    while pending:
        node, box = pending.pop()
        counts["nodes"] += 1
        x, y, X, Y = box
        if X-x <= 4 and Y-y <= 4:
            counts["leaves"] += 1
            for v in range(y+1, Y+1):
                for u in range(x+1, X+1):
                    if tables is None:
                        counts["leaf_gcds"] += 1
                        coprime = math.gcd(u, v) == 1
                    else:
                        counts["leaf_gcd_mask_reads"] += 1
                        coprime = tables[4][(v-1)*K + u-1]
                    if coprime:
                        hit = verify(form, u, v, primes, product, excluded, counts)
                        if hit is not None:
                            return hit, counts
            counts["dead_leaves"] += 1
            continue
        left, right = children(box)
        prefer_right = bool((tie_coins >> node) & 1)
        if tables is not None:
            pairs, pops, num, den, _ = tables
            delta = (box_sum(pairs, right) - box_sum(pairs, left)) * den
            delta -= (box_sum(pops, right) - box_sum(pops, left)) * num
            counts["children_scored"] += 2
            if delta:
                prefer_right = delta > 0
            else:
                counts["score_ties"] += 1
        kids = [(2*node, left), (2*node+1, right)]
        pending.extend(kids if prefer_right else reversed(kids))
    return None, counts


def matched(form, primes):
    """Archived forced_search(use_reconstruction=False) with immediate first-hit stop."""
    counts = Counter()
    product = hatch.primorial(primes)
    excluded = exclusion(form)
    seen = set()
    branches = scans = gcds = valid = 0
    hit = None
    for d, residues in hatch.composite_roots(form, primes):
        counts["moduli"] += 1
        counts["generated_slopes_in_yielded_moduli"] += len(residues)
        for r in residues:
            branches += 1
            for v in range(1, K+1):
                scans += 1
                u = r*v % d
                if 1 <= u <= K:
                    gcds += 1
                    if math.gcd(u, v) != 1:
                        continue
                    valid += 1
                    point = (u, v)
                    if point in seen:
                        continue
                    seen.add(point)
                    hit = verify(form, u, v, primes, product, excluded, counts, d)
                    if hit is not None:
                        break
            if hit is not None:
                break
        if hit is not None:
            break
    counts.update(branches=branches, failed_slopes=branches-int(hit is not None),
                  lattice_v_scans=scans, candidate_gcds=gcds, primitive_incidences=valid)
    return hit, counts
