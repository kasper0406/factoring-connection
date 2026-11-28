"""Rank candidate prime residues using possibly incorrect shared-bit observations."""

from collections import Counter
from heapq import heappop, heappush
from math import comb

if __package__:
    from .hh_residue_lists import factor_from_residue_list
else:
    from hh_residue_lists import factor_from_residue_list
from src.thmA.arith import iroot


def ranked_shared_residues(n, bits, report, *, block_bits=None, stats=None):
    """Yield (residue, inverse mass, errors) in order; ignore the known units bit.

    Every odd residue modulo 2^bits occurs exactly once if exhausted. An integer
    key combines the arity of each exact-AND block with a universal error prior.
    """
    if n % 2 == 0 or bits < 1 or not 0 <= report < 1 << bits:
        raise ValueError("expected an odd product and a report of the given width")
    if block_bits is None:
        block_bits = max(1, iroot(bits, 3))
    if block_bits < 1:
        raise ValueError("expected a positive block width")
    work = {} if stats is None else stats
    work.update(expanded=0, yielded=0, peak_queue=1)
    # Inverse mass, prefix width, first endpoint, error weight, product of arities.
    queue = [(1, 1, 1, 0, 1)]
    while queue:
        key, depth, x, errors, arities = heappop(queue)
        if depth == bits:
            work["yielded"] += 1
            yield x, key, errors
            continue
        width = min(block_bits, bits-depth)
        next_depth = depth+width
        modulus = 1 << next_depth
        report_block = (report >> depth) % (1 << width)
        children = []
        counts = Counter()
        for child in range(x, modulus, 1 << depth):
            partner = n*pow(child, -1, modulus) % modulus
            observed = (child & partner) >> depth
            children.append((child, observed))
            counts[observed] += 1
        work["expanded"] += 1
        for child, observed in children:
            new_errors = errors+(observed ^ report_block).bit_count()
            new_arities = arities*counts[observed]
            # There are next_depth-1 informative bits; a given error word of
            # weight e gets probability 1/(next_depth * C(next_depth-1,e)).
            new_key = new_arities*next_depth*comb(next_depth-1, new_errors)
            heappush(queue, (new_key, next_depth, child, new_errors, new_arities))
        work["peak_queue"] = max(work["peak_queue"], len(queue))


def factor_from_noisy_shared_bits(n, bits, report, *, block_bits=None, stats=None):
    """Split balanced n from a report, with certified completion of growing batches.

    Correctness holds for every report. Runtime depends on the rank of a true
    residue or any earlier incidental split.
    """
    work = {} if stats is None else stats
    ranking = {}
    batch, target = [], 1
    calls = 0
    for residue, key, errors in ranked_shared_residues(n, bits, report, block_bits=block_bits,
                                                      stats=ranking):
        batch.append(residue)
        if len(batch) != target and ranking["yielded"] != 1 << (bits-1):
            continue
        calls += 1
        completion = {}
        divisor = factor_from_residue_list(n, 1 << bits, batch, stats=completion)
        if divisor:
            work.update(ranking, batches=calls, last_inverse_mass=key,
                        last_errors=errors, completion=completion)
            return divisor
        batch, target = [], 2*target
    raise ValueError("no divisor compatible with the balanced semiprime promise")
