"""Factor a balanced odd semiprime given partial, possibly mistaken fixed AND readings."""

from math import comb

if __package__:
    from .hh_residue_lists import factor_from_residue_list
else:
    from hh_residue_lists import factor_from_residue_list


def _check_readings(readings):
    if any(bit not in (0, 1, None) for bit in readings):
        raise ValueError("each reading must be zero, one or an erasure (None)")


def singleton_list_bound(readings, orientations, errors, *, certified_ones=False):
    """Uniform bound on the capped endpoint list; readings start at bit one."""
    _check_readings(readings)
    if orientations < 0 or errors < 0:
        raise ValueError("the two caps must be nonnegative")
    length = len(readings) - (readings.count(1) if certified_ones else 0)
    free = min(length, readings.count(None)+orientations)
    scored = length-free
    return (1 << free)*sum(comb(scored, j) for j in range(min(scored, errors)+1))


def noisy_singleton_list(n, readings, orientations, errors, *, certified_ones=False):
    """All ordered endpoints satisfying both caps, plus the number of visited prefix states."""
    _check_readings(readings)
    if n < 1 or n % 2 == 0 or orientations < 0 or errors < 0:
        raise ValueError("expected an odd product and two nonnegative caps")
    states = [(1, 0, 0)]
    tested = 0
    for position, want in enumerate(readings, 1):
        stride = 1 << position
        following = []
        for endpoint, used_orientations, used_errors in states:
            tested += 1
            partner = n*pow(endpoint, -1, 2*stride) % (2*stride)
            unequal = bool(partner & stride)
            new_orientations = used_orientations + (unequal and want is not None)
            if new_orientations > orientations:
                continue
            # Flipping the lift adds stride to both odd factor endpoints.
            for high in (0, 1):
                shared_bit = 0 if unequal else high
                if certified_ones and want == 1 and shared_bit == 0:
                    continue
                new_errors = used_errors + (want is not None and want != shared_bit)
                if new_errors <= errors:
                    following.append((endpoint+high*stride, new_orientations, new_errors))
        states = following
        if not states:
            break
    return [endpoint for endpoint, _, _ in states], tested


def factor_from_noisy_singletons(n, readings, max_errors=None, *, certified_ones=False, stats=None):
    """Deterministic balanced completion, with optional known Hamming error cap.

    With max_errors=None, search pairs of caps in order of the proved list bound.
    If certified_ones=True, all supplied ones must be correct. Erasures are None.
    """
    readings = tuple(readings)
    _check_readings(readings)
    if n < 15 or n % 2 == 0 or (max_errors is not None and max_errors < 0):
        raise ValueError("expected an odd balanced semiprime and a nonnegative error cap")
    observed = len(readings)-readings.count(None)
    cap = observed - (readings.count(1) if certified_ones else 0)
    possible_errors = range(cap+1) if max_errors is None else (min(cap, max_errors),)
    trials = sorted((singleton_list_bound(readings, b, e, certified_ones=certified_ones), b, e)
                    for b in range(cap+1) for e in possible_errors)
    work = {} if stats is None else stats
    tested = 0
    modulus = 1 << (len(readings)+1)
    for iteration, (bound, b, e) in enumerate(trials, 1):
        candidates, nodes = noisy_singleton_list(n, readings, b, e, certified_ones=certified_ones)
        tested += nodes
        if candidates:
            divisor = factor_from_residue_list(n, modulus, candidates, stats=work)
            if divisor:
                work.update(bound=bound, orientation_cap=b, error_cap=e, trials=iteration,
                            prefix_nodes_tested=tested)
                return divisor
    return None
