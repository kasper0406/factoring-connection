"""Enumerate balanced factor candidates from a fixed menu of binomial divisibility bits."""

from math import isqrt

if __package__:
    from .adaptive_binomial_lift import exact_central_zero
else:
    from adaptive_binomial_lift import exact_central_zero


def fixed_menu(n, multiplier=3, schedule=(1,)):
    """Return (level, multiplier, half-index) queries, independent of the answers."""
    if n < 15 or n % 2 == 0:
        raise ValueError("expected a product of two distinct odd primes")
    if multiplier < 3 or multiplier % 2 == 0:
        raise ValueError("expected an odd second multiplier >= 3")
    if not schedule or any(flag not in (0, 1) for flag in schedule):
        raise ValueError("expected a nonempty periodic binary schedule")
    lower = isqrt((n - 1) // 2) + 1
    depth = lower.bit_length() - 1
    queries = []
    for level in range(1, depth):
        modulus = 1 << (level + 1)
        queries.append((level, 1, n // modulus))
        if modulus >= 4 * multiplier and schedule[level % len(schedule)]:
            queries.append((level, multiplier, multiplier * n // modulus))
    return queries


def factor_by_fixed_streams(n, zero=exact_central_zero, *, multiplier=3,
                            schedule=(1,), stats=None):
    """On odd semiprimes p<q<2p, recover p given zero(n,t) = [n divides C(2t,t)].

    Direct integer binomials are an exact small-input reference. A supplied
    evaluator is charged separately from all candidate-generation work.
    """
    queries = fixed_menu(n, multiplier, schedule)
    lower = isqrt((n - 1) // 2) + 1
    upper = isqrt(n)
    depth = lower.bit_length() - 1
    readings = [[] for _ in range(depth)]
    for level, h, t in queries:
        readings[level].append((h, bool(zero(n, t))))

    residues = [1]
    widths = [1]
    for level in range(1, depth):
        bit = 1 << level
        modulus = 2 * bit
        following = []
        for x in residues:
            y = n * pow(x, -1, modulus) % modulus
            # Adding bit to either odd coordinate toggles the highest
            # product bit. Toggling both therefore preserves the product.
            for child, partner in ((x, y), (x + bit, y ^ bit)):
                if all(bool((h * child) & (h * partner) & bit) == answer
                       for h, answer in readings[level]):
                    following.append(child)
        residues = following
        widths.append(len(residues))

    stride = 1 << depth
    trials = 0
    if stats is not None:
        stats.update(queries=len(queries), widths=widths, trials=0)
    for residue in residues:
        candidate = lower + (residue - lower) % stride
        if candidate <= upper:
            trials += 1
            if n % candidate == 0:
                if stats is not None:
                    stats["trials"] = trials
                return candidate
    if stats is not None:
        stats["trials"] = trials
    raise ValueError("no balanced factor compatible with the readings")
