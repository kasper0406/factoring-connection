"""Approximate every gcd prefix using a compressed distribution built from prime powers."""

from collections import defaultdict


def divisor_mass_jumps(factors, inverse_error=8, *, moment=1):
    """Input maps distinct primes to positive exponents; return (divisor, mass) jumps."""
    if not isinstance(inverse_error, int) or inverse_error < 1:
        raise ValueError("Expected a positive integral inverse error")
    if not isinstance(moment, int) or moment < 1:
        raise ValueError("Expected a positive integral moment")
    scale = 2 * len(factors) * inverse_error
    jumps = [(1, 1)]
    for prime, exponent in sorted(factors.items()):
        if prime < 2 or exponent < 1:
            raise ValueError("Expected prime powers with positive exponents")
        updates = defaultdict(int)
        power, jordan = 1, 1
        for _ in range(exponent + 1):
            for d, mass in jumps:
                updates[d * power] += mass * jordan
            power *= prime
            jordan = power**moment - (power // prime)**moment
        cumulative = retained = 0
        jumps = []
        for d, mass in sorted(updates.items()):
            cumulative += mass
            if scale * cumulative > (scale + 1) * retained:
                jumps.append((d, cumulative - retained))
                retained = cumulative
    return jumps


def gcd_prefix_from_jumps(jumps, endpoint):
    if endpoint < 0:
        raise ValueError("Expected nonnegative endpoint")
    return sum(mass * (endpoint // divisor) for divisor, mass in jumps)
