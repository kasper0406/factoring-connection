"""Uniform decompositions of (Z/n)^d whose labeled summands are all free."""

from math import comb
from random import SystemRandom

if __package__:
    from .idempotent_sampler import bernoulli, gl_orders
    from .projector_decompositions import (
        _factor_hit, _remove_prime, conjugate_profiles, profile_weights_from_orders,
    )
else:
    from idempotent_sampler import bernoulli, gl_orders
    from projector_decompositions import (
        _factor_hit, _remove_prime, conjugate_profiles, profile_weights_from_orders,
    )


def _global_gl_orders(n, dimension):
    """Called only on the rare correction branch. Factor n by trial division."""
    orders = [1] * (dimension + 1)
    for prime in _factor_hit(n):
        n, exponent = _remove_prime(n, prime)
        local = gl_orders(prime, exponent, dimension)
        orders = [a * b for a, b in zip(orders, local)]
    return orders


def _balanced_mass(weights, low, parts, remaining):
    """Total weight of completions using only low and low+1."""
    high_count = remaining - low * parts
    if not 0 <= high_count <= parts:
        return 0
    mass = comb(parts, high_count) * weights[low] ** (parts - high_count)
    return mass * weights[low + 1] ** high_count if high_count else mass


def _exceptional_profile(weights, dp, low, rng):
    """Sample a rank profile with some entry outside {low, low+1}."""
    parts, remaining = len(dp) - 1, len(weights) - 1
    profile = []
    balanced_prefix = True
    for left in range(parts - 1, -1, -1):
        masses = []
        for rank in range(remaining + 1):
            suffix = dp[left][remaining - rank]
            if balanced_prefix and rank in (low, low + 1):
                suffix -= _balanced_mass(weights, low, left, remaining - rank)
            masses.append(weights[rank] * suffix)
        ticket = rng.randrange(sum(masses))
        for rank, mass in enumerate(masses):
            if ticket < mass:
                break
            ticket -= mass
        profile.append(rank)
        balanced_prefix = balanced_prefix and rank in (low, low + 1)
        remaining -= rank
    return tuple(profile)


def sample_free_decomposition(modulus, dimension, parts, rng=None):
    """Draw uniformly given that the rank of each image is constant across primes.

    Expected running time is polynomial in dimension, parts and log(modulus).
    A draw can explicitly trial-factor the input on a sufficiently rare branch.
    """
    if modulus < 2 or dimension < 1 or parts < 1:
        raise ValueError("expected modulus >= 2 and positive dimension and number of parts")
    if rng is None:
        rng = SystemRandom()
    low, extra = divmod(dimension, parts)
    high_labels = set(rng.sample(range(parts), extra))
    profile = tuple(low + (index in high_labels) for index in range(parts))
    # Every reverse balancing transfer suppresses the orbit size by at least n².
    # The resulting probability of any exceptional profile is at most this gate.
    numerator = min(modulus**2, parts * (parts - 1))
    if numerator and dimension >= 2 and bernoulli(numerator, modulus**2, rng):
        orders = _global_gl_orders(modulus, dimension)
        weights, dp = profile_weights_from_orders(orders, parts)
        total = dp[parts][dimension]
        exceptional = total - _balanced_mass(weights, low, parts, dimension)
        if bernoulli(exceptional * modulus**2, total * numerator, rng):
            profile = _exceptional_profile(weights, dp, low, rng)
    return conjugate_profiles(modulus, profile, {}, rng)
