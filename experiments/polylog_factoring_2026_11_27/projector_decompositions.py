"""Uniform ordered decompositions of (Z/n)^d into s labeled summands when s divides d."""

from math import gcd
from random import SystemRandom

from sympy import diag

if __package__:
    from .idempotent_sampler import bernoulli, gl_orders, random_basis
else:
    from idempotent_sampler import bernoulli, gl_orders, random_basis


def profile_weights_from_orders(orders, parts):
    """Integer weights and suffix totals for ranks summing to d=len(orders)-1."""
    dimension = len(orders) - 1
    weights = [orders[-1] // order for order in orders]
    dp = [[1] + [0] * dimension]
    for _ in range(parts):
        previous = dp[-1]
        dp.append([sum(weights[r] * previous[t-r] for r in range(t+1))
                   for t in range(dimension+1)])
    return weights, dp


def rank_profile_weights(prime, exponent, dimension, parts):
    return profile_weights_from_orders(gl_orders(prime, exponent, dimension), parts)


def draw_rank_profile(weights, dp, rng, *, exclude=None):
    """Draw with probability proportional to the product of coordinate weights.

    Supplying one valid profile excludes it and renormalizes the probabilities.
    """
    parts, dimension = len(dp) - 1, len(weights) - 1
    matched = exclude is not None
    if matched:
        if len(exclude) != parts or sum(exclude) != dimension or min(exclude) < 0:
            raise ValueError("invalid excluded rank profile")
        excluded_suffix = [1] * (parts+1)
        for j in reversed(range(parts)):
            excluded_suffix[j] = weights[exclude[j]] * excluded_suffix[j+1]
        if excluded_suffix[0] == dp[-1][dimension]:
            raise ValueError("excluding the only possible rank profile")
    remaining = dimension
    result = []
    for j in range(parts):
        total = dp[parts-j][remaining]
        if matched:
            total -= excluded_suffix[j]
        ticket = rng.randrange(total)
        for r in range(remaining+1):
            tail = dp[parts-j-1][remaining-r]
            if matched and r == exclude[j]:
                tail -= excluded_suffix[j+1]
            mass = weights[r] * tail
            if ticket < mass:
                result.append(r)
                remaining -= r
                matched = matched and r == exclude[j]
                break
            ticket -= mass
    return tuple(result)


def _remove_prime(n, p):
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return n, exponent


def _factor_hit(hit):
    """Trial division is exact and costs O(sqrt(hit)) divisions."""
    divisor = 2
    while divisor * divisor <= hit:
        if hit % divisor == 0:
            yield divisor
            hit, _ = _remove_prime(hit, divisor)
        divisor = 3 if divisor == 2 else divisor + 2
    if hit > 1:
        yield hit


def _sparse_prime_hits(n, passes, rng):
    """Each p | n is selected independently with probability 1-(1-1/p²)^passes."""
    hits = set()
    for _ in range(passes):
        common = gcd(n, rng.randrange(n), rng.randrange(n))
        hits.update(_factor_hit(common))
    return sorted(hits)


def _slot_labels(profile):
    return [i for i, rank in enumerate(profile) for _ in range(rank)]


def conjugate_profiles(modulus, common_profile, overrides, rng):
    """Conjugate canonical block projectors by one shared uniform GL basis.

    The baseline is common_profile. Each overrides key is a pairwise coprime
    full CRT component of modulus, and its value is the local profile.
    """
    slots = _slot_labels(common_profile)
    parts, dimension = len(common_profile), len(slots)
    if modulus <= 0 or dimension < 1 or min(common_profile) < 0:
        raise ValueError("expected a positive modulus and valid dimensions")
    diagonals = [[int(slot == i) for slot in slots] for i in range(parts)]
    remaining = modulus
    for component, profile in overrides.items():
        if (component <= 1 or remaining % component or gcd(modulus//component, component) != 1
                or len(profile) != parts or min(profile) < 0 or sum(profile) != dimension):
            raise ValueError("invalid CRT profile override")
        remaining //= component
        cofactor = modulus // component
        indicator = cofactor * pow(cofactor, -1, component) % modulus
        for j, (old, new) in enumerate(zip(slots, _slot_labels(profile))):
            if old != new:
                diagonals[old][j] -= indicator
                diagonals[new][j] += indicator
    basis, inverse = random_basis(modulus, dimension, rng)
    return tuple((basis * diag(*entries) * inverse).applyfunc(lambda x: x % modulus)
                 for entries in diagonals)


def sample_decomposition(modulus, dimension, parts, rng=None):
    """Exactly sample all labeled orthogonal projectors summing to the identity.

    Local zero ranks are allowed. All positive moduli and parts dividing the
    positive dimension are supported. No factorization of the modulus is supplied.
    """
    if modulus < 1 or dimension < 1 or parts < 1 or dimension % parts:
        raise ValueError("expected a positive modulus and dimension divisible by parts")
    rng = SystemRandom() if rng is None else rng
    middle = dimension // parts
    common = (middle,) * parts
    overrides = {}
    if parts > 1 and modulus > 1:
        passes = parts * parts
        for p in _sparse_prime_hits(modulus, passes, rng):
            _, e = _remove_prime(modulus, p)
            weights, dp = rank_profile_weights(p, e, dimension, parts)
            total = dp[-1][-1]
            exceptional = total - weights[middle]**parts
            denominator = p**(2*passes)
            hit_mass = denominator - (p*p-1)**passes
            if bernoulli(exceptional * denominator, total * hit_mass, rng):
                overrides[p**e] = draw_rank_profile(weights, dp, rng, exclude=common)
    return conjugate_profiles(modulus, common, overrides, rng)


def factor_from_decomposition(modulus, projectors):
    """Return a verified proper factor if two prime components have unequal ranks."""
    if modulus < 2:
        return None
    for matrix in projectors:
        for coefficient in reversed(matrix.charpoly().all_coeffs()):
            divisor = gcd(int(coefficient), modulus)
            if divisor < modulus:
                # At the lowest nonvanishing degree the coefficient is ±1
                # at components of maximal rank and zero at all smaller ones.
                if divisor > 1:
                    return divisor
                break
    return None
