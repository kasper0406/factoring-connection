"""Uniform intermediate lattices n Z^d <= L <= Z^d in positive even dimension."""

from math import gcd, isqrt, prod
from random import SystemRandom

from sympy import diag, eye
from sympy.matrices.normalforms import hermite_normal_form

if __package__:
    from .idempotent_sampler import bernoulli, random_basis
    from .projector_decompositions import _factor_hit, _remove_prime
else:
    from idempotent_sampler import bernoulli, random_basis
    from projector_decompositions import _factor_hit, _remove_prime


def type_counts(prime, exponent, dimension):
    """Birkhoff transition weights and totals over all monotone rank profiles.

    dp[t][s] sums over t further levels, starting at a deeper level of rank s.
    The number of subgroups of (Z/p^e)^d is dp[e][0].
    """
    gaussian = [[1]]
    for n in range(1, dimension+1):
        prev = gaussian[-1]
        gaussian.append([1] + [prev[j-1] + prime**j * prev[j] for j in range(1, n)] + [1])
    transitions = [[0] * (dimension+1) for _ in range(dimension+1)]
    for s in range(dimension+1):
        for r in range(s, dimension+1):
            transitions[s][r] = prime**(s*(dimension-r)) * gaussian[dimension-s][r-s]
    dp = [[1] * (dimension+1)]
    for _ in range(exponent):
        prev = dp[-1]
        dp.append([sum(transitions[s][r] * prev[r] for r in range(s, dimension+1))
                   for s in range(dimension+1)])
    return transitions, dp


def _central_mass(transitions, exponent, k):
    return transitions[0][k] * transitions[k][k]**(exponent-1)


def draw_type(transitions, dp, rng, *, exclude_central=False):
    """Sample r1 >= ... >= re, weighted by the number of subgroups of that type."""
    dimension, exponent = len(transitions)-1, len(dp)-1
    k = dimension//2
    if exclude_central and dimension % 2:
        raise ValueError("the central profile requires even dimension")
    central_powers = [1]
    if exclude_central:
        for _ in range(exponent):
            central_powers.append(central_powers[-1] * transitions[k][k])
    matched = exclude_central
    backwards = []
    s = 0
    for t in range(exponent, 0, -1):
        total = dp[t][s]
        if matched:
            total -= transitions[s][k] * central_powers[t-1]
        ticket = rng.randrange(total)
        for r in range(s, dimension+1):
            tail = dp[t-1][r]
            if matched and r == k:
                tail -= central_powers[t-1]
            mass = transitions[s][r] * tail
            if ticket < mass:
                backwards.append(r)
                s = r
                matched = matched and r == k
                break
            ticket -= mass
    return tuple(reversed(backwards))


def lattice_from_types(modulus, dimension, default_rank, exceptions, rng):
    """Sample conditional on local types; exceptions maps p to (r1,...,re).

    The other prime powers use a free subgroup of the specified default rank.
    The output columns are the canonical HNF basis of the inverse-image lattice.
    """
    if modulus < 1 or dimension < 1 or not 0 <= default_rank <= dimension:
        raise ValueError("invalid modulus, dimension or default rank")
    diagonal = [1]*default_rank + [modulus]*(dimension-default_rank)
    for p, profile in exceptions.items():
        exponent = len(profile)
        prime_power = p**exponent
        if (p < 2 or not profile or not 0 <= profile[-1] <= profile[0] <= dimension
                or any(r < s for r, s in zip(profile, profile[1:]))
                or modulus % prime_power or modulus//prime_power % p == 0):
            raise ValueError("invalid local profile or incomplete prime power")
        for i in range(dimension):
            # A cyclic summand p^lambda occupies exactly lambda profile levels.
            length = sum(i < r for r in profile)
            if i >= default_rank:
                diagonal[i] //= prime_power
            diagonal[i] *= p**(exponent-length)
    if modulus == 1:
        return eye(dimension)
    basis, _ = random_basis(modulus, dimension, rng)
    generators = (basis * diag(*diagonal)).row_join(modulus * eye(dimension))
    return hermite_normal_form(generators, D=prod(diagonal))


def sample_sublattice(modulus, dimension, rng=None):
    """Sample exactly uniformly from all subgroups of (Z/n)^d; return their lattice HNF."""
    if modulus < 1 or dimension < 2 or dimension % 2:
        raise ValueError("expected a positive modulus and positive even dimension")
    rng = SystemRandom() if rng is None else rng
    k = dimension//2
    remainder = modulus
    exceptions = {}
    for p in (2, 3, 5, 7):
        if remainder % p == 0:
            remainder, exponent = _remove_prime(remainder, p)
            transitions, dp = type_counts(p, exponent, dimension)
            profile = draw_type(transitions, dp, rng)
            if profile != (k,)*exponent:
                exceptions[p] = profile

    if remainder > 1:
        passes = 10
        hits = set()
        for _ in range(passes):
            hits.update(_factor_hit(gcd(remainder, rng.randrange(remainder))))
        for p in sorted(hits):
            _, exponent = _remove_prime(remainder, p)
            transitions, dp = type_counts(p, exponent, dimension)
            total = dp[exponent][0]
            exceptional = total - _central_mass(transitions, exponent, k)
            denominator = p**passes
            hit_mass = denominator - (p-1)**passes
            if bernoulli(exceptional * denominator, total * hit_mass, rng):
                exceptions[p] = draw_type(transitions, dp, rng, exclude_central=True)
    return lattice_from_types(modulus, dimension, k, exceptions, rng)


def factor_from_index(modulus, index):
    """Extract a proper divisor from a lattice index supported on primes of n."""
    if modulus < 2 or index < 1:
        raise ValueError("expected a nontrivial modulus and a positive index")
    root = isqrt(modulus)
    if root*root == modulus:
        return root
    while index % modulus == 0:
        index //= modulus
    divisor = gcd(index, modulus)
    return divisor if 1 < divisor < modulus else None


def factor_from_lattice(modulus, lattice):
    return factor_from_index(modulus, abs(int(lattice.det())))
