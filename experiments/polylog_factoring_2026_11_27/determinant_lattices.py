"""Uniform integer sublattices of prescribed index, returned in column Hermite form."""

from math import gcd
from random import SystemRandom

from sympy import zeros


def _gcd_draw(modulus, count, rng):
    value = modulus
    for _ in range(count):
        value = gcd(value, rng.randrange(modulus))
        if value == 1:
            break
    return value


def _diagonal_trial(determinant, dimension, rng):
    remaining = determinant
    modulus = determinant*determinant
    trailing = []
    for count in range(1, dimension):
        entry = _gcd_draw(modulus, count, rng)
        if remaining % entry:
            return None
        trailing.append(entry)
        remaining //= entry
    return (remaining, *trailing)


def sample_determinant_lattice(determinant, dimension, rng=None):
    """Sample a uniform column HNF of positive dimension and prescribed determinant.

    Among accepted trailing diagonals the s-th entry contributes weight a^(-s).
    The first diagonal is forced. Its elimination from the HNF multiplicity gives
    exactly these weights. The proposal uses the public modulus squared so that
    even the highest allowable local valuation is strictly below the gcd cutoff.
    """
    if determinant < 1 or dimension < 1:
        raise ValueError("the determinant and dimension must both be positive")
    if rng is None:
        rng = SystemRandom()
    while True:
        diagonal = _diagonal_trial(determinant, dimension, rng)
        if diagonal is not None:
            break
    lattice = zeros(dimension)
    for i, entry in enumerate(diagonal):
        lattice[i, i] = entry
        if entry > 1:
            for j in range(i+1, dimension):
                lattice[i, j] = rng.randrange(entry)
    return lattice
