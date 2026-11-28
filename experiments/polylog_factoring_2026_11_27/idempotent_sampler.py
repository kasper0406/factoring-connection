"""Exact uniform idempotent matrices of even dimension over Z/n, without factoring n.

At p^e, rank-j idempotents form a single conjugacy class. The stabilizer is
GL_j x GL_{d-j}. For d=2k, the ratio of the size of rank k+i to rank k is
at most p^(-2e*i²): the factors F_j=prod_{a<=j}(1-p^-a) are log-convex.
The combined exceptional mass is therefore at most 2/(p²+1). Two independent
1/p² marks cover it, and exact thinning samples each exceptional local rank.
Independent uniform local changes of basis follow from CRT on GL_d(Z/n).

For the local conjugacy and counts see David Dolžan, "Bounds for the number of
idempotents in finite rings", https://doi.org/10.1080/00927872.2021.1929277.
"""

from math import gcd, isqrt
from random import SystemRandom

from sympy import Matrix, diag, zeros


def bernoulli(numerator, denominator, rng):
    """An exact rational Bernoulli draw using integer arithmetic."""
    if not 0 <= numerator <= denominator or denominator <= 0:
        raise ValueError("expected a probability between zero and one")
    return rng.randrange(denominator) < numerator


def sparse_integers(limit, rng):
    """Yield independent marked integers 2 <= m <= limit, each with probability 1/m².

    The probability of no mark in [L,H] is (L-1)(H+1)/(LH). Conditional
    on at least one, the first mark m has weight 1/(m(m-1)).
    """
    low = 2
    while low <= limit and bernoulli(limit-low+1, low*limit, rng):
        left, right = low, limit
        while left < right:
            middle = (left+right)//2
            if bernoulli(right*(middle-left+1), middle*(right-left+1), rng):
                right = middle
            else:
                left = middle+1
        yield left
        low = left+1


def _prime_by_trial(m):
    if m % 2 == 0:
        return m == 2
    return all(m % k for k in range(3, isqrt(m)+1, 2))


def candidate_prime_powers(modulus, passes, rng):
    """Yield the full exponent of each distinct marked prime divisor of modulus.

    Trial primality is exact. Its expected total division count per pass is
    bounded by sum_{m>=2} sqrt(m)/m², even without the divisibility filter.
    """
    marked = set()
    for _ in range(passes):
        marked.update(sparse_integers(modulus, rng))
    for prime in sorted(marked):
        if modulus % prime or not _prime_by_trial(prime):
            continue
        exponent, remaining = 1, modulus//prime
        while remaining % prime == 0:
            exponent += 1
            remaining //= prime
        yield prime, exponent


def gl_orders(prime, exponent, dimension):
    """Orders of GL_j(Z/(prime^exponent)) for j from 0 through dimension."""
    orders = [1]
    for j in range(1, dimension+1):
        power = (exponent-1)*(2*j-1) + j-1
        orders.append(orders[-1] * prime**power * (prime**j-1))
    return orders


def random_basis(modulus, dimension, rng):
    """Return a uniform invertible matrix and its inverse modulo modulus."""
    if modulus < 1 or dimension < 1:
        raise ValueError("expected a positive modulus and dimension")
    if modulus == 1:
        return zeros(dimension), zeros(dimension)
    while True:
        basis = Matrix(dimension, dimension, [rng.randrange(modulus) for _ in range(dimension**2)])
        integer_matrix = basis.to_DM()
        if gcd(int(integer_matrix.det()), modulus) == 1:
            break
    # Fraction-free inversion over Z avoids nonunit pivots modulo a composite.
    numerator, denominator = integer_matrix.inv_den()
    inverse = (numerator.to_Matrix()*pow(int(denominator), -1, modulus)).applyfunc(
        lambda value: value % modulus
    )
    return basis, inverse


def random_conjugate(diagonal, modulus, rng):
    """Conjugate a diagonal matrix by a uniform element of GL modulo modulus."""
    basis, inverse = random_basis(modulus, len(diagonal), rng)
    return (basis*diag(*diagonal)*inverse).applyfunc(lambda value: value % modulus)


def _rank_intervals(prime, exponent, dimension):
    """Integer lengths for each noncentral rank and the common ticket denominator."""
    orders = gl_orders(prime, exponent, dimension)
    weights = [orders[-1]//(orders[j]*orders[dimension-j]) for j in range(dimension+1)]
    # Two independent 1/p² marks put this prime in the menu with probability
    # (2p²-1)/p⁴. Divide each desired exceptional probability by that value.
    denominator = sum(weights)*(2*prime**2-1)
    weights[dimension//2] = 0
    intervals = [weight*prime**4 for weight in weights]
    if sum(intervals) > denominator:
        raise ArithmeticError("the marking probability must cover the exceptional ranks")
    return intervals, denominator


def sample_idempotent(modulus, dimension, rng=None):
    """Uniform square matrix E with E²=E modulo modulus, in positive even dimension."""
    if modulus < 1 or dimension < 2 or dimension % 2:
        raise ValueError("expected a positive modulus and positive even dimension")
    if rng is None:
        rng = SystemRandom()
    middle = dimension//2
    diagonal = [1]*middle + [0]*middle
    for prime, exponent in candidate_prime_powers(modulus, 2, rng):
        intervals, denominator = _rank_intervals(prime, exponent, dimension)
        ticket = rng.randrange(denominator)
        rank = middle
        for j, length in enumerate(intervals):
            if ticket < length:
                rank = j
                break
            ticket -= length
        if rank != middle:
            local_modulus = prime**exponent
            cofactor = modulus//local_modulus
            projector = cofactor*pow(cofactor, -1, local_modulus)
            for j in range(min(rank, middle), max(rank, middle)):
                diagonal[j] += projector if rank > middle else -projector
    return random_conjugate(diagonal, modulus, rng)
