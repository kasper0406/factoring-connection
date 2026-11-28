"""Odd-dimensional objects with uniform local laws at a balanced semiprime.

Requires n=pq for distinct primes p<q<=2p. A centered affine window couples
the two rank parities. All other local randomness is independent, including
the distance of each rank from the middle pair and its choice of basis.
"""

from math import isqrt
from random import SystemRandom

if __package__:
    from .affine_prime_window import prime_window_divisor, window_radius
    from .idempotent_sampler import gl_orders
    from .projector_decompositions import conjugate_profiles
    from .sublattice_sampler import lattice_from_types
else:
    from affine_prime_window import prime_window_divisor, window_radius
    from idempotent_sampler import gl_orders
    from projector_decompositions import conjugate_profiles
    from sublattice_sampler import lattice_from_types


def _radius_weights(p, d, projector):
    """Unnormalized weights of T=0,...,(d-1)/2, one rank per parity."""
    orders = gl_orders(p, 1, d)
    weights = []
    for r in range(d//2, -1, -1):
        size = orders[d] // (orders[r]*orders[d-r])
        if not projector:
            size //= p**(r*(d-r))  # forget the choice of complement
        weights.append(size)
    return weights


def _gated_tickets(p, q, d, interval_size, projector):
    """Conditional radius weights when a uniform candidate was p or q."""
    left, right = _radius_weights(p, d, projector), _radius_weights(q, d, projector)
    total = sum(left)*sum(right)
    exceptional = total - left[0]*right[0]
    central = 2*total - interval_size*exceptional
    if central < 0:
        raise ArithmeticError("the factor-hit probability must cover the radius correction")
    yield 0, 0, central
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i or j:
                yield i, j, interval_size*a*b


def _sample_radii(n, d, projector, rng):
    if d == 1:
        return {}
    lower = isqrt((n-1)//2) + 1
    upper = isqrt(2*n)
    candidate = rng.randrange(lower, upper+1)
    if n % candidate:
        return {}
    p, q = sorted((candidate, n//candidate))
    tickets = list(_gated_tickets(p, q, d, upper-lower+1, projector))
    draw = rng.randrange(sum(weight for _, _, weight in tickets))
    for i, j, weight in tickets:
        if draw < weight:
            return {prime: t for prime, t in ((p,i),(q,j)) if t}
        draw -= weight
    raise ArithmeticError("invalid radius ticket total")


def _rank(d, parity, radius):
    low = d//2 - radius
    return low if low % 2 == parity else d-low


def _ranks(n, d, projector, rng):
    if n < 6 or d < 1 or d % 2 != 1:
        raise ValueError("expected a balanced squarefree semiprime and a positive odd dimension")
    mark = prime_window_divisor(n, rng.randrange(n), radius=window_radius(n))
    bit = rng.randrange(2)
    radii = _sample_radii(n, d, projector, rng)
    # mark=n flips both parities. A proper mark is one known prime, even
    # when its rank happens to have the central radius.
    default = _rank(d, bit ^ (mark == n), 0)
    if 1 < mark < n:
        radii.setdefault(mark, 0)
    overrides = {}
    for p, t in radii.items():
        rank = _rank(d, bit ^ (mark % p == 0), t)
        if rank != default:
            overrides[p] = rank
    return default, overrides


def coupled_projector(n, d, rng=None):
    """Return an idempotent matrix whose two prime reductions are each uniform."""
    rng = SystemRandom() if rng is None else rng
    default, ranks = _ranks(n, d, True, rng)
    overrides = {p: (rank,d-rank) for p,rank in ranks.items()}
    return conjugate_profiles(n, (default,d-default), overrides, rng)[0]


def coupled_sublattice(n, d, rng=None):
    """Return column HNF for an intermediate lattice with uniform local subspaces."""
    rng = SystemRandom() if rng is None else rng
    default, ranks = _ranks(n, d, False, rng)
    return lattice_from_types(n, d, default, {p:(rank,) for p,rank in ranks.items()}, rng)
