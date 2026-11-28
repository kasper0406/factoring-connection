"""Recover one signed coordinate from residue decisions on a packed ternary sphere."""

from dataclasses import dataclass
from math import gcd

from sympy import jacobi_symbol

if __package__:
    from .packed_ternary_decision import _packing_parameters, _primes_until, _residue_candidates
    from .short_jacobi_grams import first_vector_orbits
else:
    from packed_ternary_decision import _packing_parameters, _primes_until, _residue_candidates
    from short_jacobi_grams import first_vector_orbits


def select_linear_x(n, limit):
    """Scan 2 through limit, returning (mixed x, None), (None, factor), or neither."""
    if n <= 1 or n % 2 == 0 or limit < 0:
        raise ValueError("expected odd n > 1 and a nonnegative limit")
    for x in range(2, limit + 1):
        a = x*x + 1
        for divisor in (gcd(x, n), gcd(a, n)):
            if 1 < divisor < n:
                return None, divisor
        if gcd(x, n) == 1 and jacobi_symbol(-a, n) == -1:
            return x, None
    return None, None


@dataclass(frozen=True)
class LinearRecovery:
    factor: int
    first_vector: tuple
    index: int
    coordinate: int
    crt_modulus: int
    queries: int
    yes: int


def _centered_crt(residues):
    value, modulus = 0, 1
    for prime, residue in residues:
        value += modulus * ((residue - value) * pow(modulus, -1, prime) % prime)
        modulus *= prime
    return (value if 2*value <= modulus else value - modulus), modulus


def recover_packed_linear_factor(n, x, decide, *, singleton_only=False):
    """Use exact existence, or (singleton_only) arbitrary bits on multiple points.

    Promise: n=pq with distinct odd primes, x>=2 coprime to n, Jacobi(-x²-1,n)=-1.
    The complete menu of decide(M, modulus, residue) calls depends only on n and x.
    The public work and menu size are polynomial in x and log(n).
    """
    if x < 2 or gcd(x, n) != 1:
        raise ValueError("expected x >= 2 coprime to n")
    scale, m = _packing_parameters(n, x)
    threshold = 4*n*n if singleton_only else 2*n
    primes = _primes_until(210*m*n*(x*x+1), threshold)
    queries = yes = 0
    candidates = []
    for u in first_vector_orbits(x):
        if gcd(*u) != 1:
            continue
        index = next(j for j, entry in enumerate(u) if entry and u.count(entry) == 1)
        axis = (0, 0, -n) if u == (0, 1, x) else None
        residues = []
        for prime in primes:
            modulus = scale * prime
            axis_class = tuple(c % prime for c in axis) if axis else None
            values = set()
            for z in _residue_candidates(u, -n*x, n*n, prime):
                if singleton_only and z == axis_class:
                    continue
                residue = tuple((a - scale*b) % modulus for a, b in zip(u, z))
                queries += 1
                if decide(m, modulus, residue):
                    yes += 1
                    values.add(z[index])
            if axis and not singleton_only:
                known = axis[index] % prime
                if known not in values:
                    raise ValueError("public coordinate frame missing from exact decisions")
                values.remove(known)
                if not values:
                    values.add(known)
            if len(values) > 1:
                raise ValueError("conflicting signed coordinates for one first vector")
            if values:
                residues.append((prime, values.pop()))
        coordinate, modulus = _centered_crt(residues)
        divisor = gcd(coordinate, n)
        if modulus > 2*n and abs(coordinate) < n and 1 < divisor < n:
            candidates.append((divisor, u, index, coordinate, modulus))
    if len(candidates) != 1:
        raise ValueError("expected exactly one first-vector orbit yielding a proper factor")
    return LinearRecovery(*candidates[0], queries, yes)
