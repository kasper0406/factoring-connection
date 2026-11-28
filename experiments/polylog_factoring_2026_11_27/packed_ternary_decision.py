"""Recover a factor from nonadaptive scalar ternary-sphere residue decisions."""

from dataclasses import dataclass
from math import gcd, isqrt, prod

from sympy import jacobi_symbol, sieve

if __package__:
    from .short_jacobi_grams import first_vector_orbits
else:
    from short_jacobi_grams import first_vector_orbits


@dataclass(frozen=True)
class PrimeReading:
    prime: int
    queries: int
    yes: int
    squares: frozenset
    symmetric_pairs: frozenset


@dataclass(frozen=True)
class PackedRecovery:
    factor: int
    coordinate_product: int
    crt_modulus: int
    packing_scale: int
    squared_radius: int
    readings: tuple


@dataclass(frozen=True)
class SymmetricRecovery:
    factor: int
    coefficients: tuple
    crt_modulus: int
    packing_scale: int
    squared_radius: int
    readings: tuple


def _primes_until(excluded, threshold):
    primes = []
    modulus = 1
    for prime in sieve:
        if gcd(prime, excluded) == 1:
            primes.append(prime)
            modulus *= prime
            if modulus > threshold:
                return primes


def _residue_candidates(u, dot, norm, prime):
    pivot = next(i for i in range(3) if u[i] % prime)
    j, k = [i for i in range(3) if i != pivot]
    inverse = pow(u[pivot], -1, prime)
    z = [0, 0, 0]
    for y in range(prime):
        z[j] = y
        for w in range(prime):
            z[k] = w
            t = (dot - u[j] * y - u[k] * w) * inverse % prime
            if (t * t + y * y + w * w - norm) % prime == 0:
                z[pivot] = t
                yield tuple(z)


def _packing_parameters(n, x):
    a = x * x + 1
    if n <= 1 or n % 2 == 0 or x < 0 or jacobi_symbol(-a, n) != -1:
        raise ValueError("expected a nonnegative x and odd n with Jacobi(−x²−1,n)=−1")
    scale = 1
    while scale < 8 * (x + 1) * (n + 1):
        scale *= 7
    return scale, (x + scale * n) ** 2 + 1


def _read_challenge(n, x, decide, power):
    scale, m = _packing_parameters(n, x)
    menu = _primes_until(210 * m * n * (x * x + 1), n**power)
    first_vectors = tuple(first_vector_orbits(x))
    readings = []
    for prime in menu:
        modulus = scale * prime
        squares = set()
        symmetric_pairs = set()
        queries = yes = 0
        for u in first_vectors:
            for z in _residue_candidates(u, -n * x, n * n, prime):
                residue = tuple((i - scale * j) % modulus for i, j in zip(u, z))
                queries += 1
                if decide(m, modulus, residue):
                    yes += 1
                    a1, a2, a3 = (c * c % prime for c in z)
                    squares.update((a1, a2, a3))
                    symmetric_pairs.add(
                        ((a1 * a2 + a2 * a3 + a3 * a1) % prime, a1 * a2 * a3 % prime)
                    )
        squares.difference_update((0, n * n % prime))
        symmetric_pairs.discard((0, 0))
        readings.append(
            PrimeReading(prime, queries, yes, frozenset(squares), frozenset(symmetric_pairs))
        )
    return scale, m, tuple(readings)


def recover_packed_factor(n, x, decide):
    """Use a predicate false on empty classes and true on singleton classes.

    Multiple-witness answers may be arbitrary, including representation parity.
    The promise is n=pq for distinct odd primes and Jacobi(−x²−1,n)=−1.
    The first-column enumeration is polynomial in x, which may exceed polylog(n).
    All queries decide(M, modulus, W_residue) are determined by n and x.
    """
    scale, m, readings = _read_challenge(n, x, decide, 24)

    count = max(len(row.squares) for row in readings)
    if not 1 <= count <= 3:
        raise ValueError("inconsistent coordinate square counts")
    value, modulus = 0, 1
    for row in readings:
        if len(row.squares) == count:
            correction = (prod(row.squares) - value) * pow(modulus, -1, row.prime) % row.prime
            value += modulus * correction
            modulus *= row.prime
    if modulus <= n**6 or not 0 < value < n**6:
        raise ValueError("coordinate product has no certified lift in the stated range")
    divisor = gcd(value, n)
    if divisor in (1, n):
        raise ValueError("the recovered coordinate product did not give a proper factor")
    return PackedRecovery(divisor, value, modulus, scale, m, readings)


def recover_packed_symmetric_factor(n, x, decide):
    """Use exact existence and CRT on two elementary symmetric coefficients."""
    scale, m, readings = _read_challenge(n, x, decide, 6)
    e2 = e3 = 0
    modulus = 1
    for row in readings:
        if len(row.symmetric_pairs) > 1:
            raise ValueError("inconsistent symmetric coefficients in one prime reading")
        c2, c3 = next(iter(row.symmetric_pairs), (0, 0))
        inverse = pow(modulus, -1, row.prime)
        e2 += modulus * ((c2 - e2) * inverse % row.prime)
        e3 += modulus * ((c3 - e3) * inverse % row.prime)
        modulus *= row.prime
    divisor = gcd(n, e2, e3)
    if not (0 < e2 < n**4 and 0 <= e3 < n**6 and 1 < divisor < n):
        raise ValueError("symmetric coefficients did not yield a valid proper factor")
    return SymmetricRecovery(divisor, (e2, e3), modulus, scale, m, readings)


def _sphere_hits(m, modulus, residue):
    """Yield the positive witness multiplicity for each successful coordinate pair."""
    if m < 0 or modulus <= 0 or len(residue) != 3:
        raise ValueError("expected a nonnegative norm, a positive modulus and three residues")
    rx, ry, rz = (v % modulus for v in residue)
    bound = isqrt(m)
    for i in range(-bound + (rx + bound) % modulus, bound + 1, modulus):
        remains = m - i * i
        ymax = isqrt(remains)
        for j in range(-ymax + (ry + ymax) % modulus, ymax + 1, modulus):
            square = remains - j * j
            k = isqrt(square)
            if k * k == square:
                hits = int(k % modulus == rz) + int(k != 0 and -k % modulus == rz)
                if hits:
                    yield hits


def brute_sphere_decision(m, modulus, residue):
    """Test one scalar sphere residue by enumeration with an exact square check."""
    return any(_sphere_hits(m, modulus, residue))


def brute_sphere_count(m, modulus, residue):
    """Count all signed ordered witnesses to the scalar sphere residue query."""
    return sum(_sphere_hits(m, modulus, residue))


def brute_sphere_parity(m, modulus, residue):
    """Return the parity of the scalar sphere witness count."""
    return bool(brute_sphere_count(m, modulus, residue) % 2)
