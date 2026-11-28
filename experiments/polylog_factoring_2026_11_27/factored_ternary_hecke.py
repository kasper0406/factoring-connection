"""Count individual residue classes on the three-square shell (x²+1)(pq)²."""

from bisect import bisect_left
from fractions import Fraction
from itertools import groupby
from math import gcd, isqrt
from operator import itemgetter

if __package__:
    from .fixed_field_walk import base_points, conjugate_vector
else:
    from fixed_field_walk import base_points, conjugate_vector


def _odd_prime(p):
    return p >= 3 and p % 2 == 1 and all(p % d for d in range(3, isqrt(p)+1, 2))


def prime_rotations(p):
    """One integral norm-p quaternion per right-unit orbit, for an odd prime p."""
    pairs = [[] for _ in range(p+1)]
    for a in range(-isqrt(p), isqrt(p)+1):
        bound = isqrt(p-a*a)
        for b in range(-bound, bound+1):
            pairs[a*a+b*b].append((a, b))
    # For p=1 mod 4 there is one odd coordinate; for p=3 mod 4 there is one
    # even coordinate. Exactly two opposite right associates put it first.
    first_parity = int(p % 4 == 1)
    for s, left in enumerate(pairs):
        for a, b in left:
            if a % 2 != first_parity:
                continue
            for c, d in pairs[p-s]:
                h = (a, b, c, d)
                if next(v for v in h if v) > 0:
                    yield h


def _scale(vector, scalar, modulus):
    return tuple(v*scalar % modulus for v in vector)


def _correction(a, prime, primitive):
    symbol = pow(-a, (prime-1)//2, prime)
    return (1 if symbol == 1 else -1) + int(primitive)


def primitive_residue_genus(x, p, q, modulus_factors, *, primitive=False):
    """Return (orbit size, mass average) for locally primitive liftable classes.

    The supplied factorization of the residue modulus is a prime->positive
    exponent mapping. Set primitive=True to count primitive global points.
    The input primes are trusted, and the formula uses no prime search.
    """
    a = x*x+1
    if x < 0 or min(p, q) <= 2 or p == q or gcd(p*q, 2*a) != 1:
        raise ValueError('require x>=0 and distinct odd primes p,q coprime to A')
    support, orbit_size = 1, 1
    for prime, exponent in modulus_factors.items():
        if prime < 2 or exponent < 1 or prime in (p, q):
            raise ValueError('require positive exponents and modulus coprime to pq')
        support *= prime
        if prime == 2:
            orbit_size *= 3*4**(exponent-1)
        else:
            local = prime**2-1 if a % prime == 0 else prime**2+prime*_correction(a, prime, False)
            orbit_size *= prime**(2*(exponent-1))*local
    count = sum(gcd(*v, support) == 1 for v in base_points(x, primitive=primitive))
    for prime in (p, q):
        count *= prime+1-_correction(a, prime, primitive)
    return orbit_size, Fraction(count, orbit_size)


class FactoredTernaryHecke:
    """Precompute the first prime's residue table for repeated exact counts.

    Inputs are x, distinct odd primes p,q coprime to x²+1, and a modulus
    coprime to pq. The modulus may be even or share primes with x²+1.
    """

    def __init__(self, x, p, q, modulus, *, primitive=False):
        if x < 0 or p == q or not _odd_prime(p) or not _odd_prime(q):
            raise ValueError('require x>=0 and distinct odd primes p,q')
        a = x*x+1
        if modulus <= 0 or gcd(p*q, a*modulus) != 1:
            raise ValueError('require positive modulus and pq coprime to A and the modulus')
        self.modulus = modulus
        self.norm_residue = a*(p*q)**2 % modulus
        points = base_points(x, primitive=primitive)
        entries = []
        for h in prime_rotations(p):
            for v in points:
                image = conjugate_vector(h, v)
                entries.append((tuple(t % modulus for t in image), 1))
        correction = _correction(a, p, primitive)
        if correction:
            entries.extend((_scale(v, p, modulus), -correction) for v in points)
        entries.sort()
        self.table = []
        for residue, group in groupby(entries, key=itemgetter(0)):
            count = sum(weight for _, weight in group)
            if count < 0:
                raise ArithmeticError('negative representation count after first Hecke step')
            if count:
                self.table.append((residue, count))
        self.inverse = pow(q, -1, modulus)
        self.inverse_rotations = [
            _scale((a0, -b, -c, -d), self.inverse, modulus)
            for a0, b, c, d in prime_rotations(q)
        ]
        self.correction = _correction(a, q, primitive)
        self.forward_images = len(points)*(p+1)
        self.base_count = len(points)
        if len(self.inverse_rotations) != q+1:
            raise ArithmeticError('incorrect number of norm-prime quaternion orbits')

    def _lookup(self, residue):
        i = bisect_left(self.table, (residue,))
        if i == len(self.table) or self.table[i][0] != residue:
            return 0
        return self.table[i][1]

    def count(self, residue):
        """Number of signed, ordered representations in the given coordinate class."""
        if len(residue) != 3:
            raise ValueError('a three-coordinate residue is required')
        m = self.modulus
        c = tuple(v % m for v in residue)
        if sum(t*t for t in c) % m != self.norm_residue:
            return 0
        total = 0
        for h in self.inverse_rotations:
            source = tuple(t % m for t in conjugate_vector(h, c))
            total += self._lookup(source)
        return total-self.correction*self._lookup(_scale(c, self.inverse, m))
