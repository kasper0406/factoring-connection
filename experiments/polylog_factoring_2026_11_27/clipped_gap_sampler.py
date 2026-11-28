"""Factored-input importance sampler for the corrected modular interval contrast."""

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, prod
from random import SystemRandom

from sympy import isprime


@dataclass(frozen=True)
class GapQuery:
    j: int
    offset: int
    width: int
    sign_multiplier: int
    sign_offset: int


def _draw_local(p, e, modulus, c, rng):
    ticket = rng.randrange(p + 2*e*(p-1))
    if ticket < p:
        return rng.randrange(modulus)
    # Each positive depth gets two equal ticket intervals, one for each sign.
    depth_index, sign_index = divmod((ticket-p)//(p-1), 2)
    period = p**(depth_index+1)
    root = c if sign_index == 0 else -c
    return root % period + period*rng.randrange(modulus//period)


class ClippedGapSampler:
    def __init__(self, factors, c0):
        """Construct the proposal from a complete factorization of an odd N > 1."""
        entries = list(factors.items() if hasattr(factors, 'items') else factors)
        if not entries or any(not isinstance(e, int) or e <= 0 or p <= 2 or not isprime(p)
                              for p, e in entries):
            raise ValueError('expected odd prime bases with positive integral exponents')
        if len({p for p, _ in entries}) != len(entries):
            raise ValueError('expected distinct prime bases')
        self.n = prod(p**e for p, e in entries)
        self.c0 = c0 % self.n
        if gcd(self.c0, self.n) != 1:
            raise ValueError('the public lower digit must be a unit')
        self.normalizer = prod(Fraction(p + 2*e*(p-1), p) for p, e in entries)
        self._local = []
        for p, e in entries:
            modulus = p**e
            other = self.n // modulus
            crt = other*pow(other, -1, modulus)
            self._local.append((p, e, modulus, crt))

    def draw(self, rng=None):
        """Draw a lower digit with mass gcd(a²-c0²,N)/(N*T), including nonunits."""
        rng = SystemRandom() if rng is None else rng
        return sum(crt*_draw_local(p, e, modulus, self.c0, rng)
                   for p, e, modulus, crt in self._local) % self.n

    def _signed_change(self, a, query, nonstationary_only):
        n, c = self.n, self.c0
        if (a+query.offset) % n >= query.width or a in (c, n-c) or gcd(a, n) != 1:
            return 0
        g = gcd(a*a-c*c, n)
        if g == 1 or (nonstationary_only and g == n):
            return 0
        inverse = pow(a, -1, n*n)
        upper = ((a+c*c*inverse) % (n*n)) // n
        velocity = (2*c*inverse) % n
        half = (n-1)//2
        # The g-fold multiplicity cancels against the sampling density.
        at_j = (half - (upper+velocity*query.j) % g) // g
        at_zero = (half - upper % g) // g
        sign = 1 if (query.sign_multiplier*a+query.sign_offset) % n <= half else -1
        return sign*(at_j-at_zero)

    def sample_value(self, a, query, *, nonstationary_only=False):
        """Contribution of a proposed digit to D(j)/N; its magnitude is at most T."""
        if not 0 <= query.width <= self.n:
            raise ValueError('window width must be between zero and the modulus')
        return self.normalizer*self._signed_change(a % self.n, query, nonstationary_only)

    def estimate(self, query, samples, *, rng=None, nonstationary_only=False):
        """Return the exact rational sample mean of independent importance samples."""
        if not isinstance(samples, int) or samples <= 0:
            raise ValueError('expected a positive sample count')
        if not 0 <= query.width <= self.n:
            raise ValueError('window width must be between zero and the modulus')
        rng = SystemRandom() if rng is None else rng
        total = sum(self._signed_change(self.draw(rng), query, nonstationary_only)
                    for _ in range(samples))
        return self.normalizer*total/samples
