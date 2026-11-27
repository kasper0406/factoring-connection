"""Simulate relative magnitudes of aligned projected cubics by padded binary comparisons."""

from math import gcd, prod
from mpmath import mp
from sympy import nextprime


class ComparisonCalibration:
    """compare(base, scale, mask1, mask2) orders aligned magnitudes at a common unit scale."""

    fractional_bits = 14

    def __init__(self, n, compare):
        if n < 5 or gcd(n, 6) != 1:
            raise ValueError("Expected an odd base coprime to six")
        self.n = n
        self.compare = compare
        self.whole_bits = 2*n.bit_length()
        self.grid = 1 << self.fractional_bits
        gains = [0]*self.whole_bits + list(range(1, self.fractional_bits+1))
        self.components = []
        prime = 3
        for fraction in gains:
            prime = int(nextprime(prime))
            while n % prime == 0:
                prime = int(nextprime(prime))
            q = prime
            while q < 10**8*len(gains):
                q *= prime
            with mp.workprec(4*q.bit_length()+80):
                gain = mp.power(2, mp.mpf(1)/(1 << fraction))
                angle = mp.acos((gain**2-5)/4)/(2*mp.pi)
                frequency = int(mp.nint(q*q*angle))
            if frequency % prime == 0:
                frequency += 1
            negative = next(x for x in range(2, prime) if pow(x, (prime-1)//2, prime) != 1)
            self.components.append((q, negative, frequency))

        pad = prod(q for q, _, _ in self.components)
        self.base = n*pad
        self.original_projector = pad*pow(pad, -1, n)
        self.negative_mask = self.original_projector
        self.flips = []
        self.scale_offset = 0
        for q, negative, frequency in self.components:
            outside = self.base//q
            projector = outside*pow(outside, -1, q)
            self.negative_mask += negative*projector
            self.flips.append((1-negative)*projector % self.base)
            self.scale_offset += frequency*outside**2
        self.negative_mask %= self.base
        self.scale_offset %= self.base**2
        self.scale_step = pad**2
        self.whole_flips = [0]
        for flip in self.flips[:self.whole_bits]:
            self.whole_flips.append((self.whole_flips[-1]+flip) % self.base)

    def magnitude(self, scale, mask):
        padded_scale = (self.scale_offset + self.scale_step*scale) % (self.base**2)
        start = self.negative_mask + (mask-1)*self.original_projector
        low, high = 0, self.whole_bits*self.grid
        while high-low > 1:
            middle = (low+high)//2
            whole, fraction = divmod(middle, self.grid)
            selected = start + self.whole_flips[whole]
            for i in range(self.fractional_bits):
                if fraction & (self.grid >> (i+1)):
                    selected += self.flips[self.whole_bits+i]
            if self.compare(self.base, padded_scale, selected % self.base, self.negative_mask):
                high = middle
            else:
                low = middle
        return mp.power(2, -mp.mpf(low+high)/(2*self.grid))


def square_unit_via_comparisons(n, u, compare):
    """Decide squareness of a unit for gcd(n,30)=1 by at most two comparisons."""
    if n < 2 or gcd(n, 30) != 1 or gcd(n, u) != 1:
        raise ValueError("Expected a unit and a modulus coprime to 30")
    # Baseline has signs (+ on n, - on 5). Candidate has (u on n, + on 5).
    baseline = 1 + n*(pow(n, -1, 5) % 5)
    candidate = u + n*(((1-u)*pow(n, -1, 5)) % 5)
    base = 5*n
    for scale in (1, 2):
        padded_scale = (25*scale + 12*n*n) % (base*base)
        if not compare(base, padded_scale, candidate % base, baseline):
            return False
    return True


def count_negative_components(n, u, calibration, confidence_bits=8, draw=None):
    """Count distinct negative Legendre signs on odd n, with failure ≤2^-confidence_bits.

    calibration(cofactor).magnitude(t,u) estimates |J(t,u)| / |J(t,1)|.
    draw(bound) must supply independent uniform integers below bound.
    """
    from secrets import randbelow
    from sympy import primerange

    if n < 1 or n % 2 == 0 or gcd(n, u) != 1 or confidence_bits < 1:
        raise ValueError("Expected an odd positive modulus, a unit mask and positive confidence")
    if draw is None:
        draw = randbelow
    bits = n.bit_length()
    known = 0
    for p in primerange(3, 16*bits + 1):
        if n % p == 0:
            known += int(pow(u, (p-1)//2, p) != 1)
            while n % p == 0:
                n //= p
    if n == 1:
        return known
    estimate = calibration(n)
    samples = 32*bits*(confidence_bits+1)
    total = mp.mpf(0)
    with mp.workprec(64 + samples.bit_length() + bits.bit_length()):
        for _ in range(samples):
            t = draw(n*n)
            while gcd(t, n) != 1:
                t = draw(n*n)
            total -= mp.log(estimate.magnitude(t, u % n), 2)
        return known + int(mp.nint(total/samples))
