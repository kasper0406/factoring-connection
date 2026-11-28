"""Encode a positive SubsetSum instance in two clipped N-squared hyperbola queries."""

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt, prod

from sympy import isprime


@dataclass(frozen=True)
class GapQuery:
    upper_digit: int
    polarity: int


@dataclass(frozen=True)
class Encoding:
    weights: tuple
    target: int
    primes: tuple
    exponents: tuple
    scale: int
    modulus: int
    coefficients: tuple
    base: int
    start: int
    width: int
    queries: tuple

    @property
    def offset(self):
        return -self.start % self.modulus

    @property
    def components(self):
        return tuple(p**e for p,e in zip(self.primes,self.exponents))


def target_scale(weights):
    total = sum(weights)
    return max(128*(3*total+1), 1024*3**(len(weights)+1))


def next_proved_prime(lower):
    candidate = max(3,(lower+1) | 1)
    while True:
        if candidate < 1 << 64:
            # SymPy uses deterministic primality tests throughout this range.
            prime = isprime(candidate)
        else:
            # The candidate is polynomial in the binary SubsetSum input length.
            prime = all(candidate % d for d in range(3,isqrt(candidate)+1,2))
        if prime:
            return candidate
        candidate += 2


def deterministic_components(weights, *, clear_predivision=True):
    """Generate prime powers; by default all bases exceed the reduction's trial bound."""
    k = len(weights)+1
    augmented_sum = 3*sum(weights)+1
    ell = 128*(k+(augmented_sum+1).bit_length()+10)**2
    scale = 1 << ell
    primes, exponents = [], []
    floor = (1 << 36)*k*(ell+100+k)
    prime = floor if clear_predivision else 2
    for _ in range(k):
        prime = next_proved_prime(prime)
        primes.append(prime)
        lower, upper = 0, scale.bit_length()+20+k.bit_length()
        while upper-lower > 1:
            middle = (upper+lower)//2
            if prime**middle > 1000*k*scale:
                upper = middle
            else:
                lower = middle
        exponents.append(upper)
    if clear_predivision:
        assert max(primes) < floor**2
        nbits = prod(p**e for p,e in zip(primes,exponents)).bit_length()
        assert min(primes) > (1 << 31)*(nbits+1)
    return scale,tuple(primes),tuple(exponents)


def nearest(x):
    return (2*x.numerator+x.denominator)//(2*x.denominator)


def nearest_unit(x, prime):
    center = nearest(x)
    if center % prime:
        return center
    return min((center-1,center+1), key=lambda value: abs(value-x))


def upper_query(base, components, quarter):
    n = prod(components)
    modulus = n*n
    coefficient = 0
    for i, q in enumerate(components):
        block = modulus//(q*q)
        residue = 2*base*pow(block,-1,q) % q
        desired = Fraction(quarter*q*q,4) if i == 0 else Fraction(0)
        h = residue + q*nearest((desired-residue)/q)
        assert abs(h-desired) <= Fraction(q,2)
        coefficient += block*h
    lift = coefficient*pow(2,-1,modulus) % modulus
    assert lift % n == base
    return GapQuery((lift-base)//n, quarter)


def encode(weights, target, primes, *, exponents=None, scale=None):
    """The supplied bases are distinct odd primes, with positive exponents."""
    weights, primes = tuple(weights), tuple(primes)
    if not weights or any(w <= 0 for w in weights) or not (0 < target < sum(weights)):
        raise ValueError("expected positive weights and an interior target")
    total = sum(weights)
    marker = 2*total+1
    augmented = (marker, *weights)
    target += marker
    scale = target_scale(weights) if scale is None else scale
    exponents = (1,)*len(primes) if exponents is None else tuple(exponents)
    if (len(primes) != len(augmented) or len(exponents) != len(primes)
            or any(p < 3 or not p % 2 for p in primes)
            or any(e < 1 for e in exponents) or scale < target_scale(weights)):
        raise ValueError("expected valid bases, exponents and scale")
    components = tuple(p**e for p,e in zip(primes,exponents))
    if (min(components) <= 1000*len(primes)*scale
            or any(gcd(p,q) != 1 for i,p in enumerate(primes) for q in primes[:i])):
        raise ValueError("expected sufficiently large, pairwise coprime prime powers")

    n = prod(components)
    coeffs = tuple((n//q)*nearest_unit(Fraction(q*w,scale),p)
                   for p,q,w in zip(primes,components,augmented))
    base = sum(coeffs) % n
    assert gcd(n,base) == 1
    limit = n//(4*scale)
    width = 1 << (limit.bit_length()-1)
    center = Fraction(n*(2*target-sum(augmented)),scale)
    start = nearest(center-Fraction(width,2))
    assert 0 < start < start+width < n//2
    assert all(not start <= public < start+width for public in (base,n-base))
    queries = tuple(upper_query(base,components,quarter) for quarter in (1,-1))
    return Encoding(augmented,target,primes,exponents,scale,n,coeffs,base,start,width,queries)
