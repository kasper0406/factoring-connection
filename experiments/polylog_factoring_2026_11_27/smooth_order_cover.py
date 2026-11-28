"""Recover a balanced prime from a large part of r²−1 supported on supplied primes.

Theorem S in notes/specialfam_results.md §5 reads shared support from n²−1.
Here the prime bases are supplied and their exponents on either factor are hidden.
"""

from fractions import Fraction
from math import gcd, isqrt

if __package__:
    from .affine_prime_window import affine_divisors, window_radius
else:
    from affine_prime_window import affine_divisors, window_radius

from src.specialfam.smoothcopp import unit_roots


def _profile_cover(n, bases, epsilon):
    """Return sufficiently large candidate moduli and their known prime powers.

    With quantum b>=1, a weight w represents the smallest ell-power >= 2^(bw).
    If b would round to zero, enumerate exact exponents. In that finite-size
    regime log(n)=O(s/epsilon), so the same composition bound applies.
    """
    size = len(bases)
    b = (epsilon.numerator*(n.bit_length()-1))//(2*size*epsilon.denominator)
    budget = ((2*n).bit_length()-1)//max(b, 1)
    upper = isqrt(2*n)
    options = []
    for prime in bases:
        # Odd prime powers divide r−1 or r+1; the dyadic part divides twice one of them.
        limit = (upper+1)*(2 if prime == 2 else 1)
        powers = [(0, 0, 1)]
        exponent, power = 0, 1
        for weight in range(1, budget+1):
            if b:
                threshold = 1 << (b*weight)
                while power < threshold:
                    power *= prime
                    exponent += 1
            else:
                power *= prime
                exponent += 1
            if power > limit:
                break
            if power != powers[-1][2]:
                powers.append((weight, exponent, power))
        options.append(powers)

    profiles = []
    active = {}
    def visit(i, remaining, modulus):
        if i == size:
            if modulus**4 >= n:
                profiles.append((modulus, active.copy()))
            return
        prime = bases[i]
        for weight, exponent, power in options[i]:
            if weight > remaining or modulus*power > 2*n:
                break
            if exponent:
                active[prime] = exponent
            visit(i+1, remaining-weight, modulus*power)
            if exponent:
                del active[prime]

    visit(0, budget, 1)
    profiles.sort(key=lambda item: item[0], reverse=True)
    return b, budget, profiles


def _decode_residue(n, modulus, residue, lower, upper, stats):
    start = -((residue-lower)//modulus)
    stop = (upper-residue)//modulus
    if start > stop:
        return None
    inverse = pow(modulus, -1, n)
    radius = window_radius(n)
    while start <= stop:
        end = min(stop, start+2*radius)
        center = (start+end)//2
        half_width = end-center
        work = {}
        found = affine_divisors(n, residue*inverse+center, half_width, lower, stats=work)
        stats['tiles'] += 1
        stats['lattices'] += int(work['method'] == 'lattice')
        stats['max_lattice_dimension'] = max(stats['max_lattice_dimension'], work['dimension'])
        for divisor in found:
            if 1 < divisor < n:
                for candidate in (divisor, n//divisor):
                    if lower <= candidate <= upper and candidate % modulus == residue:
                        stats.update(winning_method=work['method'], winning_radius=half_width)
                        return candidate
        start = end+1
    return None


def recover_smooth_order(n, bases, epsilon=Fraction(1, 8), *, stats=None):
    """Return a verified factor or None. Requires n=p*q with distinct p<q<=2p.

    The supplied bases are promised to be distinct primes. No exponents, signs
    or information about a common factor of p²−1 and q²−1 are required.
    For fixed epsilon>0, the cover is polynomial in log(n) if len(bases)
    is O(log log(n)). It guarantees a split if either hidden prime r has
    S-part of r²−1 greater than n^(1/4+epsilon).
    """
    bases = tuple(bases)
    epsilon = Fraction(epsilon)
    if n < 6 or epsilon <= 0 or any(not isinstance(p, int) or p < 2 for p in bases):
        raise ValueError("expected a semiprime, positive slack and prime bases")
    if len(set(bases)) != len(bases):
        raise ValueError("the prime bases must be distinct")
    if stats is None:
        stats = {}
    stats.update(quantum=0, budget=0, profiles=0, tried=0, residues=0, tiles=0,
                 lattices=0, max_lattice_dimension=0, winning_modulus=None, winning_residue=None,
                 winning_method=None, winning_radius=None)
    if n % 2 == 0:
        return 2
    for prime in bases:
        divisor = gcd(n, prime)
        if 1 < divisor < n:
            return divisor
    if not bases:
        return None

    lower = isqrt((n-1)//2)+1
    upper = isqrt(2*n)
    b, budget, profiles = _profile_cover(n, bases, epsilon)
    stats.update(quantum=b, budget=budget, profiles=len(profiles))
    for modulus, factors in profiles:
        stats['tried'] += 1
        for residue in unit_roots(factors):
            stats['residues'] += 1
            candidate = _decode_residue(n, modulus, residue, lower, upper, stats)
            if candidate is not None:
                stats.update(winning_modulus=modulus, winning_residue=residue)
                return candidate
    return None
