"""Read a binary factor by enumerating bounded changes in carries or their slopes.

For the outgoing carry polynomial C of a hidden multiplication, its uncarried
coefficient polynomial is L_N+(2-X)C. Enumerate J=(1-X)^order C by positions
and nonzero coefficients, for order 1 or 2. Each valid candidate requires one
integer polynomial factorization. LLL gives a polynomial-time backend in the
theorem; the runnable version uses installed FLINT or SymPy factorization.

For m-bit N each outgoing carry is at most m, by the column recurrence.
For fixed h the enumeration therefore has m^O(h) profiles with at most h
jumps. At the true profile, F=L_a L_b for some a,b>1. The digit polynomials
are monic. Some irreducible factor of L_a evaluates at 2 to an integer of
absolute value greater than 1 that divides a, yielding a proper divisor of N.
This argument applies to composite operands and allows other profiles to
return a factor earlier. Carry reversal via digit polynomials is classical
(UMass CS TR 1992-066, section 5); the search here uses sparse differences.
"""

from itertools import combinations, product
from math import gcd

try:
    from flint import fmpz_poly
except ImportError:
    fmpz_poly = None
    from sympy import Poly, Symbol


def _jumps(length, count, order, bound):
    """Enumerate sparse polynomials divisible by (1-X)^order."""
    last = length-2+order
    nonzero = tuple(range(-bound, 0)) + tuple(range(1, bound+1))
    for positions in combinations(range(1, last+1), count):
        # c[0]=0 for odd factors. The first nonzero difference is positive.
        for first in range(1, bound+1):
            for tail in product(nonzero, repeat=count-order-1):
                prefix = (first,) + tail
                total = sum(prefix)
                if order == 1:
                    if not 0 < total <= bound:
                        continue
                    values = prefix + (-total,)
                else:
                    # J(1)=J'(1)=0 determine the final two coefficients.
                    u, v = positions[-2:]
                    moment = sum(i*a for i, a in zip(positions, prefix))
                    end, rem = divmod(u*total-moment, v-u)
                    previous = -total-end
                    if rem or not 0 < end <= bound or not 0 < abs(previous) <= bound:
                        continue
                    values = prefix + (previous, end)
                yield tuple(zip(positions, values))


def _uncarry(bits, jumps, order):
    """Integrate the differences and reject impossible raw column coefficients."""
    length = len(bits)
    events = dict(jumps)
    prev = slope = 0
    coefficients = []
    for i, digit in enumerate(bits):
        delta = events.get(i, 0)
        if order == 2:
            slope += delta
            delta = slope
        carry = prev+delta
        raw = digit+2*carry-prev
        # The binary digit polynomials are monic, their degree sum is at most
        # length-1, and at most one term uses each available input degree.
        if carry < 0 or not 0 <= raw <= min(i+1, length-i):
            return None
        coefficients.append(raw)
        prev = carry
    if prev != 0:
        return None
    while coefficients[-1] == 0:
        coefficients.pop()
    if coefficients[-1] != 1:
        return None
    return coefficients


def _read_factor(n, coefficients):
    # Constant coefficient 1 makes these integer polynomials primitive.
    if fmpz_poly is not None:
        _, factors = fmpz_poly(coefficients).factor()
    else:
        poly = Poly.from_list(coefficients[::-1], Symbol('X'), domain='ZZ')
        _, factors = poly.factor_list()
    for factor, _ in factors:
        value = factor(2) if fmpz_poly is not None else factor.eval(2)
        divisor = gcd(n, int(value))
        if 1 < divisor < n:
            return divisor
    return None


def carry_factor(n, max_jumps, *, order=1, jump_bound=None, stats=None):
    """Find a proper divisor if a binary product has the promised changes.

    For order 1 count changes of outgoing carries with zero padding on both
    ends. For order 2 count changes in their first differences. jump_bound, if
    supplied, is an additional bound on the absolute nonzero differences.
    Inputs are N and public search bounds only. Return None on exhaustion.
    """
    if n < 4 or max_jumps < 0 or order not in (1, 2) or (jump_bound is not None and jump_bound < 1):
        raise ValueError("expected N>=4, nonnegative jump count, order 1 or 2, and positive jump bound")
    if stats is None:
        stats = {}
    stats.update(profiles=0, factorizations=0, winning_jumps=None)
    if n % 2 == 0:
        return 2
    length = n.bit_length()
    bits = [(n >> i) & 1 for i in range(length)]
    for count in range(min(max_jumps, length-2+order)+1):
        if 0 < count <= order:
            continue
        if order == 1:
            # Each downward change from a to b satisfies a <= 2b+1.
            # Beyond the first rise there are at most count-1 descents.
            bound = min(length, (1 << min(max(0, count-1), length.bit_length()))-1)
        else:
            bound = 2*length
        if jump_bound is not None:
            bound = min(bound, jump_bound)
        profiles = _jumps(length, count, order, bound) if count else ((),)
        for jumps in profiles:
            stats['profiles'] += 1
            coefficients = _uncarry(bits, jumps, order)
            if coefficients is not None:
                stats['factorizations'] += 1
                divisor = _read_factor(n, coefficients)
                if divisor is not None:
                    stats['winning_jumps'] = jumps
                    return divisor
    return None
