"""Decode binary carry parity templates for N=pq, p <= q <= 2p.

Integers encode polynomials over F_2 by coefficient bits. Incoming carry c_j
is encoded at bit j, so the corrected polynomial is N XOR the parity mask.
The caller supplies a finite menu; the default is the zero template.
"""

from itertools import combinations
from math import isqrt

from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_sqf_list


def binary_product(a, b):
    """Multiply coefficient bit strings over F_2."""
    product = 0
    while b:
        bit = b & -b
        product ^= a << (bit.bit_length() - 1)
        b ^= bit
    return product


def _divide(a, b):
    quotient = 0
    while a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def _gcd(a, b):
    while b:
        a, b = b, _divide(a, b)[1]
    return a


def _berlekamp(squarefree):
    """Use linear dependencies of T^(2j)-T^j modulo the squarefree input."""
    pivots = {}
    factors = [squarefree]
    for j in range(squarefree.bit_length() - 1):
        vector = _divide((1 << (2 * j)) ^ (1 << j), squarefree)[1]
        relation = 1 << j
        while vector and vector.bit_length() in pivots:
            pivot, witness = pivots[vector.bit_length()]
            vector ^= pivot
            relation ^= witness
        if vector:
            pivots[vector.bit_length()] = vector, relation
            continue
        # A fixed point of squaring reduces to 0 or 1 on each component.
        refined = []
        for factor in factors:
            common = _gcd(factor, relation)
            if 1 < common < factor:
                refined.extend((common, _divide(factor, common)[0]))
            else:
                refined.append(factor)
        factors = refined
    return factors


def polynomial_factors(bits):
    """Return (monic irreducible coefficient bits, multiplicity) pairs."""
    if bits < 1:
        raise ValueError("a nonzero polynomial is required")
    coefficients = [int(bit) for bit in bin(bits)[2:]]
    _, layers = gf_sqf_list(coefficients, 2, ZZ)
    factors = []
    for layer, multiplicity in layers:
        encoded = 0
        for bit in layer:
            encoded = (encoded << 1) | int(bit)
        factors.extend((irreducible, multiplicity) for irreducible in _berlekamp(encoded))
    return factors


def _degree_counts(factors, degree):
    """Suffix counts of divisors of each degree, including the unit."""
    counts = [[0] * (degree + 1) for _ in range(len(factors) + 1)]
    counts[-1][0] = 1
    for i in reversed(range(len(factors))):
        factor, multiplicity = factors[i]
        step = factor.bit_length() - 1
        following = counts[i + 1]
        for target in range(degree + 1):
            counts[i][target] = sum(following[target - j * step]
                                    for j in range(min(multiplicity, target // step) + 1))
    return counts


def _divisors(factors, counts, degree):
    """Enumerate just the specified degree, pruning infeasible suffixes."""
    stack = [(0, degree, 1)]
    while stack:
        i, missing, product = stack.pop()
        if i == len(factors):
            yield product
            continue
        factor, multiplicity = factors[i]
        step = factor.bit_length() - 1
        for j in range(min(multiplicity, missing // step) + 1):
            remaining = missing - j * step
            if counts[i + 1][remaining]:
                stack.append((i + 1, remaining, product))
            product = binary_product(product, factor)


def factor_from_carry_parity(n, templates=(0,), errors=0, *, cap=None, stats=None):
    """Return a verified smaller factor or None after testing the menu.

    Each template is a nonnegative integer of at most n.bit_length() bits.
    Cover its full Hamming ball of radius errors. The parity mask has bit 0
    zero, including when a template gets this public bit wrong. Each corrected
    polynomial with more than cap divisors across the eligible degrees is
    skipped. Runtime is polynomial for fixed errors, a polynomial-size menu,
    and polynomial cap. By default cap is n.bit_length()**6.
    """
    if n < 3 or n % 2 == 0:
        raise ValueError("n must be odd and at least 3")
    width = n.bit_length()
    if cap is None:
        cap = width**6
    if errors < 0 or cap < 0:
        raise ValueError("errors and cap must be nonnegative")
    lower, upper = isqrt((n - 1) // 2) + 1, isqrt(n)
    degrees = range(lower.bit_length() - 1, upper.bit_length())
    if stats is None:
        stats = {}
    stats.update(masks=0, skipped=0, divisors=0)
    seen = set()
    for template in templates:
        if template < 0 or template.bit_length() > width:
            raise ValueError("template outside the declared bit width")
        known_error = template & 1
        for weight in range(min(width - 1, errors - known_error) + 1):
            for positions in combinations(range(1, width), weight):
                mask = template ^ known_error ^ sum(1 << j for j in positions)
                if mask in seen:
                    continue
                seen.add(mask)
                stats["masks"] += 1
                factors = polynomial_factors(n ^ mask)
                counts = _degree_counts(factors, degrees.stop - 1)
                if sum(counts[0][d] for d in degrees) > cap:
                    stats["skipped"] += 1
                    continue
                for degree in degrees:
                    if counts[0][degree]:
                        for divisor in _divisors(factors, counts, degree):
                            stats["divisors"] += 1
                            if lower <= divisor <= upper and n % divisor == 0:
                                return divisor
    return None
