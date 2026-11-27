"""Exact P2 parity by simultaneous halving of T and X.

Run "python3 modrec.py T X".  Only the Python standard library is needed.
The optional --cartier path builds the reduced rational generating function
for a fixed X. It is intended for small X even when T is very large.

Write T = 2n+b. Even d give H(n, X//2). For odd d, write n = qd+r.
Then (2n+b)//d = 2q + [r >= (d+1)//2 - b]. The d=1 term is b.
Across all levels, an original d contributes exactly once, at its odd part.
"""

import argparse
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class Count:
    parity: int
    levels: int
    thresholds: int  # d >= 3
    ones: int        # thresholds that fired, before XOR


def normalize(t: int, x: int) -> tuple[int, int, int]:
    """Reduce signed inputs to a finite positive cutoff and a parity offset."""
    x = max(0, x)
    offset = 0
    if t < 0:
        # (-t'-1)//d = -1 - t'//d for t' >= 0.
        t, offset = -t - 1, x & 1
    return t, min(t, x), offset


def hyperbola(t: int, x: int) -> Count:
    t, x, parity = normalize(t, x)
    levels = thresholds = ones = 0
    while x:
        t, bit = divmod(t, 2)
        odd_rows = range(3, x + 1, 2)
        ones += sum(t % d >= (d + 1) // 2 - bit for d in odd_rows)
        thresholds += (x - 1) // 2
        parity ^= bit
        levels += 1
        x //= 2
    return Count(parity ^ (ones & 1), levels, thresholds, ones)


# Polynomials over F_2, with coefficient i stored in bit i of a Python int.
# All gcds and divisions in this section are polynomial operations.
def poly_mul(a: int, b: int) -> int:
    result = 0
    if a.bit_count() < b.bit_count():
        a, b = b, a
    while b:
        low = b & -b
        result ^= a << (low.bit_length() - 1)
        b ^= low
    return result


def poly_divmod(a: int, b: int) -> tuple[int, int]:
    if not b:
        raise ZeroDivisionError("zero polynomial")
    quotient = 0
    while a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        quotient |= 1 << shift
        a ^= b << shift
    return quotient, a


def poly_div(a: int, b: int) -> int:
    quotient, remainder = poly_divmod(a, b)
    if remainder:
        raise ArithmeticError("inexact polynomial division")
    return quotient


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_divmod(a, b)[1]
    return a


def reduce_fraction(p: int, q: int) -> tuple[int, int]:
    common = poly_gcd(p, q)
    return poly_div(p, common), poly_div(q, common)


def sum_fractions(terms: Iterable[tuple[int, int]]) -> tuple[int, int]:
    p, q = 0, 1
    for a, b in terms:
        common = poly_gcd(q, b)
        left_scale, right_scale = poly_div(b, common), poly_div(q, common)
        p = poly_mul(p, left_scale) ^ poly_mul(a, right_scale)
        q = poly_mul(q, left_scale)
    return reduce_fraction(p, q)


def floor_series(x: int) -> tuple[int, int]:
    """Sum_{t>=0} H(t,x) z^t, as a reduced P/Q over F_2.

    For x >= 1 the minimal denominator is (1+z) lcm_{d<=x}(1+z^d).
    To see this, first discard the common factor 1/(1+z) in the sum.
    At a primitive root of odd order m <= x, take the largest a with
    2^a m <= x. The unique term with pole multiplicity 2^a is d = 2^a m.
    All its competitors would require an additional odd multiplier >= 3
    or another factor 2, both exceeding x. Thus the leading pole survives,
    including for m=1. Restoring 1/(1+z) adds one degree at z=1.

    The degree contribution at order m is 2^a phi(m). This equals the
    sum of phi(2^j m) for 0 <= j <= a. Summing gives 1 + sum_{d<=x} phi(d).
    """
    return sum_fractions(
        (1 << d, poly_mul(3, (1 << d) | 1)) for d in range(1, x + 1)
    )


def boundary_series(x: int, bit: int) -> tuple[int, int]:
    """Series in n of the odd-row correction to H(2n+bit,x).

    Proof of the obstruction. Each odd d contributes the period polynomial
    z^h + ... + z^(d-1) over 1+z^d, with h=(d+1)//2-bit. All poles are
    simple. At a primitive root w of odd order m>1, a row can have a pole
    only for d an odd multiple of m. For all such rows, 2h = 1-2bit mod m,
    so the period polynomial evaluates to the same nonzero (w^h+1)/(w+1).
    The derivative of 1+z^d at w is always w^(-1). The equal residues
    survive summation exactly when the number of odd multiples is odd,
    equivalently when floor(x/m) mod 4 is 1 or 2.

    At 1, row d contributes (d-1)//2+bit mod 2. Summing over d gives
    epsilon = K(K-1)/2 + K*bit mod 2 where K=(x+1)//2. Thus the minimal
    denominator degree D is epsilon plus phi(m) summed over the surviving
    odd orders m>1. In particular every such order in (x/2,x] survives.

    More precisely D ~ (2G/pi^2) x^2, with G = sum_{j>=0} (-1)^j/(2j+1)^2.
    The survival indicator is sum_{k<=x/m} chi_4(k) over the integers.
    The sum of phi(m) over odd m <= y is asymptotic to 2 y^2/pi^2 by
    Mobius inversion: sum the odd k <= y/d weighted by k for each odd
    divisor d with weight mu(d). The leading constant is
    (1/4) product_{p odd}(1-p^(-2)). Interchange the sums over m and k
    and divide by x^2. The k-th term is bounded in magnitude by O(1/k^2),
    giving the stated limit. The m=1 adjustment is bounded.

    Over a splitting field the sequence is a sum of nonzero multiples of
    w^n at distinct odd-order roots. Taking n -> 2n+b squares these roots
    and rescales their coefficients by units. Squaring permutes each
    surviving set of primitive roots. Every iterated binary section has
    exactly the same reduced denominator. The smaller H term has no pole
    of odd order m>x/2. Its child denominators divide its denominator by
    the Cartier identity below. Thus even the complete H series retains
    every odd order m>1 in (x/2,x] after any binary digit path.

    Translation of n preserves the boundary denominator too.
    Even a block of 2D coefficients forces linear recurrence order D:
    two proper fractions of orders D and L<D have cross-difference of
    degree <D+L. Agreement through 2D coefficients forces equality.

    For an explicit sum of fractions with binomial denominators 1+z^e,
    every prime x/2 < p <= x, p>2, must divide some e. Otherwise the
    expression is regular at every primitive p-th root. The product of the
    listed absolute exponents is therefore divisible by the product of the
    primes in this interval. By the prime number theorem the exponents
    alone require Omega(x) binary digits. This applies to representations
    of the whole correction sequence for a fixed cutoff, including every
    section. A single coefficient can in principle exploit its index.
    """
    return sum_fractions(
        ((1 << d) - (1 << ((d + 1) // 2 - bit)), (1 << d) | 1)
        for d in range(1, x + 1, 2)
    )


def cartier(p: int, q: int, bit: int) -> tuple[int, int]:
    """Reduced series whose coefficient n is [z^(2n+bit)] P/Q."""
    product = poly_mul(p, q) >> bit
    section = 0
    index = 0
    while product:
        section |= (product & 1) << index
        product >>= 2
        index += 1
    # Q(z)^2 = Q(z^2) in characteristic two.
    return reduce_fraction(section, q)


def coefficient(t: int, series: tuple[int, int]) -> int:
    """Extract a coefficient of a rational series with Q(0)=1, t >= 0."""
    if t < 0:
        raise ValueError("negative power series index")
    while t:
        t, bit = divmod(t, 2)
        series = cartier(*series, bit)
    return series[0] & 1


def hyperbola_cartier(t: int, x: int) -> int:
    t, x, offset = normalize(t, x)
    return coefficient(t, floor_series(x)) ^ offset


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("t", type=int)
    parser.add_argument("x", type=int)
    parser.add_argument("--cartier", action="store_true")
    args = parser.parse_args()
    if args.cartier:
        print(hyperbola_cartier(args.t, args.x))
    else:
        print(hyperbola(args.t, args.x).parity)
