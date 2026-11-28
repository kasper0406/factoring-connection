"""Locate the smaller balanced prime with an oracle for a complex Mellin magnitude.

For N = p*q with odd primes p < q <= 2p, use the normalized divisor moment

    D_N(s) = product_{r=p,q} (1 + (1+1/r) r^(2-s)).

The oracle estimates V(a) = N b² |D_N(2 + iπ/log a)|, b = bit_length(N),
to a specified multiplicative distortion. This module supplies the decoder only.

To bound V for N >= 2^16, write d = |r-a| and m_r for the magnitude of
the individual r factor. For ceil(sqrt(N/2)) <= a <= floor(sqrt(N)),
r/a lies in [1/sqrt(2), 2]. For d > 0, (br m_r)² = b² + (c d)² with

    c = π b/log(a) * sqrt(1+1/r) * r|log(r/a)|/d * sinc(π log(r/a)/(2 log a)).

Here sinc(x) = sin(x)/x. The ratio t log(t)/(t-1) increases across
this range and lies between 0.8 and 1.4. The sinc argument is at most
π/15 in absolute value, so sinc exceeds 0.99. Also 2/log(2) <= b/log(a)
<= 34/(15 log(2)) and r >= 181. These estimates give 7 <= c <= 15.
At d = 0 we have br m_r = b exactly. By Cauchy-Schwarz and the triangle
inequality, W_r/2 <= br m_r <= 3 W_r/2, where W_r = b+10d. Multiplying
gives W_p W_q/4 <= V <= 9 W_p W_q/4.
"""

from fractions import Fraction
from math import ceil, gcd, isqrt
from numbers import Rational


def _trial(n, start, stop):
    for candidate in range(start, stop + 1):
        divisor = gcd(n, candidate)
        if 1 < divisor < n:
            return divisor
    raise ValueError("no balanced prime found under the input and oracle promises")


def factor_from_mellin_magnitude(n, oracle, distortion=2):
    """Return p using rational answers V(a)/K <= oracle(n, a) <= K V(a).

    Requires distinct odd primes p < q <= 2p and n = p*q. K is the rational
    distortion, at least one and two by default. For n >= 2^16 the decoder
    takes O(K² log n) oracle calls and O(K² log n) final gcds. The oracle
    may choose either error direction on each call. All decisions use exact
    arithmetic, with polynomial cost for K polynomial in the bit length
    and oracle answers of polynomial bit length.
    """
    if n < 3 or n % 2 == 0:
        raise ValueError("expected an odd balanced semiprime")
    if not isinstance(distortion, Rational) or distortion < 1:
        raise ValueError("expected rational distortion at least one")
    k = Fraction(distortion)
    root = isqrt(n)
    if root * root == n:
        raise ValueError("expected distinct prime factors")
    a = isqrt(n // 2) + 1
    if n < 2**16:
        return _trial(n, a, root)

    b = n.bit_length()
    while a <= root:
        answer = oracle(n, a)
        if not isinstance(answer, Rational) or answer <= 0:
            raise ValueError("expected a positive rational oracle answer")
        # Maintain a <= p. Then ceil(sqrt(n)) <= q <= n/a, and
        # n/a-a = (sqrt(n)-a)(sqrt(n)/a+1) <= 2.5 (sqrt(n)-a).
        # Hence S <= W_q <= 2.5 S and W_p/(4k) <= z <= 6k W_p.
        scale = b + 10 * (root + 1 - a)
        z = Fraction(answer, scale)
        if z <= 12 * k * b:
            return _trial(n, a, min(root, a + ceil(48 * k**2 * b / 10)))
        # Since z > 12kb, this advances by at least z/(120k), which is
        # at least (p-a)/(48k²). Its upper bound is the integer gap p-a.
        a += ceil((z - 6 * k * b) / (60 * k))
    raise ValueError("oracle advanced beyond the balanced factor interval")
