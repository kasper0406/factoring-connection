"""Check which information survives the last Euclidean edge of b/N."""

from fractions import Fraction
from math import gcd, isqrt

from manin59 import LEVEL, word


def terminal(b, n):
    """Return signed previous denominator and positive reduced final denominator."""
    a, d = b, n
    previous, current, sign = 1, 0, -1
    while d:
        q, rem = divmod(a, d)
        previous, current = current, q * current + previous
        if rem == 0:
            return sign * previous, current
        a, d, sign = d, rem, -sign


def check(last=511):
    endpoints = semiprimes = 0
    for n in range(3, last + 1, 2):
        terms = [terminal(b, n) for b in range(n)]
        assert sorted(Fraction(r, d) for r, d in terms) == [
            Fraction(k, n) for k in range(-(n // 2), n // 2 + 1)
        ]
        if n % LEVEL:
            # The repository charts a projective row by its first / second entry.
            expected = sorted(
                n * pow(k, -1, LEVEL) % LEVEL if k % LEVEL else LEVEL
                for k in range(-(n // 2), n // 2 + 1)
            )
            assert sorted(word(b, n)[-1] for b in range(n)) == expected
        modulus = 1 << (n.bit_length() + 2)
        value = sum(n * pow(d, -1, modulus) for _, d in terms) % modulus
        assert value == sum(gcd(b, n) for b in range(n)) % modulus
        p = next((h for h in range(2, isqrt(n) + 1) if n % h == 0), None)
        if p is not None:
            first = next(b for b in range(1, isqrt(n) + 1) if n // terms[b][1] > 1)
            assert first == p
            q = n // p
            if p != q and all(q % h for h in range(2, isqrt(q) + 1)):
                assert (4 * n + 1 - value) // 2 == p + q
                semiprimes += 1
        endpoints += n
    return {
        "odd_moduli": (last - 1) // 2,
        "endpoints": endpoints,
        "distinct_prime_semiprimes": semiprimes,
    }


if __name__ == "__main__":
    print(check())
