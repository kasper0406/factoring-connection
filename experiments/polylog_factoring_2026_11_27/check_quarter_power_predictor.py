"""Evaluate planted power exponents and verify the direct resonance gcd."""

from fractions import Fraction
from math import gcd, lcm

from check_quarter_fourier import quarter


def check(p, q):
    n = p*q
    hidden = 1 + p*(-2*pow(p, -1, q) % q)
    units = [c for c in range(1, n) if gcd(c, n) == 1]
    for label, exponent in (
        ("identity", 1),
        ("planted", 1 + lcm((p-1)//2, (q-1)//2)),
        ("arbitrary", 65_537),
    ):
        correlation = resonances = 0
        for c in units:
            v = pow(c, exponent-1, n)
            correlation += quarter(n, c*v) * quarter(n, c*hidden)
            if v not in (hidden, n-hidden):
                continue
            assert 1 < gcd(v-1, n) < n
            resonances += 1
        if label == "planted":
            assert 2*resonances == len(units)
        print(n, label, Fraction(correlation, len(units)), Fraction(resonances, len(units)))


if __name__ == "__main__":
    for p, q in [(31, 43), (103, 107), (211, 239)]:
        check(p, q)
