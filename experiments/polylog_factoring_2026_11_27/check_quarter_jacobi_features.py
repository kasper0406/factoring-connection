"""Compute the best fixture classifier within a fixed finite feature family."""

from collections import defaultdict
from fractions import Fraction
from math import gcd

from sympy import jacobi_symbol

from check_quarter_fourier import quarter


def check(p, q):
    n = p*q
    hidden = 1+p*(-2*pow(p, -1, q) % q)
    totals = defaultdict(int)
    units = 0
    for c in range(1, n):
        if gcd(c, n) != 1:
            continue
        symbols = tuple(int(jacobi_symbol(value % n, n)) for value in
                        (c, c+1, c*c+1, c*c+c+1))
        features = (4*c//n, c % 3, symbols)
        totals[features] += quarter(n, c*hidden)
        units += 1
    advantage = Fraction(sum(map(abs, totals.values())), 2*units)
    print(n, len(totals), advantage, f"{float(advantage):.6f}")


if __name__ == "__main__":
    for p, q in [(31, 37), (101, 109), (211, 239)]:
        check(p, q)
