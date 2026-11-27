"""Evaluate the fourth angular coefficient by summing one coordinate analytically."""

from math import isqrt

from mpmath import mp


def row_resolvent(s):
    root = mp.sqrt(s)
    return mp.pi * mp.coth(mp.pi * root) / root


def pole_estimate(n):
    if n < 1:
        raise ValueError("expected a positive integer")
    with mp.workprec(80 + 12 * n.bit_length()):
        eta = 1 / mp.mpf(8*n)
        total = mp.mpc(0)
        for x in range(isqrt(n) + 1):
            u = x*x
            s = u - n + 1j*eta
            _, d1, d2, d3 = mp.taylor(row_resolvent, s, 3)
            row = -d1 - (2*s + 6*u)*d2 - (s*s + 6*u*s + u*u)*d3
            total += row if x == 0 else 2*row
        return eta**4 * total / 4


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    print(int(mp.nint(pole_estimate(parser.parse_args().N).real)))
