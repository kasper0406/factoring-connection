"""Compute the fourth angular coefficient from an unfactored positive index."""

from math import isqrt

from mpmath import mp


def fourth_coefficient(n):
    if n < 1:
        raise ValueError("expected a positive index")
    width = n.bit_length()
    length = n + 1
    with mp.workprec(64 + 12 * width):
        total = mp.mpf(0)
        for k in range(length // 2 + 1):
            theta0, theta2, theta4 = 1, 0, 0
            for x in range(1, isqrt(n) + 1):
                square = x * x
                residue = (k * square) % length
                term = 2 * mp.expjpi(mp.mpf(2 * residue) / length)
                theta0 += term
                theta2 += square * term
                theta4 += square * square * term
            sample = (theta0 * theta4 - 3 * theta2**2) / 2
            # At index n = length - 1 the Fourier demodulator is z itself.
            weight = 1 if k == 0 or 2 * k == length else 2
            total += weight * (sample * mp.expjpi(mp.mpf(2 * k) / length)).real
        return int(mp.nint(total / length))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    print(fourth_coefficient(parser.parse_args().N))
