"""Check spectral sparsification on modular affine branch masks."""

from fractions import Fraction
from math import gcd
from random import Random

from check_quarter_fourier import exact_spectrum
from quarter_fourier_recovery import estimator, heavy_frequencies


def report(p, q, rng):
    n = p*q
    hidden = 1+p*(-2*pow(p,-1,q) % q)

    def left(c):
        return int(hidden*c % n < n//2)

    def right(c):
        return int((7*c+13) % n < 3*n//4)

    values = [left(c)*right(c) for c in range(n)]
    coeff = exact_spectrum(values)
    l1 = sum(map(abs, coeff))
    bound = prod_l1(left, right, n)
    assert l1 <= bound+1e-9
    threshold = Fraction(1, 6)
    search = heavy_frequencies(
        n, threshold, estimator(n, values.__getitem__, rng, 16_384)
    )
    assert hidden in search or n-hidden in search
    assert any(1 < gcd(frequency-1, n) < n for frequency in search)
    dropped_energy = sum(abs(value)**2 for value in coeff if abs(value) < threshold)
    assert dropped_energy <= float(threshold)*l1+1e-9
    print(n, len(search), f"{l1:.3f}", f"{bound:.3f}")


def prod_l1(left, right, n):
    return sum(map(abs,exact_spectrum([left(c) for c in range(n)])))*sum(
        map(abs,exact_spectrum([right(c) for c in range(n)]))
    )


if __name__ == "__main__":
    rng = Random(261131)
    for p,q in [(11,13),(17,23),(29,37)]:
        report(p,q,rng)
