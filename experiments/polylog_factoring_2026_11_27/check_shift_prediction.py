"""Check the public translated-window baseline and sampled Fourier recovery."""

from fractions import Fraction
from random import Random

from quarter_fourier_recovery import recover


def sign(n, shift, value):
    half = (n-1)//2
    left = 1 <= (value-shift) % n <= half
    right = 1 <= (-value-shift) % n <= half
    return -1 if left ^ right else 1


def mean_numerator(n, shift):
    half = (n-1)//2
    summed = (-2*shift) % n
    overlap = max(0, half-abs(summed-(half+1)))
    return n-4*(half-overlap)


def check_means():
    for n in range(5, 160, 2):
        total = 0
        for a in range(n):
            numerator = mean_numerator(n, a)
            assert sum(sign(n,a,x) for x in range(n)) == numerator
            total += abs(numerator)
        assert 2*total == n*n+1


def check_recovery(p, q, rng):
    n = p*q
    u = 1+p*(-2*pow(p,-1,q) % q)
    shift = n//8
    calls = [0]

    def predictor(scale):
        calls[0] += 1
        exact = sign(n, shift, scale*u)
        return -exact if rng.random() < .12 else exact

    d, frequencies = recover(n, predictor, rng, Fraction(1, 4), 4, 16_384, probe_samples=False)
    assert d in (p, q)
    assert any((h*u) % n in frequencies for h in (-4,-3,-2,-1,1,2,3,4))
    print(n, d, calls[0], len(frequencies))


if __name__ == "__main__":
    check_means()
    print("Exact shift-only baseline checked on 78 odd moduli.")
    rng = Random(261129)
    for p, q in [(11,13), (101,109), (503,937)]:
        check_recovery(p, q, rng)
