"""Recover factoring frequencies from noisy signed phase interval values."""

from fractions import Fraction
from math import isqrt, sqrt
from random import Random

import numpy as np

from check_salie_phase_bins import interval_values
from quarter_fourier_recovery import recover


def check(p, q, rng):
    n = p*q
    u = 1+p*(-2*pow(p,-1,q)%q)
    values, _, _, _ = interval_values(n,u)
    assert np.mean(values**2) <= 4
    integer_values = np.rint(values*sqrt(n)).astype(np.int64)
    assert np.max(np.abs(values*sqrt(n)-integer_values)) < 1e-8
    calls = [0]

    def predictor(scale):
        calls[0] += 1
        bound = isqrt(n)//8
        approximate = int(integer_values[scale])+rng.randint(-bound,bound)
        return max(-4,min(4,approximate/sqrt(n)))/4

    divisor, frequencies = recover(
        n,predictor,rng,Fraction(1,11),1,65_536,probe_samples=False
    )
    assert divisor in (p,q)
    assert any(f in (u,n-u) for f in frequencies)
    print(n, divisor, calls[0], len(frequencies))


if __name__ == "__main__":
    rng = Random(261132)
    for p,q in [(11,13),(13,17),(31,37),(101,109)]:
        check(p,q,rng)
