"""Verify the Fourier support bounds for literal cubic and quartic terms."""

from math import gcd, prod, tau

import numpy as np


def polynomial(grid, modulus, degree):
    square = grid * grid % modulus
    if degree == 3:
        return grid * (3 - square) % modulus
    return square * (square - 2) % modulus


def derivative(grid, modulus, degree):
    square = grid * grid % modulus
    if degree == 3:
        return 3 * (1 - square) % modulus
    return 4 * grid * (square - 1) % modulus


def check(factors, scale, degree):
    n = prod(factors)
    modulus = n * n
    y = np.arange(modulus, dtype=np.int64)
    terms = np.exp((1j * tau / modulus) * (scale * polynomial(y, modulus, degree) % modulus))
    observed = np.fft.fft(terms) / modulus

    a = np.arange(n, dtype=np.int64)
    fibers = scale * derivative(a, n, degree) % n
    predicted = np.zeros_like(observed)
    largest_fiber = 0
    for lower_frequency in np.unique(fibers):
        rows = a[fibers == lower_frequency]
        largest_fiber = max(largest_fiber, len(rows))
        weights = np.zeros(n, dtype=np.complex128)
        base = scale * polynomial(rows, modulus, degree) - lower_frequency * rows
        weights[rows] = np.exp((1j * tau / modulus) * (base % modulus))
        predicted[lower_frequency::n] = np.fft.fft(weights) / n

    error = float(np.max(np.abs(observed - predicted)))
    assert error < 1e-11
    bound = (degree - 1)**2
    assert largest_fiber <= bound
    assert float(np.max(np.abs(observed))) <= largest_fiber / n + 1e-11
    energy = np.abs(observed)**2
    support = int(np.count_nonzero(energy > 1e-18))
    assert support * largest_fiber**2 >= modulus
    best_n_energy = float(np.partition(energy, -n)[-n:].sum())
    assert best_n_energy <= n * (largest_fiber / n)**2 + 1e-11
    return support / modulus, best_n_energy, error


if __name__ == '__main__':
    checks = []
    for factors in ((5, 7), (7, 11), (11, 13), (17, 19), (19, 23)):
        n = prod(factors)
        for scale in (1, n - 1, n + 1, 2 * n - 1):
            assert gcd(scale, n) == 1
            for degree in (3, 4):
                checks.append(check(factors, scale, degree))
    print(f'{len(checks)} full-grid transforms passed')
    print(f'Minimum observed support density {min(row[0] for row in checks):.6f}')
    print(f'Maximum energy retained by the strongest N frequencies {max(row[1] for row in checks):.6f}')
    print(f'Maximum transform discrepancy {max(row[2] for row in checks):.3g}')
