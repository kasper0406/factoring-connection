"""Measure an a priori square-root phase against upper binomial bits on all prime pairs."""

from math import cos, isqrt, pi, sqrt, tau

import numpy as np


CONSTANT = 4 * (2 * sqrt(2) - 1) / (3 * pi**2 * sqrt(2))


def public_score(n, h, period, precision=48):
    """Evaluate the score from n alone with fixed-point phase reduction."""
    denominator = 1 << precision
    scale = 2 * h * denominator
    phase = isqrt(n * scale**2) // period % denominator
    return -cos(tau * (phase / denominator) + pi / 4)


def interval_primes(left):
    high = 2 * left
    sieve = np.ones(high + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, isqrt(high) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    return np.flatnonzero(sieve[left:]) + left


def cohort(left, h, period):
    primes = interval_primes(left)
    lam = h * left / period
    bits = h * primes % period >= period // 2
    correlation = gain = positive = count = 0
    for start in range(0, len(primes), 96):
        x = primes[start : start + 96]
        products = x[:, None] * primes
        score = -np.cos(tau * (2 * h / period) * np.sqrt(products) + pi / 4)
        centered = (bits[start : start + len(x), None] & bits).astype(float) - 0.25
        distinct = x[:, None] != primes
        correlation += np.sum((score * centered)[distinct])
        adjustment = 0.25 * sqrt(period / h) * score / np.sqrt(np.sqrt(products))
        gain += np.sum((2 * centered * adjustment - adjustment**2)[distinct])
        positive += np.sum((centered[distinct] + 0.25))
        count += np.count_nonzero(distinct)
    cross = correlation / count
    n = int(primes[1]) * int(primes[-1])
    direct = -cos(tau * (2 * h / period) * sqrt(n) + pi / 4)
    assert abs(public_score(n, h, period) - direct) < 1e-10
    prevalence = positive / count
    best_constant_gain = gain / count - (prevalence - 0.25)**2
    return count, lam, prevalence, cross * sqrt(lam), best_constant_gain


if __name__ == "__main__":
    print(f"Limiting scaled correlation {CONSTANT:.6f}")
    print("R h R/D lambda ordered_pairs positive scaled_correlation gain_vs_best_constant")
    for left in (1 << 14, 1 << 15):
        for h in (1, 3):
            for periods in (16, 32, 64):
                count, lam, positive, cross, gain = cohort(left, h, left // periods)
                assert gain > 0 and cross > 0, (left, h, periods)
                print(f"{left} {h} {periods} {lam:.0f} {count} "
                      f"{positive:.4f} {cross:.6f} {gain:.7f}")
