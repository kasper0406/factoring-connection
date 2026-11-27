"""Score a public square-bin sign estimate against the public cosine baseline."""

from math import pi
from random import Random

import numpy as np

from check_square_only_score import examples


def train(n, bins=32):
    counts = np.zeros(bins)
    votes = np.zeros(bins)
    for square, sign in examples(n, 100_000, Random(3101)):
        index = bins*square//n
        counts[index] += 1
        votes[index] += sign
    return votes / np.maximum(counts, 1)


def predictions(n, scales, bias):
    visible = 2*np.cos(2*pi*scales/n)
    baseline = 2/pi*(np.sqrt(np.maximum(0, 4-visible**2))
                     + visible*np.arcsin(visible/2))
    # Difference between sign-conditioned and unconditional torus means.
    response = np.sign(visible)*(abs(visible) + 4/pi - baseline)
    square = scales**2 % n
    corrected = baseline + response*bias[len(bias)*square//n]
    return baseline, corrected


def check(p, q):
    n = p*q
    bias = train(n)
    # Known factors enter only the scoring calculation.
    u = 1 + p*((-2*pow(p, -1, q)) % q)
    held_out = []
    rng = Random(8429)
    while len(held_out) < 100_000:
        c = rng.randrange(1, n)
        if np.gcd(c, n) == 1:
            held_out.append(c)
    all_scales = np.arange(1, n, dtype=np.int64)
    units = all_scales[np.gcd(all_scales, n) == 1]
    for label, scales in [('independent draws', np.array(held_out)), ('all units', units)]:
        baseline, corrected = predictions(n, scales, bias)
        truth = abs(2*np.cos(2*pi*scales/n) + 2*np.cos(2*pi*((u*scales) % n)/n))
        gain = np.mean((truth-baseline)**2 - (truth-corrected)**2)
        print(f'{n} {label}: MSE gain {gain:+.8f}')


if __name__ == '__main__':
    check(101, 103)
    check(1009, 1013)
