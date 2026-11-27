"""Recover a factor from coarse disagreement rates; simulate the oracle separately."""

from fractions import Fraction
from itertools import combinations
from math import gcd
from random import Random


def folded(x):
    return abs(x - round(x))


def recover_factor(n, rate):
    shifts = [1 << j for j in range(n.bit_length() + 1)]

    def distance(h):
        return min(Fraction(1, 2), max(Fraction(0), Fraction(rate(h % n)) / 2))

    distances = [distance(h) for h in shifts]
    t, anchor = next((h, r) for h, r in zip(shifts, distances)
                     if Fraction(1, 10) <= r <= Fraction(2, 5))
    phases = []
    for h, r in zip(shifts, distances):
        observation = distance(t + h)
        phase = min((r, -r), key=lambda s: abs(folded(anchor + s) - observation))
        phases.append(phase % 1)

    phase = phases.pop()
    for estimate in reversed(phases):
        phase = min((phase / 2, (phase + 1) / 2),
                    key=lambda s: folded(s - estimate))
    root = round(phase * n) % n
    factor = gcd(root - 1, n)
    if root * root % n != 1 or factor in (1, n):
        raise ValueError("Oracle responses did not decode to a nonpublic root")
    return factor


if __name__ == "__main__":
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    pairs = list(combinations(primes, 2)) + [(2**61 - 1, 2**127 - 1)]
    rng = Random(485)
    for p, q in pairs:
        n = p * q
        u = (1 - 2 * p * pow(p, -1, q)) % n
        for mode in ("rounded", "low", "high", "mixed"):
            def oracle(h):
                exact = 2 * folded(Fraction(h * u, n))
                if mode == "rounded":
                    return Fraction(round(50 * exact), 50)
                sign = {"low": -1, "high": 1}.get(mode)
                if sign is None:
                    sign = rng.choice((-1, 1))
                return exact + Fraction(sign, 100)

            assert recover_factor(n, oracle) in (p, q)
    print(f"Recovered factors in {4 * len(pairs)} oracle simulations.")
