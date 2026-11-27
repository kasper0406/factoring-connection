"""Train and score a square-only interval predictor using public random roots."""

from math import gcd, log, sqrt
from random import Random


def examples(n, count, rng):
    for _ in range(count):
        h = rng.randrange(1, n)
        while gcd(h, n) != 1:
            h = rng.randrange(1, n)
        yield h*h % n, 1 if 4*min(h, n-h) < n else -1


def fit_and_score(n, bins=32, count=100000):
    votes = [0] * bins
    for square, label in examples(n, count, Random(3101)):
        votes[bins*square//n] += label
    labels = [1 if vote >= 0 else -1 for vote in votes]
    hits = sum(labels[bins*square//n] == label
               for square, label in examples(n, count, Random(8429)))
    return hits / count, sqrt(log(200) / (2*count))


if __name__ == "__main__":
    for modulus in (10403, 1022117):
        accuracy, radius = fit_and_score(modulus)
        print(f"{modulus}: {accuracy:.5f} (99% Hoeffding radius {radius:.5f})")
