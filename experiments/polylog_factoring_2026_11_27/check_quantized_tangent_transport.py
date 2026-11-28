"""Exact scalar phase cells and optimistic fixed-product AND forecasts on toy grids."""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from math import comb, isqrt


def menu(scale, height, bins):
    return [(a, b, d, bins) for a in range(1, height + 1)
            for b in range(a, height + 1) for d in (scale // 2, scale, 2 * scale)]


def phase_bin(product, probe):
    a, b, denominator, bins = probe
    return (isqrt(4 * bins**2 * a * b * product) // denominator) % bins


def real_cuts(scale, probes):
    """Interior breaks in the variable N, with exact comparison to the interval ends."""
    breaks = set()
    for a, b, denominator, bins in probes:
        coefficient = 4 * bins**2 * a * b
        first = isqrt(coefficient * scale**2) // denominator + 1
        last = isqrt(coefficient * (2 * scale)**2 - 1) // denominator
        for i in range(first, last + 1):
            breaks.add(Fraction((i * denominator)**2, coefficient))
    return sorted(breaks)


def integer_starts(breaks):
    return sorted({-(-point.numerator // point.denominator) for point in breaks})


def check_boundaries():
    scale = 32
    probes = menu(scale, 3, 4)
    starts = set(integer_starts(real_cuts(scale, probes)))
    previous = tuple(phase_bin(scale**2, probe) for probe in probes)
    for product in range(scale**2 + 1, 4 * scale**2):
        current = tuple(phase_bin(product, probe) for probe in probes)
        assert (current != previous) == (product in starts)
        previous = current


def check_fixed_product(scale, low_product, start, length):
    """Allow an oracle to fit the Hamming-optimal prediction separately in each cell."""
    probes = menu(scale, 3, 4)
    cuts = real_cuts(scale, probes)
    starts = integer_starts(cuts)
    histograms = [Counter() for _ in range(len(starts) + 1)]
    mask = (1 << length) - 1
    for low in range(1, scale, 2):
        x = scale + low
        y = scale + low_product * pow(low, -1, scale) % scale
        cell = bisect_right(starts, x*y)
        word = ((x & y) >> start) & mask
        histograms[cell][word] += 1

    population = scale // 2
    total = Counter()
    hamming = 0
    predictions = set()
    close = 0
    radius = length // 8
    for histogram in histograms:
        total.update(histogram)
        count = sum(histogram.values())
        if not count:
            continue
        prediction = 0
        for bit in range(length):
            ones = sum(c for word, c in histogram.items() if word & (1 << bit))
            hamming += min(ones, count - ones)
            if 2 * ones > count:
                prediction |= 1 << bit
        predictions.add(prediction)
        close += sum(c for word, c in histogram.items()
                     if (word ^ prediction).bit_count() <= radius)

    ideal = [Fraction(3**(length - word.bit_count()), 4**length)
             for word in range(1 << length)]
    variation = sum(abs(Fraction(total[w], population) - prob)
                    for w, prob in enumerate(ideal)) / 2
    ideal_ball = sum(Fraction(comb(length, d) * 3**(length-d), 4**length)
                     for d in range(radius + 1))
    assert Fraction(close, population) <= len(predictions) * ideal_ball + variation
    assert len(predictions) <= len(cuts) + 1
    return len(cuts) + 1, len(predictions), hamming / (population * length), float(variation)


if __name__ == '__main__':
    check_boundaries()
    for exponent in (12, 15, 18):
        scale = 1 << exponent
        start, length = exponent // 2, exponent // 3
        for low_product in (1, 82791 % scale):
            cells, words, loss, variation = check_fixed_product(scale, low_product, start, length)
            print(f'k={exponent} n={low_product} j={start} b={length} cells={cells} '
                  f'forecast_words={words} Hamming={loss:.6f} block_TV={variation:.6f}')
