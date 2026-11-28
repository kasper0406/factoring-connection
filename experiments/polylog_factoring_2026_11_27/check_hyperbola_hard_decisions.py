"""Check the public hard-decision windows, exact prime counts, and contour triangle."""

from bisect import bisect_left, bisect_right
from fractions import Fraction
from math import ceil, floor, isqrt, log, log1p


def windows(lam):
    """Closed intervals for K=λ(4−N/R²)/2. The corner interval is always positive."""
    yield Fraction(0), Fraction(1, 2)
    j = 1
    while j**3 <= lam:
        center = Fraction(2 * j + 1, 2)
        radius = Fraction(1, 4 * j)
        yield center - radius, center + radius
        j += 1


def public_decision(n, band, h, period):
    """O(1) integer comparisons for the common dyadic band model."""
    lam = h * band // period
    if n < 2 * band**2 or n > 4 * band**2:
        return False
    numerator = (4 * band**2 - n) * h
    denominator = 2 * band * period
    if 2 * numerator <= denominator:
        return True
    j = numerator // denominator
    return j >= 1 and j**3 <= lam and (
        2 * j * abs(2 * numerator - (2 * j + 1) * denominator) <= denominator
    )


def primes_in_band(band):
    high = 2 * band
    sieve = bytearray([1]) * (high + 1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(high) + 1):
        if sieve[p]:
            sieve[p * p :: p] = bytes((high - p * p) // p + 1)
    return [p for p in range(band, high + 1) if sieve[p]]


def count_prime_windows(primes, band, h, period):
    """Count all distinct ordered pairs in the predicted windows by binary search."""
    lam = h * band // period
    positive = [p for p in primes if h * p % period >= period // 2]
    source = 2 * band**2
    selected = successes = corner = 0
    for low_k, high_k in windows(lam):
        low = max(source, ceil(2 * source - source * high_k / lam))
        high = floor(2 * source - source * low_k / lam)
        index = bisect_left(primes, (low + 2 * band - 1) // (2 * band))
        for p in primes[index:]:
            a, b = (low + p - 1) // p, high // p
            diagonal = int(a <= p <= b)
            total = bisect_right(primes, b) - bisect_left(primes, a) - diagonal
            selected += total
            if low_k == 0:
                corner += total
            if h * p % period >= period // 2:
                successes += bisect_right(positive, b) - bisect_left(positive, a) - diagonal
    return selected, successes, corner


def contour_height(lam, length):
    """2λ times the uniform continuous contour excess in the upper corner."""
    def partner(u):
        return 2 * lam * (length - u) / (2 * lam - u)

    cuts = {0.0, length}
    for i in range(1, ceil(2 * length)):
        u = i / 2
        cuts.update((u, partner(u)))
    points = sorted(cuts)
    result = 0.0
    for left, right in zip(points, points[1:]):
        mid = (left + right) / 2
        yes = floor(2 * mid) % 2 == 0 and floor(2 * partner(mid)) % 2 == 0
        weight = 2 * lam * log1p((right - left) / (2 * lam - right))
        result += weight if yes else -weight
    return result


def verify():
    checks = 0
    for lam in (2**24, 2**32, 2**48):
        for j in (3, 9, 20, 64, 128):
            for shift in (-0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5):
                delta = shift / j
                length = j + 0.5 + delta
                height = contour_height(lam, length)
                error = 6 * abs(delta) + 20 * (length + 1)**3 / lam + 1e-10
                assert abs(height - (0.5 - 2 * j * abs(delta))) <= error
                checks += 1
    for band in (256, 512, 1024):
        primes = primes_in_band(band)
        for h in (1, 3):
            period = band // 16
            selected, successes, corner = count_prime_windows(primes, band, h, period)
            expected = expected_hits = expected_corner = 0
            for p in primes:
                for q in primes:
                    if p == q or not public_decision(p * q, band, h, period):
                        continue
                    expected += 1
                    expected_hits += min(h * p % period, h * q % period) >= period // 2
                    expected_corner += (4 * band**2 - p * q) * h <= band * period
            assert (selected, successes, corner) == (expected, expected_hits, expected_corner)
    print(f"Verified {checks} contour heights and six direct prime-pair censuses.")


if __name__ == "__main__":
    verify()
    print("R h lambda selected actual_yes corner_only_yes net_gain scaled_gain")
    for band in (2**17, 2**19, 2**21):
        primes = primes_in_band(band)
        denominator = len(primes) * (len(primes) - 1)
        for h in (1, 3):
            for scale in (32, 128):
                lam = h * scale
                selected, successes, corner = count_prime_windows(primes, band, h, band // scale)
                gain = (2 * successes - selected) / denominator
                print(band, h, lam, selected, successes, corner,
                      f"{gain:.8g}", f"{gain * lam**2 / log(lam):.6g}")
