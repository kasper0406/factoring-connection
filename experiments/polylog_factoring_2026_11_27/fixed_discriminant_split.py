"""Split odd inputs at a common discriminant and cutoff."""

from functools import cache
from math import gcd

from fixed_window_split import mixed_root, reference_count, split_from_windows


def query_form(n, start):
    half = (n - 1) // 2
    radius = half // 2
    origin = start + radius + 1
    scale = 1 << (16 * n * n).bit_length()
    leading = scale * scale
    base = leading * (radius + 1)**2
    diagonal = base + (-leading - base) % n
    diagonal += n * (diagonal % 2)
    shear = scale * origin % n
    shear += n * (1 - shear % 2)
    cutoff = (base + (scale * radius + 2*n)**2 + 3*n - 1) // n
    return (leading, 2 * scale * shear, diagonal + shear**2), cutoff


def corrected_residue(n, start, disk_count):
    half = (n - 1) // 2
    form, cutoff = query_form(n, start)
    residue = disk_count(n, form, cutoff)
    if half % 2 == 0:
        endpoint = start + half + 1
        residue -= 2 * int((endpoint*endpoint - 1) % n == 0)
    for root in (1, n - 1):
        residue -= 2 * int(1 <= (root - start) % n <= half)
    return residue % 4


def find_factor(n, disk_count):
    half = (n - 1) // 2
    center = half // 2
    if half % 2 and mixed_root(n, center + 1):
        return gcd(center, n)

    @cache
    def residue(start):
        return corrected_residue(n, start, disk_count)

    baseline = residue(center)

    def parity(start):
        return ((residue(start) - baseline) % 4) // 2

    return split_from_windows(n, parity)


if __name__ == "__main__":
    windows = splits = prime_powers = 0
    for n in range(5, 321, 2):
        half = (n - 1) // 2
        center = half // 2
        signatures, baselines = set(), set()
        for start in range(center + 1):
            (a, b, c), cutoff = query_form(n, start)
            discriminant = 4*a*c - b*b
            signatures.add((cutoff, discriminant))
            assert gcd(gcd(a, b), c) == 1
            assert 3*cutoff**2 < discriminant
            expected = sum(mixed_root(n, x) for x in range(start + 1, start + half + 1))
            baselines.add((corrected_residue(n, start, reference_count) - 2*expected) % 4)
            windows += 1
        assert len(signatures) == len(baselines) == 1
        divisor = find_factor(n, reference_count)
        has_mixed = any(mixed_root(n, x) for x in range(n))
        if has_mixed:
            assert divisor and 1 < divisor < n and n % divisor == 0
            splits += 1
        else:
            assert divisor is None
            prime_powers += 1
    print(f"{windows} windows, {splits} splits, {prime_powers} primes or prime powers")
