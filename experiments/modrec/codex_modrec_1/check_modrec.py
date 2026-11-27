"""Run the independent checks and reproduce the P2 cost table with Python 3."""

from itertools import accumulate
from math import isqrt
from random import Random
from statistics import median
from time import perf_counter

from modrec import (
    boundary_series, cartier, coefficient, floor_series, hyperbola, hyperbola_cartier,
)


def totients(limit):
    # Invert sum_{d|n} phi(d) = n by subtraction. Diagnostic only.
    phi = list(range(limit + 1))
    for d in range(1, limit // 2 + 1):
        value = phi[d]
        for multiple in range(2 * d, limit + 1, d):
            phi[multiple] -= value
    return phi


def boundary_degree(x, bit, phi):
    odd_count = (x + 1) // 2
    at_one = (odd_count * (odd_count - 1) // 2 + bit * odd_count) & 1
    return at_one + sum(
        phi[m] for m in range(3, x + 1, 2) if x // m % 4 in (1, 2)
    )


def check_small():
    max_t, max_x = 512, 128
    shells = [0] * (max_t + 1)
    direct_checks = dense_checks = 0
    for x in range(max_x + 1):
        # Add the row of points (x,y) explicitly; prefix in the product xy.
        if x:
            for y in range(1, max_t // x + 1):
                shells[x * y] += 1
        series = floor_series(x) if x <= 32 else None
        for t, count in enumerate(accumulate(shells)):
            expected = count & 1
            result = hyperbola(t, x)
            assert result.parity == expected, (t, x, count)
            cutoff = min(t, x)
            assert result.thresholds == cutoff - cutoff.bit_length()
            direct_checks += 1
            if series is not None and t <= 128:
                assert coefficient(t, series) == expected, (t, x, "series")
                dense_checks += 1

    signed = 0
    for t in range(-64, 1):
        for x in range(-2, 49):
            expected = sum(t // d for d in range(1, x + 1)) & 1
            assert hyperbola(t, x).parity == expected
            if x <= 16:
                assert hyperbola_cartier(t, x) == expected
            signed += 1
    print(f"P2 lattice enumerations {direct_checks}, rational path {dense_checks}, signed {signed}")


def check_denominators():
    limit = 96
    phi = totients(limit)
    full = boundaries = children = values = 0
    for x in range(limit + 1):
        _, q = floor_series(x)
        assert q.bit_length() - 1 == (1 + sum(phi[:x + 1]) if x else 0)
        full += 1
        for bit in (0, 1):
            p, q = boundary_series(x, bit)
            assert q.bit_length() - 1 == boundary_degree(x, bit, phi)
            boundaries += 1
            for digit in (0, 1):
                assert cartier(p, q, digit)[1] == q
                children += 1
            # Compare this correction itself with the actual odd floor rows.
            for n in (0, 1, 10, x * x + 5):
                expected = sum((2 * n + bit) // d for d in range(1, x + 1, 2)) & 1
                assert coefficient(n, (p, q)) == expected
                values += 1
    print(f"reduced denominators H {full}, boundary {boundaries}, section children {children}, boundary values {values}")


def linear_complexity(bits):
    """Berlekamp–Massey over F_2, with packed coefficient vectors."""
    history, order = 0, 0
    connection, previous = 1, 1
    last_update = -1
    for i, bit in enumerate(bits):
        history = (history << 1) | bit
        if (connection & history).bit_count() & 1:
            corrected = connection ^ (previous << (i - last_update))
            if 2 * order <= i:
                previous, last_update = connection, i
                order = i + 1 - order
            connection = corrected
    return order


def check_finite_blocks():
    phi = totients(96)
    print("X bit measured_boundary_order_on_2D_bits")
    for x in (8, 16, 32, 64, 96):
        for bit in (0, 1):
            degree = boundary_degree(x, bit, phi)
            start = x * x + 7
            values = (
                sum((2 * n + bit) // d for d in range(1, x + 1, 2)) & 1
                for n in range(start, start + 2 * degree)
            )
            order = linear_complexity(values)
            assert order == degree, (x, bit, order, degree)
            print(x, bit, order)


def check_large():
    rng = Random(20261127)
    for i in range(128):
        bits = rng.randrange(35, 41)
        t = (1 << (bits - 1)) | rng.getrandbits(bits - 1)
        x = rng.randrange(1, 49) if i < 32 else rng.randrange(isqrt(t) + 1)
        expected = sum(t // d for d in range(1, x + 1)) & 1
        assert hyperbola(t, x).parity == expected
        if i < 32:
            assert hyperbola_cartier(t, x) == expected
    print("35-40-bit floor-sum comparisons 128, including 32 on the rational path")


def check_circle_claim():
    prime_limit = 200
    primes = [p for p in range(3, prime_limit, 2)
              if all(p % d for d in range(2, isqrt(p) + 1))]
    shells = [0] * (prime_limit ** 2 + 1)
    for y in range(2, prime_limit + 1):
        for x in range(1, y):
            norm = x * x + y * y
            if norm < len(shells):
                shells[norm] += 1
    distinct = squares = 0
    for i, p in enumerate(primes):
        for q in primes[i:]:
            if p * q % 4 != 1:
                continue
            if p == q:
                assert shells[p * q] == int(p % 4 == 1)
                squares += 1
            else:
                assert shells[p * q] == 2 * int(p % 4 == 1)
                distinct += 1
    print(f"P1 strict octant shell checks {distinct} distinct-prime products, {squares} prime squares")


def benchmark():
    cases = [(1 << k) for k in (16, 24, 32, 36)] + [(1 << 39) + 12345]
    phi = totients(isqrt(max(cases)))
    print("T X H levels thresholds ones full_degree boundary_degree median_ms")
    for t in cases:
        x = isqrt(t)
        expected = sum(t // d for d in range(1, x + 1)) & 1
        times = []
        for _ in range(5):
            start = perf_counter()
            result = hyperbola(t, x)
            times.append(perf_counter() - start)
        assert result.parity == expected
        degree = boundary_degree(x, t & 1, phi)
        print(t, x, result.parity, result.levels, result.thresholds, result.ones,
              1 + sum(phi[:x + 1]), degree, f"{1000 * median(times):.3f}")


if __name__ == "__main__":
    check_small()
    check_denominators()
    check_finite_blocks()
    check_large()
    check_circle_claim()
    benchmark()
