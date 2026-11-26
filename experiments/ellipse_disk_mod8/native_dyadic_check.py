"""Independent checks for the tangent-block evaluator, diagonal disks and binary rank."""

from math import gcd, isqrt

from native_dyadic import block_count, block_size, one_form_split, prefix


def poly_gcd(a, b):
    while b:
        while a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def prime_powers(n):
    powers = []
    for p in range(2, isqrt(n) + 1):
        power = 1
        while n % p == 0:
            power *= p
            n //= p
        if power > 1:
            powers.append(power)
    if n > 1:
        powers.append(n)
    return powers


def disk(n, a, d, bound):
    total = 0
    for y in range(-isqrt(n * bound // d), isqrt(n * bound // d) + 1):
        radius = isqrt((n * bound - d * y * y) // a)
        for x in range(-radius, radius + 1):
            val = a * x * x + d * y * y
            total += val > 0 and val % n == 0
    return total


def main():
    ranks = disks = prefixes = blocks = 0
    splits = primes = probes = 0
    for n in range(3, 512, 2):
        roots = [x for x in range(n) if x * x % n == 1]
        bits = sum(1 << x for x in roots)
        period = (1 << n) | 1
        powers = prime_powers(n)
        expected_rank = n
        for power in powers:
            expected_rank = expected_rank // power * (power - 1)
        assert n + 1 - poly_gcd(bits, period).bit_length() == expected_rank
        if len(powers) == 2:
            mixed = bits ^ (1 << 1) ^ (1 << (n - 1))
            assert n + 1 - poly_gcd(mixed, period).bit_length() == n - 1
        ranks += 1
        midpoint = (n - 1) // 2
        a = n + 1
        d = a * n * ((midpoint**2 + n - 1) // n) - 1
        baseline = d + a * n * n

        def single_oracle(modulus, leading, trailing, cutoff, expected=(n, a, d)):
            nonlocal probes
            assert (modulus, leading, trailing) == expected
            assert 3 * cutoff**2 < 4 * leading * trailing
            probes += 1
            return disk(modulus, leading, trailing, cutoff) % 8

        divisor = one_form_split(n, single_oracle)
        if n >= 5 and any(n % p == 0 for p in range(2, isqrt(n) + 1)):
            assert divisor and 1 < divisor < n and n % divisor == 0
            splits += 1
        else:
            assert divisor is None
            primes += 1
        for cutoff in range(midpoint + 1):
            exact = sum(x <= cutoff for x in roots)
            if n <= 127:
                bound = (d + a * cutoff**2 + n - 1) // n
                assert gcd(a, d) == 1 and gcd(a, baseline) == 1
                assert 3 * bound**2 < 4 * a * d
                assert disk(n, a, d, bound) - disk(n, a, baseline, bound) == 4 * exact
                disks += 1
            if n <= 127 or cutoff in (midpoint, midpoint // 2):
                assert prefix(n, cutoff) == exact % 2
                prefixes += 1
    for n in (11009, 16789, 50707, 1022117):
        length = block_size(n)
        for start in range(1, n // 2 + 1, length):
            size = min(length, n // 2 + 1 - start)
            got = block_count(n, start, size)
            want = sum(x * x % n == 1 for x in range(start, start + size))
            assert got == want
            blocks += 1
    print(
        f"{ranks} rank checks, {disks} disk differences, {prefixes} prefixes, {blocks} larger blocks"
    )
    print(
        f"one fixed form: {splits} composites split, {primes} prime leaves, {probes} queries"
    )


if __name__ == "__main__":
    main()
