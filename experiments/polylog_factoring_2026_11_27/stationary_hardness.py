"""Check a counting-subset-sum encoding in clipped stationary magnitudes."""

from itertools import product
from math import cos, fsum, isqrt, prod, tau
from random import Random


def odd_primes(count):
    primes = []
    candidate = 3
    while len(primes) < count:
        if all(candidate % p for p in primes if p <= isqrt(candidate)):
            primes.append(candidate)
        candidate += 2
    return primes


def encode(weights, target):
    """Construct the known factorization, unit scale, and cyclic window."""
    size = len(weights)
    signed_target = 2 * target - sum(weights)
    denominator = 16 * (sum(weights) + abs(signed_target) + 1)
    minimum_power = (1 << (size + 8)) * size * denominator
    components = []
    for p, weight in zip(odd_primes(size), weights):
        power = p
        while power < minimum_power:
            power *= p
        numerator = power * weight // denominator
        if numerator % p == 0:
            numerator += 1
        components.append((power, numerator))

    n = prod(power for power, _ in components)
    lifted_parts = []
    for power, numerator in components:
        complement = n // power
        lift_numerator = numerator * pow(complement, -1, power) % power
        lifted_parts.append(complement**2 * lift_numerator)

    divisor = 4 * denominator
    low = ((4 * signed_target - 1) * n + divisor - 1) // divisor
    high = ((4 * signed_target + 1) * n + divisor - 1) // divisor
    return n, lifted_parts, low % n, high - low


def check(weights, target):
    n, lifted_parts, start, width = encode(weights, target)
    modulus = n * n
    scale = sum(lifted_parts) % modulus
    roots = []
    wanted = 0
    for signs in product((-1, 1), repeat=len(weights)):
        root = sum(sign * value for sign, value in zip(signs, lifted_parts)) % modulus
        assert root * root % modulus == scale * scale % modulus
        in_window = (root - start) % n < width
        is_solution = sum(weight for sign, weight in zip(signs, weights) if sign == 1) == target
        assert in_window == is_solution
        wanted += is_solution
        roots.append(root)
    assert len(set(roots)) == 2**len(weights)

    phases = [cos(tau * ((2 * root) % modulus) / modulus) for root in roots]

    def stationary_sum(low, span):
        return fsum(phase for root, phase in zip(roots, phases) if (root - low) % n < span)

    clipped = stationary_sum(start, width)
    half = (n - 1) // 2
    half_difference = stationary_sum(start, half) - stationary_sum(start + width, half)
    assert abs(clipped - wanted) < 1 / 8
    assert abs(half_difference - clipped) < 1e-8
    # Add a worst-direction allowed error to each of the two evaluations.
    for error in (-1 / 4, 1 / 4):
        assert round(half_difference + error) == wanted
    return n.bit_length(), wanted, abs(clipped - wanted)


if __name__ == '__main__':
    rng = Random(261127)
    fixtures = [([1, 1, 2, 3, 5], 4), ([3, 6, 9, 12], 5), ([1] * 10, 5), ([9, 8], 30)]
    for size in range(2, 11):
        for _ in range(5):
            weights = [rng.randrange(1, 251) for _ in range(size)]
            fixtures.append((weights, rng.randrange(sum(weights) + 1)))
    results = [check(weights, target) for weights, target in fixtures]
    print(f'{len(results)} encodings verified on moduli up to {max(r[0] for r in results)} bits')
    print(f'Largest solution count {max(r[1] for r in results)}')
    print(f'Maximum real stationary discrepancy {max(r[2] for r in results):.3g}')
