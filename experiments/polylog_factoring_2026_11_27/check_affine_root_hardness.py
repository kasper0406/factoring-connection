"""Check the parsimonious affine-root and half-window encodings."""

from itertools import product
from math import gcd, prod
from random import Random


def primes():
    candidate = 3
    while True:
        if all(candidate % divisor for divisor in range(3, int(candidate**0.5) + 1, 2)):
            yield candidate
        candidate += 2


def instance(weights):
    """Return a known-factor modulus and a unit scale for the signed sums."""
    count = len(weights)
    mass = sum(weights)
    denominator = 4 * mass + 4
    moduli = []
    stream = primes()
    for _ in weights:
        prime = next(stream)
        modulus = prime
        while modulus <= 8 * count * denominator:
            modulus *= prime
        moduli.append(modulus)
    modulus = prod(moduli)
    terms = []
    for weight, local_modulus in zip(weights, moduli, strict=True):
        local = (2 * local_modulus * weight + denominator) // (2 * denominator)
        prime = next(
            divisor
            for divisor in range(3, local_modulus + 1, 2)
            if local_modulus % divisor == 0
        )
        if local % prime == 0:
            local += 1
        terms.append((modulus // local_modulus) * local)
    scale = sum(terms) % modulus
    assert gcd(scale, modulus) == 1
    return modulus, moduli, terms, scale, denominator


def roots(terms, scale, modulus):
    answer = []
    inverse = pow(scale, -1, modulus)
    for signs in product((-1, 1), repeat=len(terms)):
        coordinate = (
            sum(sign * term for sign, term in zip(signs, terms, strict=True))
            % modulus
        )
        root = coordinate * inverse % modulus
        assert root * root % modulus == 1
        answer.append((signs, coordinate, root))
    assert len({root for _, _, root in answer}) == len(answer)
    return answer


def cyclic_interval(value, end, length, modulus):
    return (value - (end - length + 1)) % modulus < length


def cumulative_count(weights, target):
    if target < 0:
        return 0, 0, 0, 0
    if target >= sum(weights):
        count = 1 << len(weights)
        return count, count, 0, 0
    modulus, _, terms, scale, denominator = instance(weights)
    encoded = roots(terms, scale, modulus)
    signed_target = 2 * target - sum(weights)
    # This is floor(modulus * (signed_target + 1) / denominator), including
    # Python's floor convention when the numerator is negative.
    end = modulus * (signed_target + 1) // denominator
    length = (modulus - 1) // 2
    count = sum(
        cyclic_interval(coordinate, end, length, modulus)
        for _, coordinate, _ in encoded
    )
    truth = sum(
        sum(weight for bit, weight in zip(bits, weights, strict=True) if bit) <= target
        for bits in product((0, 1), repeat=len(weights))
    )
    return count, truth, modulus.bit_length(), len(encoded)


def check():
    rng = Random(271828)
    cases = max_bits = roots_checked = 0
    for count in range(1, 11):
        for _ in range(20):
            weights = [rng.randrange(1, 18) for _ in range(count)]
            for target in range(-1, sum(weights) + 2):
                got, expected, bits, root_count = cumulative_count(weights, target)
                assert got == expected
                cases += 1
                max_bits = max(max_bits, bits)
                roots_checked += root_count
    return {
        "cumulative_instances": cases,
        "encoded_roots_checked": roots_checked,
        "maximum_modulus_bits": max_bits,
    }


if __name__ == "__main__":
    print(check())
