"""Check product-residue fibers and measure certified decisions on balanced prime pairs."""

from collections import Counter
from itertools import combinations
from math import isqrt
from random import Random

from sympy import nextprime, primerange

if __package__:
    from .cofactor_arch_bits import arch_outcomes, floor_lift
    from .cofactor_local_bits import local_outcomes, local_witness, screened_outcomes
else:
    from cofactor_arch_bits import arch_outcomes, floor_lift
    from cofactor_local_bits import local_outcomes, local_witness, screened_outcomes


def check_local():
    masks_checked = 0
    for j in range(1, 10):
        modulus = 2 << j
        masks = [1 << j] + [(1 << j) | (1 << k) for k in range(1, j)]
        units = range(1, modulus, 2)
        for n in units:
            pairs = [(u, n*pow(u, -1, modulus) % modulus) for u in units]
            for mask in masks:
                yes = sum(bool((u & mask) and (v & mask)) for u, v in pairs)
                actual = {s for s, count in ((True, yes), (False, len(units)-yes)) if count}
                assert actual == local_outcomes(n, mask)
                assert 2*yes <= len(units) if mask.bit_count() == 1 else 2*yes >= len(units)
                masks_checked += 1
                for answer in (True, False):
                    witness = local_witness(n, mask, answer)
                    if witness is None:
                        assert answer not in actual
                    else:
                        u, v = witness
                        assert (u*v-n) % modulus == 0
                        assert bool((u & mask) and (v & mask)) == answer
    return masks_checked


def check_large_fibers():
    rng = Random(2187)
    lifts = 0
    for j in (5, 20, 60, 128, 256):
        modulus = 2 << j
        for _ in range(10):
            n = (4096*modulus**6 + rng.getrandbits(5*j)) | 1
            mask = (1 << j) | (1 << rng.randrange(1, j+1))
            for answer in (True, False):
                a, b = local_witness(n, mask, answer)
                x, y = floor_lift(n, modulus, a, b)
                assert x <= y <= 2*x and y == n//x
                assert x % modulus == a and y % modulus == b
                assert bool((x & mask) and (y & mask)) == answer
                lifts += 1
    return lifts


def cohort(pairs):
    counts = {1: Counter(), 2: Counter()}
    for p, q in pairs:
        n = p*q
        low_floor = isqrt(n//2)
        top_weight = 1 << low_floor.bit_length()
        distance = top_weight-low_floor
        # Exclude bit positions that already exceed both possible factor sizes.
        positions = [1 << j for j in range(1, isqrt(2*n).bit_length())]
        for size in (1, 2):
            work = counts[size]
            for bits in combinations(positions, size):
                mask = sum(bits)
                answer = bool((p & mask) and (q & mask))
                local = local_outcomes(n, mask)
                real = arch_outcomes(n, mask)
                joint = screened_outcomes(n, mask)
                assert answer in joint and joint == local & real
                forced_positive = size == 2 and bits[1] == top_weight and bits[0] >= distance
                assert (real == {True}) == forced_positive
                if real == {False}:
                    assert (16*bits[0])**2 > n
                work['total'] += 1
                work['local'] += len(local) == 1
                work['real'] += len(real) == 1
                work['either'] += len(joint) == 1
                work['certified_yes'] += joint == {True}
                work['certified_no'] += joint == {False}
    return counts


def verify():
    print("Complete local residue and mask checks:", check_local())
    print("Constructive floor lifts through 256-bit indices:", check_large_fibers())
    fixed = [(p, q) for p, q in combinations(primerange(17, 150), 2) if q < 2*p]
    rng = Random(5731)
    held = []
    for bits in (24, 48, 96):
        for _ in range(4):
            p, q = sorted(int(nextprime((1 << (bits-1))+rng.randrange(1 << (bits-1))))
                          for _ in range(2))
            assert p < q < 2*p
            held.append((p, q))
    for name, pairs in (("small", fixed), ("held-out", held)):
        print(name, len(pairs), "prime pairs")
        for size, counts in cohort(pairs).items():
            print(f"  {size} selected positions: {dict(counts)}")


if __name__ == "__main__":
    verify()
