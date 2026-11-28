"""Run public subset-product targets with and without their dictionary words.

Given a b-bit benchmark modulus, sample b+64 units before sampling any
targets. Draw each target as a uniform binary word in those generators.
Both treatments get the same Euclidean rows and the same cofactor groups.
The control gives each target an independent formal variable. The labelled
treatment uses its public word. A dependency with even sign and prime
columns and even word exponents supplies two square roots by public
multiplication. Check at most 40 independent dependencies per modulus.

By default use four 64-bit semiprimes with 4,096 targets each. To change
the scale pass e.g. --bits 80 --pool 8192. Search counts exclude nonunits
encountered and rejected in preparing the dictionary; report those
separately. Timings cover setup, both treatments, and gcd verification.
"""

import argparse
from collections import Counter
from math import gcd, isqrt, prod
from random import Random
from time import perf_counter

if __package__:
    from . import check_cross_target_dependencies as control
else:
    import check_cross_target_dependencies as control


def set_bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def sample_targets(n, bits, seed, pool):
    rng = Random(f"relation-dictionary:{bits}:{seed}")
    generators = []
    rejected = 0
    while len(generators) < bits + 64:
        unit = rng.randrange(1, n)
        if gcd(unit, n) == 1:
            generators.append(unit)
        else:
            rejected += 1
    words = [rng.getrandbits(len(generators)) for _ in range(pool)]
    targets = []
    for word in words:
        t = 1
        for index in set_bits(word):
            t = t * generators[index] % n
        targets.append(t)
    return generators, rejected, words, targets


def extract(n, generators, words, records, witness):
    signed_product = 1
    denominator = 1
    target_counts = Counter()
    for i in set_bits(witness):
        target, u, v = records[i]
        signed_product *= u * v
        denominator = denominator * v % n
        target_counts[target] += 1
    assert signed_product > 0
    left = isqrt(signed_product)
    assert left * left == signed_product

    exponents = [0] * len(generators)
    for target, count in target_counts.items():
        for index in set_bits(words[target]):
            exponents[index] += count
    right = denominator
    for g, exponent in zip(generators, exponents):
        assert exponent % 2 == 0
        right = right * pow(g, exponent // 2, n) % n
    assert (left * left - right * right) % n == 0
    factor = gcd(left - right, n)
    if 1 < factor < n:
        assert n % factor == 0
        return factor
    return None


def dictionary_dependencies(n, generators, words, groups, primes):
    pivots = {}
    records = []
    rows = tested = 0
    factors = []
    for parity, labels, source in control.projected_rows(groups, primes, True, records):
        word = 0
        for target in set_bits(labels):
            word ^= words[target]
        row = parity | (word << (len(primes) + 1))
        rows += 1
        dependency = control.insert(row, source, pivots)
        if dependency is not None and tested < 40:
            assert dependency
            tested += 1
            factor = extract(n, generators, words, records, dependency)
            if factor is not None:
                factors.append(factor)
    nullity = rows - len(pivots)
    assert tested == min(40, nullity)
    return rows, nullity, tested, factors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=64)
    parser.add_argument("--pool", type=int, default=4096)
    args = parser.parse_args()
    if (args.bits < 5 or args.bits in (6, 8)) or args.pool < 1:
        parser.error("use a positive pool and a bit size with two distinct benchmark primes (5, 7 or ≥9)")
    primes = control.primes_upto(args.bits ** 2)
    primorial = prod(primes)
    print("bits,pool,seed,N,rejected,rows,base_nullity,control_even,control_tested,"
          "control_splits,word_even,word_tested,word_splits,factors,seconds")
    for seed in range(4):
        n = control.balanced_semiprime(args.bits, seed).N
        start = perf_counter()
        generators, rejected, words, targets = sample_targets(n, args.bits, seed, args.pool)
        groups = control.collect_groups(n, targets, primorial)
        rows, base_nullity, _, even, tested, splits = control.measure(n, targets, groups, primes, True)
        word_rows, word_even, word_tested, factors = dictionary_dependencies(
            n, generators, words, groups, primes)
        assert word_rows == rows and even <= word_even <= base_nullity
        print(args.bits, args.pool, seed, n, rejected, rows, base_nullity, even, tested,
              splits, word_even, word_tested, len(factors), "|".join(map(str, sorted(set(factors)))),
              f"{perf_counter() - start:.3f}", sep=",", flush=True)


if __name__ == "__main__":
    main()
