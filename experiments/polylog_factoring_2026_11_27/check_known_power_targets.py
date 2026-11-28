"""Measure the value of known power exponents on the same Euclidean rows.

For t=2^e and u=tv modulo N, the signed integer uv·2^(e mod 2)
has the public modular square root v·2^ceil(e/2). Compare the same
targets with and without these exponents. Reuse the one-opaque-cofactor
policy of check_cross_target_dependencies. Count a factor only after
checking both the integer square and the modular congruence.
"""

import argparse
from math import gcd, isqrt, prod
from random import Random

if __package__:
    from . import check_cross_target_dependencies as control
else:
    import check_cross_target_dependencies as control


def extract(n, records, witness, extra_twos, roots):
    product = 1
    modular_root = 1
    while witness:
        bit = witness & -witness
        witness ^= bit
        label, u, v = records[bit.bit_length() - 1]
        product *= (u * v) << extra_twos[label]
        modular_root = modular_root * v * roots[label] % n
    assert product > 0
    integer_root = isqrt(product)
    assert integer_root * integer_root == product
    assert (integer_root * integer_root - modular_root * modular_root) % n == 0
    factor = gcd(integer_root - modular_root, n)
    if 1 < factor < n:
        assert n % factor == 0
        return factor
    return None


def measure_known(n, extra_twos, roots, groups, primes):
    """Eliminate atom parities with public roots attached to each original row."""
    pivots = {}
    records = []
    projected = tested = splits = 0
    two_mask = sum(bit << i for i, bit in enumerate(extra_twos))
    assert primes[0] == 2  # Sign is column zero, so 2 is column one.
    for base, labels, witness in control.projected_rows(groups, primes, True, records):
        if (labels & two_mask).bit_count() % 2:
            base ^= 1 << 1
        projected += 1
        dependency = control.insert(base, witness, pivots)
        if dependency is not None and tested < 40:
            tested += 1
            splits += extract(n, records, dependency, extra_twos, roots) is not None
    return projected - len(pivots), tested, splits


def make_targets(n, bits, seed, count, cohort):
    if cohort == 'powers':
        rng = Random(f'known-power-targets:{bits}:{seed}')
        exponents = [rng.getrandbits(bits + 64) for _ in range(count)]
        return ([pow(2, e, n) for e in exponents], [e % 2 for e in exponents],
                [pow(2, (e + 1) // 2, n) for e in exponents])
    rng = Random(f'relation-baseline:{bits}:{seed}')
    roots = []
    while len(roots) < count:
        x = rng.randrange(1, n)
        if gcd(x, n) == 1:
            roots.append(x)
    return [x * x % n for x in roots], [0] * count, roots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pools', type=int, nargs='*', default=[1024, 2048, 4096, 8192])
    parser.add_argument('--bits', type=int, default=64)
    parser.add_argument('--cohort', choices=['powers', 'squares'], default='powers')
    args = parser.parse_args()
    if min(args.pools) < 1 or (args.bits < 5 or args.bits in (6, 8)):
        parser.error('use a positive pool and a bit size with two distinct benchmark primes (5, 7 or ≥9)')
    primes = control.primes_upto(args.bits ** 2)
    primorial = prod(primes)
    print('bits,pool,cohort,seed,convergents,projected,no_witness_deps,no_witness_tested,'
          'no_witness_splits,known_deps,known_tested,known_splits')
    for count in args.pools:
        totals = [0] * 8
        for seed in range(4):
            n = control.balanced_semiprime(args.bits, seed).N
            targets, extra_twos, roots = make_targets(n, args.bits, seed, count, args.cohort)
            groups = control.collect_groups(n, targets, primorial)
            row_count = sum(len(group) for group in groups.values())
            projected, _, _, zero_labels, tested, splits = control.measure(
                n, targets, groups, primes, True)
            known = measure_known(n, extra_twos, roots, groups, primes)
            counts = (row_count, projected, zero_labels, tested, splits, *known)
            print(args.bits, count, args.cohort, seed, *counts, sep=',', flush=True)
            totals = [a + b for a, b in zip(totals, counts)]
        print(args.bits, count, args.cohort, 'sum', *totals, sep=',', flush=True)


if __name__ == '__main__':
    main()
