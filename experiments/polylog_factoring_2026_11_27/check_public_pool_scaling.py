"""Scale the public-unit, one-opaque-atom control on matched 64-bit inputs.

Each pool is an independent replay of a prefix of the same seeded stream.
Timings include generation of the targets, all Euclidean rows, grouping,
elimination, and verification of up to 40 zero-label dependencies per
modulus and pool. The prime table and benchmark modulus are prepared once.
Run directly or with -m from the repository root.
"""

import argparse
from math import gcd, prod
from random import Random
from time import perf_counter

if __package__:
    from . import check_cross_target_dependencies as control
else:
    import check_cross_target_dependencies as control


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pools", type=int, nargs="*", default=[1024, 2048, 4096, 8192])
    parser.add_argument("--bits", type=int, default=64)
    args = parser.parse_args()
    if min(args.pools) < 1 or args.bits < 5 or args.bits in (6, 8):
        parser.error("use a positive pool and a bit size with two distinct benchmark primes (5, 7 or ≥9)")
    primes = control.primes_upto(args.bits ** 2)
    primorial = prod(primes)
    print("bits,pool,seed,convergents,projected_rows,nullity,label_rank,zero_labels,tested,splits,seconds")
    for pool in args.pools:
        totals = [0] * 7
        total_seconds = 0.0
        for seed in range(4):
            n = control.balanced_semiprime(args.bits, seed).N
            start = perf_counter()
            rng = Random(f"relation-baseline:{args.bits}:{seed}")
            targets = []
            while len(targets) < pool:
                t = rng.randrange(1, n)
                if gcd(t, n) == 1:
                    targets.append(t)
            groups = control.collect_groups(n, targets, primorial)
            convergents = sum(len(group) for group in groups.values())
            counts = (convergents, *control.measure(n, targets, groups, primes, True))
            seconds = perf_counter() - start
            print(args.bits, pool, seed, *counts, f"{seconds:.3f}", sep=",", flush=True)
            totals = [a + b for a, b in zip(totals, counts)]
            total_seconds += seconds
            del groups
        print(args.bits, pool, "sum", *totals, f"{total_seconds:.3f}", sep=",", flush=True)


if __name__ == "__main__":
    main()
