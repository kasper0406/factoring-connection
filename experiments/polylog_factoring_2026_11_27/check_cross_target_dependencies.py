"""Measure which Euclidean relations retain hidden-root randomness.

Draw 1,024 targets for each of four balanced semiprimes per bit size. Both
cohorts use the same pseudorandom units x. One supplies t=x² mod N to the
collector and the control supplies t=x. Keep the private x out of the
relation collector and out of every gcd.

For u=tv mod N, strip primes up to bits² from |uv|. A square remainder
is free. In the opaque-atom treatment, identical nonsquare remainders can
also cancel in pairs. Let A contain the remaining sign and base-prime
parities, and B the target-use parities. After cofactor cancellation the
dimension of B(ker A) is rank([A;B])−rank(A). A dependency in ker A with
zero label B uses each target an even number of times, and thus permits a
gcd using only public values. Verify the first 40 such basis vectors.

Run directly or as a module from the repository root, optionally passing
bit sizes, for example 48 64 80 96. All ranks are computed per modulus.
"""

import argparse
from collections import Counter, defaultdict
from math import gcd, isqrt, prod
from pathlib import Path
from random import Random
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.instances import balanced_semiprime
from src.smoothbuild.core import primes_upto


def collect_groups(n, targets, primorial):
    groups = defaultdict(list)
    for label, t in enumerate(targets):
        old_u, old_v = n, 0
        u, v = t, 1
        while u:
            assert (u - v * t) % n == 0
            if gcd(v, n) == 1:
                magnitude = u * abs(v)
                rough = magnitude
                while (common := gcd(rough, primorial)) > 1:
                    rough //= common
                atom = 1 if isqrt(rough) ** 2 == rough else rough
                groups[atom].append((magnitude // rough, label, u, v))
            quotient, remainder = divmod(old_u, u)
            old_u, u = u, remainder
            old_v, v = v, old_v - quotient * v
    return groups


def sign_and_base_parity(smooth, negative, primes):
    row = int(negative)
    for column, prime in enumerate(primes, 1):
        while smooth % prime == 0:
            smooth //= prime
            row ^= 1 << column
        if smooth == 1:
            break
    assert smooth == 1
    return row


def insert(row, witness, pivots):
    """Return the reduced witness for a dependency, or None for a new pivot."""
    while row:
        lead = row.bit_length() - 1
        if lead not in pivots:
            pivots[lead] = row, witness
            return None
        pivot, source = pivots[lead]
        row ^= pivot
        witness ^= source
    return witness


def public_split(n, targets, records, witness):
    magnitude = 1
    denominator = 1
    uses = Counter()
    while witness:
        bit = witness & -witness
        witness ^= bit
        label, u, v = records[bit.bit_length() - 1]
        magnitude *= u * v
        denominator = denominator * v % n
        uses[label] += 1
    assert magnitude > 0 and all(count % 2 == 0 for count in uses.values())
    left = isqrt(magnitude)
    assert left * left == magnitude
    right = denominator
    for label, count in uses.items():
        right = right * pow(targets[label], count // 2, n) % n
    assert (left * left - right * right) % n == 0
    factor = gcd(left - right, n)
    if 1 < factor < n:
        assert n % factor == 0
        return factor
    return None


def projected_rows(groups, primes, allow_opaque, records):
    """Yield parity, target labels, and source mask after cancelling whole cofactors."""
    for atom, group in groups.items():
        if atom != 1 and (not allow_opaque or len(group) < 2):
            continue
        anchor = None
        for smooth, label, u, v in group:
            base = sign_and_base_parity(smooth, v < 0, primes)
            labels = 1 << label
            witness = 1 << len(records)
            records.append((label, u, v))
            if atom != 1:
                if anchor is None:
                    anchor = base, labels, witness
                    continue
                base ^= anchor[0]
                labels ^= anchor[1]
                witness ^= anchor[2]
            yield base, labels, witness


def measure(n, targets, groups, primes, allow_opaque):
    base_pivots = {}
    labelled_pivots = {}
    records = []
    rows = tested = splits = 0
    for base, labels, witness in projected_rows(groups, primes, allow_opaque, records):
        rows += 1
        insert(base, 0, base_pivots)
        row = base | (labels << (len(primes) + 1))
        dependency = insert(row, witness, labelled_pivots)
        if dependency is not None and tested < 40:
            assert dependency
            tested += 1
            splits += public_split(n, targets, records, dependency) is not None
    nullity = rows - len(base_pivots)
    label_rank = len(labelled_pivots) - len(base_pivots)
    zero_labels = nullity - label_rank
    assert tested == min(40, zero_labels)
    return (rows, nullity, label_rank, zero_labels, tested, splits)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bits', type=int, nargs='*', default=[48, 64])
    args = parser.parse_args()
    print('bits,seed,cohort,atoms,rows,nullity,label_rank,zero_labels,tested,splits')
    for bits in args.bits:
        primes = primes_upto(bits ** 2)
        primorial = prod(primes)
        totals = {}
        for seed in range(4):
            n = balanced_semiprime(bits, seed).N
            rng = Random(f'relation-baseline:{bits}:{seed}')
            units = []
            while len(units) < 1024:
                unit = rng.randrange(1, n)
                if gcd(unit, n) == 1:
                    units.append(unit)
            for cohort in ('squares', 'units'):
                targets = [x * x % n for x in units] if cohort == 'squares' else units
                groups = collect_groups(n, targets, primorial)
                for allow_opaque in (False, True):
                    key = cohort, 'opaque' if allow_opaque else 'square'
                    counts = measure(n, targets, groups, primes, allow_opaque)
                    print(bits, seed, *key, *counts, sep=',')
                    previous = totals.get(key, (0,) * len(counts))
                    totals[key] = tuple(a + b for a, b in zip(previous, counts))
        for key, counts in totals.items():
            print(bits, 'sum', *key, *counts, sep=',')


if __name__ == '__main__':
    main()
