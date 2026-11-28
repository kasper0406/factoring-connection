"""Compare the noisy singleton decoder to every endpoint in small dyadic fibers."""

from itertools import combinations, product
from math import comb
from random import Random

from sympy import primerange

if __package__:
    from .check_cofactor_zero_queries import direct_zero
    from .noisy_singleton_lists import factor_from_noisy_singletons, noisy_singleton_list, singleton_list_bound
else:
    from check_cofactor_zero_queries import direct_zero
    from noisy_singleton_lists import factor_from_noisy_singletons, noisy_singleton_list, singleton_list_bound


def check_direct():
    total = 0
    for w in range(5):
        m = 2 << w
        for n in range(1, m, 2):
            pairs = [(a, n*pow(a, -1, m) % m) for a in range(1, m, 2)]
            for readings in product((None, 0, 1), repeat=w):
                observed = sum(1 << j for j, bit in enumerate(readings, 1) if bit is not None)
                ones = sum(1 << j for j, bit in enumerate(readings, 1) if bit == 1)
                costs = [(a, ((a^b) & observed).bit_count(), (((a&b)^ones) & observed).bit_count(),
                          (a & b & ones) == ones) for a, b in pairs]
                for cap in range(observed.bit_count()+1):
                    for errors in range(observed.bit_count()+1):
                        for certified in (False, True):
                            got, visits = noisy_singleton_list(n, readings, cap, errors, certified_ones=certified)
                            expected = {a for a, used, mistakes, valid in costs
                                        if used <= cap and mistakes <= errors and (not certified or valid)}
                            bound = singleton_list_bound(readings, cap, errors, certified_ones=certified)
                            assert set(got) == expected and len(got) <= bound
                            assert visits <= w*bound
                            total += 1
    return total


def check_native_inputs():
    total = mistakes = erasures = 0
    rng = Random(30102)
    for p, q in combinations(primerange(11, 110), 2):
        if q >= 2*p:
            continue
        n = p*q
        bits = [int(direct_zero(n, 1 << j)) for j in range(1, 5)]
        assert bits == [(p & q) >> j & 1 for j in range(1, 5)]
        for certified in (False, True):
            readings = bits.copy()
            erased = rng.randrange(4)
            readings[erased] = None
            erasures += 1
            choices = [j for j in range(4) if j != erased and (not certified or readings[j] == 1)]
            if choices:
                j = rng.choice(choices)
                readings[j] ^= 1
                mistakes += 1
            # Half the runs omit both budgets, allowing the decoder to choose them.
            result = factor_from_noisy_singletons(n, readings, 1 if (total+total//2) % 2 else None,
                                                  certified_ones=certified)
            assert result in (p, q)
            total += 1
    return total, mistakes, erasures


def check_planted_list():
    p, q = 2639587, 3890503
    n = p*q
    w = 20
    readings = [(p & q) >> j & 1 for j in range(1, w+1)]
    for j in (1, 11):
        readings[j] = None
    for j in (5, 15):
        readings[j] ^= 1
    cap = sum(((p ^ q) >> j & 1) for j, bit in enumerate(readings, 1) if bit is not None)
    residues, tested = noisy_singleton_list(n, readings, cap, 2)
    erased = readings.count(None)
    naive = (1 << (erased+cap))*sum(comb(w-erased, i) for i in range(3))
    improved = singleton_list_bound(readings, cap, 2)
    assert p % (2 << w) in residues and improved < naive
    return len(residues), tested, cap, improved, naive


def check_zero_baseline():
    # All ordered endpoints at 48 fixed uniformly drawn product residues. Grant
    # the exact two caps of each endpoint, then compute its compatible list size.
    w = 12
    m = 2 << w
    rng = Random(74319)
    outcomes = {}
    wrong = 0
    for n in rng.sample(range(1, m, 2), 48):
        cells = [[0]*(w+1) for _ in range(w+1)]
        for a in range(1, m, 2):
            b = n*pow(a, -1, m) % m
            e = (a & b).bit_count()-1
            cells[(a^b).bit_count()][e] += 1
            wrong += e
        for b in range(w+1):
            for e in range(w+1):
                if cells[b][e]:
                    size = sum(sum(row[:e+1]) for row in cells[:b+1])
                    outcomes[size] = outcomes.get(size, 0)+cells[b][e]
    total = sum(outcomes.values())
    cumulative = 0
    for size in sorted(outcomes):
        cumulative += outcomes[size]
        if 2*cumulative >= total:
            return total, 1-wrong/(w*total), size, 1 << w


if __name__ == "__main__":
    print("Exhaustive capped lists:", check_direct())
    print("Recovery from corrupted native readings (runs, mistakes, erasures):", check_native_inputs())
    print("Planted (residues, prefix visits, orientation cap, bound, naive bound):", check_planted_list())
    print("All-zero baseline (endpoints, accuracy, median capped list, whole fiber):", check_zero_baseline())
