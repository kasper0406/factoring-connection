"""Measure sparse signed differences at the certified lattice precision."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, product

from check_lehman_window_baseline import cohorts
from check_shared_suffix_split import suffix_precision
from check_sum_low_bits import odd_square_roots
from check_sum_residue_lattice import query_classes


@lru_cache(None)
def weight(x, bits):
    """Minimum number of signed powers of two summing to x modulo 2**bits."""
    if bits == 0:
        return 0
    if x % 2 == 0:
        return weight(x//2, bits-1)
    return 1 + min(weight(x//2, bits-1),
                   weight(((x+1)//2) % (1 << (bits-1)), bits-1))


def sparse_split(n, budget):
    low, high, bits = suffix_precision(n)
    modulus = 1 << bits
    orientation = 1 if n % 4 == 1 else -1
    seen = set()
    for size in range(budget+1):
        for positions in combinations(range(2,bits),size):
            for signs in product((-1,1),repeat=size):
                gap = sum(sign*(1 << j) for sign,j in zip(signs,positions)) % modulus
                if gap in seen:
                    continue
                seen.add(gap)
                half = gap//2
                roots = odd_square_roots(half*half+orientation*n,bits)
                classes = [(half+r) % modulus for r in roots]
                factor = query_classes(n,low,high,modulus,classes,1,0)
                if factor:
                    return factor,len(seen)
    return None,len(seen)


if __name__ == '__main__':
    for size, planted, cases in cohorts():
        if planted or size == 16:
            continue
        histogram = Counter()
        for p,q in cases:
            n = p*q
            _,_,bits = suffix_precision(n)
            orientation = 1 if n % 4 == 1 else -1
            level = weight((q-orientation*p) % (1 << bits), bits)
            histogram[level] += 1
            if level <= 2:
                divisor, attempts = sparse_split(n,2)
                assert divisor in (p,q)
                print(n,divisor,attempts)
        print(size,sorted(histogram.items()))
