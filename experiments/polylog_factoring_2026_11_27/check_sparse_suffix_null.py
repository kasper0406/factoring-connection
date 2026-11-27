"""Conditional local coverage of the two-term signed-gap search."""

from functools import lru_cache
from itertools import combinations, product
from math import prod

from check_lehman_window_baseline import cohorts
from check_shared_suffix_split import suffix_precision


@lru_cache(None)
def proposals(bits, forced):
    mod = 1 << bits
    gaps = {0}
    for size in (1,2):
        for positions in combinations(range(2,bits),size):
            for signs in product((-1,1),repeat=size):
                gaps.add(sum(sign*(1 << j) for sign,j in zip(signs,positions))%mod)
    return {gap for gap in gaps if gap%8 == forced}


if __name__ == '__main__':
    for bits in range(3,11):
        m = 1 << bits
        for n in range(1,16,2):
            en = n if n%4 == 1 else -n
            counts = [0]*m
            for x in range(1,m,2):
                counts[(x-en*pow(x,-1,m))%m] += 1
            assert all(c == (4 if (d+en-1)%8 == 0 else 0) for d,c in enumerate(counts))
    for size,planted,cases in cohorts():
        if planted or size == 16:
            continue
        probabilities = []
        hits = 0
        for p,q in cases:
            n = p*q
            _,_,k = suffix_precision(n)
            e = 1 if n%4 == 1 else -1
            gaps = proposals(k,(1-e*n)%8)
            hits += (q-e*p)%(1 << k) in gaps
            probabilities.append(len(gaps)/(1 << (k-3)))
        print(size,hits,sum(probabilities),1-prod(1-v for v in probabilities))
