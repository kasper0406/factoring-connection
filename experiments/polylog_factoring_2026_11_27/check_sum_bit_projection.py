"""Certify that a fixed odd-prime sieve leaves every admissible sum ending."""

from math import isqrt, prod

import numpy as np


def ending_count(n, bits):
    """Size of the image u + N/u on odd residues modulo 2^bits."""
    def at_least(j):
        if j == 1:
            roots = 1
        elif j == 2:
            roots = 2 if n % 4 == 1 else 0
        else:
            roots = 4 if n % 8 == 1 else 0
        return roots << (bits-j)

    cutoff = (bits+1)//2
    return (sum((at_least(j)-at_least(j+1))//(2 << j) for j in range(1,cutoff))
            + at_least(cutoff)//(1 << (bits//2)))


def certify(n, bits, primes):
    step = 1 << bits
    period = prod(primes)
    width = isqrt(9*n//2) - isqrt(4*n-1)
    block_limit = (width//step)//2
    blocks = (period+block_limit-1)//block_limit
    # Balanced integer partition of the period. Two adjacent blocks fit
    # inside the minimum number of possible lifts of any ending.
    edges = np.arange(blocks+1, dtype=np.int64)*period//blocks
    assert 2*int(np.diff(edges).max()) <= width//step
    candidate = edges[:-1].copy()
    ends = edges[1:]
    masks = []
    for p in primes:
        squares = {x*x % p for x in range(p)}
        masks.append(np.array([((step*t)**2-4*n)%p in squares for t in range(p)]))
    prefix = prod(primes[:5])
    residues = np.arange(prefix)
    keep = np.logical_and.reduce([mask[residues%p] for mask,p in zip(masks[:5],primes[:5])])
    allowed = residues[keep]
    cyclic = np.append(allowed, allowed[0]+prefix)
    jump = cyclic[np.searchsorted(cyclic,residues)]-residues
    while len(candidate):
        candidate += jump[candidate%prefix]
        assert np.all(candidate < ends), 'No certificate in this block'
        good = np.logical_and.reduce([mask[candidate%p] for mask,p in zip(masks[5:],primes[5:])])
        candidate = candidate[~good]+1
        ends = ends[~good]
    return blocks, width//step


if __name__ == '__main__':
    for k in range(1,13):
        for n in range(1,32,2):
            mod = 1 << k
            assert ending_count(n,k)==len({(x+n*pow(x,-1,mod))%mod for x in range(1,mod,2)})
    n = 118464083331556819902753234409
    primes = (3,5,7,11,13,17,19,23,29,31)
    print(ending_count(n,26), certify(n,26,primes))
