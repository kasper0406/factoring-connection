"""Held-out prediction from small-prime valuations of degree-five curve entries."""

from math import comb
from random import Random

import numpy as np
from sympy import isprime


def valuation(n, k, prime):
    other = n-k
    value = 0
    while n:
        n //= prime
        k //= prime
        other //= prime
        value += n-k-other
    return value


def features(n):
    m, t = (n-1)//2, (n-1)//5
    return [valuation(m,k,r) for r in (2,3,5,7,11,13) for k in (t,2*t)]


def sample_prime(rng, lo, hi, residue):
    start = lo + (residue-lo)%5
    while True:
        n = start + 5*rng.randrange((hi-1-start)//5+1)
        if isprime(n):
            return n


def sample(seed, per_class):
    rng = Random(seed)
    scale = 1 << 39
    rows, labels = [], []
    for label in (-1,1):
        for _ in range(per_class):
            residue = rng.choice((2,3) if label == 1 else (1,4))
            p = sample_prime(rng,scale,5*scale//4,residue)
            q = sample_prime(rng,3*scale//2,7*scale//4,pow(residue,-1,5))
            rows.append(features(p*q))
            labels.append(label)
    return np.asarray(rows), np.asarray(labels)


if __name__ == '__main__':
    for m in range(1,40):
        for k in range(m+1):
            for r in (2,3,5,7,11,13):
                v = valuation(m,k,r)
                assert comb(m,k) % r**v == 0 and comb(m,k) % r**(v+1) != 0
    train, y = sample(20261127,1024)
    center = train.mean(axis=0)
    weights = np.linalg.lstsq(train-center,y,rcond=None)[0]
    test, target = sample(20261128,2048)
    correct = ((test-center) @ weights >= 0) == (target == 1)
    print([int(correct[target==c].sum()) for c in (-1,1)], int(correct.sum()), len(test))
