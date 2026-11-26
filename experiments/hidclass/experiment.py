"""Reproduce the frozen design in prereg.md with Python's standard library."""
import csv
from math import gcd, isqrt, sqrt
from random import Random
from statistics import NormalDist


def prime(n):
    # The bases 2,3,5,7,11,13 are deterministic below 3,474,749,660,383.
    small = (2, 3, 5, 7, 11, 13)
    if any(n % a == 0 for a in small):
        return n in small
    d = n - 1
    e = (d & -d).bit_length() - 1
    d >>= e
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(e - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sample_prime(rng, bits, m, residue):
    low, high = 1 << (bits - 1), 1 << bits
    first = low + (residue - low) % m
    while True:
        n = rng.randrange(first, high, m)
        if prime(n):
            return n


def cf_features(n):
    disc = 7 * n
    a0 = isqrt(disc)
    a, h, d = a0, 0, 1
    prev, denom = 0, 1  # q[-1], q[0]
    parity = 0
    splits = set()
    for _ in range(64):
        h = a * d - h
        d = (disc - h * h) // d
        a = (a0 + h) // d
        parity ^= a & 1
        prev, denom = denom, a * denom + prev
        for value in (denom, d):
            g = gcd(n, value)
            if 1 < g < n:
                splits.add(g)
    return (parity, denom % 2), splits


def jacobi(a, n):
    a %= n
    s = 1
    while a:
        k = (a & -a).bit_length() - 1
        if k % 2 and n % 8 in (3, 5):
            s = -s
        a >>= k
        if a % 4 == n % 4 == 3:
            s = -s
        a, n = n % a, a
    return s if n == 1 else 0


def cubic_feature(rng, n):
    sign = 1
    splits = set()
    zeros = 0
    for _ in range(4):
        _, b, c = (rng.randrange(n) for _ in range(3))
        u, v, w = 2*b - 3*c, b + c, 2*c - b  # a - sigma(a)
        # Determinant of multiplication in the basis 1,x,x^2
        norm = (u * ((u+2*w)*(u-v+3*w) - (2*v-w)*(v-w))
                - w * (v*(u-v+3*w) - w*(2*v-w))
                + (v-w) * (v*(v-w) - w*(u+2*w)))
        g = gcd(n, norm)
        if g > 1:
            zeros += 1
            if g < n:
                splits.add(g)
        sign *= jacobi(norm, n)
    return ((int(sign == 1),) if sign else (None,)), splits, zeros


def describe(cells, level):
    z = NormalDist().inv_cdf((1 + level) / 2)
    # Independent strata and labels, equally weighted strata
    delta = variance = 0.0
    for c in cells:
        for label in (0, 1):
            hits, total = c[label]
            rate = hits / total
            delta += (2*label - 1) * rate
            variance += rate * (1-rate) / total
    scale = 50 / len(cells)  # percentage points of accuracy over chance
    return tuple(round(scale * v, 4) for v in
                 (delta, delta - z*sqrt(variance), delta + z*sqrt(variance)))


def main():
    rng = Random(20261126)
    sizes = range(40, 81, 8)
    features = ('cf_sum', 'cf_denom', 'cubic')
    counts = {(b,f): [[0,0], [0,0]] for b in sizes for f in features}
    split_rows = 0
    zero_symbols = 0
    with open('samples.csv', 'w', newline='') as out:
        rows = csv.writer(out)
        rows.writerow(('bits','target','label','p','q','N',*features,'proper_gcds','zero_symbols'))
        for bits in sizes:
            jobs = []
            for label in (0, 1):
                pairs = [(1,1),(6,6)] if label == 0 else [(2,4),(4,2),(3,5),(5,3)]
                for i in range(1024):
                    jobs.append((1, label, 16, (1,3) if label == 0 else (9,11)))
                    jobs.append((2, label, 7, pairs[i % len(pairs)]))
            rng.shuffle(jobs)
            for target, label, m, (s,t) in jobs:
                while True:
                    p = sample_prime(rng,bits//2,m,s)
                    q = sample_prime(rng,bits//2,m,t)
                    n = p*q
                    if p != q and n.bit_length() == bits:
                        break
                if target == 1:
                    values, splits = cf_features(n)
                    zeros = 0
                    values = (*values, None)
                else:
                    vals, splits, zeros = cubic_feature(rng,n)
                    values = (None, None, *vals)
                split_rows += bool(splits)
                zero_symbols += zeros
                rows.writerow((bits,target,label,p,q,n,*values,' '.join(map(str,sorted(splits))),zeros))
                for f, value in zip(features, values):
                    if value is not None:
                        cell = counts[bits,f][label]
                        cell[0] += value
                        cell[1] += 1
    print('Valid rows by design: 1,024 per target, label and bit size.')
    print('Rows with a proper gcd:', split_rows, 'Zero cubic symbols:', zero_symbols)
    print('Bits | CF quotient-sum parity | CF denominator parity | cubic sign product')
    print('Each tuple gives signed advantage in percentage points and its 95% normal interval.')
    for bits in sizes:
        print(bits, *(describe([counts[bits,f]], .95) for f in features))
    for f in features:
        cells = [counts[b,f] for b in sizes]
        print(f, 'pooled', describe(cells,.95), 'adjusted', describe(cells,1-.05/3),
              'halves', describe(cells[:3],.95), describe(cells[3:],.95))


if __name__ == '__main__':
    main()
