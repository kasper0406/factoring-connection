"""Run from scratch with ../../../.venv/bin/python -B pilot.py > trials.csv."""
import csv
from fractions import Fraction
from math import isqrt
from pathlib import Path
from random import Random
import sys

import gmpy2
from flint import fmpz

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from src.floorsum.lattice_block import best_params, lattice_block


def prime(rng, bits):
    while True:
        x = rng.getrandbits(bits) | 1 | (1 << (bits - 1))
        if fmpz(x).is_prime():
            return x


def feasible(N, X):
    bits = N.bit_length() // 2
    maximum = (1 << bits) - 1
    lo = max(1 << (bits - 1), (N + maximum - 1) // maximum)
    return lo + 2 * X, isqrt(N) - 2 * X


def instance(rng, n, exponent):
    num, den = exponent
    while True:
        p, q = sorted((prime(rng, n // 2), prime(rng, n // 2)))
        N = p * q
        if p == q or N.bit_length() != n:
            continue
        X = int(gmpy2.iroot(N**num, den)[0])
        left, right = feasible(N, X)
        if left <= p <= right:
            return N, p - rng.randrange(2 * X), X, p, q


def slopes(N, c, cap):
    out = {(1, 1)}
    a, b = N, c*c
    u0, u1, v0, v1 = 0, 1, 1, 0
    while b:
        h, r = divmod(a, b)
        u0, u1 = u1, h*u1 + u0
        v0, v1 = v1, h*v1 + v0
        if max(u1, v1) > cap:
            break
        if min(u1, v1) > 0:
            out.add((u1, v1))
        a, b = b, r
    return sorted(out)


def tangent(N, lo, X):
    c = lo + X // 2
    best = None
    witness = 0
    for u, v in slopes(N, c, 2 * N.bit_length()):
        target = u*c + v*N // c
        for t in range(target - 4, target + 5):
            D = t*t - 4*u*v*N
            if D < 0:
                continue
            k = isqrt(D)
            if D - k*k > (k+1)**2 - D:
                k += 1
            gap = abs(D-k*k)
            for numerator in (t-k, t+k):
                if not 2*u*lo <= numerator < 2*u*(lo+X):
                    continue
                if gap == 0:
                    d, rem = divmod(numerator, 2*u)
                    if rem == 0 and 1 < d < N and N % d == 0:
                        witness = d
                else:
                    if best is None or gap*best.denominator < (2*k+1)*best.numerator:
                        best = Fraction(gap, 2*k+1)
    return witness, best


def score(N, lo, X):
    m, dim, W = best_params(N, isqrt(N//2), nmax=10)
    W = min(X, max(1, W))
    left, right = feasible(N, X)
    children = []
    for start in (lo, lo+X):
        block = start + (X-W)//2
        if W == 1:
            factors = [block] if N % block == 0 else []
        else:
            factors = lattice_block(N, block, m, dim, W)
        assert all(block <= d < block+W and N % d == 0 for d in factors)
        t, soft = tangent(N, start, X)
        overlap = max(0, min(start+X, right+1) - max(start, left))
        children.append((next(iter(factors), 0), t, soft, overlap, block))
    return children, W


def main():
    holdout = sys.argv[1:] == ['--holdout']
    if sys.argv[1:] and not holdout:
        raise SystemExit('usage: pilot.py [--holdout]')
    out = csv.writer(sys.stdout)
    out.writerow("bits exponent trial N lo X p q truth W block0 block1 copp tangent soft0 soft1 soft_pred support_pred coin pred".split())
    for i, n in enumerate((48, 72, 96)):
        for j, exponent in enumerate(((13, 50), (3, 10), (2, 5))):
            if holdout and (i,j) != (2,2):
                continue
            rng = Random(731906541 if holdout else 2026112500 + 10*i+j)
            for row in range(3072 if holdout else 384):
                N, lo, X, p, q = instance(rng, n, exponent)
                children, W = score(N, lo, X)
                (c0,t0,s0,g0,b0),(c1,t1,s1,g1,b1) = children
                truth = (p-lo)//X
                assert truth in (0, 1) and q >= lo+2*X and p*q == N
                for c,t,_,_,b in children:
                    assert (c != 0) == (b <= p < b+W)
                    assert c in (0,p) and t in (0,p)
                soft_tie, support_tie, coin = (rng.randrange(2) for _ in range(3))
                soft = soft_tie if s0 == s1 else int(s0 is None or (s1 is not None and s1 < s0))
                support = support_tie if g0 == g1 else int(g1 > g0)
                copp = 0 if c0 else 1 if c1 else -1
                tan = 0 if t0 else 1 if t1 else -1
                pred = copp if copp >= 0 else tan if tan >= 0 else soft
                assert all(k == truth for k in (copp,tan) if k >= 0)
                out.writerow((n, str(Fraction(*exponent)), row, N, lo, X, p, q, truth,
                              W,b0,b1,copp,tan,s0,s1,soft,support,coin,pred))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
