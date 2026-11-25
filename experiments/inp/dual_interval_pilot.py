"""Fixed diagnostic on paired public interval centres. Factors are used in the harness only."""
from itertools import combinations
from math import comb, gcd
from random import Random
from fpylll import IntegerMatrix, LLL
from gmpy2 import is_prime

W, M = 8, 4
SEED = 114307
PER_SIDE = 36

def prime(rng, bits):
    while True:
        t = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(t):
            return t

def score(N, a, X):
    # Integer scaling of the dual to the power-ideal coefficient lattice.
    scale = (N**M) * X**(W - 1)
    rows = [[N**min(i, M) * comb(j, i) * (-a)**(j-i) * X**(W-1-j)
             if i <= j else 0 for j in range(W)] for i in range(W)]
    red = IntegerMatrix.from_matrix(rows)
    LLL.reduction(red)
    basis = [tuple(int(z) for z in row) for row in red]
    candidates = list(basis)
    for u, v in combinations(basis, 2):
        for sign in (1, -1):
            candidates.append(tuple(x + sign*y for x,y in zip(u,v)))
    best, hit = float("inf"), False
    for v in candidates:
        if not v[0]:
            continue
        numer = abs(v[0]) // X**(W-1)
        while numer and numer % N == 0:
            numer //= N
        if 1 < gcd(numer,N) < N:
            hit = True
        t = max(-1., min(1., v[1]/v[0]))
        signed_scale = (a + t*X)**M / scale * (1 if v[0] > 0 else -1)
        error = sum((signed_scale * c - t**j)**2 for j,c in enumerate(v))
        best = min(best, error)
    return best, hit

def run():
    rng = Random(SEED)
    for n in (48,64,80):
        for fraction in (0.12, 0.28, 0.34):
            X = 1 << round(n*fraction)
            taken = [0,0]
            wins = hits = siblings = 0
            while min(taken) < PER_SIDE:
                p,q = sorted((prime(rng,n//2),prime(rng,n//2)))
                N = p*q
                if p==q or N.bit_length()!=n:
                    continue
                lo = (p//(4*X))*4*X
                truth = (p-lo)//(2*X)
                if taken[truth] >= PER_SIDE or q < lo+4*X:
                    continue
                taken[truth] += 1
                scores = [score(N, lo + i*X, X) for i in (1,3)]
                (s,h), (t,k) = scores[truth], scores[1-truth]
                wins += (s<t) + 0.5*(s==t)
                hits += h
                siblings += k
            print(n,fraction,sum(taken),wins,hits,siblings,flush=True)

if __name__ == '__main__':
    run()
