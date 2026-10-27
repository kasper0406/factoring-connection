# Tiny sanity checks for conn_top (survey). Exact arithmetic only.
import numpy as np
from math import gcd, comb
from sympy import jacobi_symbol, legendre_symbol, primerange

def unitary_cayley_spectrum(N):
    A = np.array([[1 if gcd((i-j)%N, N)==1 else 0 for j in range(N)] for i in range(N)])
    ev = np.round(np.linalg.eigvalsh(A)).astype(int)
    vals, cnt = np.unique(ev, return_counts=True)
    return dict(zip(vals.tolist(), cnt.tolist()))

for p,q in [(3,5),(5,7),(7,11)]:
    N=p*q; s=unitary_cayley_spectrum(N)
    exp={ (p-1)*(q-1):1, -(q-1):p-1, -(p-1):q-1, 1:(p-1)*(q-1)}
    print("unitary Cayley", N, s, "expected", exp, s==exp)

# Jacobi-sum over Z/N of a cubic = a_p * a_q
def ap(f,p): return -sum(legendre_symbol(f(x)%p,p) if f(x)%p else 0 for x in range(p))
f=lambda x: x**3+2*x+3
for p,q in [(101,103),(211,307)]:
    N=p*q
    S=sum(jacobi_symbol(f(x)%N,N) if gcd(f(x),N)==1 else 0 for x in range(N))
    # terms with gcd>1 contribute 0 in Jacobi product too
    print("sum Jacobi", N, S, "a_p a_q", ap(f,p)*ap(f,q), S==ap(f,p)*ap(f,q))

# lozenge count of hexagon (1,x,x) = C(2x,x): gcd with N reveals p for x in [p/2, q/2)
p,q=1009,1201; N=p*q
hits=[x for x in range(1,700) if 1<gcd(comb(2*x,x)%N,N)<N]
print("C(2x,x) splits N for x in", (min(hits),max(hits)) if hits else None, "p/2,q/2 =", p/2, q/2)

# Sum-free density (Diananda-Yap/Green-Ruzsa) brute force for tiny cyclic groups
from itertools import combinations
def maxsumfree(N):
    best=0
    # greedy over arithmetic-progression-in-dilate candidates: exact search too costly; use known extremal shapes check for N<=15
    for r in range(N//2+1, 0, -1):
        for A in combinations(range(1,N), r):
            S=set(A)
            if all(((a+b)%N) not in S for a in A for b in A):
                return r
    return 0
for N in [10,14,15]:
    print("max sum-free Z/%d ="%N, maxsumfree(N))
