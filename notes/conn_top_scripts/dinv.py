# d-invariants of lens spaces L(N,a) (Ozsvath-Szabo recursion), exact rationals.
# Check codex's "collision" claim: d_i = d_j with j != a-1-i (mod N) => gcd split.
from fractions import Fraction as F
from math import gcd, isqrt
from functools import lru_cache
import random, sys
sys.setrecursionlimit(10000)
def d(p,q,i):
    if p==1: return F(0)
    return F(-1,4) + F((2*i+1-p-q)**2, 4*p*q) - d(q, p%q, i%q)
from sympy import randprime
random.seed(1)
rows=[]
for bits in (12,14,16):
    for t in range(6):
        p=randprime(2**(bits//2-1),2**(bits//2)); q=randprime(2**(bits//2-1),2**(bits//2))
        if p==q: continue
        N=p*q; a=isqrt(N)
        while gcd(a,N)!=1: a+=1
        vals={}
        for i in range(N):
            vals.setdefault(d(N,a,i),[]).append(i)
        nontriv=0; splits=0; frac_only=0
        for v,idx in vals.items():
            for x in range(len(idx)):
                for y in range(x+1,len(idx)):
                    i,j=idx[x],idx[y]
                    if (i+j-(a-1))%N!=0:
                        nontriv+=1
                        g=gcd((i-j)*(i+j-a+1),N)
                        if 1<gcd(i-j,N)<N or 1<gcd(i+j-a+1,N)<N: splits+=1
        # fractional-part collisions (mod 1) that are non-conjugate: count i with partner
        rows.append((bits,N,a,len(vals),nontriv,splits))
        print(bits,N,"distinct d:",len(vals),"nontrivial collisions:",nontriv,"of which split:",splits, flush=True)

# structure of collisions for one N
p,q=101,307; N=p*q; a=isqrt(N)
while gcd(a,N)!=1: a+=1
vals={}
for i in range(N): vals.setdefault(d(N,a,i),[]).append(i)
ex=[]
for v,idx in vals.items():
    for x in range(len(idx)):
        for y in range(x+1,len(idx)):
            i,j=idx[x],idx[y]
            if (i+j-(a-1))%N!=0: ex.append((i,j,gcd(i-j,N),gcd(i+j-a+1,N)))
print("N",N,"a",a,"collisions",len(ex))
from collections import Counter
print("gcd(i-j) pattern", Counter((g1,g2) for _,_,g1,g2 in ex))
print("sample", ex[:12])
print("|i-j| dist (min, median, max)", min(abs(i-j) for i,j,_,_ in ex), sorted(abs(i-j) for i,j,_,_ in ex)[len(ex)//2], max(abs(i-j) for i,j,_,_ in ex))
