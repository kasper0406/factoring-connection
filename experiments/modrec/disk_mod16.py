"""E2b: whole-disk count mod 16 = 1 + 4*(Sq(T) + 2*G(T)) mod 16, where
G(T) = #{n<=T: f(n) = 2 mod 4} = #{n = p^e s <= T: p=1 mod 4, e=1 mod 4, s in sq or 2sq, p not | s}
Sq(T) = sum_{n<=T, n in sq u 2sq} f(n) mod 4, f = r_2/4.  So the mod-16 disk count is a parity of a
split-prime count (Polymath4/TCH class) plus a square-indexed sum."""
import numpy as np
from math import isqrt
T = 10**6
f = np.zeros(T+1, dtype=np.int64)
# f(n) = sum_{d|n} chi(d)
for d in range(1, T+1, 2):
    f[d::d] += 1 if d % 4 == 1 else -1
# direct disk counts
C = np.zeros(T+1, dtype=np.int64)
r2 = np.zeros(T+1, dtype=np.int64)
for x in range(-isqrt(T), isqrt(T)+1):
    ys = np.arange(-isqrt(T - x*x), isqrt(T - x*x)+1)
    np.add.at(r2, x*x + ys*ys, 1)
C = np.cumsum(r2)
# classification of f mod 4 via factorization-free sieve of p^e (p=1 mod 4) structure is f itself;
# check: f odd <=> n in sq u 2sq
n = np.arange(T+1)
sq = np.zeros(T+1, bool)
for m in range(0, isqrt(T)+1):
    sq[m*m] = True
    if 2*m*m <= T: sq[2*m*m] = True
assert np.all(((f[1:] % 2) == 1) == sq[1:])
G = np.cumsum((f % 4) == 2) ; G[0] = 0
Sq = np.cumsum(np.where(sq, f, 0)); 
lhs = C % 16
rhs = (1 + 4*(Sq + 2*G)) % 16
print('identity holds for all T<=1e6:', bool(np.all(lhs[1:] == rhs[1:])))
# classification of f = 2 mod 4: exactly one p=1(4) with odd exponent, which is 1 mod 4; checked on a sample with sympy
import sympy as sp, random
rng = random.Random(3); bad = 0
for _ in range(20000):
    m = rng.randrange(1, T+1)
    fa = sp.factorint(m)
    odd1 = [(p, e) for p, e in fa.items() if p % 4 == 1 and e % 2]
    odd3 = [(p, e) for p, e in fa.items() if p % 4 == 3 and e % 2]
    pred = (not odd3) and len(odd1) == 1 and odd1[0][1] % 4 == 1
    bad += pred != (f[m] % 4 == 2)
print('f=2 mod 4 classification mismatches:', bad)
