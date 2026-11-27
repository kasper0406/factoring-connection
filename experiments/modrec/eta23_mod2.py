"""E3: a(N) for f = eta(z)eta(23z) = sum a(n) q^n. Mod 2, Euler's pentagonal theorem gives
a(N) = #{(u,v): u^2 + 23 v^2 = 24N, u = v = -1 mod 6} (mod 2): a representation parity at index
24N > N (no descent). Check vs q-expansion and vs the cubic-splitting formula on N=pq."""
import numpy as np, sympy as sp
from math import isqrt
L = 100000
def euler(L, s):
    c = np.zeros(L+1, dtype=np.int64)
    k = 0
    while True:
        done = True
        for kk in ([k] if k == 0 else [k, -k]):
            e = s*kk*(3*kk-1)//2
            if e <= L:
                c[e] += (-1)**abs(kk); done = False
        if done and k > 0: break
        k += 1
    return c
A = euler(L, 1); B = euler(L, 23)
prod = np.zeros(L+1, dtype=np.int64)
nzB = np.nonzero(B)[0]
for j in nzB:
    prod[j:] += B[j]*A[:L+1-j]
a = np.zeros(L+2, dtype=np.int64); a[1:] = prod  # a(n) = coeff of q^(n-1) in product
def pent(N):
    c = 0
    M = 24*N
    for v in range(-isqrt(M//23), isqrt(M//23)+1):
        if v % 6 != 5: continue
        w = M - 23*v*v; u = isqrt(w)
        if u*u == w:
            for uu in {u, -u}:
                if uu % 6 == 5: c += 1
    return c
bad = sum((pent(N) - a[N]) % 2 != 0 for N in range(1, 3000))
print('pentagonal parity mismatches N<3000:', bad)
cub = sp.Poly(sp.Symbol('x')**3 - sp.Symbol('x') + 1)
def ap(p):
    if p == 23: return 1
    if sp.jacobi_symbol(-23, p) == -1: return 0
    r = len(sp.polys.galoistools.gf_factor_sqf(cub.all_coeffs()[::-1][::-1], p, sp.polys.domains.ZZ)[1])
    return 2 if r == 3 else -1
bad2 = n = 0
for N in range(15, L):
    fa = sp.factorint(N)
    if len(fa) == 2 and all(e == 1 for e in fa.values()) and 23 not in fa and 2 not in fa:
        p, q = fa
        n += 1; bad2 += (ap(p)*ap(q) != a[N])
print('a(pq) = a(p)a(q) via cubic splitting: semiprimes', n, 'mismatches', bad2)
