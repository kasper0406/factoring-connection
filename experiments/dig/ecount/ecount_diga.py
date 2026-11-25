"""Run from scratch with the repo virtualenv and Python -B."""
import json
import math
from pathlib import Path

import gmpy2
import numpy as np
from cypari2 import Pari
from flint import fmpz_mod_poly_ctx, fmpz_poly

SEED = 370231
BITS = (40, 60, 80)
SAMPLES = 2400
MODEL = [0, 0, 0, -1, 1]

def public_features(N):
    R = fmpz_mod_poly_ctx(N)
    x = R.gen()
    f2 = x**3 - x + 1
    f3 = x**4 - 2*x**2 + 4*x - pow(3, -1, N)
    signs = []
    split = False
    for f in (f2, f3):
        u = pow(x, N, f) - x
        f_int = fmpz_poly([int(a) for a in f.coeffs()])
        for c in (0, 1):
            lifted = fmpz_poly([int(a) for a in (u+c).coeffs()])
            v = int(f_int.resultant(lifted)) % N
            d = math.gcd(v, N)
            split |= 1 < d < N
            signs.append(int(gmpy2.jacobi(v, N)))
    return signs, split

def draw_prime(rng, bits, used):
    while True:
        t = int(rng.integers(1 << (bits-1), 1 << bits)) | 1
        if t not in used and t > 23 and gmpy2.is_prime(t):
            return t

def direct_count(p):
    return p + 1 + sum(gmpy2.jacobi((x*x*x-x+1) % p, p) for x in range(p))

def center(a, groups):
    a = np.asarray(a, dtype=float).copy()
    for g in groups:
        a[g] -= a[g].mean(axis=0)
    norm = np.sqrt((a*a).sum(axis=0))
    return a/norm

def main():
    rng = np.random.default_rng(SEED)
    pari = Pari()
    curve = pari.ellinit(MODEL)
    records = []
    splits = {}
    for n in BITS:
        used = set()
        splits[n] = 0
        while len(used) < 2*SAMPLES:
            p,q = (draw_prime(rng,n//2,used) for _ in range(2))
            N = p*q
            if p==q or N.bit_length()!=n:
                continue
            used.update((p,q))
            features, split = public_features(N)
            M = (p+1-int(pari.ellap(curve,p))) * (q+1-int(pari.ellap(curve,q)))
            character = int(gmpy2.jacobi(-23,N))
            assert character==1 or M%2==0
            if n==40 and len(used)==2:
                assert M == direct_count(p)*direct_count(q)
            if split:
                splits[n] += 1
            else:
                records.append(((n,N%24,character),features,M%2,M%3))
    strata = [r[0] for r in records]
    groups = [np.array([i for i,t in enumerate(strata) if t==s]) for s in sorted(set(strata))]
    feat = np.array([r[1] for r in records])
    Y = np.array([[p,int(m==0),int(m==1)] for _,_,p,m in records])
    A,B = center(feat,groups),center(Y,groups)
    corr = A.T @ B
    observed = np.max(np.abs(corr))
    null = []
    for _ in range(999):
        idx = np.arange(len(Y))
        for g in groups:
            idx[g] = rng.permutation(g)
        null.append(np.abs(A.T @ B[idx]).max())
    planted = Y[:,0].copy()
    mask = rng.random(len(Y)) < .2
    planted[mask] = feat[mask,0]>0
    plant = float(A[:,0] @ center(planted,groups))
    by_size = {}
    for n in BITS:
        rows = [r for r in records if r[0][0]==n]
        by_size[n] = dict(clean=len(rows),odd=sum(r[2] for r in rows),
                          mod3=[sum(r[3]==k for r in rows) for k in range(3)])
    result = dict(samples=SAMPLES,splits=splits,by_size=by_size,
                  corr=corr.tolist(),max_r=float(observed),p=(1+sum(float(v)>=observed for v in null))/1000,
                  null95=float(np.quantile(null,.95)),control_r=plant,
                  odd_with_negative_sign=sum(r[2] for r in records if r[0][2]==-1))
    out = json.dumps(result,indent=2,default=float)
    Path(__file__).with_name("ecount_diga_result.json").write_text(out+"\n")
    print(out)
if __name__ == "__main__":
    main()
