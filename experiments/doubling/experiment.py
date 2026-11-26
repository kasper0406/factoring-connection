import argparse
import ctypes
import json
import math
import multiprocessing
import random
import time
from functools import lru_cache
from pathlib import Path

import flint
import gmpy2
from flint import fmpz

HERE = Path(__file__).resolve().parent
LIB = ctypes.CDLL(str(next((Path(flint.__file__).parent / ".dylibs").glob("libflint.*.dylib"))))
W = ctypes.c_ulong

def native(name, ret, *args):
    f = getattr(LIB, name)
    f.restype, f.argtypes = ret, list(args)
    return f

preinv = native("n_preinvert_limb", W, W)
factorial_native = native("n_factorial_mod2_preinv", W, W, W, W)

@lru_cache(maxsize=128)
def factorial(n, p):
    return factorial_native(n, p, preinv(p))

def central(n, p):
    value = 1
    while n:
        n, d = divmod(n, p)
        if 2*d >= p:
            return 0
        value = value * factorial(2*d, p) * pow(factorial(d, p), -2, p) % p
    return value

def exact(N):
    (p, e), (q, f) = sorted((int(s), e) for s, e in fmpz(N).factor())
    assert e == f == 1
    m = N//4
    u, v = central(m, p), central(m, q)
    return u + p*((v-u)*pow(p, -1, q) % q), p, q

def gauss_a(p):
    x = int(fmpz(p-1).sqrtmod(p))
    r, t = p, min(x, p-x)
    while t*t > p:
        r, t = t, r % t
    b = math.isqrt(p-t*t)
    assert b*b+t*t == p
    a = t if t % 2 else b
    return a if a % 4 == 1 else -a

def lucas(P, Q, k, N):
    # Consecutive Lucas terms U_j and U_(j+1).
    u, v = 0, 1
    for d in bin(k)[2:]:
        a = u*(2*v-P*u) % N
        b = (v*v-Q*u*u) % N
        if d == '0':
            u, v = a, b
        else:
            u, v = b, (P*b-Q*a) % N
    return u

SMALL = (3, 5, 7, 11, 13, 17, 31)
MODULI = (2, 3, 4, 5, 7, 8, 16)

def residue_features(x, N):
    g = math.gcd(x, N)
    j = int(gmpy2.jacobi(x, N))
    return dict(zero=x == 0, nonunit=g>1, split=1<g<N,
                **{f'jac{v}': j == v for v in (-1, 0, 1)})

def features(N, B):
    m = N//4
    pre = {f'bin{s}':central(m, s)==0 for s in SMALL}
    pre.update({f'root{d}': 2*(m % d)>=d for d in (2,4,8,16,32)})
    splits = {}
    for P,Q in ((3,1),(1,-1),(4,1)):
        for j in (1,2):
            x = lucas(P,Q,j*m,N)
            name = f'U{P}_{Q}_{j}'
            pre.update({f'{name}.{k}':v for k,v in residue_features(x,N).items()})
            g = math.gcd(x,N)
            if 1 < g < N:
                splits[name] = g
    post = residue_features(B, N)
    g = math.gcd(B,N)
    if 1<g<N:
        splits['B'] = g
    post['upper'] = 2*B>N
    for d in MODULI:
        for r in range(d):
            post[f'mod{d}={r}'] = B % d == r
    return pre, post, splits

def prime_draw(rng, bits, c):
    while True:
        s = rng.randrange(1<<(bits-3), 1<<(bits-2))*4 + c
        if fmpz(s).is_prime():
            return s

def samples(bits):
    rng = random.Random(2026112600+bits)
    seen = set()
    for c in (1,3):
        i = 0
        while i < 128:
            p,q = sorted((prime_draw(rng,bits//2,c),prime_draw(rng,bits//2,c)))
            N = p*q
            if p == q or N in seen or N.bit_length() != bits:
                continue
            seen.add(N)
            yield bits,i,N,p,q,c
            i += 1

def compute(data):
    bits,i,N,p,q,c = data
    start = time.perf_counter()
    B, rp, rq = exact(N)
    elapsed = time.perf_counter()-start
    assert (p,q) == (rp,rq) and (B == 0) == (c == 3)
    a = b = None
    if c == 1:
        a,b = gauss_a(p),gauss_a(q)
        for s,t,aa in ((p,q,a),(q,p,b)):
            assert central((s-1)//4,s) == (2*aa)%s
            assert B % s == 2*aa*central((t-1)//4,s)%s
    pre,post,splits = features(N,B)
    return dict(bits=bits,index=i,N=N,p=p,q=q,c=c,B=B,ap=a,aq=b,
                pre=pre,post=post,splits=splits,baseline_s=elapsed)

if __name__ == '__main__':
    pars = argparse.ArgumentParser()
    pars.add_argument('--workers', type=int, default=4)
    args = pars.parse_args()
    primes = [p for p in range(3,500,2) if fmpz(p).is_prime()]
    checks = 0
    for p in primes:
        for q in primes:
            N = p*q
            if p<q and N<1500 and p%4==q%4:
                assert exact(N)[0] == math.comb(N//2,N//4)%N
                checks += 1
    print('Exact small pairs',checks,flush=True)
    with multiprocessing.Pool(args.workers) as pool:
        for bits in (32,48,64):
            with (HERE/f'rows{bits}.jsonl').open('w') as out:
                for row in pool.imap(compute, samples(bits),chunksize=2):
                    out.write(json.dumps(row)+'\n')
            print('Finished',bits,flush=True)
