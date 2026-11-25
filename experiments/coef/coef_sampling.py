"""Use the factoring-connection virtualenv. All outputs live next to this script."""
import json
import math
import random
import sys
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

import gmpy2
import numpy as np
from cypari2 import Pari

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from src.ecount.zn_curve import ec_mul, Split

WIDTHS = (40,60,80)
ELLS = (2,3)
def strip_l(l):
    m = math.lcm(*range(1,44))
    while m % l == 0:
        m //= l
    return m

def sample(N, draws):
    read = [0,0]
    witness = [False,False]
    try:
        for k in draws:
            P = ec_mul(k, (0,1), -1, N)
            for col,l in enumerate(ELLS):
                Q = ec_mul(strip_l(l), P, -1, N)
                T = ec_mul(l,Q,-1,N)
                if Q is not None:
                    witness[col] |= T is None
                    v = Q[0]+1
                    g = math.gcd(v,N)
                    if 1 < g < N:
                        raise Split(g,"readout")
                    read[col] += int(gmpy2.jacobi(v,N))
    except Split as e:
        assert 1 < e.factor < N and N % e.factor == 0
        return None
    return tuple((v>0)-(v<0) for v in read),witness

def prime(bits,rng,used,pari):
    while True:
        x = rng.getrandbits(bits) | (1<<(bits-1)) | 1
        if x not in used and gmpy2.is_prime(x) and pari.isprime(x):
            return x

def census(l):
    mats = []
    for a,b,c,d in product(range(l),repeat=4):
        det = (a*d-b*c)%l
        if det:
            mats.append((det,(1+det-a-d)%l))
    return {str(D):[sum(1 for d,v in mats if d==D and v==out) for out in range(l)]
            for D in range(1,l)}, {str(D):[sum(1 for d,v in mats for e,w in mats
                        if d*e%l==D and v*w%l==out) for out in range(l)] for D in range(1,l)}

def main():
    prime_rng = random.Random(252323)
    scalar_rng = random.Random(317221)
    pari = Pari()
    E = pari.ellinit([0,0,0,-1,1])
    splits = Counter()
    data = []
    for size in WIDTHS:
        used = set()
        for i in range(800):
            while True:
                p,q = (prime(size//2,prime_rng,used,pari) for _ in range(2))
                N = p*q
                if p!=q and N.bit_length()==size:
                    used.update((p,q))
                    break
            draws = [scalar_rng.randrange(1,2**(2*size)+1) for _ in range(6)]
            res = sample(N,draws)
            if res is None:
                splits[size] += 1
                continue
            M = (p+1-int(pari.ellap(E,p)))*(q+1-int(pari.ellap(E,q)))
            read,wit=res
            assert all(not w or M%l==0 for l,w in zip(ELLS,wit))
            data.append(((size,N%3,int(gmpy2.jacobi(-23,N))),i<400,read,wit,M))
    output = {"splits":splits,"census":{l:census(l) for l in ELLS},"targets":{}}
    for pos,l in enumerate(ELLS):
        base = defaultdict(lambda:np.ones(l))
        bins = defaultdict(lambda:np.zeros(l))
        for key,train,read,wit,M in data:
            if train:
                base[key][M%l] += 1
                bins[(key,read[pos])][M%l] += 1
        out = []
        for key,train,read,wit,M in data:
            if train:
                continue
            b=base[key]/base[key].sum()
            a=bins[(key,read[pos])] + 20*b
            a=a/a.sum()
            if wit[pos]:
                a=np.zeros(l)
                a[0]=1
            out.append((key[0], math.log(a[M%l]/b[M%l]),int(b.argmax()==M%l),int(a.argmax()==M%l)))
        A=np.array(out)
        per=[A[A[:,0]==k,1] for k in WIDTHS]
        rng=np.random.default_rng(91337)
        draws=np.zeros(4000)
        for sub in per:
            draws += rng.choice(sub,(4000,len(sub))).sum(axis=1)/len(A)
        output["targets"][l]=dict(clean_test=[len(v) for v in per],gains=[v.mean() for v in per],
              overall=A[:,1].mean(),ci99=np.quantile(draws,[.005,.995]).tolist(),
              acc=A[:,2:].mean(axis=0).tolist(),witnesses=sum(wit[pos] for _,_,_,wit,_ in data))
    txt=json.dumps(output,indent=2)
    Path(__file__).with_name("coef_sampling_result.json").write_text(txt+"\n")
    print(txt)
if __name__ == '__main__':
    main()
