import ctypes as ct
import json
import math
import random
import time
from pathlib import Path

from fpylll import CVP, IntegerMatrix, LLL

SEED = 20261126
PRIMES = [i for i in range(2, 10000) if all(i % d for d in range(2,math.isqrt(i)+1))]
ODD = PRIMES[1:]

class Result(ct.Structure):
    _fields_ = [('a', ct.c_int), ('b', ct.c_int), ('prefix', ct.c_int)]

search = ct.CDLL(str(Path('wheel.dylib').resolve())).exact
search.argtypes = [ct.POINTER(ct.c_int), ct.c_int, ct.c_uint32, ct.c_int, ct.POINTER(Result)]
search.restype = None

def controls(P, rng):
    k = len(P)
    bits = math.prod(P).bit_length()
    while True:
        R = P[:]
        for _ in range(100):
            i, new = rng.randrange(k), rng.choice(ODD[:3*k])
            if new not in R and (math.prod(R)//R[i]*new).bit_length() == bits:
                R[i] = new
        if sorted(R) != P: break
    m = (k+1)//2
    g = 2**((bits-0.5)/m)
    pool = [p for p in ODD if g*.5 <= p <= g*2]
    while True:
        S = rng.sample(pool,m)
        if math.prod(S).bit_length() == bits: break
    return sorted(R), sorted(S)

def signature(x, wheel):
    sig = 0
    for i,p in enumerate(wheel):
        r = pow(x,(p-1)//2,p)
        if r == 0: return None
        sig |= int(r != 1)<<i
    return sig

def kernel(n):
    d, result = 2, 1
    while d*d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                n //= d
                continue
            result *= d
        d = 3 if d == 2 else d+2
    return result*n

def exhaustive(wheel,t,rat):
    p = (ct.c_int*len(wheel))(*wheel)
    r = Result()
    start = time.process_time()
    search(p, len(wheel), t, rat, ct.byref(r))
    a,b = r.a,r.b
    h = kernel(a)*kernel(b) if a and rat else a
    if a and rat:
        h //= math.gcd(kernel(a),kernel(b))**2
        d = math.gcd(a,b)
        a,b = a//d,b//d
    return dict(a=a,b=b,h=h,prefix=r.prefix,secs=time.process_time()-start)

def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit

def lattice(wheel, target, seed, sqrt_weights=False):
    start = time.process_time()
    rng = random.Random(seed)
    k = len(wheel)
    primes = [p for p in PRIMES if p not in wheel][:2*k+8]
    cols = {p:signature(p,wheel) for p in primes}
    best, support = None, []
    rounds = 0
    for j in range(12):
        if j: rng.shuffle(primes)
        n = len(primes)
        if sqrt_weights:
            wt = [round(1024*math.sqrt(math.log2(p))) for p in primes]
        else:
            wt = [round(32*math.log2(p)) for p in primes]
        piv = {}
        basis = []
        def eliminate(v, mask):
            while v:
                pos = v.bit_length()-1
                if pos not in piv: break
                w,e = piv[pos]
                v,mask = v^w,mask^e
            return v,mask
        for i,p in enumerate(primes):
            v,mask = eliminate(cols[p], 1<<i)
            row = [0]*n
            if v:
                piv[v.bit_length()-1] = v,mask
                row[i] = 2*wt[i]
            else:
                for l in bits(mask): row[l] = wt[l]
            basis.append(row)
        rem,mask = eliminate(target,0)
        if rem: break
        v = tuple(wt[i] if (mask>>i)&1 else 0 for i in range(n))
        M = IntegerMatrix.from_matrix(basis)
        LLL.reduction(M,delta=.99)
        near = CVP.babai(M,v)
        for sol in (list(bits(mask)), [i for i in range(n) if ((v[i]-near[i])//wt[i])%2]):
            h = math.prod(primes[i] for i in sol)
            if best is None or h<best:
                best,support=h,[primes[i] for i in sol]
        rounds += 1
        if time.process_time()-start>=.25: break
    a,b=1,1
    for p in sorted(support,reverse=True):
        if a<=b: a *= p
        else: b *= p
    return dict(a=a if best else 0,b=b if best else 0,h=best or 0,
                rounds=rounds,secs=time.process_time()-start)

def run():
    rng = random.Random(SEED)
    with open('trials.jsonl','w',buffering=1) as f:
        for k in [11,18,24,31]:
            P=ODD[:k]
            for j in range(32):
                target = rng.randrange(1<<k)
                R,S = controls(P,rng)
                for typ, wheel in zip('PRS',(P,R,S)):
                    t=target% (1<<len(wheel))
                    base=dict(k=len(wheel),bits=math.prod(wheel).bit_length(),wheel=wheel,
                              typ=typ,trial=j,t=t,overlap=len(set(wheel)&set(P)))
                    for name in ('scan','lll','xor'):
                        if name=='lll': result=lattice(wheel,t,SEED+1000*k+j)
                        else: result=exhaustive(wheel,t,int(name=='xor'))
                        if result['h']:
                            assert signature(result['h'],wheel)==t
                            assert signature(result['a']*result['b'],wheel)==t
                            assert math.gcd(result['a'],result['b'])==1
                        f.write(json.dumps(dict(**base,method=name,**result))+'\n')
            print('completed k',k,flush=True)

if __name__=='__main__': run()
