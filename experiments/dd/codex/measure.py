import argparse
import csv
import ctypes as ct
import itertools
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path('/Volumes/git/factoring-connection')
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'experiments/subthird/codex'), str(ROOT)]
import experiment as ex
from flint import acb_poly, arb, ctx, fmpq_mat
from sympy import integer_nthroot

BITS = (60,72,84,96,112)
DEGREE = dict(zip(BITS,(7,8,9,10,12)))
T, ONE = ex.T, ex.one
GRID = 2**56
LIB = ct.CDLL(str(HERE/'nearest.so'))
I64PTR, INTPTR = ct.POINTER(ct.c_int64), ct.POINTER(ct.c_int)
LIB.create.argtypes = [ct.c_int,I64PTR]
LIB.create.restype = ct.c_void_p
LIB.solve.argtypes = [ct.c_void_p,I64PTR,INTPTR]
LIB.solve.restype = None
LIB.destroy.argtypes = [ct.c_void_p]
LIB.destroy.restype = None
ctx.prec = 384

def signs():
    for u in itertools.count(1):
        yield u
        yield -u

def poles(m,p):
    q,r=divmod(p,m-2)
    rr=[]
    for i in range(m-2):
        d=q+(i<r)
        for c in signs():
            candidate=T**d+c
            if all(candidate.gcd(v)==ONE for v in rr):
                rr.append(candidate)
                break
    return rr

def select(A,m,p,n):
    rr=poles(m,p)
    den=math.prod(rr)
    for b in range(int(integer_nthroot(A,n)[0]),1,-1):
        rem=A-b**n
        digits=[]
        for _ in range(n-1):
            d=(rem+b//2)%b-b//2
            digits.append(d)
            rem=(rem-d)//b
        base=ex.Poly(digits+[rem,1])
        for u in [0,*itertools.islice(signs(),32)]:
            P=base+u*(T-b)
            if not ex.irreducible(P*P-2):
                continue
            for j in itertools.islice(signs(),16):
                Q=P+j*(T-b)*den
                if not ex.irreducible(Q*Q-2):
                    continue
                assert P(b)==Q(b)==A and P*P!=Q*Q
                gens=[[ONE,T]+[X//R for R in rr] for X in (P,Q)]
                assert all(v(b)==w(b) for v,w in zip(*gens))
                for X,gs in zip((P,Q),gens):
                    cols=ex.expand(X,X*X-2,gs)
                    assert len(ex.pivots([ex.coordinate_columns(v,2*n) for v in cols]))==2*m
                return b,rr,P,Q,gens
    raise ValueError('no registered pair')

def normal(rows,m):
    # Two rational equations per k-row for its annihilator. Clear rational denominators.
    if not rows:
        return [1,0]+[0]*(2*m-2)
    eq=[]
    for row in rows:
        eq.append([v for a,b in zip(row[::2],row[1::2]) for v in (a,2*b)])
        eq.append([v for a,b in zip(row[::2],row[1::2]) for v in (b,a)])
    red,rank=fmpq_mat(eq).rref()
    assert rank==2*len(rows)
    piv=[next(j for j in range(2*m) if red[i,j]) for i in range(rank)]
    free=next(j for j in range(2*m) if j not in piv)
    den=math.lcm(*(int(red[i,free].q) for i in range(rank)))
    c=[0]*(2*m)
    c[free]=den
    for i,j in enumerate(piv):
        c[j]=int(-red[i,free]*den)
    assert all(sum(x*y for x,y in zip(r,c))==0 for r in eq)
    return c

def gridpoint(v):
    k=int((v.mid()+arb(1)/2).floor().unique_fmpz())
    assert abs(v-k)<1
    return k

def fixed_threshold(positions,m):
    assert len(positions)==4*m and max(map(abs,positions)) < 2**57
    ar=ct.c_int64*(4*m)
    result=(ct.c_int*(2*m))()
    ptr=LIB.create(m,ar(*positions))
    rows=[]
    sizes=[]
    try:
        while len(rows)<m:
            c=normal(rows,m)
            # Multiplication by sqrt(2) takes (a,b) to (2b,a).
            tags=[v for a,b in zip(c[::2],c[1::2]) for v in (a,b,2*b,a)]
            assert sum(map(abs,tags))<2**62
            LIB.solve(ptr,ar(*tags),result)
            w=tuple(result)
            assert any(sum(x*y for x,y in zip(w,tags[k::2])) for k in (0,1))
            rows.append(w)
            sizes.append(sum(sum(v*a for v,a in zip(positions[k::2],w))**2 for k in (0,1)))
        assert ex.kranks(rows)==m
    finally:
        LIB.destroy(ptr)
    return max(sizes)

def endpoints(v,places):
    scale=10**places
    return (int((v*scale).lower().floor().unique_fmpz()),
            int((v*scale).upper().ceil().unique_fmpz()))

def measure_root(z,s,gs):
    vals=[acb_poly(ex.coeffs(g))(z) for g in gs]
    m=len(vals)
    M=arb(1)
    for v in vals:
        M=M.max(abs(v))
    positions=[gridpoint(u) for v in vals for w in (v,v*s) for u in ((w/M*GRID).real,(w/M*GRID).imag)]
    best=fixed_threshold(positions,m)
    # Each real/imag component is within one grid unit, hence all rows within 2m sqrt(2).
    tgrid=arb(best).sqrt()
    tn=(tgrid+arb(0,2*m*arb(2).sqrt()))/GRID
    loss=-tn.log()
    ratio=loss/(m*arb(m).log())
    lo,hi=endpoints(loss,9)
    rl,rh=endpoints(ratio,9)
    assert hi-lo<1000 and rh-rl<1000
    return lo,hi,rl,rh

def cell(arg):
    bits,N,m,p=arg
    n=max(DEGREE[bits],p+2)
    b,rr,P,Q,gens=select(math.isqrt(N+2),m,p,n)
    manifest=dict(bits=bits,m=m,p=p,n=n,b=b,R=list(map(ex.coeffs,rr)),P=ex.coeffs(P),Q=ex.coeffs(Q))
    data=[]
    for side,X,gs in zip('PQ',(P,Q),gens):
        pol=acb_poly(ex.coeffs(X))
        roots=(X*X-2).complex_roots()
        assert len(roots)==2*n and all(mult==1 for z,mult in roots)
        cache={}
        for i,(z,_) in enumerate(roots):
            value=pol(z)
            assert value.real>0 or value.real<0
            sign=1 if value.real>0 else -1
            s=sign*arb(2).sqrt()
            assert value.overlaps(s)
            kind='real' if z.imag==0 else 'complex'
            if i not in cache:
                val=measure_root(z,s,gs)
                if kind=='complex':
                    j,=[j for j,(w,_) in enumerate(roots) if z.conjugate().overlaps(w)]
                    assert j!=i
                    cache[j]=val
            else:
                val=cache[i]
            # Rounded root coordinates identify each isolated ball; root index is FLINT's output index.
            x,y=endpoints(z.real,12)[0],endpoints(z.imag,12)[0]
            data.append((bits,m,p,n,side,i,kind,sign,x,y,*val))
    return manifest,data

def fixtures():
    out={}
    for bits in BITS:
        N=ex.fixture(bits,19)
        fac=ex.pari.factorint(N)
        primes=[int(a) for a in fac[0]]
        assert len(primes)==2 and math.prod(primes)==N and N.bit_length()==bits
        assert all(q%8 in (1,7) and ex.pari.isprime(q) for q in primes)
        A=math.isqrt(N+2)
        assert A*A-2==N
        out[bits]=dict(N=N,A=A,primes=primes)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--workers',type=int,default=3)
    args=ap.parse_args()
    fx=fixtures()
    (HERE/'fixtures.json').write_text(json.dumps(fx,indent=1)+'\n')
    jobs=[(bits,fx[bits]['N'],m,p) for m in range(3,9) for p in range(m,2*m+1) for bits in BITS]
    with ProcessPoolExecutor(args.workers) as pool, (HERE/'cells.jsonl').open('w') as cf, (HERE/'embeddings.csv').open('w') as ef:
        writer=csv.writer(ef)
        writer.writerow('bits m p n side root kind sqrt_sign floor_re_1e12 floor_im_1e12 loss_lo_1e9 loss_hi_1e9 ratio_lo_1e9 ratio_hi_1e9'.split())
        for manifest,rows in pool.map(cell,jobs):
            cf.write(json.dumps(manifest)+'\n')
            writer.writerows(rows)
            ef.flush()
            cf.flush()
            print(manifest['bits'],manifest['m'],manifest['p'],len(rows),flush=True)

if __name__=='__main__':
    main()
