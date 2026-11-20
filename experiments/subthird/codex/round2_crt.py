"""Registered exploratory split-slope decoder. Run with the repo's .venv/bin/python -B."""
import json
import random
from math import ceil,gcd,isqrt,prod
from pathlib import Path
from cypari2 import Pari
from sympy import nextprime,primerange,sqrt_mod
from flint import fmpz

pari=Pari()
def ceil_fourth(n):
    x=isqrt(isqrt(n))
    return x+(x**4<n)

def fixture(bits,rng):
    lo=3*(1<<(bits//2))//4
    hi=(1<<(bits//2))-1000
    p,q=(int(nextprime(rng.randrange(lo,hi))) for _ in range(2))
    assert p!=q and (p*q).bit_length()==bits
    return p*q

def form(N,primes,rng):
    generators=[]
    for p in primes:
        if p!=2 and pow(-N%p,(p-1)//2,p)==1:
            w=int(sqrt_mod(-N,p))
            generators.append(pari.Qfb(p,2*w,(w*w+N)//p))
        if len(generators)>=32: break
    f=pari.Qfb(1,0,N)
    for step in range(1000):
        f=f.qfbcomp(rng.choice(generators))
        if step>=20 and 16*int(f[0])**2>=N and int(f[0])**2<=4*N:
            a,b,c=map(int,(f[0],f[1],f[2]))
            assert b*b-4*a*c==-4*N
            return a,b,c
    raise ValueError("desired shape not reached")

def shortest(res,mod):
    a,b = (mod,0),(res,1)
    while True:
        aa=sum(x*x for x in a)
        bb=sum(x*x for x in b)
        if aa<bb: a,b=b,a;aa,bb=bb,aa
        inner=a[0]*b[0]+a[1]*b[1]
        if 2*abs(inner)<=bb:
            return b if b[1]>=0 else (-b[0],-b[1])
        t=(2*inner+bb)//(2*bb)
        a=(a[0]-t*b[0],a[1]-t*b[1])

def is_smooth(v,P):
    h=gcd(v,int(P % v))
    while h>1:
        v//=h
        h=gcd(v,h)
    return v==1

def ground(a,b,c,K,P):
    hits=0
    total=0
    for v in range(1,K+1):
        vals=[a*u*u+b*u*v+c*v*v for u in range(1,K+1) if gcd(u,v)==1]
        total+=len(vals)
        # A product tree shares the expensive reductions of the primorial.
        levels=[[fmpz(x) for x in vals]]
        while len(levels[-1])>1:
            prev=levels[-1]
            levels.append([prev[i]*prev[i+1] if i+1<len(prev) else prev[i] for i in range(0,len(prev),2)])
        residues=[P%levels[-1][0]]
        for sub in reversed(levels[:-1]):
            residues=[residues[i//2]%v for i,v in enumerate(sub)]
        hits+=sum(is_smooth(v,h) for v,h in zip(vals,residues))
    return total,hits

def run(bits,seed):
    rng=random.Random(8841703*bits+seed)
    N=fixture(bits,rng)
    lg=(bits-1).bit_length()
    B=1<<ceil_fourth(bits*lg**3)
    K=1
    while (16*(K.bit_length()-1))**4<3**4*bits**3*lg: K*=2
    assert K*K<=2_000_000
    primes=list(primerange(B+1))
    P=fmpz(prod(primes))
    a,b,c=form(N,primes,rng)
    def hit(u,v): return is_smooth(a*u*u+b*u*v+c*v*v,P)
    roots=[]
    for p in primes:
        if p==2 or a%p==0 or N%p==0 or pow((-N)%p,(p-1)//2,p)!=1: continue
        h=int(sqrt_mod(-N,p))
        inv=pow(a,-1,p)
        roots.append((p,((-b//2+h)*inv%p,(-b//2-h)*inv%p)))
    seen=set()
    target=2*K*K
    lattice_ops=0
    terminals=0
    for _ in range(512):
        batch=rng.sample(roots,min(30,len(roots)))
        mod=1
        branches=[0]
        for p,zz in batch:
            inv=pow(mod,-1,p)
            nmod=mod*p
            ext=[]
            for z in zz:
                for r in branches:
                    x=r+mod*((z-r)*inv%p)
                    u,v=shortest(x,nmod)
                    assert (u-x*v)%nmod==0
                    ext.append((max(abs(u),abs(v)),x,u,v))
            lattice_ops+=len(ext)
            best=sorted(ext)[:4]
            branches=[x for _,x,_,_ in best]
            mod=nmod
            if mod>target:
                terminals+=len(best)
                for _,_,u,v in best:
                    if 1<=u<=K and 1<=v<=K and gcd(u,v)==1: seen.add((u,v))
                break
    decoded=sum(hit(u,v) for u,v in sorted(seen))
    control=set()
    while len(control)<len(seen):
        u,v=rng.randint(1,K),rng.randint(1,K)
        if gcd(u,v)==1: control.add((u,v))
    direct=sum(hit(u,v) for u,v in control)
    total,count=ground(a,b,c,K,P)
    return dict(bits=bits,seed=seed,N=N,form=[a,b,c],B=B,K=K,batches=512,reductions=lattice_ops,
                terminal_branches=terminals,distinct_tested=len(seen),decoded=decoded,control=direct,
                primitive_points=total,smooth_points=count)

if __name__=='__main__':
    dst=Path(__file__).with_name('round2_crt.jsonl')
    with dst.open('w') as out:
        for bits in (48,64,80,96):
            for seed in range(4):
                r=run(bits,seed)
                print(json.dumps(r),file=out,flush=True)
                print(bits,seed,r['decoded'],r['control'],flush=True)
