"""Exact arithmetic, admissible relation graph and global target ownership."""
import collections
import concurrent.futures
import functools
import hashlib
import json
import math
import os
import random
import site
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
site.addsitedir("/Volumes/git/factoring-connection/.venv/lib/python3.12/site-packages")
import gmpy2 as gm
import sympy as sp
from cypari2 import Pari
sys.path.append(str(Path(__file__).resolve().parent.parent))
from generate import polynomial

B=262139
YS=(2000,32768)
SEED=7328119
PR=np.array(list(sp.primerange(2,B+1)),dtype=np.int32)
SMALL=PR[PR<=YS[-1]]
PI={int(p):j for j,p in enumerate(SMALL)}
FORBID=gm.primorial(B)//gm.primorial(YS[0])
ULTRA=2**127-1
PARI=Pari()
SPLITS=["fit"]*8+["val"]*6+["test"]*10
GRID=[60,68,76,84,92,98,104,110,62,72,82,92,102,110,60,66,72,78,84,90,96,100,106,110]


def js(path,v):
    Path(path).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")


def ev(c,a,b):
    return sum(t*a**j*b**(len(c)-1-j) for j,t in enumerate(c))


def factors(n):
    v=PARI(abs(n)).factor()
    return [(int(p),int(e)) for p,e in zip(v[0],v[1])]


@functools.lru_cache(120000)
def safe(n):
    n=abs(int(n))
    return n!=0 and (n<=2000 or gm.gcd(n,FORBID)==1)


def is_donor(a,b):
    z=(a*0x45d9f3b+b*0x27d4eb2d)&0xffffffff
    z=((z^(z>>16))*0x45d9f3b)&0xffffffff
    return ((z^(z>>16))&3)==0


def reduce_basis(q,r):
    u,v=(q,0),(r,1)
    norm=lambda x: x[0]**2+x[1]**2
    while True:
        if norm(v)<norm(u): u,v=v,u
        k=round((u[0]*v[0]+u[1]*v[1])/norm(u))
        if not k: return u,v
        v=(v[0]-k*u[0],v[1]-k*u[1])


def parents(c,m,rng):
    pools=[[[],[]] for _ in range(4)]
    used=set()
    attempts=0
    while min(len(l) for v in pools for l in v)<12:
        a,b=rng.randrange(-160,161),rng.randrange(1,161)
        attempts+=1
        if attempts>1500000: raise ValueError("parent pool exhausted")
        if math.gcd(a,b)!=1: continue
        F=abs(ev(c,a,b))
        fs=factors(F)
        for q,e in fs:
            if not B<q<=8*B or q in used: continue
            used.add(q)
            r=a*pow(b,-1,q)%q
            seeded=int(e==1 and all(p<=B or p==q for p,_ in fs) and all(p<=B for p,_ in factors(a-m*b)))
            band=min(3,int(math.log(q/B,8)*4))
            if len(pools[band][seeded])>=48: continue
            u,v=reduce_basis(q,r)
            signatures=[math.log(q),math.log(math.hypot(*v)/math.hypot(*u))]
            for si in (0,1):
                vals=[]
                for i,j in [(1,0),(0,1),(1,1),(1,-1)]:
                    A=i*u[0]+j*v[0]; D=i*u[1]+j*v[1]
                    norm=abs(A-m*D) if si==0 else abs(ev(c,A,D))//q
                    vals.append(math.log(max(norm,1)))
                signatures.append(np.median(vals))
            pools[band][seeded].append(dict(q=q,r=r,parent=[a,b],seed=seeded,band=band,
                                             basis=[u,v],signature=signatures))
    chosen=[]
    for band,(controls,seeds) in enumerate(pools):
        x=np.array([r["signature"] for r in seeds]); y=np.array([r["signature"] for r in controls])
        cost=np.sum(((x[:,None]-y)/np.array([.22,.50,1.,2.]))**2,axis=2)
        # Select four distinct minimum-cost pairs per stratum from the candidate inventory.
        for j in range(4):
            i,h=np.unravel_index(np.argmin(cost),cost.shape)
            for rec in [seeds[i],controls[h]]:
                rec.update(pair=4*band+j,match_cost=float(cost[i,h]))
                chosen.append(rec)
            cost[i,:]=np.inf; cost[:,h]=np.inf
    rng.shuffle(chosen)
    return chosen,attempts


def cf(a,b):
    ds=[]
    while b:
        k,r=divmod(a,b); ds.append(k); a,b=b,r
    return ds


def convergent(seq):
    h0,h,k0,k=0,1,1,0
    for d in seq: h0,h=h,d*h+h0; k0,k=k,d*k+k0
    return h,k


def canonical(a,b):
    return (-a,-b) if b<0 else (a,b)


def prime_alpha(c):
    """Uniform affine p-adic integral, with exact recursive treatment of multiple roots."""
    pol=PARI.Polrev(c)
    disc=abs(int(pol.poldisc()))
    arr=np.zeros((2,2))
    counts=[]
    def integral(co,p):
        content=0
        while all(x%p==0 for x in co):
            co=[x//p for x in co]; content+=1
        roots=PARI.polrootsmod(PARI.Polrev(co),p)
        ans=float(content)
        for rt in roots:
            r=int(rt)
            if sum(j*v*r**(j-1) for j,v in enumerate(co) if j)%p:
                ans+=1/(p-1)
            else:
                new=[sum(co[j]*math.comb(j,k)*r**(j-k)*p**k for j in range(k,len(co))) for k in range(len(co))]
                ans+=integral(new,p)/p
        return ans
    for pv in PR:
        p=int(pv)
        nr=len(PARI.polrootsmod(pol,p))
        counts.append(nr)
        if p<=YS[0]: continue
        integ=nr/(p-1) if disc%p else integral(c,p)
        da=[math.log(p)*(1/(p-1)-p/((p+1)*(p-1))),
            math.log(p)*(1/(p-1)-p*integ/(p+1))]
        for k,y in enumerate(YS):
            if p>y: arr[k]+=da
    return arr,np.array(counts,dtype=np.uint8)


def make(k):
    t=time.monotonic()
    rng=random.Random(SEED+17027*k)
    deg=3+k//24
    p,q,m,c=polynomial(deg,GRID[k%24],rng,PARI)
    latt,attempts=parents(c,m,rng)
    forbidden_parents={tuple(r["parent"]) for r in latt}
    coords=[]; blocks=[]; used=set(); priors=[]
    for h,r in enumerate(latt):
        u,v=r["basis"]; qp=r["q"]; A,D=r["parent"]
        history=[]
        n=0
        while n<352 or len(history)<12:
            i,j=rng.randrange(-42,43),rng.randrange(-42,43)
            a,b=canonical(i*u[0]+j*v[0],i*u[1]+j*v[1])
            if b<=0 or math.gcd(a,b)!=1 or (a,b) in forbidden_parents or (a,b) in used: continue
            if not safe(a*D-A*b): continue
            if is_donor(a,b):
                if len(history)<12: history.append((a,b,qp)); used.add((a,b))
            elif n<352:
                coords.append((a,b,qp)); blocks.append(h); used.add((a,b)); n+=1
        priors.append(history)
    n=0
    while n<11264:
        R,H=(4800,2100) if n%2==0 else (100000,48000)
        a,b=rng.randrange(-R,R+1),rng.randrange(1,H+1)
        if (a,b) in used or (a,b) in forbidden_parents or is_donor(a,b) or math.gcd(a,b)!=1 or not a: continue
        coords.append((a,b,1)); blocks.append(32+n%16); used.add((a,b)); n+=1
    nt=len(coords)
    reg={v:i for i,v in enumerate(coords)}
    def node(a,b,q=1):
        key=(a,b,q)
        if key not in reg: reg[key]=len(coords); coords.append(key)
        return reg[key]
    # Pre-approve fixed projective maps on both norm sides via integer resultants.
    poly=PARI.Polrev(c); x=PARI('x'); maps=[]; transforms=[]
    for f,h in [(2,1),(3,1),(-2,1),(1,2),(1,3),(2,3),(3,2),(-1,2),(-2,3)]:
        for typ in [0,1,-1]:
            A,L,D=(f,0,h) if typ==0 else (h,typ*f,h)
            transformed=sum(v*(A*x+L)**j*D**(deg-j) for j,v in enumerate(c))
            res1=int(poly.polresultant(transformed)); res0=A*m+L-D*m
            if safe(res1) and safe(res0): maps.append((typ,f,h,res0,res1))
    for A,L,C,D in [(0,1,1,0),(-1,0,0,1),(1,1,0,1),(0,-1,1,0)]:
        transformed=sum(v*(A*x+L)**j*(C*x+D)**(deg-j) for j,v in enumerate(c))
        rr=int(poly.polresultant(transformed)); rg=A*m+L-m*(C*m+D)
        transforms.append([A,L,C,D,rg,rr])
    shifts=[(i,j) for i,j in [(1,0),(-1,0),(0,1),(1,1),(-1,1),(2,1),(1,2),(-2,1),(2,3)]
            if safe(i-m*j) and safe(ev(c,i,j))]
    # columns = shifts, Farey, admissible SL2 image, product pair, sum/difference pair, four earlier q donors.
    edges=np.full((nt,12),-1,dtype=np.int32)
    for j,(a,b,qp) in enumerate(coords[:nt]):
        def usable(A,D):
            return D>0 and (A,D) not in forbidden_parents and is_donor(A,D) and math.gcd(A,D)==1 and safe(a*D-A*b)
        sh=[]
        for u,v in shifts:
            if usable(a+u,b+v): sh.append(node(a+u,b+v))
            if len(sh)==2: break
        edges[j,:len(sh)]=sh
        A,D=convergent(cf(a,b)[:-1])
        if usable(A,D): edges[j,2]=node(A,D)
        for A,L,C,D,rg,rr in transforms:
            if not (safe(rg) and safe(rr)): continue
            aa,bb=canonical(A*a+L*b,C*a+D*b)
            if usable(aa,bb): edges[j,3]=node(aa,bb); break
        for typ,f,h,_,_ in maps:
            col=4 if typ==0 else 6
            if edges[j,col]>=0 or b%h or not safe(a*h-f*b): continue
            bb=b//h
            if typ==0:
                if a%f: continue
                aa=a//f
            else:
                if (a-typ*f*bb)%h: continue
                aa=(a-typ*f*bb)//h
            if usable(aa,bb) and (f,h) not in forbidden_parents and (f,h) not in used:
                edges[j,col:col+2]=[node(aa,bb),node(f,h)]
        stop=min(blocks[j],32); count=0
        for prev in reversed(range(stop)):
            for A,D,Q in priors[prev]:
                if safe(a*D-A*b):
                    edges[j,8+count]=node(A,D,Q); count+=1; break
            if count==4: break
    nr=len(coords)
    sizes=np.zeros((nr,2,2),dtype=float); w=np.zeros_like(sizes); local=np.zeros_like(sizes)
    raw=np.zeros((nr,2)); signs=np.zeros((nr,2),np.int8); labels=np.zeros((nr,2),np.uint8)
    flags=np.zeros((nr,2,2),np.int8); leak=np.zeros_like(flags); digits=np.zeros((nr,2,6),np.float32)
    vi=[]; ve=[]; ptr=[0]
    maxnorm=0
    for j,(a,b,qp) in enumerate(coords):
        F=ev(c,a,b); assert F%qp==0
        for s,v in enumerate([a-m*b,F//qp]):
            n=abs(v); assert 0<n<ULTRA
            maxnorm=max(maxnorm,n* (qp if s else 1)); signs[j,s]=1 if v>0 else -1
            raw[j,s]=math.log(n); fac=factors(n)
            assert math.prod(t**e for t,e in fac)==n
            labels[j,s]=int(all(t<=B for t,e in fac))
            for t,e in fac:
                if t<=YS[-1]: vi.append(s*len(SMALL)+PI[t]); ve.append(e)
            for l,y in enumerate(YS):
                S=math.prod(t**e for t,e in fac if t<=y); C=n//S
                sm=math.prod(t**e for t,e in fac if y<t<=B)
                sizes[j,l,s]=math.log(C); local[j,l,s]=math.log(S); w[j,l,s]=math.log(sm)/math.log(C) if C>1 else 1.
                flags[j,l,s]=1 if C<=B else (-1 if C<(2003 if l==0 else 32771)**2 else 0)
                assert not flags[j,l,s] or labels[j,s]==int(flags[j,l,s]>0)
                leak[j,l,s]=any(y<t<=2*y for t,e in fac)
                bits=bin(C)[2:]; runs=[len(z) for z in __import__('re').findall(r'0+|1+',bits)]
                counts=collections.Counter(runs)
                entropy=-sum((f/len(runs))*math.log(f/len(runs)) for f in counts.values())
                digits[j,l,3*s:3*s+3]=[C.bit_count()/len(bits),entropy,(C/2**(len(bits)-1))-1]
        ptr.append(len(vi))
    assert maxnorm<ULTRA
    alpha,rootcounts=prime_alpha(c)
    path=Path(f"data/n{k:02}.npz")
    np.savez_compressed(path,coords=np.array(coords,dtype=np.int64),blocks=np.array(blocks,dtype=np.int16),
        edges=edges,sizes=sizes,local=local,raw=raw,signs=signs,y=labels,w=w,flags=flags,leak=leak,digits=digits,
        ptr=np.array(ptr,np.int32),vi=np.array(vi,np.int16),ve=np.array(ve,np.uint8),alpha=alpha,rootcounts=rootcounts)
    info=dict(id=k,N=str(p*q),factors=[str(p),str(q)],degree=deg,m=m,c=c,split=SPLITS[k%24],bits=GRID[k%24],
              parents=latt,parent_search=attempts,shifts=shifts,maps=maps,transforms=transforms,nt=nt,nr=nr,
              maxnorm=str(maxnorm),seconds=time.monotonic()-t)
    js(f"data/n{k:02}.json",info)
    return k


if __name__=="__main__":
    start=time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as pool:
        for k in pool.map(make,range(48)): print(k,flush=True)
    infos=[json.loads(Path(f"data/n{k:02}.json").read_text()) for k in range(48)]
    old={v["N"] for f in ["../data/manifest.json","../power_cohort/data/manifest.json"] for v in json.loads(Path(f).read_text())}
    assert len({v["N"] for v in infos})==48 and not old.intersection(v["N"] for v in infos)
    seen=set()
    for it in sorted(infos,key=lambda v:(SPLITS.index(v["split"]),v["id"])):
        k=it["id"]
        with np.load(f"data/n{k:02}.npz") as z: xy=z["coords"][:it["nt"],:2]
        keys=((xy[:,0]+2**30)<<32)+xy[:,1]
        mask=np.array([int(key) not in seen for key in keys]); seen.update(map(int,keys))
        np.save(f"data/keep{k:02}.npy",mask)
    js("data/manifest.json",infos); js("logs/generation.json",dict(seconds=time.monotonic()-start,unique=len(seen)))
