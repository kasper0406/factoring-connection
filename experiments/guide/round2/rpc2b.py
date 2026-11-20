"""Matched exact-root arithmetic control on the square locus of a derived irreducible quartic."""
import concurrent.futures
import json
import math
import random
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from rdata import B,SEED,YS,PR,SMALL,PI,ev,gm,PARI,prime_alpha,factors,js,polynomial,GRID
from rbase import B1,impose,saddle_grid
from rfeatures import geometric
from rcontrols import model,measured
from rstats import block_arrays,hierarchical,write,holm
sys.path.insert(0,"/Volumes/git/factoring-connection")
from src.cadobench.smoothcount import smooth_part


POWER_PC=os.environ.get('ROUND2_PC2_POWER')=='1'
SUF='_power' if POWER_PC else ''
NFORMS=4 if POWER_PC else 1

def generate(index):
    k=24+index//NFORMS; fid=index%NFORMS
    it=json.loads(Path(f'data/n{k:02}.json').read_text()); c=it['c']; m=it['m']; rng=random.Random(SEED+909*k)
    if POWER_PC:
        nrng=random.Random(97328119+1103*k)
        P,Q,_,_=polynomial(4,GRID[k%24],nrng,PARI)
        rng=random.Random(14900121+8777*k+2521*fid)
        while True:
            c=[rng.randrange(-9,10) for _ in range(4)]+[1]
            if c[0]==0 or (c[1]==c[3]==0):continue
            f=PARI.Polrev(c)
            if not f.polisirreducible():continue
            rp=PARI.polrootsmod(f,P);rq=PARI.polrootsmod(f,Q)
            if len(rp)==0 or len(rq)==0:continue
            m=int(rp[0])+P*((int(rq[0])-int(rp[0]))*pow(P,-1,Q)%Q)
            break
        it=dict(N=str(P*Q))
    E=[c[0],c[2],1]; O=[c[1],c[3]]
    co=[int(v) for v in np.polynomial.polynomial.polymul(np.array(E,dtype=object),np.array(E,dtype=object))]
    for j,v in enumerate(np.polynomial.polynomial.polymul(np.array(O,dtype=object),np.array(O,dtype=object))):co[j+1]-=int(v)
    assert PARI.Polrev(co).polisirreducible() and ev(co,m*m,1)%int(it['N'])==0
    xy=[]; seen=set(); arm=[]
    while len(xy)<(10000 if POWER_PC else 20000):
        a,b=(rng.randrange(2,116),rng.randrange(2,116)) if POWER_PC else (rng.randrange(2,360),rng.randrange(2,240))
        if math.gcd(a,b)!=1 or (a*a,b*b) in seen:continue
        X,Y=a*a,b*b; assert ev(co,X,Y)==ev(c,a,b)*ev(c,-a,b)
        xa,ya=X+rng.randrange(-9,10),Y+rng.randrange(-9,10)
        if xa<=0 or ya<=0 or (xa,ya) in seen or math.gcd(xa,ya)!=1 or (math.isqrt(xa)**2==xa and math.isqrt(ya)**2==ya):continue
        seen.update([(X,Y),(xa,ya)]); xy.extend([(X,Y),(xa,ya)]);arm.extend([1,0])
    alpha,_=prime_alpha(co); xy=np.array(xy,dtype=np.int64)
    g2,g3,g8,_=geometric(co,m*m,xy); n=len(xy)
    raw=np.zeros((n,2));sizes=np.zeros((n,2,2));local=sizes.copy();w=sizes.copy();flag=np.zeros_like(w,dtype='b')
    signs=np.zeros((n,2),np.int8);y=np.empty((n,2),np.uint8)
    localkey=np.empty((n,2),dtype='U50'); vis=[[] for _ in range(n)];ves=[[] for _ in range(n)]
    for s in (0,1):
        ns=[int(X)-m*m*int(Y) if s==0 else ev(co,int(X),int(Y)) for X,Y in xy]
        signs[:,s]=np.sign(np.array(ns,dtype=object)).astype(np.int8);ns=list(map(abs,ns))
        sm=list(map(int,smooth_part(ns,gm.primorial(B)))); y[:,s]=[x==v for x,v in zip(ns,sm)]
        raw[:,s]=[math.log(x) for x in ns]
        for l,Yc in enumerate(YS):
            sy=list(map(int,smooth_part(ns,gm.primorial(Yc))))
            if s==1:localkey[:,l]=[str(v) for v in sy]
            for j,(x,t,d) in enumerate(zip(ns,sy,sm)):
                C=x//t;ss=math.log(C);sizes[j,l,s]=ss;local[j,l,s]=math.log(t)
                w[j,l,s]=math.log(d//t)/ss if ss else 1.
                flag[j,l,s]=1 if C<=B else -1 if C<(2003 if l==0 else 32771)**2 else 0
                if l==1:
                    for prime,e in factors(t):vis[j].append(s*len(SMALL)+PI[prime]);ves[j].append(e)
    ptr=np.r_[0,np.cumsum([len(v) for v in vis])]; vi=np.concatenate([np.asarray(v,'i') for v in vis])
    np.savez_compressed(f'data/pc2b{SUF}_{index if POWER_PC else k}.npz',N=str(it['N']),f=c,m=str(m),F=co,coords=xy,raw=raw,sizes=sizes,local=local,w=w,flag=flag,sign=signs,y=y,
          arm=np.array(arm,dtype='b'),keys=localkey,alpha=alpha,ptr=ptr,vi=vi,ve=np.concatenate([np.asarray(v,'B') for v in ves]),
          g2=g2,g3=g3,g8=g8)


def main(l):
    parts={k:[] for k in ('base','coord','ok','resp','fl','p0','n','arm','g2','g3','g8','pair')}
    vals=[];grid,pv=saddle_grid(YS[l])
    for index in range(24*NFORMS):
        k=24+index//NFORMS
        src=Path(f'data/pc2b{SUF}_{index}.npz') if POWER_PC else Path(f'data/pc2b_{k}.npz')
        z=np.load(src);s=z['sizes'][:,l];loc=z['local'][:,l];raw=z['raw'];flags=z['flag'][:,l];sign=z['sign'];arm=z['arm']; keys=z['keys'][:,l]
        # Maximal within-cell matching without replacement, random deterministic order within each cell.
        pool={};rng=np.random.default_rng(SEED+23*k)
        for j in rng.permutation(len(arm)):
            if flags[j,1]!=0:continue
            key=(keys[j],int(sign[j,1]),int(np.floor(s[j,1]/.2)))
            pool.setdefault(key,[[],[]])[arm[j]].append(j)
        pairs=[(a,b) for grp in pool.values() for a,b in zip(*grp)]
        ix=np.array(pairs).reshape(-1)
        assert np.all(keys[ix[::2]]==keys[ix[1::2]]) and np.all(abs(s[ix[::2],1]-s[ix[1::2],1])<=.2)
        nf=len(s);alpha=np.broadcast_to(z['alpha'][l],(nf,2));pa=np.column_stack([np.interp(s[:,i]+alpha[:,i],grid,pv) for i in range(2)])
        p0=impose(np.column_stack((pa.prod(1),pa,pa)),flags)[:,:3]; yy=z['y']
        allv=[np.column_stack((s,loc,raw,sign,flags,alpha,s[:,0]*s[:,1],p0)),z['coords'],
              np.column_stack(((flags==0).all(1),flags==0,flags==0)),np.column_stack((yy.prod(1),yy,z['w'][:,l])),
              flags,p0,np.full(nf,k),arm,z['g2'],z['g3'],z['g8']]
        for key,v in zip(list(parts)[:-1],allv):parts[key].append(v[ix])
        parts['pair'].append(np.repeat(np.arange(len(pairs)),2))
        vals.append(sparse.csr_matrix((z['ve'],z['vi'],z['ptr']),shape=(nf,2*len(SMALL)))[ix]);z.close()
    d=SimpleNamespace(**{k:np.concatenate(v) for k,v in parts.items()});d.val=sparse.vstack(vals,format='csr');d.cut=YS[l];d.flags=d.fl
    info=json.loads(Path('data/manifest.json').read_text());role=np.array([v['split'] for v in info])[d.n]
    fit=role=='fit';test=role=='test';p=np.zeros(len(role))
    for fold in (0,1,2):
        train=fit&((d.n%2!=fold) if fold<2 else True)
        ids=np.unique(d.n[train]);calids=ids[1::3];tr=np.flatnonzero(train&~np.isin(d.n,calids));ca=np.flatnonzero(train&np.isin(d.n,calids))
        # Keep a side-only PC2b design using its identical outcome for ancillary heads where deterministic.
        base=B1().fit(d,tr,ca,d.resp)
        ix=np.flatnonzero((fit&(d.n%2==fold)) if fold<2 else ~fit)
        p[ix]=base.predict(d,ix)[:,2]
    tr=np.flatnonzero(fit);te=np.flatnonzero(test); y=d.resp[:,2];rows=[];scaling=[]
    for name,x in [('G2',d.g2),('G3',d.g3),('G8',d.g8),('oracle_square_locus',d.arm[:,None])]:
        new=model(d.base[:,:4],x.astype('f'),y,p,tr,te,name=='G8')
        r,_,_=measured(y[te],p[te],new,d.n[te],d.pair[te]%16);rows.append(dict(ycut=YS[l],model=name,**r))
        u=d.base[te,1]/np.log(B);bins=np.digitize(u,np.quantile(u,[1/3,2/3]))
        for bn in range(3):
            mask=bins==bn;ii=te[mask]
            rr,mu,_=measured(y[ii],p[ii],new[mask],d.n[ii],d.pair[ii]%16)
            sr=y[ii][d.arm[ii]==1];orate=y[ii][d.arm[ii]==0]
            scaling.append(dict(ycut=YS[l],model=name,bin=bn,mean_u=float(u[mask].mean()),split_hits=int(sr.sum()),
             ordinary_hits=int(orate.sum()),split_n=len(sr),ordinary_n=len(orate),gain=rr['gain'],theoretical_u_log2=float(u[mask].mean()*np.log(2))))
    write(f'results/pc2b{SUF}_scaling_{l}.csv',scaling)
    a=te[::2];b=te[1::2];assert np.all(d.arm[a]==0) and np.all(d.arm[b]==1)
    v=np.column_stack((y[b]-y[a],(y[b]-p[b])-(y[a]-p[a]),y[a],y[b]))
    av,bt,pt=hierarchical(block_arrays(v,d.n[b],d.pair[b]%16))
    desc=dict(ycut=YS[l],exact_local_matched_pairs=len(b),ordinary_rate=av[2],split_rate=av[3],raw_RD=av[0],
        raw_CI=np.quantile(bt[:,0],[.005,.995]).tolist(),adjusted_RD=av[1],adjusted_CI=np.quantile(bt[:,1],[.005,.995]).tolist())
    return rows,desc


if __name__=='__main__':
    if len(sys.argv)<2:
        with concurrent.futures.ProcessPoolExecutor(4) as pool:list(pool.map(generate,range(24*NFORMS)))
    with threadpool_limits(limits=1):
        o=[main(l) for l in (0,1)]
    rows=o[0][0]+o[1][0]
    for r,p in zip(rows,holm([r['p'] for r in rows])):r['Holm']=p
    write(f'results/pc2b{SUF}.csv',rows);js(f'results/pc2b{SUF}_match.json',[v[1] for v in o])
