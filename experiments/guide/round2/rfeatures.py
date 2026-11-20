"""Global feature allowlist. Own prime factors and own cofactor digits never enter these designs."""
import concurrent.futures
import functools
import json
import math
import sys
import time
from pathlib import Path
import numpy as np
from scipy.special import softmax
from rdata import YS,PR,SEED,ULTRA,cf,js
import gmpy2


def cf_float(x):
    seq=[]
    for _ in range(28):
        k=math.floor(x); seq.append(k)
        if x-k<1e-13: break
        x=1/(x-k)
    return seq


def prefix(x,y):
    for j,(a,b) in enumerate(zip(x,y)):
        if a!=b: return j
    return min(len(x),len(y))


def coord_sieve(xy,y):
    a,b=xy[:,:2].T
    numbers=[abs(a),b,abs(a+b),abs(a-b)]
    lim=max(int(t.max()) for t in numbers)+1
    prod=np.ones(lim,dtype=np.int64); rad=prod.copy()
    for p in PR[PR<=min(y,lim-1)]:
        p=int(p); rad[p::p]*=p; k=p
        while k<lim: prod[k::k]*=p; k*=p
    def fraction(x):
        return np.log(prod[x])/np.log(np.maximum(2,x))
    return np.column_stack([fraction(t) for t in numbers]+[(prod[abs(a)]==abs(a)),(prod[b]==b),
             (np.log(rad[abs(a)])+np.log(rad[b]))/np.maximum(1,np.log(prod[abs(a)])+np.log(prod[b]))])


def naf(n):
    n=abs(n); cnt=0
    while n:
        if n&1: n-=2-(n%4);cnt+=1
        n//=2
    return cnt


def digital(n):
    n=abs(n); s=bin(n)[2:]; seen=set(); buf=''
    for ch in s:
        buf+=ch
        if buf not in seen: seen.add(buf);buf=''
    t=n; tern=0
    while t: t,r=divmod(t,3);tern+=int(r!=0)
    return [n.bit_count(),tern,naf(n),len(seen)+bool(buf)]


def geometric(c,m,xy):
    a,b=xy[:,:2].T; x=a/b; n=len(x); deg=len(c)-1
    roots=sorted(np.roots(c[::-1]),key=lambda v:(round(v.real,10),v.imag))
    real=[r.real for r in roots if abs(r.imag)<1e-7]
    rr=real+[float(m)]
    seqs=[cf_float(v) for v in rr]
    g2=np.zeros((n,15),dtype=float); lb=np.log(np.maximum(b,2))
    h0=np.array([math.log2(1+1/(j*(j+2))) for j in range(1,6)])
    h0=np.r_[h0,math.log2(1+1/6)-math.log2(1+1/11),math.log2(1+1/11)]
    for j,(aa,bb) in enumerate(zip(a.tolist(),b.tolist())):
        ds=cf(aa,bb); ds1=ds[1:] or [1]
        hh=np.bincount([min(t-1,4) if t<=5 else (5 if t<=10 else 6) for t in ds1],minlength=7)/len(ds1)
        mask=hh>0
        g2[j,:4]=[len(ds)/lb[j],math.log(max(ds1))/lb[j],sum(map(math.log,ds1))/len(ds1),np.sum(hh[mask]*np.log(hh[mask]/h0[mask]))]
        for i,seq in enumerate(seqs): g2[j,4+i]=prefix(ds,seq)
    for i,r in enumerate(rr): g2[:,9+i]=-np.log(np.maximum(1e-30,abs(x-r)))/lb
    g2[:,-1]=-np.log(np.maximum(1e-15,np.min([abs(x-np.round(x*q)/q) for q in range(1,17)],axis=0)))
    logs=np.column_stack([np.log(np.maximum(1e-100,np.abs(a-r*b))) for r in roots])
    logs-=logs.mean(1,keepdims=True)
    emb=np.zeros((n,4)); emb[:,:deg]=logs
    chamber=np.zeros((n,4)); chamber[:,:len(real)]=np.column_stack([np.sign(x-r) for r in real]) if real else np.empty((n,0))
    sm=softmax(-logs,axis=1); poly=np.polynomial.Polynomial(c); f=poly(x); fp=poly.deriv()(x); fpp=poly.deriv(2)(x)
    hes=(deg-1)*(deg*f*fpp-(deg-1)*fp*fp)
    g3=np.column_stack((emb,np.linalg.norm(logs,axis=1),logs.min(1),logs.max(1),-(sm*np.log(sm)).sum(1),chamber,
                        np.log1p(abs(fp)/np.maximum(abs(f),1e-100)),np.log(np.maximum(abs(hes),1e-100))-2*(deg-2)/deg*np.log(np.maximum(abs(f),1e-100))))
    uvals=np.unique(np.r_[abs(a),b]); dig={int(v):digital(int(v)) for v in uvals}
    bits=np.arange(26,dtype=np.int64)
    mask=(1<<max(96,int(m).bit_length()+int(b.max()).bit_length()+2))-1
    def carry(i,j):
        i&=mask; j&=mask
        return i.bit_count()+j.bit_count()-(i+j).bit_count()
    g8=np.column_stack(((abs(a[:,None])>>bits)&1,(b[:,None]>>bits)&1,a<0,
             [dig[int(abs(v))] for v in a],[dig[int(v)] for v in b],
             [(carry(int(t),int(v)),carry(int(t),-m*int(v))) for t,v in zip(a,b)]))
    assert max(max(abs(a)),max(b))<2**len(bits)
    iv={int(v):pow(int(v),-1,ULTRA) for v in np.unique(b)}
    cos=[]
    for aa,bb in zip(a.tolist(),b.tolist()):
        r=aa*iv[bb]%ULTRA
        cos.append([int(gmpy2.legendre(r,ULTRA)),((r^(r>>43)^(r>>83))%1021)/1021])
    return [g2.astype('f'),g3.astype('f'),g8.astype('f'),np.array(cos,dtype='f')]


def geomfile(k):
    it=json.loads(Path(f'data/n{k:02}.json').read_text())
    with np.load(f'data/n{k:02}.npz') as z: xy=z['coords'][:it['nt']]
    g2,g3,g8,cos=geometric(it['c'],it['m'],xy)
    np.savez_compressed(f'data/g{k:02}.npz',a0=coord_sieve(xy,YS[0]).astype('f'),a1=coord_sieve(xy,YS[1]).astype('f'),g2=g2,g3=g3,g8=g8,cos=cos)
    return k


def load_global(d,pred,tag='A',null_labels=None):
    out=[[] for _ in range(8)]; dom=[[] for _ in range(8)]
    tgt=[]; anchors=[[] for _ in range(8)]
    for it,full,ts in zip(d.info,d.idx,d.tidx):
        k=it['id']; use=d.target[ts]; t=ts[use]; tgt.append(t)
        with np.load(f'data/n{k:02}.npz') as z:
            ed=z['edges'][use]; digit=z['digits'][:,d.l]
        with np.load(f'data/g{k:02}.npz') as z:
            g1,g2,g3,g8,cos=[z[key][use] for key in [f'a{d.l}','g2','g3','g8','cos']]
        if tag in ('B','Bnull'):
            # Paid references for the within-N split must themselves belong to the same buffered partition.
            reg=d.region[full]
            same=(ed>=0)&(reg[np.maximum(ed,0)]==d.region[t,None])&(d.region[t,None]>=0)
            for start in (4,6):
                same[:,start:start+2]&=same[:,start:start+2].all(1,keepdims=True)
            ed=np.where(same,ed,-1)
        ys=(d.resp if null_labels is None else null_labels)[full]
        res=ys-pred[full]; valid=ed>=0
        edges=np.maximum(ed,0); features=res[edges]*valid[:,:,None]
        def avg(cols):
            return features[:,cols].sum(1)/np.maximum(1,valid[:,cols].sum(1)[:,None])
        g4=np.column_stack((avg(slice(0,4)),features[:,:4].reshape(len(t),-1),valid[:,:4]))
        g5=np.column_stack((avg(slice(4,8)),features[:,4:8].reshape(len(t),-1),valid[:,4:8]))
        bs=d.bid[t]; qmask=bs<32
        seed=np.array([r['seed'] for r in it['parents']]+[.5]*16)[bs]
        wt=valid[:,8:]*np.array([1,.65,.65**2,.65**3])
        wt/=np.maximum(1e-9,wt.sum(1,keepdims=True))
        ew=(features[:,8:]*wt[:,:,None]).sum(1)
        dg=(digit[edges[:,8:]]*wt[:,:,None]).sum(1)
        g7=np.column_stack((ew,dg,cos,valid[:,8:].sum(1)))
        allm=np.ones(len(t),bool)
        for j,(z,mask,a) in enumerate(zip([g1,g2,g3,g4,g5,seed[:,None],g7,g8],
             [allm,allm,allm,valid[:,:4].any(1),valid[:,4:8].any(1),qmask,allm,allm],
             [g1[:,0]+g1[:,1]+g1[:,2]+g1[:,3],g2[:,0],g3[:,4],g4[:,0],g5[:,0],seed,ew[:,0],g8[:,1]])):
            out[j].append(z.astype('f'));dom[j].append(mask);anchors[j].append(a.astype('f'))
    return np.concatenate(tgt),[np.concatenate(v) for v in out],[np.concatenate(v) for v in dom],[np.concatenate(v) for v in anchors]


if __name__=='__main__':
    with concurrent.futures.ProcessPoolExecutor(4) as pool:
        for k in pool.map(geomfile,range(48)): print(k,flush=True)
