"""All within-cell permutations and exact-stratum conditional-logit checks on frozen scores."""
import json
from pathlib import Path

import numpy as np
from scipy.special import logit
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from rbase import clip,SEED,YS
from rfit import NAMES,HEADS
from rstats import data,write


def sets(keys):
    order=np.lexsort(tuple(keys.T));ordered=keys[order]
    chunks=np.split(order,np.flatnonzero(np.any(ordered[1:]!=ordered[:-1],axis=1))+1)
    return [g for g in chunks if len(g)>1]


def exact_cl(groups,yy,eta,delta):
    informative=[g for g in groups if 0<yy[g].sum()<len(g)]
    if not informative:return None,0,0
    def gh(beta):
        gradient,information=0.,0.
        for g in informative:
            y=yy[g];d=delta[g];v=eta[g]+beta*d; m=int(y.sum())
            if m>len(g)//2:y=1-y;d=-d;v=-v;m=len(g)-m
            weights=np.exp(np.clip(v-v.max(),-200,0))
            z=np.zeros(m+1);z[0]=1;z1=np.zeros(m+1);z2=np.zeros(m+1)
            for w,x in zip(weights,d):
                z2[1:]+=w*(z2[:-1]+2*x*z1[:-1]+x*x*z[:-1])
                z1[1:]+=w*(z1[:-1]+x*z[:-1]);z[1:]+=w*z[:-1]
            if z[-1]==0:continue
            u=z1[-1]/z[-1];gradient+=y@d-u;information+=max(0,z2[-1]/z[-1]-u*u)
        return gradient,information
    beta=0.
    for _ in range(12):
        grad,inf=gh(beta)
        if inf<1e-12:break
        step=np.clip(grad/inf,-5,5);beta=np.clip(beta+step,-40,40)
        if abs(step)<1e-5:break
    return float(beta),len(informative),sum(len(g) for g in informative)


def full(l,split):
    p=Path(f'results/A_{l}');d=data(p/f'{split}_data.npz');fs=[data(p/f'{split}_{nm}.npz') for nm in NAMES]
    basekey=np.column_stack((d['n'],d['sign'],np.floor(d['s'][:,:2]/.02))).astype(int)
    dom=np.sum(np.column_stack([f['domain'] for f in fs])*2**np.arange(9),axis=1)
    flags=np.sum(d['ok'][:,:3]*2**np.arange(3),axis=1)
    groups=sets(np.column_stack((basekey,dom,flags)))
    groups.sort(key=len);ix=np.concatenate(groups); lens=np.array(list(map(len,groups)));starts=np.r_[0,lens.cumsum()[:-1]]
    run=[(int(starts[np.flatnonzero(lens==k)[0]]),int(np.sum(lens==k)),int(k)) for k in np.unique(lens)]
    res=[];deltas=[];entries=[]
    # Conditional logit uses each family's full eligible domain, without unnecessary common-support refinement.
    cache={}
    for nm,f in zip(NAMES,fs):
        for h,head in enumerate(HEADS):
            sel=f['domain']&d['ok'][:,h];y=d['y'][:,h];bp=d['p'][:,h];fp=f['pred'][:,h]
            delta=logit(clip(fp))-logit(clip(bp)) if h<3 else fp-bp
            D=delta[ix].astype(float);D-=np.repeat(np.add.reduceat(D,starts)/lens,lens)
            deltas.append(D);res.append((y[ix]-bp[ix])*sel[ix])
            beta=None;ncell=None;nrow=None
            if h<3:
                key=sel.tobytes()
                if key not in cache:
                    eligible=np.flatnonzero(sel)
                    cache[key]=[eligible[g] for g in sets(basekey[eligible])]
                beta,ncell,nrow=exact_cl(cache[key],y.astype(float),logit(clip(bp)).astype(float),delta.astype(float))
            entries.append(dict(ycut=YS[l],split=split,family=nm,head=head,full_conditional_beta=beta,informative_cells=ncell,informative_rows=nrow))
    D=np.stack(deltas,axis=1);R=np.stack(res,axis=1);stats=[(D*R).sum(axis=0)]
    rng=np.random.default_rng(SEED+1194);orig=np.arange(len(ix))
    for _ in range(1000):
        perm=orig.copy()
        for st,n,k in run:
            order=rng.random((n,k)).argsort(axis=1)
            perm[st:st+n*k]=(st+np.arange(n)[:,None]*k+order).ravel()
        stats.append((D[perm]*R).sum(axis=0))
    stats=np.asarray(stats)
    pv=np.column_stack([rankdata(-col,method='max')/len(col) for col in stats.T]);order=np.argsort(pv[0]);prev=0.
    for j,h in enumerate(order):
        prev=max(prev,np.mean(pv[:,order[j:]].min(1)<=pv[0,h]))
        entries[h].update(full_perm_p=pv[0,h],full_WY=prev)
    return entries


if __name__=='__main__':
    with threadpool_limits(limits=1):out=sum([full(l,s) for l in [0,1] for s in ['val','test']],[])
    write('results/all_cell_permutations.csv',out)
