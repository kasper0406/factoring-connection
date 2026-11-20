"""L1 positive fixture and hit-normalized empirical scaling of the intentionally leaky PC2a."""
import numpy as np
from scipy.special import logit
from threadpoolctl import threadpool_limits
from rbase import Data,YS,B,SEED
from rdata import gm,js
from rcontrols import model,measured
from rstats import write


def run(l):
    d=Data(l); t=np.flatnonzero(d.target);probe=np.random.default_rng(SEED).choice(t,2048,replace=False)
    D=gm.primorial(2*YS[l])//gm.primorial(YS[l]); forbidden=gm.primorial(B)//gm.primorial(YS[l])
    assert gm.gcd(forbidden,D)==D and D>1
    for i in probe:
        a,b=map(int,d.coord[i,:2]);m=d.info[int(d.n[i])]['m'];g=a-m*b
        assert gm.gcd(g+D,D)==gm.gcd(g,D)
        assert (gm.gcd(g+D,D)>1)==bool(d.leak[i,0])
    p=np.load(f'models/A_{l}_pred.npy',mmap_mode='r')[t]; rows=[]
    for h in (0,1,2):
        y=d.resp[t,h];s=d.base[t,:4];tr=np.flatnonzero((d.role[t]==0)&d.ok[t,h]);te=np.flatnonzero((d.role[t]==2)&d.ok[t,h])
        pred=model(s,d.leak[t].astype('f'),y,p[:,h],tr,te)
        u=s[te,:2]/np.log(B);U=(u*np.log(np.maximum(u,1e-7))).sum(1) if h==0 else u[:,h-1]*np.log(u[:,h-1])
        bins=np.digitize(U,np.quantile(U,[1/3,2/3]))
        for b in range(3):
            ix=bins==b; ii=te[ix]
            r,mu,bt=measured(y[ii],p[ii,h],pred[ix],d.n[t[ii]],d.bid[t[ii]],y[ii]*U[ix])
            rows.append(dict(ycut=YS[l],head=h,bin=b,**r,mean_U=U[ix].mean(),I_hit=mu[0]/max(1e-12,mu[1]),theta=mu[0]/max(1e-12,mu[2])))
    return rows


if __name__=='__main__':
    with threadpool_limits(limits=1): rows=run(0)+run(1)
    write('results/pc2a_scaling.csv',rows);js('results/L1_positive_fixture.json',dict(cutoffs=list(YS),checks_per_cutoff=2048,rejected_parallel_shift=True))
