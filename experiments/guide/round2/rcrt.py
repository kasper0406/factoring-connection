"""Recompute the CRT regional null with the final buffered angular and q partition."""
import numpy as np
from threadpoolctl import threadpool_limits
from rbase import Data,YS,PR,SEED,project
from rcontrols import model,measured
from rstats import write

def run(l):
    d=Data(l);t=np.flatnonzero(d.target);n=d.n[t];region=d.region[t];keep=d.role[t]==0
    pr=PR[(PR>YS[l])&(PR<=2*YS[l])][:160];labels=np.zeros(len(t),'f')
    raw=[]
    for it,ts in zip(d.info,d.tidx):
        with np.load("data/g%02d.npz"%it['id']) as z:raw.append(z['g8'][d.target[ts]])
    for k in np.unique(n[keep]):
        for reg in (0,1):
            ix=np.flatnonzero((n==k)&(region==reg));code,_=project(d.coord[t[ix]],pr,int(pr[-1]))
            rng=np.random.default_rng(SEED+29011+int(k)*23+reg)
            roots=np.array([rng.integers(int(pp)+1) for pp in pr]);labels[ix]=(code==roots).any(1)
    p=np.full(len(t),1-np.prod(1-1/(pr.astype(float)+1)));tr=np.flatnonzero(keep&(region==0));te=np.flatnonzero(keep&(region==1))
    new=model(d.base[t,:4],np.concatenate(raw),labels,p,tr,te,True)
    r,_,_=measured(labels[te],p[te],new,n[te],d.bid[t[te]])
    write(f'results/crt_{l}.csv',[dict(ycut=YS[l],**r)])

if __name__=='__main__':
    with threadpool_limits(limits=1):
        for l in [0,1]:run(l)
