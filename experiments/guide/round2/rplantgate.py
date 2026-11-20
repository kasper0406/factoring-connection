"""Apply the registered validation and microcell gates to each native plant as a power diagnostic."""
import csv
import sys
from pathlib import Path
import numpy as np
from scipy.special import expit,logit
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

from rbase import Data,SEED,YS,clip
from rfeatures import load_global
from rcontrols import model,intercept,measured
from rstats import data,write,paired_indices
from rfit import NAMES,threshold


def work(l):
    d=Data(l);p=np.load(f'models/A_{l}_pred.npy',mmap_mode='r')
    t,g,doms,anchors=load_global(d,p); phase=d.role[t];b=p[t];s=d.base[t,:4];ok=d.ok[t]
    dest=Path(f'results/A_{l}');bundles={split:{h:[] for h in (0,2)} for split in ('val','test')}
    for j,name in enumerate(NAMES[:8]):
        tr=np.flatnonzero((phase==0)&doms[j]&ok[:,2]);th=threshold(anchors[j][tr]);rng=np.random.default_rng(SEED+1311+33*j+l)
        qq=expit(logit(clip(b[:,2]))+.35*np.where(anchors[j]>th,1.,-1.));yy=(rng.random(len(t))<qq).astype('f')
        other=np.flatnonzero(phase>0);ph=model(s,g[j],yy,b[:,2],tr,other,j==7)
        for code,split in [(1,'val'),(2,'test')]:
            select=phase==code; z=data(dest/f'{split}_{name}.npz')
            for h,lab,oracle,fit in [(0,z['plant_y'],z['plant_oracle'],z['plant_p']),(2,yy[select],qq[select],ph[phase[other]==code])]:
                eligible=doms[j][select]&ok[select,h];bp=b[select,h]
                interceptshift=float(np.mean(logit(clip(intercept(bp[eligible],oracle[eligible])))-logit(clip(bp[eligible]))))
                bundles[split][h].append((lab,fit,expit(logit(clip(bp))+interceptshift),eligible))
    rows=[]
    for split in ('val','test'):
        common=data(dest/f'{split}_data.npz'); ia,ib=paired_indices(common,[v[3] for v in bundles[split][0]])
        for h,parts in bundles[split].items():
            terms=[];counts=[];lower=[]
            for yy,pr,bp,use in parts:
                l1=logit(clip(pr))-logit(clip(bp));x=l1[ia]-l1[ib]
                disc=use[ia]&(yy[ia]!=yy[ib])&use[ib]
                gradient=np.where(disc,x*(yy[ia]-expit(logit(clip(bp[ia]))-logit(clip(bp[ib])))),0)
                terms.append(gradient);counts.append(int(disc.sum()))
                _,_,bt=measured(yy[use],bp[use],pr[use],common['n'][use],common['block'][use])
                lower.append(float(np.quantile(bt[:,0]/np.maximum(bt[:,1],1e-12),.005)))
            mat=np.stack(terms,axis=1);rng=np.random.default_rng(SEED+804)
            stats=np.vstack([mat.sum(0),rng.choice([-1.,1.],(1000,len(ia)))@mat])
            pv=np.column_stack([rankdata(-v,method='max')/len(v) for v in stats.T]);ord=np.argsort(pv[0]); adj=np.zeros(8);prev=0.
            for j,k in enumerate(ord):
                prev=max(prev,np.mean(np.min(pv[:,ord[j:]],1)<=pv[0,k]));adj[k]=prev
            for k,name in enumerate(NAMES[:8]):
                rows.append(dict(ycut=YS[l],family=name,head='joint' if h==0 else 'algebraic',split=split,
                              cell_p=pv[0,k],val_or_test_WY=adj[k],discordant_pairs=counts[k],per_hit_lo=lower[k]))
    write(f'results/native_gates_{l}.csv',rows)


if __name__=='__main__':
    with threadpool_limits(limits=1):work(int(sys.argv[1]))
