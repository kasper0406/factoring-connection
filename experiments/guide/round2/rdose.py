"""Descriptive native-feature power curve on the same physical candidates."""
import sys
from pathlib import Path
import numpy as np
from scipy.special import expit,logit
from threadpoolctl import threadpool_limits
from rbase import Data,SEED,YS,clip
from rfeatures import load_global
from rfit import Project,Correction,SmallNet,corrected,threshold,NAMES
from rcontrols import intercept,measured,pow_reps
from rstats import write,holm


def run(l):
    d=Data(l);p=np.load(f'models/A_{l}_pred.npy',mmap_mode='r')
    t,glob,doms,anchors=load_global(d,p);s=d.base[t,:4];role=d.role[t];rows=[]; stored={}
    for j,name in enumerate(NAMES[:8]):
        for h in (0,2):
            tr=np.flatnonzero((role==0)&doms[j]&d.ok[t,h]);te=np.flatnonzero((role==2)&doms[j]&d.ok[t,h])
            ti=np.random.default_rng(SEED+11).choice(tr,min(100000,len(tr)),replace=False)
            prj=Project().fit(s,glob[j],ti)
            xt=prj.transform(s[ti],glob[j][ti]);xe=prj.transform(s[te],glob[j][te])
            anchor=anchors[j];th=threshold(anchor[tr]); ss=np.where(anchor>th,1.,-1.)
            pp=p[t,h];uniform=np.random.default_rng(SEED+7381+169*j+11*h+l).random(len(t))
            for odds in (4,8):
                q=expit(logit(clip(pp))+.5*np.log(odds)*ss);yy=(uniform<q).astype('f')
                ok=np.ones((len(ti),1),bool)
                mod=Correction().fit(xt,yy[ti,None],pp[ti,None],ok,[True])
                net=SmallNet().fit(xt,yy[ti,None],pp[ti,None],ok,[True]) if j==7 else None
                ph=corrected(xe,mod,net,pp[te,None],[True]).ravel()
                base=intercept(pp[te],q[te]);r,_,_=measured(yy[te],base,ph,d.n[t[te]],d.bid[t[te]])
                pv,lo=pow_reps(q[te],base,ph,d.n[t[te]],d.bid[t[te]],SEED+413*j+73*h+l)
                key=f'{name}_{h}_{odds}';stored['p_'+key]=pv;stored['l_'+key]=lo
                rows.append(dict(ycut=YS[l],family=name,head=h,odds=odds,**r))
    write(f'results/dose_{l}.csv',rows);np.savez(f'results/dose_reps_{l}.npz',**stored)


def aggregate():
    import csv
    rows=sum([list(csv.DictReader(open(f'results/dose_{l}.csv'))) for l in range(2)],[])
    for h in ('0','2'):
        for odds in ('4','8'):
            sub=[r for r in rows if r['head']==h and r['odds']==odds]; ps=[];ls=[]
            for r in sub:
                l=YS.index(int(r['ycut']));key=f"{r['family']}_{h}_{odds}"
                z=np.load(f'results/dose_reps_{l}.npz');ps.append(z['p_'+key]);ls.append(z['l_'+key])
            adj=np.stack([holm(v) for v in np.array(ps).T]).T;pw=((adj<.01)&(np.array(ls)>0)).mean(1)
            for r,a,p in zip(sub,holm([float(r['p']) for r in sub]),pw):r.update(Holm=a,detected=bool(a<.01 and float(r['lo'])>0),conditional_test_power=p)
    write('results/native_dose_response.csv',rows)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        if sys.argv[1]=='aggregate':aggregate()
        else:run(int(sys.argv[1]))
