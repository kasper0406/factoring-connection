"""Locked validation screen and independent, paired held-out inference."""
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import expit,logit
from scipy.stats import rankdata,t as student
from threadpoolctl import threadpool_limits
from rbase import clip,js,YS,SEED
from rfit import HEADS,NAMES


def write(name,rows):
    with open(name,'w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)


def data(path):
    with np.load(path) as v:return dict(v)


def members(a):
    key,idx,cnt=np.unique(a,return_inverse=True,return_counts=True)
    return key,idx,cnt


def block_arrays(values,n,b):
    if values.ndim==1:values=values[:,None]
    key,ix,cnt=members(n.astype(int)*128+b.astype(int))
    mean=np.column_stack([np.bincount(ix,weights=x)/cnt for x in values.T])
    return [mean[key//128==k] for k in np.unique(key//128)]


def hierarchical(groups,seed=SEED):
    rng=np.random.default_rng(seed); nr=1999
    mu=np.stack([r.mean(0) for r in groups])
    # Repeated draws of an N get independent inner resamples of its q/box blocks.
    outer=rng.multinomial(len(groups),np.full(len(groups),1/len(groups)),size=nr)
    boot=np.zeros((nr,mu.shape[1]))
    for j,r in enumerate(groups):
        for copies in np.unique(outer[:,j]):
            if copies==0:continue
            rows=outer[:,j]==copies
            counts=rng.multinomial(int(copies)*len(r),np.full(len(r),1/len(r)),size=rows.sum())
            boot[rows]+=(counts@r)/(len(groups)*len(r))
    mean=mu.mean(0); se=mu.std(0,ddof=1)/np.sqrt(len(mu))
    p=np.ones_like(mean); nz=se>0;p[nz]=student.sf(mean[nz]/se[nz],len(mu)-1);p[(~nz)&(mean>0)]=0
    return mean,boot,p


def holm(p):
    p=np.array(p);o=np.argsort(p);r=p.copy()
    r[o]=np.minimum(1,np.maximum.accumulate(p[o]*np.arange(len(p),0,-1)));return r


def loss(y,p,b):
    p=p.astype(float)
    if b:return -np.where(y,np.log(clip(p)),np.log1p(-clip(p)))
    return (y-p)**2


def paired_indices(d,doms):
    # Random, disjoint pairs inside every common-support, within-N/sign/size cell.
    n=len(d['n']); rng=np.random.default_rng(SEED+17)
    mask=np.zeros(n,np.int32)
    for j,v in enumerate(doms): mask|=v.astype(np.int32)<<j
    for j in range(3):mask|=d['ok'][:,j].astype(np.int32)<<(10+j)
    s=np.floor(d['s'][:,:2]/.02).astype(int)
    keys=np.column_stack((d['n'],d['sign'],s,mask)).astype(int)
    order=np.lexsort(tuple([rng.random(n)]+list(keys.T)))
    k=keys[order]; start=np.r_[True,np.any(k[1:]!=k[:-1],axis=1)]
    offset=np.arange(n)-np.maximum.accumulate(np.where(start,np.arange(n),0))
    left=np.flatnonzero((offset[:-1]%2==0)&~start[1:])
    return order[left],order[left+1]


def cell_tests(path,split,forms,screen=False):
    d=data(path/f'{split}_data.npz'); fs=[data(path/f'{split}_{f}.npz') for f in forms]
    ia,ib=paired_indices(d,[z['domain'] for z in fs])
    terms=[]; keys=[]; info=[]
    for nm,z in zip(forms,fs):
        for h in range(6):
            use=z['domain'][ia]&d['ok'][ia,h]; binary=h<3
            y1=d['y'][ia,h];y2=d['y'][ib,h]
            pa,pb=d['p'][ia,h],d['p'][ib,h]
            if binary:
                dp=logit(clip(pa))-logit(clip(pb))
                c1=logit(clip(z['pred'][ia,h]))-logit(clip(pa));c2=logit(clip(z['pred'][ib,h]))-logit(clip(pb))
                discord=use&(y1!=y2);dx=(c1-c2).astype(float); e=expit(dp)
                # Conditional logit score at zero retains the fixed baseline odds offset.
                score=np.where(discord,dx*(y1-e),0.)
                beta=0.
                for _ in range(12):
                    pr=expit(dp[discord]+beta*dx[discord]); dd=dx[discord]
                    step=np.sum(dd*(y1[discord]-pr))/(1e-8+np.sum(dd*dd*pr*(1-pr)))
                    beta=float(np.clip(beta+np.clip(step,-5,5),-40,40))
                a0=(y1-y2)*(pa-pb); a1=(y1-y2)*(z['pred'][ia,h]-z['pred'][ib,h])
                concord=lambda v: (v>0)+.5*(v==0)
                aucdiff=float(np.mean(concord(a1[discord])-concord(a0[discord]))) if discord.any() else None
                auc_by_n=[np.mean(concord(a1[discord&(d['n'][ia]==n)])-concord(a0[discord&(d['n'][ia]==n)])) for n in np.unique(d['n'][ia][discord])]
                if auc_by_n:aucdiff=float(np.mean(auc_by_n))
            else:
                dx=(z['pred'][ia,h]-pa)-(z['pred'][ib,h]-pb)
                score=np.where(use,dx*((y1-pa)-(y2-pb)),0)
                discord=use;beta=None;aucdiff=None
            keys.append(f'{nm}/{HEADS[h]}');terms.append(score)
            info.append(dict(pairs=int(use.sum()),informative_pairs=int(discord.sum()),conditional_beta=beta,delta_size_stratified_auc=aucdiff))
    scores=np.column_stack(terms).astype(float); dist=[scores.sum(0)]
    rng=np.random.default_rng(SEED+804)
    for _ in range(20):dist.extend(rng.choice([-1.,1.],(50,len(ia)))@scores)
    dist=np.asarray(dist)
    pr=np.column_stack([rankdata(-q,method='max')/len(q) for q in dist.T])
    observed=pr[0]; order=np.argsort(observed);adjusted=observed.copy();prev=0.
    for j,h in enumerate(order):
        prev=max(prev,np.mean(np.min(pr[:,order[j:]],axis=1)<=observed[h]));adjusted[h]=prev
    return {k:dict(**e,cell_p=float(p),WY_p=float(a)) for k,e,p,a in zip(keys,info,observed,adjusted)}


def screen_all():
    out={}
    for l in range(2):
        p=Path(f'results/A_{l}')
        cc=cell_tests(p,'val',NAMES,True)
        for v in cc.values():v['survives']=bool(v['WY_p']<=.025)
        out[str(YS[l])]=cc
    js('results/validation_lock.json',out)
    js('results/lock_provenance.json',dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),prereg=hashlib.sha256(Path('PREREGISTRATION.md').read_bytes()).hexdigest(),
       validation_sha256=hashlib.sha256(Path('results/validation_lock.json').read_bytes()).hexdigest()))


def test(l,tag,forms,screened):
    path=Path(f'results/{tag}_{l}');split='test' if tag in ('A','Afull') else 'region'
    d=data(path/f'{split}_data.npz')
    cells=cell_tests(path,split,forms)
    rows=[]; plants=[]; extras=[]; bins_out=[]
    den=lambda v:(float(v[0]),float(v[1]))
    for name in forms:
        z=data(path/f'{split}_{name}.npz')
        for h,head in enumerate(HEADS):
            use=z['domain']&d['ok'][:,h]
            y=d['y'][use,h]; p=d['p'][use,h]; new=z['pred'][use,h];n=d['n'][use];bb=d['block'][use]
            binary=h<3
            g=loss(y,p,binary)-loss(y,new,binary)
            av,boot,pval=hierarchical(block_arrays(np.column_stack((g,y)),n,bb))
            ci=np.quantile(boot[:,0],[.005,.995])
            norm_ci=np.quantile(boot[:,0]/np.maximum(boot[:,1],1e-15),[.005,.995]) if binary else [None,None]
            cell=cells[f'{name}/{head}'];val=screened.get(f'{name}/{head}',{})
            allix=z['domain']; allg=loss(d['y'][allix,h],d['p'][allix,h],binary)-loss(d['y'][allix,h],z['pred'][allix,h],binary)
            full=np.mean([v.mean() for v in block_arrays(allg,d['n'][allix],d['block'][allix])])
            rows.append(dict(ycut=YS[l],family=name,head=head,rows=int(use.sum()),hits=int(y.sum()) if binary else None,
                gain=float(av[0]),ci_lo=ci[0],ci_hi=ci[1],p=float(pval[0]),per_hit=av[0]/max(av[1],1e-15) if binary else None,
                per_hit_lo=norm_ci[0],per_hit_hi=norm_ci[1],full_gain=full,screened=val.get('survives',False),
                val_WY=val.get('WY_p'),**cell))
            if h==0:
                u=d['s'][use,:2]/np.log(262139);U=np.sum(u*np.log(np.maximum(u,1e-9)),axis=1)
                bins=np.digitize(U,np.quantile(U,[1/3,2/3]))
                for b in range(3):
                    ss=bins==b
                    bins_out.append(dict(ycut=YS[l],family=name,bin=b,U_mean=U[ss].mean(),rows=int(ss.sum()),hits=int(y[ss].sum()),gain=g[ss].mean(),
                          I_over_U=float(g[ss].sum()/max(1e-15,np.sum(y[ss]*U[ss])))))
        use=z['domain']&d['ok'][:,0]; y=z['plant_y'][use];p=d['p'][use,0]
        # Best common offset given known plant probabilities, eliminating a prevalence-only improvement.
        q=z['plant_oracle'][use]; eta=logit(clip(p)); b=0.
        for _ in range(8):
            pr=expit(eta+b);b+=np.sum(q-pr)/max(1e-8,np.sum(pr*(1-pr)))
        g=loss(y,expit(eta+b),True)-loss(y,z['plant_p'][use],True)
        mu,bt,pt=hierarchical(block_arrays(g,d['n'][use],d['block'][use]))
        lo,hi=np.quantile(bt[:,0],[.005,.995])
        plants.append(dict(ycut=YS[l],family=name,gain=mu[0],ci_lo=lo,ci_hi=hi,p=pt[0],hits=int(y.sum()),rows=len(y),
                exposure=z['plant_group'][use].mean(),oracle_intercept=b))
    for h in range(3):
        for domain in ('uncertain','full'):
            use=d['ok'][:,h] if domain=='uncertain' else np.ones(len(d['n']),bool)
            g=np.column_stack([loss(d['y'][use,h],v[use,h],True) for v in [d['b0'],d['p']]])
            mu,bt,pt=hierarchical(block_arrays(g,d['n'][use],d['block'][use]))
            extras.append(dict(ycut=YS[l],head=HEADS[h],domain=domain,rows=int(use.sum()),hits=int(d['y'][use,h].sum()),B0_loss=mu[0],B1_loss=mu[1]))
    # Actual within-block covariance of side residuals and its two-level CI.
    use=d['ok'][:,0];r=d['y'][use,1:3]-d['p'][use,1:3]
    groups=block_arrays(np.column_stack((r,r.prod(1))),d['n'][use],d['block'][use])
    groups=[(v[:,2]-v[:,0]*v[:,1])[:,None] for v in groups]
    mu,bt,_=hierarchical(groups)
    cov=dict(ycut=YS[l],centered_cov=mu[0],ci99=np.quantile(bt[:,0],[.005,.995]).tolist())
    baseinfo=json.loads(Path('data/manifest.json').read_text())
    parent=[]
    for h,head in enumerate(HEADS[:5]):
        for domain in ['uncertain','full']:
            clusters=[];npairs=0; counts=[0,0];events=[0.,0.]
            for k in np.unique(d['n']):
                pairs={}
                for j,rec in enumerate(baseinfo[k]['parents']):
                    use=(d['n']==k)&(d['block']==j)
                    if domain=='uncertain':use&=d['ok'][:,h]
                    if not use.any():continue
                    yy=d['y'][use,h];p=d['p'][use,h]
                    seed=rec['seed'];counts[seed]+=int(use.sum());events[seed]+=float(yy.sum())
                    pairs.setdefault(rec['pair'],{})[seed]=np.array([yy.mean(),(yy-p).mean()])
                ds=[v[1]-v[0] for v in pairs.values() if len(v)==2]
                if ds:clusters.append(np.stack(ds));npairs+=len(ds)
            mu,bt,_=hierarchical(clusters)
            parent.append(dict(ycut=YS[l],head=head,domain=domain,pairs=npairs,moduli=len(clusters),
                seed_n=counts[1],rough_n=counts[0],seed_events=events[1],rough_events=events[0],raw_RD=mu[0],
                raw_lo=np.quantile(bt[:,0],.005),raw_hi=np.quantile(bt[:,0],.995),B1_residual_RD=mu[1],
                adjusted_lo=np.quantile(bt[:,1],.005),adjusted_hi=np.quantile(bt[:,1],.995)))
    return rows,plants,extras,bins_out,cov,parent


def evaluate(tag='A'):
    lock=json.loads(Path('results/validation_lock.json').read_text())
    forms=NAMES if tag in ('A','Afull') else NAMES[:7]
    out=[test(l,tag,forms,lock[str(YS[l])]) for l in range(2)]
    cols=[sum([v[j] for v in out],[]) for j in [0,1,2,3,5]]
    rows,plants,extras,bins_out,parent=cols
    # Screened-out hypotheses remain in the family; also report unconditional Holm values.
    for table in (rows,plants):
        raw=holm([r['p'] for r in table])
        for r,v in zip(table,raw):r.update(Holm=float(v))
    for r in rows:
        r['declared']=bool(tag=='A' and r['screened'] and r['Holm']<.01 and r['ci_lo']>0 and r['cell_p']<.01 and
                           (r['head'] not in HEADS[:3] or (r['per_hit_lo']>.02 and r['conditional_beta']>0)))
    for r in plants:r['detected']=bool(r['Holm']<.01 and r['ci_lo']>0)
    for suffix,a in zip(['families','plants','baselines','scaling','parents'],[rows,plants,extras,bins_out,parent]):
        write(f'results/{tag}_{suffix}.csv',a)
    js(f'results/{tag}_covariance.json',[v[4] for v in out])


if __name__=='__main__':
    with threadpool_limits(limits=1):
        if sys.argv[1]=='screen':screen_all()
        else:evaluate(sys.argv[1])
