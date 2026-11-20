"""Oracle-feature, deliberately unsafe arithmetic, native-feature power, and CRT controls."""
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.special import expit,logit
from threadpoolctl import threadpool_limits
from rbase import Data,SEED,YS,B,PR,clip,project,js
from rfit import Project,Correction,SmallNet,corrected,threshold,NAMES
from rfeatures import load_global
from rstats import block_arrays,hierarchical,holm,write,loss

def model(s,x,y,p,train,test,mlp=False,cap=100000):
    rng=np.random.default_rng(SEED+11); tr=rng.choice(train,min(cap,len(train)),replace=False)
    proj=Project().fit(s,x,tr); xx=proj.transform(s,x)
    ok=np.ones((len(tr),1),bool)
    c=Correction().fit(xx[tr],y[tr,None],p[tr,None],ok,[True])
    nn=SmallNet().fit(xx[tr],y[tr,None],p[tr,None],ok,[True]) if mlp else None
    return corrected(xx[test],c,nn,p[test,None],[True]).ravel()


def intercept(p,q):
    eta=logit(clip(p)); b=0.
    for _ in range(12):
        pp=expit(eta+b);b+=np.sum(q-pp)/max(1e-9,np.sum(pp*(1-pp)))
    return expit(eta+b)


def measured(y,p0,p1,n,b,extra=None):
    g=loss(y,p0,True)-loss(y,p1,True)
    vals=np.column_stack((g,y)) if extra is None else np.column_stack((g,y,extra))
    mu,bt,pt=hierarchical(block_arrays(vals,n,b))
    return dict(gain=float(mu[0]),lo=float(np.quantile(bt[:,0],.005)),hi=float(np.quantile(bt[:,0],.995)),p=float(pt[0]),hits=int(y.sum()),n=len(y)),mu,bt


def pow_reps(oracle,p0,new,n,blocks,seed):
    rng=np.random.default_rng(seed)
    key,ix,cnt=np.unique(n.astype(int)*128+blocks,return_inverse=True,return_counts=True)
    delta=logit(clip(new)).astype(float)-logit(clip(p0))
    c=np.log1p(-clip(new)).astype(float)-np.log1p(-clip(p0))
    rates=[]
    for _ in range(100):
        y=rng.random(len(oracle))<oracle
        rates.append(np.bincount(ix,weights=c+y*delta)/cnt)
    rates=np.array(rates).T;groups=[rates[key//128==k] for k in np.unique(n)]
    mu,bt,pr=hierarchical(groups)
    return pr,np.quantile(bt,.005,axis=0)


def run(l):
    d=Data(l); pred=np.load(f'models/A_{l}_pred.npy',mmap_mode='r')
    t,glob,doms,anchors=load_global(d,pred)
    s=d.base[t,:4]; n=d.n[t];b=d.bid[t];phase=d.role[t];y=d.resp[t];p=pred[t]
    ok=d.ok[t];U=np.sum((s[:,:2]/np.log(B))*np.log(np.maximum(1e-8,s[:,:2]/np.log(B))),axis=1)
    pc1=[];pc2=[];native=[];crt=[];pvals={};lows={}
    for j,theta in enumerate([.01,.02,.05,.10,.25]):
        rng=np.random.default_rng(SEED+83);sig=np.sqrt(2*theta*np.maximum(0,U));G=rng.normal(size=len(t))+y[:,0]*sig
        warp=np.arcsinh(G)
        x=np.column_stack((warp,warp*np.sqrt(np.maximum(0,U)),np.tanh(G),rng.normal(size=(len(t),3)))).astype('f')
        tr=np.flatnonzero((phase==0)&ok[:,0]); te=np.flatnonzero((phase==2)&ok[:,0])
        pr=model(s,x,y[:,0],p[:,0],tr,te)
        z,mu,bt=measured(y[te,0],p[te,0],pr,n[te],b[te],U[te]*y[te,0])
        oracle=expit(logit(clip(p[te,0]))+sig[te]*G[te]-.5*sig[te]**2)
        oo,om,_=measured(y[te,0],p[te,0],oracle,n[te],b[te],U[te]*y[te,0])
        ratio=mu[0]/(theta*mu[2])
        pc1.append(dict(ycut=YS[l],theta=theta,**z,recovery=ratio,oracle_recovery=om[0]/(theta*om[2]),gate=.8<=ratio<=1.15))
    for h in (0,1,2):
        # The rational bit equals divisibility of the parallel shift a+P_(y,2y], b.
        x=d.leak[t].astype('f');tr=np.flatnonzero((phase==0)&ok[:,h]);te=np.flatnonzero((phase==2)&ok[:,h])
        new=model(s,x,y[:,h],p[:,h],tr,te)
        z,mu,_=measured(y[te,h],p[te,h],new,n[te],b[te])
        pc2.append(dict(ycut=YS[l],head=['joint','rational','algebraic'][h],**z,nats_per_hit=mu[0]/max(mu[1],1e-12)))
    for j,name in enumerate(NAMES):
        if j==8:continue
        domain=doms[j]; x=glob[j]; a=anchors[j]; rng=np.random.default_rng(SEED+1311+33*j+l)
        for h in (0,2):
            tr=np.flatnonzero((phase==0)&ok[:,h]&domain);te=np.flatnonzero((phase==2)&ok[:,h]&domain)
            if h==0:
                root=Path(f'results/A_{l}'); z=np.load(root/f'test_{name}.npz')
                use=ok[phase==2,0]&domain[phase==2];q=z['plant_oracle'][use]; phat=z['plant_p'][use];yy=z['plant_y'][use]
            else:
                th=threshold(a[tr]);qall=expit(logit(clip(p[:,h]))+.35*np.where(a>th,1.,-1.))
                yyall=(rng.random(len(a))<qall).astype('f')
                phat=model(s,x,yyall,p[:,h],tr,te,j==7)
                q=qall[te];yy=yyall[te]
            bp=intercept(p[te,h],q)
            r,_,_=measured(yy,bp,phat,n[te],b[te])
            pv,lo=pow_reps(q,bp,phat,n[te],b[te],SEED+202*j+94*h+l)
            key=name+'_'+str(h);pvals[key]=pv;lows[key]=lo
            native.append(dict(ycut=YS[l],family=name,head='joint' if h==0 else 'algebraic',**r))
    # Random CRT challenge on the within-N split. Each held-out sector has fresh hidden roots.
    pr=PR[(PR>YS[l])&(PR<=2*YS[l])][:160]
    c=np.zeros(len(t),dtype='f')
    for k in np.unique(n[phase==0]):
        for region in (0,1):
            ix=np.flatnonzero((n==k)&(d.region[t]==region))
            codes,_=project(d.coord[t[ix]],pr,int(pr[-1]))
            rng=np.random.default_rng(SEED+29011+int(k)*23+region)
            roots=np.array([rng.integers(int(pp)+1) for pp in pr])
            c[ix]=(codes==roots).any(1)
    bp=np.full(len(t),1-np.prod(1-1/(pr.astype(float)+1)))
    tr=np.flatnonzero((phase==0)&(d.region[t]==0));te=np.flatnonzero((phase==0)&(d.region[t]==1))
    new=model(s,glob[7],c,bp,tr,te,True)
    r,_,_=measured(c[te],bp[te],new,n[te],b[te])
    crt.append(dict(ycut=YS[l],**r))
    for name,rows in [('pc1',pc1),('pc2a',pc2),('native',native),('crt',crt)]:write(f'results/{name}_{l}.csv',rows)
    np.savez(f'results/native_reps_{l}.npz',**{f'p_{k}':v for k,v in pvals.items()},**{f'l_{k}':v for k,v in lows.items()})


def aggregate():
    import csv
    rs=sum([list(csv.DictReader(open(f'results/native_{l}.csv'))) for l in range(2)],[])
    for head,h in [('joint',0),('algebraic',2)]:
        parts=[r for r in rs if r['head']==head];ps=[];ls=[]
        for r in parts:
            l=YS.index(int(r['ycut'])); z=np.load(f'results/native_reps_{l}.npz')
            key=r['family']+'_'+str(h);ps.append(z['p_'+key]);ls.append(z['l_'+key])
        adj=np.stack([holm(v) for v in np.array(ps).T]).T
        power=((adj<.01)&(np.array(ls)>0)).mean(1)
        for r,fp,pw in zip(parts,holm([float(v['p']) for v in parts]),power):
            r.update(Holm=fp,detected=(fp<.01 and float(r['lo'])>0),conditional_test_power=pw)
    write('results/native_sensitivity.csv',rs)


if __name__=='__main__':
    with threadpool_limits(limits=1):
        if sys.argv[1]=='aggregate':aggregate()
        else:run(int(sys.argv[1]))
