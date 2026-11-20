"""Buchstab/saddle model and cross-modulus, saturated local-state comparator."""
import ctypes
import functools
import json
import math
import sys
import time
import warnings
from pathlib import Path

import joblib
import numpy as np
from scipy import sparse
from scipy.interpolate import BSpline
from scipy.special import expit, logit
from sklearn.ensemble import HistGradientBoostingClassifier,HistGradientBoostingRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression,Ridge
from sklearn.preprocessing import SplineTransformer
from threadpoolctl import threadpool_limits
from rdata import B,YS,SMALL,PR,SEED,js

EP=1e-6
clip=lambda x:np.clip(x,EP,1-EP)
lib=ctypes.CDLL(str(Path('local.dylib').resolve())).project
lib.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int,ctypes.c_void_p]


def project(coords,pr,y):
    coords=np.ascontiguousarray(coords[:,:2],dtype=np.int64); pr=np.ascontiguousarray(pr,dtype=np.int32)
    mods=pr.copy()
    while np.any(mods<=y//pr):
        sel=mods<=y//pr; mods[sel]*=pr[sel]
    out=np.empty((len(coords),len(pr)),dtype=np.uint16)
    lib(coords.ctypes.data,len(coords),pr.ctypes.data,mods.ctypes.data,len(pr),out.ctypes.data)
    return out,mods


@functools.lru_cache(2)
def saddle_grid(y):
    # Density per unit x from the finite Euler product. Buchstab solves (u*w)'=w(u-1).
    ps=PR[(PR>y)&(PR<=B)].astype(float); lp=np.log(ps)
    a=np.linspace(.68,1.22,500)
    grid=[]; vals=[]
    du=.001; us=np.arange(0,20+du,du); om=np.zeros_like(us)
    om[(us>=1)&(us<=2)]=1/us[(us>=1)&(us<=2)]
    for j in range(2001,len(us)): om[j]=(us[j-1]*om[j-1]+du*om[j-1000])/us[j]
    for beta in a:
        v=np.exp(-beta*lp); s=np.sum(lp*v/(1-v))
        var=np.sum(lp**2*v/(1-v)**2)
        ld=(beta-1)*s - np.log1p(-v).sum()-.5*np.log(2*np.pi*var)
        rough=np.interp(s/math.log(y),us,om)/math.log(y)
        grid.append(s); vals.append(np.clip(np.exp(min(0.,ld))/max(rough,1e-12),EP,1-EP))
    return np.array(grid[::-1]),np.array(vals[::-1])


def impose(p,flags):
    p=p.copy()
    for s in range(2):
        p[:,1+s]=np.where(flags[:,s]==0,p[:,1+s],flags[:,s]>0)
        p[:,3+s]=np.where(flags[:,s]==0,p[:,3+s],flags[:,s]>0)
    p[:,0]=np.where((flags<0).any(1),0,np.where(flags[:,0]>0,p[:,2],np.where(flags[:,1]>0,p[:,1],p[:,0])))
    p[:,:3]=clip(p[:,:3])
    return p


class Data:
    def __init__(self,l):
        self.l=l; self.cut=YS[l]; self.info=json.loads(Path("data/manifest.json").read_text())
        parts={name:[] for name in ["coord","base","flags","resp","n","target","bid","region","p0","leak"]}
        vals=[]; self.idx=[]; self.tidx=[]; offset=0
        grid,pv=saddle_grid(self.cut)
        for it in self.info:
            k=it["id"]; nt=it["nt"]
            with np.load(f"data/n{k:02}.npz") as z:
                xy=z['coords']; s=z['sizes'][:,l,:]; loc=z['local'][:,l,:]; raw=z['raw']; f=z['flags'][:,l]
                signs=z['signs']; nf=len(xy); alpha=np.broadcast_to(z['alpha'][l],(nf,2)); side=z['y']; w=z['w'][:,l]
                t=np.zeros(nf,bool); t[:nt]=np.load(f"data/keep{k:02}.npy")
                bid=np.full(nf,-1,np.int16); bid[:nt]=z['blocks']
                ang=(np.arctan2(xy[:,0],2.2*xy[:,1])+np.pi/2)/np.pi
                tile=np.minimum(3,(ang*4).astype(int)); gap=np.minimum(ang*4-tile,1-(ang*4-tile))>=.10
                bands=np.array([r['band'] for r in it['parents']]+[0]*16)
                qband={r['q']:r['band']%2 for r in it['parents']}
                band=np.array([qband.get(int(q),-1) for q in xy[:,2]],dtype=np.int8)
                region=np.where(gap & ((band==-1)|(band==tile%2)),tile%2,-1)
                vv=sparse.csr_matrix((z['ve'],z['vi'],z['ptr']),shape=(nf,2*len(SMALL)),dtype=np.float32)
                pa=np.column_stack([np.interp(s[:,i]+alpha[:,i],grid,pv) for i in range(2)])
                p0=impose(np.column_stack((pa.prod(1),pa,pa)),f)[:,:3]
                base=np.column_stack((s,loc,raw,signs,f,alpha,s[:,0]*s[:,1],p0))
                for name,v in zip(parts,[xy,base.astype(np.float32),f,np.column_stack((side.prod(1),side,w)).astype(np.float32),
                             np.full(nf,k,dtype=np.int16),t,bid,region.astype(np.int8),p0.astype(np.float32),z['leak'][:,l]]):
                    parts[name].append(v)
                vals.append(vv); self.idx.append(np.arange(offset,offset+nf)); self.tidx.append(np.arange(offset,offset+nt))
                offset+=nf
        for name,v in parts.items(): setattr(self,name,np.concatenate(v))
        self.val=sparse.vstack(vals,format='csr')
        self.role=np.array([['fit','val','test'].index(it['split']) for it in self.info])[self.n]
        self.ok=np.column_stack(((self.flags==0).all(1),self.flags==0,self.flags==0))


def screen(d,idx,resp):
    rng=np.random.default_rng(SEED); sel=rng.choice(idx,min(6000,len(idx)),replace=False)
    pr=SMALL[SMALL<=d.cut]
    # All maximal prime powers are examined in nested fit data. 8-bin score is multi-outcome.
    scores=[]; residual=resp[sel].copy()
    residual[:,:3]-=d.p0[sel]; residual-=residual.mean(0)
    for j in range(0,len(pr),64):
        C,M=project(d.coord[sel],pr[j:j+64],d.cut)
        bins=np.minimum(7,C.astype(float)*8/(M+M/pr[j:j+64])).astype(int)
        out=np.zeros(bins.shape[1])
        for h in range(bins.shape[1]):
            cnt=np.bincount(bins[:,h],minlength=8)+30.
            score=sum(np.sum(np.bincount(bins[:,h],weights=r,minlength=8)**2/cnt) for r in residual.T)
            out[h]=score
        scores.extend(out)
    top=np.unique(np.r_[np.arange(8),np.argsort(scores)[-40:]])
    available=np.r_[np.arange(len(pr)),len(SMALL)+np.arange(len(pr))]
    v=d.val[sel][:,available]
    strength=np.sum((v.T@residual)**2,axis=1)/(np.asarray(v.power(2).sum(0)).ravel()+12)
    vi=available[np.unique(np.r_[np.arange(8),len(pr)+np.arange(8),np.argsort(strength)[-56:]])]
    return pr[top],vi


class B1:
    def xd(self,d,idx):
        raw=d.base[idx]; code,mods=project(d.coord[idx],self.pr,d.cut)
        return np.column_stack((raw,code.astype(np.float32),d.val[idx][:,self.vi].toarray()))

    def linear(self,x):
        return sparse.hstack((self.spl.transform(x[:,:6]),sparse.csr_matrix((x[:,6:16]-self.mu)/self.sd)),format='csr')

    def fit(self,d,train,cal,resp):
        rng=np.random.default_rng(SEED+93)
        tr=rng.choice(train,min(65000,len(train)),replace=False); ca=rng.choice(cal,min(40000,len(cal)),replace=False)
        self.pr,self.vi=screen(d,tr,resp)
        xd=self.xd(d,tr); xc=self.xd(d,ca)
        self.spl=SplineTransformer(n_knots=30,degree=3,knots="uniform",include_bias=False,sparse_output=True).fit(xd[:,:6])
        self.mu=xd[:,6:16].mean(0); self.sd=np.maximum(.1,xd[:,6:16].std(0))
        xl=self.linear(xd); xlc=self.linear(xc)
        self.models=[]
        for j in range(5):
            ok=d.ok[tr,j]; oc=d.ok[ca,j]; y=resp[tr,j][ok]
            if j<3 and (len(np.unique(y))<2 or not oc.any()):
                self.models.append((None,None,None))
                continue
            if j<3:
                gam=LogisticRegression(C=1.2,max_iter=200,solver='lbfgs')
                tree=HistGradientBoostingClassifier(max_iter=50,max_depth=4,max_leaf_nodes=16,min_samples_leaf=60,l2_regularization=4,early_stopping=False,random_state=SEED)
            else:
                gam=Ridge(alpha=12,solver='lsqr'); tree=HistGradientBoostingRegressor(max_iter=50,max_depth=4,max_leaf_nodes=16,min_samples_leaf=60,l2_regularization=4,early_stopping=False,random_state=SEED)
            gam.fit(xl[ok],y); tree.fit(xd[ok],y)
            fn=lambda m,a:m.predict_proba(a)[:,1] if j<3 else m.predict(a)
            pred=.65*fn(gam,xlc[oc])+.35*fn(tree,xc[oc])
            iso=IsotonicRegression(out_of_bounds="clip",y_min=1e-5,y_max=1-1e-5).fit(pred,resp[ca,j][oc])
            self.models.append((gam,tree,iso))
        return self

    def predict(self,d,idx):
        out=np.empty((len(idx),5),dtype=np.float32)
        for j in range(0,len(idx),35000):
            ix=idx[j:j+35000]; xd=self.xd(d,ix); xl=self.linear(xd); parts=[]
            for h,(gam,tree,iso) in enumerate(self.models):
                if gam is None:
                    parts.append(d.p0[ix,h])
                    continue
                fn=lambda m,a:m.predict_proba(a)[:,1] if h<3 else m.predict(a)
                p=.65*fn(gam,xl)+.35*fn(tree,xd)
                parts.append(np.clip(.7*iso.predict(p)+.3*p,EP,1-EP))
            out[j:j+len(ix)]=impose(np.column_stack(parts),d.flags[ix])
        return out


def fitbase(l,region=False,null=False):
    t=time.monotonic(); d=Data(l); lab=d.resp.copy(); tag=('Bnull' if null else 'B') if region else 'A'
    dest=Path('models')/f'{tag}_{l}'
    if null:
        rng=np.random.default_rng(SEED+943+l)
        lab[:,1:3]=rng.random((len(lab),2))<d.p0[:,1:]
        lab[:,0]=lab[:,1]*lab[:,2]
    eligible=d.target&(d.role==0)
    if region: eligible &= d.region==0
    fits=sorted(set(d.n[eligible]))
    folds={k:(k%24)%2 for k in fits}
    out=np.zeros((len(d.n),5),dtype=np.float32)
    for exc in [0,1,2]:
        pool=[k for k in fits if exc==2 or folds[k]!=exc]
        caids=pool[1::4] or pool[-1:]
        train=np.flatnonzero(eligible&np.isin(d.n,[k for k in pool if k not in caids]))
        cal=np.flatnonzero(eligible&np.isin(d.n,caids))
        model=B1().fit(d,train,cal,lab)
        joblib.dump(model,f'{dest}_{exc}.joblib',compress=2)
        if exc!=2: test=np.flatnonzero(np.isin(d.n,[k for k in fits if folds[k]==exc]))
        else: test=np.flatnonzero(d.role!=0)
        out[test]=model.predict(d,test)
        if exc==2 and region:
            ix=np.flatnonzero((d.role==0)&(d.region==1)&d.target)
            out[ix]=model.predict(d,ix)
    np.save(f'{dest}_pred.npy',out)
    js(f'logs/base_{tag}_{l}.json',dict(seconds=time.monotonic()-t,nodes=len(out)))


if __name__=='__main__':
    with threadpool_limits(limits=1):
        fitbase(int(sys.argv[1]),len(sys.argv)>2,len(sys.argv)>3)
