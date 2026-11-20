"""Cross-fitted offsets with train-only size projection and family-specific alternatives."""
import gc
import json
import sys
import time
from pathlib import Path

import joblib
import numpy as np
from scipy.special import expit,logit
from sklearn.preprocessing import SplineTransformer
from sklearn.tree import DecisionTreeRegressor
from threadpoolctl import threadpool_limits
from rbase import Data,clip,SEED,js
from rfeatures import load_global

HEADS=['joint','rational','algebraic','W0','W1','cross']
NAMES=[f'G{i}' for i in range(1,9)]+['all']


def targets(d,pred,t,lab):
    y=lab[t]; p=pred[t]
    return (np.column_stack((y,(y[:,1]-p[:,1])*(y[:,2]-p[:,2]))),
            np.column_stack((p,p[:,0]-p[:,1]*p[:,2])),
            np.column_stack((d.ok[t],d.ok[t,0])))


class Project:
    def fit(self,s,z,i):
        self.spl=SplineTransformer(n_knots=20,degree=3,include_bias=False).fit(s[i])
        x=self.spl.transform(s[i]).astype(float)
        x=np.column_stack((np.ones(len(x)),x))
        self.co=np.linalg.solve(x.T@x+np.eye(x.shape[1])*.1,x.T@z[i])
        val=z[i]-x@self.co; self.scale=np.maximum(val.std(0),.015)
        return self
    def transform(self,s,z):
        out=np.empty_like(z,dtype=np.float32)
        for i in range(0,len(z),40000):
            sl=slice(i,i+40000); x=self.spl.transform(s[sl])
            out[sl]=np.clip((z[sl]-self.co[0]-x@self.co[1:])/self.scale,-6,6)
        return out


def pred_delta(base,delta,binary):
    return expit(logit(clip(base))+np.clip(delta,-5,5)) if binary else base+delta


class Correction:
    def fit(self,z,y,base,ok,binary):
        self.binary=binary; self.co=[]; self.trees=[]
        for j,b in enumerate(binary):
            ix=ok[:,j]; x=z[ix].astype(float); yy=y[ix,j]; p=base[ix,j]
            beta=np.zeros(z.shape[1])
            if len(x)>10:
                for _ in range(5 if b else 1):
                    pp=pred_delta(p,x@beta,b); h=np.maximum(pp*(1-pp),.0003) if b else np.ones(len(p))
                    penalty=(.00035*len(p)+1)
                    step=np.linalg.solve((x.T*h)@x+penalty*np.eye(len(beta)),x.T@(yy-pp)-penalty*beta)
                    beta+=np.clip(step,-1,1)
            self.co.append(beta)
            heads=[]; delta=np.zeros(len(x))
            for i in range(5 if len(x)>60 and np.var(yy)>0 else 0):
                pp=pred_delta(p,delta,b)
                h=np.maximum(pp*(1-pp),.002) if b else np.ones(len(p))
                grad=(yy-pp)/h
                tree=DecisionTreeRegressor(max_depth=2,min_samples_leaf=max(80,len(x)//80),random_state=SEED)
                tree.fit(x,grad,sample_weight=h)
                delta+=.18*np.clip(tree.predict(x),-2,2); heads.append(tree)
            self.trees.append(heads)
        return self

    def delta(self,z):
        out=z@np.array(self.co,dtype='f').T*.65
        for j,trees in enumerate(self.trees):
            for tr in trees: out[:,j]+=.35*.18*np.clip(tr.predict(z),-2,2)
        return out


class SmallNet:
    """One hidden tanh layer with the frozen logit/identity offsets in the loss."""
    def fit(self,z,y,p,ok,binary):
        rng=np.random.default_rng(SEED+37); nf=z.shape[1]; nh=24
        self.w=rng.normal(0,1/np.sqrt(nf),(nf,nh)).astype('f')
        self.c=np.zeros(nh,'f'); self.v=np.zeros((nh,y.shape[1]),'f')
        pars=[self.w,self.c,self.v]; ms=[np.zeros_like(a) for a in pars]; vs=[a.copy() for a in ms]
        by=np.array(binary); norm=np.maximum(.025,np.mean(ok*np.where(by,p*(1-p),1),axis=0))
        step=0
        for _ in range(10):
            perm=rng.permutation(len(z))
            for start in range(0,len(z),512):
                ii=perm[start:start+512]; x=z[ii]; h=np.tanh(x@self.w+self.c); d=h@self.v
                pp=np.where(by,expit(logit(clip(p[ii]))+np.clip(d,-5,5)),p[ii]+d)
                er=(pp-y[ii])*ok[ii]/(len(ii)*norm)
                dh=(er@self.v.T)*(1-h*h)
                grads=[x.T@dh+.003*self.w,dh.sum(0),h.T@er+.01*self.v];step+=1
                for par,mm,vv,gr in zip(pars,ms,vs,grads):
                    mm*=.9;mm+=.1*gr; vv*=.999;vv+=.001*gr**2
                    par-=.0007*(mm/(1-.9**step))/(np.sqrt(vv/(1-.999**step))+1e-7)
        return self
    def delta(self,x): return np.tanh(x@self.w+self.c)@self.v


def corrected(z,co,nn,base,binary):
    out=np.empty_like(base,dtype='f')
    for st in range(0,len(z),24000):
        sl=slice(st,st+24000); d=co.delta(z[sl])
        if nn is not None: d=.65*d+.35*nn.delta(z[sl])
        for j,b in enumerate(binary):
            out[sl,j]=pred_delta(base[sl,j],d[:,j],b)
    if out.shape[1]>=5: out[:,3:5]=np.clip(out[:,3:5],0,1)
    return out


def threshold(a):
    v,n=np.unique(a,return_counts=True)
    if len(v)<2:return float(v[0])
    i=np.argmin(abs(n.cumsum()[:-1]-.5*len(a)))
    return (float(v[i])+float(v[i+1]))/2


def run(l,tag):
    start=time.monotonic()
    d=Data(l); pred=np.load(f'models/{tag}_{l}_pred.npy',mmap_mode='r')
    lab=d.resp
    if tag=='Bnull':
        lab=d.resp.copy();rng=np.random.default_rng(SEED+943+l)
        lab[:,1:3]=(rng.random((len(lab),2))<d.p0[:,1:]);lab[:,0]=lab[:,1]*lab[:,2]
    t,glob,domains,anchors=load_global(d,pred,tag,lab)
    y,p,ok=targets(d,pred,t,lab)
    phase=d.role[t].copy()
    if tag in ('B','Bnull'):phase=np.where(d.role[t]!=0,-1,np.where(d.region[t]==0,0,np.where(d.region[t]==1,2,-1)))
    s=d.base[t,:4]; bflag=np.array([1,1,1,0,0,0],bool)
    fold=d.n[t];bid=d.bid[t]
    seal=Path('results')/f'{tag}_{l}';seal.mkdir(exist_ok=True)
    splitnames={1:'val',2:'test'} if tag in ('A','Afull') else {2:'region'}
    for code,name in splitnames.items():
        sel=phase==code
        np.savez_compressed(seal/f'{name}_data.npz',p=p[sel],y=y[sel],ok=ok[sel],n=fold[sel],block=bid[sel],
                            s=s[sel],sign=d.base[t[sel],6:8],b0=d.p0[t[sel]],flags=d.flags[t[sel]],leak=d.leak[t[sel]],
                            xy=d.coord[t[sel],:2],region=d.region[t[sel]])
    origlen=len(glob)
    for j,name in enumerate(NAMES):
        if tag in ('B','Bnull') and j in (7,8):continue
        if j==origlen:
            z=np.column_stack(glob); dom=np.ones(len(t),bool);anchor=anchors[5]
        else:z=glob[j];dom=domains[j];anchor=anchors[j]
        train=np.flatnonzero((phase==0)&dom)
        rng=np.random.default_rng(SEED+l+51)
        ti=rng.choice(train,min(100000,len(train)),replace=False)
        proj=Project().fit(s,z,ti)
        zz=proj.transform(s,z)
        co=Correction().fit(zz[ti],y[ti],p[ti],ok[ti],bflag)
        nn=SmallNet().fit(zz[ti],y[ti],p[ti],ok[ti],bflag) if j in (7,8) else None
        th=threshold(anchor[ti][ok[ti,0]])
        sign=np.where(anchor>th,1.,-1.)
        basep=clip(p[:,0]); planted_p=expit(logit(basep)+.35*sign)
        plant_y=(np.random.default_rng(SEED+733*j+31*l).random(len(t))<planted_p).astype('f')
        cplant=Correction().fit(zz[ti],plant_y[ti,None],basep[ti,None],ok[ti,:1],[True])
        nplant=SmallNet().fit(zz[ti],plant_y[ti,None],basep[ti,None],ok[ti,:1],[True]) if j in (7,8) else None
        for code,sn in splitnames.items():
            keep=phase==code
            phat=corrected(zz[keep],co,nn,p[keep],bflag)
            ff=d.flags[t[keep]]
            # Use corrected side estimates where joint status is reduced to one uncertain side.
            from rbase import impose
            phat=impose(phat,ff)
            pp=corrected(zz[keep],cplant,nplant,basep[keep,None],[True]).ravel()
            np.savez_compressed(seal/f'{sn}_{name}.npz',pred=phat,domain=dom[keep],
                plant_y=plant_y[keep],plant_p=pp,plant_oracle=planted_p[keep].astype('f'),plant_group=(sign[keep]>0))
        joblib.dump((proj,co,nn,cplant,nplant,th),f'models/{tag}_{l}_{name}.joblib',compress=2)
        print(tag,l,name,len(ti),flush=True)
        del zz;gc.collect()
    js(f'logs/families_{tag}_{l}.json',dict(seconds=time.monotonic()-start))


if __name__=='__main__':
    with threadpool_limits(limits=1):run(int(sys.argv[1]),sys.argv[2] if len(sys.argv)>2 else 'A')
