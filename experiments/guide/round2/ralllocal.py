"""Unscreened all-prime-power B1 robustness, with sparse prediction on exactly used columns."""
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path('.extra').resolve()))
import lightgbm as lgb
import joblib
import numpy as np
from scipy import sparse
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits
from rbase import Data,B1,SMALL,lib,SEED,impose,clip,js
from rstats import loss


class Local:
    def __init__(self,d):
        self.d=d;pr=SMALL[SMALL<=d.cut]
        self.mods=np.array([(int(p),m) for p in pr for m in self.powers(int(p),d.cut)],dtype=np.int32)
        self.ids=np.r_[np.arange(len(pr)),len(SMALL)+np.arange(len(pr))]
        self.nb=d.base.shape[1];self.ne=self.nb+len(self.mods);self.ncol=self.ne+len(self.ids)

    @staticmethod
    def powers(p,y):
        m=p
        while m<=y:
            yield m
            m*=p

    def matrix(self,ix,columns=None):
        d=self.d
        cc=np.arange(self.ncol) if columns is None else np.asarray(columns)
        nc=cc[(cc>=self.nb)&(cc<self.ne)]-self.nb;pr,mod=[np.ascontiguousarray(v,dtype='i') for v in self.mods[nc].T]
        xy=np.ascontiguousarray(d.coord[ix,:2],dtype=np.int64);code=np.empty((len(ix),len(pr)),np.uint16)
        lib(xy.ctypes.data,len(ix),pr.ctypes.data,mod.ctypes.data,len(pr),code.ctypes.data)
        vi=self.ids[cc[cc>=self.ne]-self.ne]
        arr=np.column_stack((d.base[ix][:,cc[cc<self.nb]],code,d.val[ix][:,vi].toarray())).astype('f',copy=False)
        if columns is None:return arr
        mat=sparse.csr_matrix(arr)
        mat.indices=cc[mat.indices].astype('i');mat._shape=(len(ix),self.ncol)
        return mat


def work(l):
    t=time.monotonic(); d=Data(l); state=Local(d)
    old=np.load(f'models/A_{l}_pred.npy',mmap_mode='r'); new=np.array(old); report=[]
    fitids=np.unique(d.n[d.target&(d.role==0)])
    for fold in (0,1,2):
        ids=[k for k in fitids if fold==2 or (k%24)%2!=fold];calids=ids[1::4];trids=[k for k in ids if k not in calids]
        train=np.flatnonzero(d.target&np.isin(d.n,trids));cal=np.flatnonzero(d.target&np.isin(d.n,calids))
        rng=np.random.default_rng(SEED+fold)
        tr=rng.choice(train,min(16000,len(train)),replace=False);ca=rng.choice(cal,min(20000,len(cal)),replace=False)
        parent=joblib.load(f'models/A_{l}_{fold}.joblib');oldcal=parent.predict(d,ca)
        test=np.flatnonzero((d.role!=0) if fold==2 else np.isin(d.n,[k for k in fitids if (k%24)%2==fold]))
        cfg=dict(max_bin=31,min_data_in_leaf=45,feature_pre_filter=False,max_depth=4,num_leaves=16,
                 learning_rate=.08,lambda_l2=5.,num_threads=2,verbosity=-1,force_col_wise=True,seed=SEED)
        X=state.matrix(tr)
        ds=lgb.Dataset(X,label=d.resp[tr,0],params=cfg,free_raw_data=True).construct()
        for h in range(5):
            use=d.ok[tr,h]; oc=d.ok[ca,h];yt=d.resp[tr,h]
            if len(np.unique(yt[use]))<2 or not oc.any():
                report.append(dict(fold=fold,head=h,weight=0.,used=0));continue
            ds.set_label(yt);ds.set_weight(use.astype('f'))
            model=lgb.train(dict(cfg,objective='binary' if h<3 else 'regression'),ds,num_boost_round=40)
            used=np.flatnonzero(model.feature_importance()>0)
            if not len(used):
                report.append(dict(fold=fold,head=h,weight=0.,used=0));continue
            # Compare sparse deployment against full-column predictions on fit rows.
            np.testing.assert_allclose(model.predict(X[:24],num_threads=1),model.predict(state.matrix(tr[:24],used),num_threads=1),atol=1e-8)
            tree=model.predict(state.matrix(ca,used),num_threads=1);p0=oldcal[oc,h];pp=tree[oc];yc=d.resp[ca,h][oc]
            objective=lambda w:loss(yc,(1-w)*p0+w*pp,h<3).mean()
            w=min([0.,1.,float(minimize_scalar(objective,bounds=(0,1),method='bounded').x)],key=objective)
            report.append(dict(fold=fold,head=h,weight=w,used=int(len(used))))
            model.save_model(f'models/alllocal_{l}_{fold}_{h}.txt')
            if w:
                for st in range(0,len(test),18000):
                    ix=test[st:st+18000];val=model.predict(state.matrix(ix,used),num_threads=1)
                    new[ix,h]=(1-w)*old[ix,h]+w*val
        print(l,fold,len(tr),state.ncol,flush=True)
        del ds,X
    new=impose(new,d.flags);np.save(f'models/Afull_{l}_pred.npy',new)
    js(f'logs/alllocal_{l}.json',dict(columns=state.ncol,models=report,seconds=time.monotonic()-t))


if __name__=='__main__':
    with threadpool_limits(limits=1):work(int(sys.argv[1]))
