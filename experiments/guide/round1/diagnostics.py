"""Descriptive frozen-model sensitivity without blending, and data provenance counts."""
import csv
import hashlib
import json
import time
from pathlib import Path

import joblib
import numpy as np
from threadpoolctl import threadpool_limits

from experiment import *
from settings import REFS
from summarize import interval, table


def run():
    start,cpu=time.monotonic(),time.process_time()
    root=Path("results/true")
    size,params,models=joblib.load(root/"models.joblib")
    b=load_split("test")
    raw=size.scores(b)
    p0=calibrated(raw,params)
    x=make_x(b,raw)
    groups=[(it["degree"],it["id"],g,np.flatnonzero((b.n==it["id"])&(b.geo[:,0]==g)))
            for it in manifest() if it["split"]=="test" for g in (0,1)]
    p_log=np.array([loss(b.y[ix,0],p0[ix,0]).mean() for _,_,_,ix in groups])
    by_form={}
    for name,m in models.items():
        pc=calibrated(m.raw(x),m.cal)[:,0]
        by_form[name]=np.array([loss(b.y[ix,0],pc[ix]).mean() for _,_,_,ix in groups])
    rows=[]
    for d in (3,4):
        keep=np.array([dg==d for dg,_,_,_ in groups])
        for name in FORMS:
            scores=by_form[name]
            for ref,s in [("absolute",scores),("logistic",p_log-scores),
                          ("same_learner_base",by_form[name.split("_")[0]+"_base"]-scores)]:
                estimates=s[keep].reshape(-1,2).mean(axis=1)
                rows.append([d,name,ref,*interval(estimates)])
    table(root/"unblended_joint.csv",["degree","model","comparator","moduli","estimate","SE","CI95_low","CI95_high"],rows)
    rows=[]
    for it in manifest():
        k=it["id"]
        with np.load(f"data/n{k:02}.npz") as z:
            y,coords=z["sides"],z["coords"]
        for g in (0,1):
            for stream,centers in [("target",np.arange(len(y))>=REFS),("reference",np.arange(len(y))<REFS)]:
                yy=y[(coords[:,3]==g)&centers]
                rows.append([k,it["degree"],it["split"],it["bits"],it["N"],["box","special_q"][g],stream,len(yy),
                             int(yy[:,0].sum()),int(yy[:,1].sum()),int(yy.prod(axis=1).sum())])
    table(Path("results/dataset_counts.csv"),["N_id","degree","split","bits","N","geometry","stream","rows",
          "rational","algebraic","joint"],rows)
    for mode in MODES:
        with np.load(Path("results")/mode/"heldout_scores.npz") as z:
            y=z["y"]
            assert np.array_equal(z["n"],b.n) and np.array_equal(z["geometry"],b.geo[:,0])
            if mode=="true":
                assert np.array_equal(y,b.y)
            if mode=="shuffled":
                assert np.array_equal(y[:,0],y[:,1]*y[:,2])
                for _,_,_,ix in groups:
                    assert np.array_equal(y[ix].sum(axis=0),b.y[ix].sum(axis=0))
            if mode=="planted":
                assert np.array_equal(y[:,1:],b.y[:,1:])
                assert np.all((z["oracle"]>0)&(z["oracle"]<1))
    dump("logs/diagnostics_time.json",dict(wall=time.monotonic()-start,cpu=time.process_time()-cpu))


if __name__=="__main__":
    with threadpool_limits(limits=1,user_api="blas"),threadpool_limits(limits=2,user_api="openmp"):
        run()
