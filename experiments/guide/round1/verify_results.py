"""Recompute the registered loss contrasts and operating counts from the saved test predictions."""
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import t

from experiment import manifest, dump
from settings import PER_FAMILY, POWER


def verify():
    started,cpu=time.monotonic(),time.process_time()
    records=manifest()
    ids=[it["id"] for it in records if it["split"]=="test"]
    checked=[]
    modes=["true","shuffled","planted"] if POWER else ["true","shuffled","planted","planted_strong"]
    if (Path("results/global_shuffle")/"heldout_summary.json").exists():
        modes.append("global_shuffle")
    for mode in modes:
        directory=Path("results")/mode
        sel=json.loads((directory/"selection_lock.json").read_text())
        out=json.loads((directory/"heldout_summary.json").read_text())
        assert sel["prereg_sha256"]==hashlib.sha256(Path("PREREGISTRATION.md").read_bytes()).hexdigest()
        assert sel["model_sha256"]==hashlib.sha256((directory/"models.joblib").read_bytes()).hexdigest()
        with np.load(directory/"heldout_scores.npz") as z:
            n,geo,y=z["n"],z["geometry"],z["y"]
            if mode=="global_shuffle":
                with np.load("results/true/heldout_scores.npz") as true:
                    np.testing.assert_array_equal(y.sum(axis=0),true["y"].sum(axis=0))
                assert np.array_equal(y[:,0], y[:,1]*y[:,2])
            best=z[out["locked"]]
            assert len(y)==sum(v["targets"] for v in records if v["split"]=="test")
            assert set(np.unique(n))==set(ids)
            assert np.isfinite(best).all() and (best>0).all() and (best<1).all()
            for ref in ("logistic",out["matched"]):
                baseline=z[ref]
                gain=np.where(y,np.log(best/baseline),np.log((1-best)/(1-baseline)))
                for d in (3,4):
                    by_n=[]
                    for k in ids:
                        if (k<PER_FAMILY)!=(d==3): continue
                        by_n.append(np.mean([gain[(n==k)&(geo==g)].mean(axis=0) for g in (0,1)],axis=0))
                    ds=np.asarray(by_n)
                    p_sign=np.zeros(3)
                    for mask in range(2**len(ds)):
                        signs=1-2*((mask>>np.arange(len(ds)))&1)
                        p_sign+=(np.mean(ds*signs[:,None],axis=0)>=ds.mean(axis=0)-1e-15)
                    p_sign/=2**len(ds)
                    for j,h in enumerate(("joint","rational","algebraic")):
                        stat=out["by_degree"][str(d)]["comparisons"][h][ref]["loss"]
                        np.testing.assert_allclose([ds[:,j].mean(),ds[:,j].std(ddof=1)/np.sqrt(len(ds))],
                                                   [stat["mean"],stat["se"]],atol=1e-12,rtol=1e-7)
                        assert p_sign[j]==stat["sign_p"]
            with open(directory/"yield_ledger.csv",newline="") as f:
                for row in csv.DictReader(f):
                    if row["u_bin"]!="-1": continue
                    idx=np.flatnonzero((n==int(row["N_id"]))&(geo==int(row["geometry"]=="special_q")))
                    cutoff=len(idx)//5
                    kept=np.argsort(-best[idx,0],kind="stable")[:cutoff]
                    assert int(row["budget"])==cutoff and int(row[out["locked"]])==int(y[idx[kept],0].sum())
        checked.append(mode)
    dump("results/verification.json",dict(verified_modes=checked,unique_test_N=len(ids),
         wall=time.monotonic()-started,cpu=time.process_time()-cpu))


if __name__=="__main__":
    verify()
