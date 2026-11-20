"""Frozen-model held-out evaluation. Outputs row/cluster-level data underlying the report."""
import csv
import hashlib
import itertools
import sys
import time
from pathlib import Path

import joblib
import numpy as np
from scipy.stats import t
from sklearn.metrics import roc_auc_score
from threadpoolctl import threadpool_limits
from experiment import *
from settings import PER_FAMILY, PERM_PER_GEOMETRY


def records(path, rows):
    with open(path,"w",newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def inference(d):
    d = np.asarray([v for v in d if v is not None],dtype=float)
    if len(d) < 2:
        return dict(n=len(d), mean=None, se=None, ci95=None, sign_p=None)
    mean = np.mean(d)
    se = d.std(ddof=1)/np.sqrt(len(d))
    half = t.ppf(.975,len(d)-1)*se
    signs = np.array(list(itertools.product((1,-1), repeat=len(d))))
    p = float(np.mean((signs*d).mean(axis=1) >= mean-1e-15))
    return dict(n=len(d),mean=float(mean),se=float(se),ci95=[float(mean-half),float(mean+half)],sign_p=p)


def permutation(model, batch, x, p0, sel_p, root, degree):
    rng = np.random.default_rng(SEED+119*degree)
    ids, strata = [], []
    offset = 0
    for k in np.unique(batch.n):
        if (k<PER_FAMILY) != (degree==3):
            continue
        for g in (0,1):
            ix = np.flatnonzero((batch.n==k)&(batch.geo[:,0]==g))
            cuts = np.quantile(p0[ix],np.arange(1,10)/10)
            sampled = rng.choice(ix,PERM_PER_GEOMETRY,replace=False)
            ids.extend(sampled)
            bins = np.searchsorted(cuts,p0[sampled],side="left")
            for j in np.unique(bins):
                strata.append(offset+np.flatnonzero(bins==j))
            offset += len(sampled)
    ids = np.array(ids)
    y = batch.y[ids,0]
    x0 = np.ascontiguousarray(x[ids,:model.dim])
    buf = x0.copy()
    original = float(loss(y,sel_p[ids]).mean())
    assert np.allclose(sel_p[ids],model.joint(x0,p0[ids]),atol=2e-6)
    perm_losses=[]
    ix = np.arange(len(ids))
    for _ in range(99):
        for group in strata:
            ix[group] = rng.permutation(group)
        buf[:,27:] = x0[ix,27:]
        pred = model.joint(buf,p0[ids])
        perm_losses.append(float(loss(y,pred).mean()))
    np.savez_compressed(root/f"permutation_d{degree}.npz", rows=ids, loss=perm_losses)
    # Degenerate models return exactly their input base score; numerical ties count against detection.
    p = (1+np.sum(np.array(perm_losses) <= original+1e-13))/(1+len(perm_losses))
    return dict(rows=len(ids),draws=len(perm_losses),original_loss=original,
                mean_permuted_loss=np.mean(perm_losses), p=float(p))


def top(y, p):
    n = len(y)//5
    ix = np.argsort(-p,kind="stable")[:n]
    return n,int(y[ix].sum())


def evaluate(mode, teacher):
    root = Path("results")/mode
    start = time.monotonic()
    lock = json.loads((root/"selection_lock.json").read_text())
    assert hashlib.sha256((root/"models.joblib").read_bytes()).hexdigest()==lock["model_sha256"]
    selected = lock["selected"]
    matched = selected.split("_")[0]+"_base"
    size, cal, models = joblib.load(root/"models.joblib")
    batch = load_split("test",mode,teacher)
    t0 = time.monotonic()
    raw = size.scores(batch)
    p0 = calibrated(raw,cal)
    timers = {"logistic":time.monotonic()-t0}
    x = make_x(batch,raw)
    predictions = {"logistic":p0}
    for name, model in models.items():
        t0 = time.monotonic()
        predictions[name] = model.predict(x,p0)
        timers[name] = time.monotonic()-t0
    np.savez_compressed(root/"heldout_scores.npz", n=batch.n, geometry=batch.geo[:,0].astype(np.uint8),
                        y=batch.y, oracle=batch.oracle, **predictions)
    by_n = {it["id"]:it for it in manifest()}
    rows, ledger = [], []
    for k in np.unique(batch.n):
        it = by_n[k]
        for g in (0,1):
            take = (batch.n==k)&(batch.geo[:,0]==g)
            geom = ["box","special_q"][g]
            y = batch.y[take]
            for name,probs in predictions.items():
                p = probs[take]
                for h,head in enumerate(HEADS):
                    yh = y[:,h]
                    both = len(np.unique(yh))==2
                    budget,found = top(yh,p[:,h])
                    rows.append(dict(N_id=int(k),degree=it["degree"],bits=it["bits"],geometry=geom,head=head,
                        model=name,rows=len(yh),positives=int(yh.sum()),loss=float(loss(yh,p[:,h]).mean()),
                        auc=float(roc_auc_score(yh,p[:,h])) if both else None,budget=budget,found=found))
            u = np.max(batch.base[take,:2],axis=1)
            bins = np.searchsorted([2.5,3.5,4.5],u,side="right")
            for b in [-1,0,1,2,3]:
                idx = np.arange(len(u)) if b==-1 else np.flatnonzero(bins==b)
                entry = dict(N_id=int(k),degree=it["degree"],bits=it["bits"],geometry=geom,u_bin=b,rows=len(idx),
                             positives=int(y[idx,0].sum()))
                for name in ["logistic",matched,selected]:
                    entry["budget"],entry[name] = top(y[idx,0],predictions[name][take,0][idx])
                ledger.append(entry)

    records(root/"by_modulus_geometry.csv",rows)
    records(root/"yield_ledger.csv",ledger)
    lookup = {(r["N_id"],r["geometry"],r["head"],r["model"]):r for r in rows}
    summary = {}
    for d in (3,4):
        keys = [k for k in np.unique(batch.n) if by_n[k]["degree"]==d]
        comps = {}
        for head in HEADS:
            comps[head] = {}
            for ref in ["logistic", matched]:
                comps[head][ref] = {}
                for metric in ["loss","auc"]:
                    ds=[]
                    for k in keys:
                        delta=[]
                        for g in ["box","special_q"]:
                            v=lookup[(k,g,head,selected)][metric]
                            b=lookup[(k,g,head,ref)][metric]
                            if v is not None and b is not None:
                                delta.append(b-v if metric=="loss" else v-b)
                        # Both geometries must define AUC for an N-level paired difference.
                        if len(delta)==2:
                            ds.append(np.mean(delta))
                    comps[head][ref][metric] = inference(ds)
        perm = permutation(models[selected],batch,x,p0[:,0],predictions[selected][:,0],root,d)
        joint = [v["loss"] for v in comps["joint"].values()]
        passer = all(z["mean"] > 3*z["se"] and z["sign_p"]<=.025 for z in joint)
        summary[str(d)] = dict(comparisons=comps,permutation=perm,passes_main_comparisons=bool(passer),
                               passes_family=bool(passer and perm["p"]<=.05))
        if mode.startswith("planted"):
            take = np.isin(batch.n,keys)
            q = clip(batch.oracle[take].astype(np.float64))
            base = clip(p0[take,0])
            oracle_gain = q*np.log(q/base)+(1-q)*np.log((1-q)/(1-base))
            # Equal weight each of the sixteen N/geometry cells, as for the primary metric.
            n = batch.n[take]
            geo = batch.geo[take,0]
            summary[str(d)]["oracle_expected_gain_vs_logistic"] = float(np.mean([
                oracle_gain[(n==k)&(geo==g)].mean() for k in keys for g in (0,1)]))
    out = dict(mode=mode,locked=selected,matched=matched,by_degree=summary,
               crosses_rule=all(v["passes_family"] for v in summary.values()),seconds=time.monotonic()-start,
               score_seconds=timers)
    dump(root/"heldout_summary.json",out)
    print(mode,selected,out["crosses_rule"],[(d,v["comparisons"]["joint"]["logistic"]["loss"],
                                               v["permutation"]["p"]) for d,v in summary.items()], flush=True)


if __name__=="__main__":
    wall,cpu=time.monotonic(),time.process_time()
    with threadpool_limits(limits=1,user_api="blas"),threadpool_limits(limits=2,user_api="openmp"):
        modes=sys.argv[1:] or MODES
        assert all((Path("results")/m/"selection_lock.json").exists() for m in MODES)
        size,cal,_=joblib.load("results/true/models.joblib")
        for m in modes:
            evaluate(m,(size,cal))
    dump("logs/testing_time_"+ "_".join(modes)+".json",dict(wall=time.monotonic()-wall,cpu=time.process_time()-cpu))
