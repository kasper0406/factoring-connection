"""Decompose the group-preserving placebo and run a globally exchangeable placebo."""
import json
import time
from pathlib import Path
from scipy.special import logit

from evaluate import evaluate
from fit_models import fit
from experiment import *
from threadpoolctl import threadpool_limits


def decompose():
    root=Path("results/shuffled")
    summary=json.loads((root/"heldout_summary.json").read_text())
    rng=np.random.default_rng(SEED+771)
    out={}
    with np.load(root/"heldout_scores.npz") as z:
        y=z["y"][:,0]; n=z["n"]; geo=z["geometry"]
        selected=z[summary["locked"]][:,0]
        for comparator in ["logistic",summary["matched"]]:
            p0=z[comparator][:,0]
            out[comparator]={}
            for d in [3,4]:
                observed,expectation,replicas=[],[],[]
                for it in manifest():
                    if it["split"]!="test" or it["degree"]!=d: continue
                    for g in [0,1]:
                        ix=(n==it["id"])&(geo==g)
                        ys=y[ix]; pa=selected[ix]; pb=p0[ix]
                        delta_logit=logit(pa)-logit(pb)
                        offset=np.mean(np.log1p(-pa)-np.log1p(-pb))
                        delta=offset+np.dot(ys,delta_logit)/len(ys)
                        expect=offset+ys.mean()*delta_logit.mean()
                        events=int(ys.sum())
                        sims=np.array([offset+delta_logit[rng.choice(len(ys),events,replace=False)].sum()/len(ys)
                                       for _ in range(499)])
                        observed.append(delta);expectation.append(expect);replicas.append(sims)
                obs,ex=float(np.mean(observed)),float(np.mean(expectation))
                sims=np.mean(replicas,axis=0)
                out[comparator][str(d)]=dict(observed_gain=obs,exchangeable_expected_gain=ex,
                    within_cell_excess=obs-ex,randomization_se=float(sims.std(ddof=1)),
                    p_upper=float((1+(sims>=obs-1e-15).sum())/500))
    dump(root/"within_cell_placebo.json",out)


if __name__=="__main__":
    start,cpu=time.monotonic(),time.process_time()
    with threadpool_limits(limits=1,user_api="blas"),threadpool_limits(limits=2,user_api="openmp"):
        decompose()
        fit("global_shuffle")
        evaluate("global_shuffle",None)
    dump("logs/additional_negative_time.json",dict(wall=time.monotonic()-start,cpu=time.process_time()-cpu))
