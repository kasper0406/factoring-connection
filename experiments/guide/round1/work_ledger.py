"""Replay the three real-data rankings with paid candidate generation, scoring and exact confirmation."""
import csv
import json
import random
import site
import sys
import time
from pathlib import Path

import joblib
import numpy as np
from scipy.special import logit
from threadpoolctl import threadpool_limits

from experiment import BOUND, REFS, SEED, calibrated, clip, dump, geometry, guard_table, history_edges, manifest
from settings import POWER

# Exact arithmetic modules supplied by the read-only repository environment.
site.addsitedir("/Volumes/git/factoring-connection/.venv/lib/python3.12/site-packages")
from cypari2 import Pari
from generate import LOCAL, polynomial

clock = time.perf_counter


def local_inputs(rows, m, coeff):
    base=np.zeros((len(rows),24), dtype=np.float32)
    cof=[]
    degree=len(coeff)-1
    for j,(a,b,_,_) in enumerate(rows.tolist()):
        norms=(abs(a-m*b), abs(sum(c*a**k*b**(degree-k) for k,c in enumerate(coeff))))
        pair=[]
        for s,n in enumerate(norms):
            for k,p in enumerate(LOCAL):
                if n%p==0:
                    e=0
                    while n%p==0:
                        e+=1
                        n//=p
                    base[j,2+s*11+k]=e
            pair.append(n)
        cof.append(pair)
    base[:,:2]=np.log(np.asarray(cof,dtype=np.float64))/np.log(BOUND)
    return base,cof


def confirm(c, pari):
    fac=pari(c).factor() if POWER else pari(c).factor(BOUND)
    return all(int(p)<=BOUND for p in fac[0])


def main():
    import math
    start,cpu=clock(),time.process_time()
    pari=Pari()
    root=Path("results/true")
    known=json.loads((root/"heldout_summary.json").read_text())
    selected,matched=known["locked"],known["matched"]
    size,params,models=joblib.load(root/"models.joblib")
    with open(root/"yield_ledger.csv",newline="") as f:
        expected={(int(r["N_id"]),r["geometry"]):r for r in csv.DictReader(f) if r["u_bin"]=="-1"}
    has_history=selected.endswith("history")
    table_cost=0
    if has_history:
        t0=clock()
        table=guard_table((31*224+15*512)*768 + (31*224+15*768)*512)
        table_cost=clock()-t0
    rows=[]
    for item in manifest():
        if item["split"]!="test": continue
        k,m,coeff=item["id"],item["m"],item["coefficients_ascending"]
        t0=clock()
        residue=(m+15)%31-15
        coords=[]
        for g in (0,1):
            for b in range(1+256*g,257+256*g):
                if ((b-1)//8)%4!=3: continue
                for i in range(-224,224):
                    a=31*i+residue*b if g else i
                    if math.gcd(a,b)==1:
                        coords.append((a,b,i,g))
        coords=np.array(coords,dtype=np.int32)
        base,cof=local_inputs(coords,m,coeff)
        common_local=clock()-t0
        t0=clock()
        raw=size.predict(base)
        p0=calibrated(raw,params)[:,0]
        common_score=clock()-t0
        t0=clock()
        xb=np.column_stack((base,logit(clip(raw)).astype(np.float32)))
        flexible_setup=clock()-t0
        t0=clock()
        xg=np.column_stack((xb,geometry(coords)))
        geometry_setup=clock()-t0
        history_setup=0
        if has_history:
            # Recover the seeded reference stream's initial state; the polynomial is already public.
            rng=random.Random(SEED+1009*k)
            p,q,m2,c2=polynomial(item["degree"],item["bits"],rng,pari)
            assert m2==m and c2==coeff
            t0=clock()
            references, seen=[],set()
            block_size=1024 if POWER else REFS
            for _ in range(REFS//block_size):
                block=[]
                for j in range(block_size):
                    g=j%2
                    while True:
                        b,i=rng.randrange(513,769),rng.randrange(-224,224)
                        a=31*i+residue*b if g else i
                        if math.gcd(a,b)==1 and (a,b) not in seen:
                            seen.add((a,b))
                            block.append((a,b,i,g))
                            break
                rng.shuffle(block)
                references.extend(block)
            references=np.array(references,dtype=np.int32)
            _,co_ref=local_inputs(references,m,coeff)
            ys=np.array([[confirm(c,pari) for c in pair] for pair in co_ref],dtype=np.uint8)
            yref=np.column_stack((ys.prod(axis=1),ys))
            edge=history_edges(np.vstack([references,coords]),table)
            vals=yref[edge]
            xg=np.column_stack((xg,vals[:,:4].mean(axis=1).astype(np.float32),vals.mean(axis=1).astype(np.float32)))
            history_setup=clock()-t0
        names=["logistic",matched,selected]
        rotation=k%3
        for name in names[rotation:]+names[:rotation]:
            t0=clock()
            pred=p0 if name=="logistic" else models[name].joint(xb if name==matched else xg,p0)
            ix=[]
            for g in (0,1):
                take=np.flatnonzero(coords[:,3]==g)
                ix.append(take[np.argsort(-pred[take],kind="stable")[:len(take)//5]])
            scoring=clock()-t0
            t0=clock()
            found=[sum(confirm(cof[j][0],pari) and confirm(cof[j][1],pari) for j in group) for group in ix]
            confirmation=clock()-t0
            for j,g in enumerate(("box","special_q")):
                assert found[j]==int(expected[(k,g)][name]), (k,name,g,found[j],expected[(k,g)][name])
            extra=(flexible_setup if name!="logistic" else 0)+(geometry_setup+history_setup if name==selected else 0)
            rows.append(dict(N_id=k,degree=item["degree"],policy=name,candidates=len(coords),checks=sum(map(len,ix)),
                relations=sum(found),norms_local_s=common_local,logistic_s=common_score,extra_inputs_s=extra,
                rank_score_s=scoring,confirm_s=confirmation,seconds=common_local+common_score+extra+scoring+confirmation))
    with open(root/"work_seconds.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    totals={}
    for name in ["logistic",matched,selected]:
        subset=[r for r in rows if r["policy"]==name]
        secs=sum(r["seconds"] for r in subset)+(table_cost if name==selected else 0)
        found=sum(r["relations"] for r in subset)
        totals[name]=dict(seconds=secs,relations=found,relations_per_second=found/secs)
    dump(root/"work_summary.json",dict(policies=totals,guard_table_seconds=table_cost,
          elapsed=clock()-start,cpu=time.process_time()-cpu,includes_training=False))


if __name__=="__main__":
    with threadpool_limits(limits=1,user_api="blas"),threadpool_limits(limits=1,user_api="openmp"):
        main()
