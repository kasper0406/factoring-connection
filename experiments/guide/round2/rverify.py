"""Independent recomputation of contrasts, global control checks and match balance."""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t as T
from rdata import js,YS,ev
from rstats import holm


def main():
    info=json.loads(Path('data/manifest.json').read_text())
    a=json.loads(Path('results/audit.json').read_text())
    locks=json.loads(Path('results/lock_provenance.json').read_text())
    assert hashlib.sha256(Path('PREREGISTRATION.md').read_bytes()).hexdigest()==locks['prereg']
    assert hashlib.sha256(Path('results/validation_lock.json').read_bytes()).hexdigest()==locks['validation_sha256']
    rows=list(csv.DictReader(open('results/A_families.csv')))
    recalced=[]
    for l,cut in enumerate(YS):
        path=Path(f'results/A_{l}')
        with np.load(path/'test_data.npz') as z:d=dict(z)
        assert len(np.unique(d['n']))==20
        assert all(info[k]['split']=='test' for k in np.unique(d['n']))
        for fam in [f'G{i}' for i in range(1,9)]+['all']:
            with np.load(path/f'test_{fam}.npz') as z:pr=z['pred'];dom=z['domain']
            assert np.isfinite(pr).all()
            for h,head in enumerate(['joint','rational','algebraic','W0','W1','cross']):
                rr=next(r for r in rows if r['ycut']==str(cut) and r['family']==fam and r['head']==head)
                ix=dom&d['ok'][:,h];y=d['y'][ix,h].astype(float);p0=d['p'][ix,h].astype(float);pa=pr[ix,h].astype(float)
                if h<3:
                    p0=np.clip(p0,1e-6,1-1e-6);pa=np.clip(pa,1e-6,1-1e-6)
                    v=np.where(y,np.log(pa/p0),np.log((1-pa)/(1-p0)))
                else:v=(y-p0)**2-(y-pa)**2
                nn=d['n'][ix];bb=d['block'][ix];ds=[]
                for k in np.unique(nn):
                    sel=nn==k
                    ds.append(np.mean([v[sel&(bb==j)].mean() for j in np.unique(bb[sel])]))
                mean=np.mean(ds);se=np.std(ds,ddof=1)/np.sqrt(len(ds))
                np.testing.assert_allclose(mean,float(rr['gain']),rtol=3e-5,atol=1e-10)
                assert int(rr['rows'])==len(y)
                recalced.append((rr,T.sf(mean/se,len(ds)-1) if se else float(mean<=0)))
    adjusted=holm([p for _,p in recalced]);np.testing.assert_allclose(adjusted,[float(r['Holm']) for r,_ in recalced],atol=2e-6)
    # Exact global coordinate ownership and fresh modulus identities for the additional arithmetic control.
    original={it['N'] for f in ['../data/manifest.json','../power_cohort/data/manifest.json','data/manifest.json'] for it in json.loads(Path(f).read_text())}
    extra={}
    for i in range(96):
        with np.load(f'data/pc2b_power_{i}.npz') as z:
            n=str(z['N']);m=int(str(z['m']));co=list(map(int,z['F'])); assert ev(co,m*m,1)%int(n)==0
            assert n not in original
            if i//4 in extra:assert n==extra[i//4]
            extra[i//4]=n
    assert len(set(extra.values()))==24
    summaries={}
    for split in ['fit','val','test']:
        res=[]
        for it in info:
            if it['split']!=split:continue
            by={}
            for r in it['parents']:by.setdefault(r['pair'],{})[r['seed']]=r
            for pair in by.values():
                s,c=pair[1],pair[0]
                assert s['band']==c['band'] and s['q']!=c['q']
                res.append(np.array(s['signature'])-np.array(c['signature']))
        v=np.array(res);summaries[split]=dict(pairs=len(v),mean_difference=v.mean(0).tolist(),
            median_abs_difference=np.median(abs(v),axis=0).tolist(),p90_abs_difference=np.quantile(abs(v),.90,axis=0).tolist())
    js('results/match_balance.json',summaries)
    js('results/verification.json',dict(verified_contrasts=len(recalced),audited_nodes=a['nodes'],unique_primary_N=len(info),
       unique_additional_control_N=len(set(extra.values())),primary_targets=a['targets'],test_N=20,validation_lock_digest=locks['validation_sha256']))


if __name__=='__main__':main()
