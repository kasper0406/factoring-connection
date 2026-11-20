"""Independent audit of saved integer records, forbidden classes, and split membership."""
import hashlib
import json
import math
import time
from pathlib import Path

import gmpy2
import numpy as np
from cypari2 import Pari

from settings import BOUND as B, POWER, REFS
Y = 31
PROD_LOCAL=math.prod([2,3,5,7,11,13,17,19,23,29,31])
SMALL=[2,3,5,7,11,13,17,19,23,29,31]
FORBIDDEN=gmpy2.primorial(B)//PROD_LOCAL
pari=Pari()


def check():
    start,cpu=time.monotonic(),time.process_time()
    infos=json.loads(Path("data/manifest.json").read_text())
    counts=dict(targets=0,references=0,side_norms=0,edges=0,independent_full_factorizations=0)
    assert len(set(it["N"] for it in infos))==len(infos)
    if POWER:
        old=json.loads(Path("../data/manifest.json").read_text())
        assert not {it["N"] for it in infos}.intersection(it["N"] for it in old)
    for it in infos:
        k=it["id"]
        path=Path(f"data/n{k:02}.npz")
        assert hashlib.sha256(path.read_bytes()).hexdigest()==it["sha256"]
        p,q,n,m = (int(it[j]) for j in ("p","q","N","m"))
        deg=it["degree"]
        coeff=it["coefficients_ascending"]
        assert p*q==n and p!=q and pari(p).isprime() and pari(q).isprime()
        assert n.bit_length()==it["bits"] and p.bit_length()==q.bit_length()
        assert sum(c*m**j for j,c in enumerate(coeff))==n
        assert pari.Polrev(coeff).polisirreducible()
        with np.load(path) as z:
            coords=z["coords"]
            word=z["cofactor_words"]
            val=z["valuation"]
            y=z["sides"]
        with np.load(f"data/edges{k:02}.npz") as z:
            edges=z["indices"]
        a,b,i,g=coords.T.astype(np.int64)
        assert np.all(np.gcd(a,b)==1) and np.all((i >= -224) & (i <= 223))
        assert np.all(a==np.where(g,31*i+((m+15)%31-15)*b,i))
        assert np.all((b[:REFS]>=513)&(b[:REFS]<=768))
        assert np.all((b[REFS:]>=1+256*g[REFS:])&(b[REFS:]<=256+256*g[REFS:]))
        assert np.all(((b[REFS:]-1)//8)%4==["fit","cal","select","test"].index(it["split"]))
        assert len(np.unique(coords[:,:2],axis=0))==len(coords)
        assert (edges[:,1:]>edges[:,:-1]).all() and edges.shape==(len(a)-REFS,16)
        det=np.abs(a[REFS:,None]*b[edges]-b[REFS:,None]*a[edges])
        assert det.min()>0
        assert all(gmpy2.gcd(FORBIDDEN,int(d))==1 for d in np.unique(det))
        sample=np.random.default_rng(2307+k)
        ix=set(sample.choice(REFS,16,replace=False).tolist())
        for geo in (0,1):
            ix.update(sample.choice(np.flatnonzero((np.arange(len(g))>=REFS)&(g==geo)),64,replace=False).tolist())
        # Independent homogeneous Horner evaluation, local divisibility, and reference records too.
        for j,(av,bv,_,_) in enumerate(coords.tolist()):
            algebraic=coeff[-1]
            bp=bv
            for c in reversed(coeff[:-1]):
                algebraic=algebraic*av+c*bp
                bp*=bv
            for s,norm in enumerate((abs(m*bv-av),abs(algebraic))):
                co=int(word[j,s,0])+(int(word[j,s,1])<<64)
                assert math.gcd(co,PROD_LOCAL)==1
                assert co*math.prod(pr**int(e) for pr,e in zip(SMALL,val[j,s]) if e)==norm
                if j in ix:
                    factors=pari(co).factor()
                    largest=max((int(t) for t in factors[0]),default=1)
                    assert (largest<=B)==bool(y[j,s])
                    counts["independent_full_factorizations"]+=1
        counts["targets"]+=len(edges)
        counts["references"]+=REFS
        counts["side_norms"]+=2*len(y)
        counts["edges"]+=int(edges.size)
    counts.update(wall=time.monotonic()-start,cpu=time.process_time()-cpu)
    Path("results/integer_audit.json").write_text(json.dumps(counts,indent=2)+"\n")
    print(counts)


if __name__=="__main__":
    check()
