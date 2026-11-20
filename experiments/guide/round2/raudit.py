"""Independent exact checks on the records, candidate ownership and every relation edge."""
import concurrent.futures
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from rdata import B,YS,ULTRA,FORBID,gm,ev,js,PARI,canonical
sys.path.insert(0,"/Volumes/git/factoring-connection")
from src.cadobench.smoothcount import smooth_part


def audit(k):
    it=json.loads(Path(f'data/n{k:02}.json').read_text());m=it['m'];c=it['c'];nt=it['nt']
    keep=np.load(f'data/keep{k:02}.npy')
    with np.load(f'data/n{k:02}.npz') as z:
        xy=z['coords'];ed=z['edges']; bid=z['blocks'];s=z['sizes'];w=z['w'];ys=z['y'];flags=z['flags']
    F=lambda a,b:ev(c,a,b)
    assert PARI.Polrev(c).polisirreducible() and ev(c,m,1)==int(it['N'])
    assert math.prod(map(int,it['factors']))==int(it['N']) and all(PARI(int(x)).isprime() for x in it['factors'])
    a,b,q=xy.T;assert np.all(np.gcd(a,b)==1) and np.all(b>0)
    norms=[[],[]]
    for aa,bb,qq in xy.tolist():
        f=F(aa,bb);assert f%qq==0 and abs(f)<ULTRA
        norms[0].append(abs(aa-m*bb));norms[1].append(abs(f)//qq)
    for side,ns in enumerate(norms):
        sb=list(map(int,smooth_part(ns,gm.primorial(B))))
        assert all((v==n)==bool(y) for v,n,y in zip(sb,ns,ys[:,side]))
        for l,Y in enumerate(YS):
            loc=list(map(int,smooth_part(ns,gm.primorial(Y))))
            co=np.array([math.log(n//t) for n,t in zip(ns,loc)])
            want_w=np.array([math.log(v//t)/ss if ss>0 else 1. for v,t,ss in zip(sb,loc,co)])
            np.testing.assert_allclose(s[:,l,side],co,atol=1e-12);np.testing.assert_allclose(w[:,l,side],want_w,atol=1e-12)
    assert int(it['maxnorm'])<ULTRA
    pn=PARI.Polrev(c);x=PARI('x')
    def safe_num(n):return n!=0 and gm.gcd(abs(n),FORBID)==1
    for A,D in it['shifts']:assert safe_num(A-m*D) and safe_num(F(A,D))
    asserted=set()
    for typ,f,h,rg,rf in it['maps']:
        A,L,D=(f,0,h) if not typ else (h,typ*f,h)
        assert rg==A*m+L-D*m
        assert rf==int(pn.polresultant(sum(v*(A*x+L)**j*D**(len(c)-1-j) for j,v in enumerate(c))))
        assert safe_num(rf) and safe_num(rg);asserted.add((typ,f,h))
    sl=[]
    for A,L,C,D,rg,rf in it['transforms']:
        assert rg==A*m+L-m*(C*m+D)
        assert rf==int(pn.polresultant(sum(v*(A*x+L)**j*(C*x+D)**(len(c)-1-j) for j,v in enumerate(c))))
        if safe_num(rf) and safe_num(rg):sl.append((A,L,C,D))
    for j in range(ed.shape[1]):
        ix=np.flatnonzero(ed[:,j]>=0);other=ed[ix,j];assert np.all(other>=nt)
        det=a[ix]*b[other]-a[other]*b[ix]
        assert all(safe_num(int(dd)) for dd in np.unique(np.abs(det)))
        if j in (0,1):assert set(zip((a[other]-a[ix]).tolist(),(b[other]-b[ix]).tolist()))<=set(map(tuple,it['shifts']))
        if j==2:assert np.all(np.abs(det)==1)
        if j==3:
            for aa,bb,A,D in zip(a[ix],b[ix],a[other],b[other]):
                assert any(canonical(t*aa+l*bb,r*aa+h*bb)==(A,D) for t,l,r,h in sl)
        if j in (4,6):
            aa,bb=a[other],b[other];ff,hh=a[ed[ix,j+1]],b[ed[ix,j+1]]
            assert np.all(b[ix]==bb*hh)
            assert np.all(a[ix]==aa*ff) if j==4 else np.all((a[ix]==aa*hh+ff*bb)|(a[ix]==aa*hh-ff*bb))
        if j>=8:
            qidx={r['q']:i for i,r in enumerate(it['parents'])}
            assert all(qidx[int(qq)]<min(int(bl),32) for qq,bl in zip(q[other],bid[ix]))
    for i,r in enumerate(it['parents']):
        qp=r['q'];A,D=r['parent']; assert qp>B and PARI(qp).isprime() and A%qp==r['r']*D%qp
        f=abs(F(A,D));assert f%qp==0
        vals=[abs(A-m*D),f//qp];parts=smooth_part(vals,gm.primorial(B))
        assert int(all(a==b for a,b in zip(vals,parts)))==r['seed']
        ix=np.flatnonzero(bid==i);assert np.all(q[ix]==qp) and np.all((a[ix]-r['r']*b[ix])%qp==0)
        assert all(safe_num(int(v)) for v in np.unique(np.abs(a[ix]*D-A*b[ix])))
        assert not np.any((a[:nt]==A)&(b[:nt]==D))
    return dict(k=k,nodes=len(xy),targets=int(keep.sum()),edge_counts=(ed[keep]>=0).sum(0).tolist(),parents=32)


if __name__=='__main__':
    t=time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(4) as pool: audits=list(pool.map(audit,range(48)))
    xy=[]
    for k in range(48):
        with np.load(f'data/n{k:02}.npz') as z:
            a,b=z['coords'][:len(z['blocks']),:2][np.load(f'data/keep{k:02}.npy')].T
        xy.append(((a+2**30)<<32)+b)
    xy=np.concatenate(xy);assert len(np.unique(xy))==len(xy)
    js('results/audit.json',dict(nodes=sum(r['nodes'] for r in audits),targets=len(xy),parent_links=48*32*352,
       ordinary_links=np.sum([r['edge_counts'] for r in audits],axis=0).tolist(),per_N=audits,seconds=time.monotonic()-t))
