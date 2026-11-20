"""Independent all-row check at the first real and upper-half-plane root on each side."""
import csv
import itertools
import json
from pathlib import Path
from measure import ex, acb_poly, arb, endpoints, HERE

def direct(X,gs,z):
    m=len(gs)
    vals=[acb_poly(ex.coeffs(g))(z) for g in gs]
    base=arb(2).sqrt() if acb_poly(ex.coeffs(X))(z).real>0 else -arb(2).sqrt()
    M=arb(1)
    for v in vals:
        M=M.max(abs(v))
    def half(values):
        ev=[t for v in values for t in (v,base*v)]
        return [(u,sum(v*a for a,v in zip(u,ev))) for u in itertools.product((-1,0,1),repeat=len(ev))]
    L,R=half(vals[:m//2]),half(vals[m//2:])
    ordered=[]
    for u,v in L:
        for w,t in R:
            r=u+w
            if next((i for i in r if i),0)>0:
                d=abs(v+t)
                ordered.append((d.lower(),d.upper(),r))
    ordered.sort(key=lambda e:e[0])
    basis=[]
    upper=arb(0)
    for lo,hi,r in ordered:
        if ex.kranks(basis+[r])>len(basis):
            basis.append(r)
            upper=upper.max(hi)
            if len(basis)==m:
                return endpoints((M/upper).log(),9)[0],endpoints((M/lo).log(),9)[1]
    raise AssertionError('no basis')

def main():
    cells=[json.loads(l) for l in (HERE/'cells.jsonl').read_text().splitlines()]
    with (HERE/'embeddings.csv').open() as f:
        rows={(int(r['bits']),int(r['m']),int(r['p']),r['side'],int(r['root'])):r for r in csv.DictReader(f)}
    for m in (3,4,5):
        cell=next(c for c in cells if (c['m'],c['p'],c['bits'])==(m,m,60))
        poles=list(map(ex.Poly,cell['R']))
        for side in 'PQ':
            X=ex.Poly(cell[side])
            gs=[ex.one,ex.T,*[X//v for v in poles]]
            seen=set()
            for i,(z,_) in enumerate((X*X-2).complex_roots()):
                kind='real' if z.imag==0 else 'complex'
                if kind in seen or z.imag<0:
                    continue
                seen.add(kind)
                lo,hi=direct(X,gs,z)
                r=rows[60,m,m,side,i]
                assert lo<=int(r['loss_hi_1e9']) and hi>=int(r['loss_lo_1e9'])
                print(m,side,i,kind,lo,hi,flush=True)
            assert seen=={'real','complex'}
if __name__=='__main__':
    main()
