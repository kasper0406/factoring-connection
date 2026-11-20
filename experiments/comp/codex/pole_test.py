import json
import math
from concurrent.futures import ProcessPoolExecutor
from flint import acb_poly,arb,fmpz_poly
from test_comp import DATA,HERE,DEC,endpoint,decide

def cell(job):
    c,old,N=job
    m,n=c['m'],c['n']
    rr=list(map(fmpz_poly,c['R']))
    P,Q=[fmpz_poly(c[k]) for k in 'PQ']
    f=P*P-2
    cross=abs(int(f.resultant(Q-P)))
    assert cross and cross%N==0
    den=math.prod(abs(int(f.resultant(R)))**2 for R in rr)
    Z=arb(den).log()
    H=arb(0)
    B=arb(0)
    check=arb(0)
    tags={}
    for side,X in enumerate([P,Q]):
        generators=[fmpz_poly([1]),fmpz_poly([0,1]),*[X//R for R in rr]]
        gx=[acb_poly(g.coeffs()) for g in generators]
        px=acb_poly(X.coeffs())
        poles=[acb_poly(R.coeffs()) for R in rr]
        for z,mult in (X*X-2).complex_roots():
            assert mult==1 and (px(z).real>0 or px(z).real<0)
            sign=int(px(z).real>0)
            lm=max(abs(g(z)).upper() for g in gx)  # only used as an upper sanity check
            maximum=arb(1)
            for g in gx[1:]: maximum=maximum.max(abs(g(z)))
            assert maximum.lower()<=lm
            lm=maximum.log()
            H+=lm
            logs=[abs(R(z)).log() for R in poles]
            plus=sum(v.max(arb(0)) for v in logs)
            neg=sum((-v).max(arb(0)) for v in logs)
            B+=(neg-lm).max(arb(0))
            check+=plus-neg
            tags.setdefault((side,sign),[]).append((z,plus))
    assert check.overlaps(Z)
    lows=[]
    for sign in (0,1):
        values=[u+v-(m-2)*abs(z-w).log() for z,u in tags[0,sign] for w,v in tags[1,sign]]
        lower=values[0]
        for v in values[1:]: lower=lower.min(v)
        lows.append(lower)
    L=m*arb(m).log()+(2*m-1)*(1+arb(2).sqrt()).log()
    bound=(m-2)*arb(cross).log()/n+n*sum(lows)-Z-B-4*n*L
    naive=H-4*n*L
    X=arb(N).log()/(2*((2*n+m-1)//m)+1)
    assert naive>=bound or (naive-bound).contains(0)
    assert naive<arb(old['Q'][1])/DEC and bound<arb(old['Q'][1])/DEC
    res={k:c[k] for k in ('bits','m','p','n')}
    res.update(split=int(all(sum(v for _,v in R.roots())==R.degree() for R in rr)),
               constant=int(all((P%R).degree()<=0 for R in rr)),
               geometric=decide(bound,X/4),direct=decide(naive,X/4),cross_norm=cross)
    for k,v in [('bound',bound),('direct_bound',naive),('Z',Z),('overlap',B)]:
        res[k]=[endpoint(DEC*v),endpoint(DEC*v,True)]
    return res

if __name__=='__main__':
    cells=[json.loads(l) for l in (DATA/'cells.jsonl').read_text().splitlines()]
    old=[json.loads(l) for l in (HERE/'cells_comp.jsonl').read_text().splitlines()]
    fx=json.loads((DATA/'fixtures.json').read_text())
    assert all((a['bits'],a['m'],a['p'])==(b['bits'],b['m'],b['p']) for a,b in zip(cells,old,strict=True))
    with ProcessPoolExecutor(4) as pool, (HERE/'pole_results.jsonl').open('w') as f:
        results=list(pool.map(cell,[(c,r,fx[str(c['bits'])]['N']) for c,r in zip(cells,old)]))
        for r in results: f.write(json.dumps(r)+'\n')
    print('cells',len(results),'geometric',sum(r['geometric'] for r in results),'direct',sum(r['direct'] for r in results),
          'split',sum(r['split'] for r in results),'constant',sum(r['constant'] for r in results),
          'both',sum(r['constant'] and r['split'] for r in results))
