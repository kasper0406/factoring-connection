import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from flint import arb, fmpq_poly, fmpz_poly

from test_comp import DATA, HERE, DEC, endpoint, decide

T=fmpz_poly([0,1])
ONE=fmpz_poly([1])

def norm_at(X,rr,a):
    g=a[0]+a[1]*T+sum((c*(X//v) for c,v in zip(a[2:],rr)),fmpz_poly([]))
    return int((X*X-2).resultant(g))

def norm_pencil(P,E,rr,a,degree_bound):
    values=[norm_at(P+j*E,rr,a) for j in range(degree_bound+1)]
    poly=fmpq_poly([])
    choose=fmpq_poly([1])
    for j in range(degree_bound+1):
        poly+=choose*values[0]
        values=[y-x for x,y in zip(values,values[1:])]
        choose=choose*fmpq_poly([-j,1])/(j+1)
    assert all(c.q==1 for c in poly)
    poly=fmpz_poly([int(c) for c in poly])
    assert all(poly(j)==norm_at(P+j*E,rr,a) for j in [-1,0,1,degree_bound+1])
    return poly

def analyze(job):
    cell,N=job
    P,Q=[fmpz_poly(cell[k]) for k in 'PQ']
    rr=list(map(fmpz_poly,cell['R']))
    m,n=cell['m'],cell['n']
    X=arb(N).log()/(2*((2*n+m-1)//m)+1)
    rows=[]
    for name,a in [('plus',[1]*m),('alternate',[(-1)**j for j in range(m)])]:
        C=norm_pencil(P,Q-P,rr,a,2*(cell['p']+1))
        b=cell['b']
        ev=a[0]+a[1]*b+sum(t*int((P//v)(b)) for t,v in zip(a[2:],rr))
        content,fac=C.factor()
        loss=arb(0)
        mixed=0
        all_roots=0
        for h,e in fac:
            roots=h.complex_roots()
            assert sum(k for _,k in roots)==h.degree()
            costs=[]
            for j in (0,1):
                costs.append(sum((k*arb(0).max(-(2*abs(z-j)).log()) for z,k in roots),arb(0)))
            loss+=e*costs[0].min(costs[1])
            mixed+=int(costs[0]>0 and costs[1]>0)
            all_roots+=e*h.degree()
        assert all_roots==C.degree()
        lc=abs(int(C.leading_coefficient()))
        lead=arb(lc).log()
        bound=arb(lc*abs(int(content))).log()-C.degree()*arb(2).log()-loss
        obs=abs(int(C(0)*C(1)))
        assert (arb(obs).log()-bound).lower()>=0 or (arb(obs).log()-bound).contains(0)
        r={k:cell[k] for k in ['bits','m','p','n']}
        r.update(row=name,unit=int(math.gcd(ev,N)==1),norm=obs,degree=C.degree(),
                 coefficients=list(map(int,C)),mixed_factors=mixed,lead_X=decide(lead,X),bound_X=decide(bound,X))
        for k,v in [('Gamma',loss),('bound',bound),('lead',lead)]:
            r[k]=[endpoint(v*DEC),endpoint(v*DEC,True)]
        rows.append(r)
    return rows

if __name__=='__main__':
    cells=[json.loads(l) for l in (DATA/'cells.jsonl').read_text().splitlines()]
    fx=json.loads((DATA/'fixtures.json').read_text())
    with ProcessPoolExecutor(4) as pool, (HERE/'pencils.jsonl').open('w') as f:
        for rows in pool.map(analyze,[(c,fx[str(c['bits'])]['N']) for c in cells]):
            for r in rows:
                f.write(json.dumps(r)+'\n')
            f.flush()
