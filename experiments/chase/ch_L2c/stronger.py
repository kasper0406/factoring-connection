"""The grouped lower-norm comparison specified in stronger_prereg.txt."""
from itertools import groupby
from sympy import factorint
from sweep import ex, itertools, json, math, orbit, pair_units, row_prefix


def compare(bits, seed, n, E):
    ex.n, ex.d = n, 2*n
    N = ex.fixture(bits,seed)
    _, b, P, Q = ex.select(N)
    fs = [X*X-2 for X in (P,Q)]
    gs = [ex.expand(X,f,ex.generators(X,"divided",3)) for X,f in zip((P,Q),fs)]
    residues = [int(g(b))%N for g in gs[0]]
    points = [z for z in itertools.product(range(-E,E+1),repeat=6)
              if next((a for a in z if a),0)>0 and
              math.gcd(sum(a*v for a,v in zip(z,residues)),N)==1]
    rows,_=ex.rows_at(ex.coeffs(fs[0]),list(map(ex.coeffs,gs[0])),
                      ex.coeffs(fs[1]),list(map(ex.coeffs,gs[1])),192)
    spaces=[row_prefix(row,points,3,E)[0] for row in rows]
    cap=3**(4*n)*2**((bits+2*n-1)//(2*n))
    B={60:2048,72:8192}[bits]
    units=pair_units()
    records=[]
    allprimes=set()
    for z in sorted({orbit(z,units,E) for z in points}):
        vals=[ex.element(g,z) for g in gs]
        norms=[abs(int(f.resultant(a))) for f,a in zip(fs,vals)]
        F=math.prod(norms)
        if 0<F<=cap:
            primes=sorted(set(factorint(norms[0])) | set(factorint(norms[1])))
            allprimes.update(primes)
            inside=any(ex.kranks(s+[z])==2 for s in spaces)
            records.append((F,inside,max(primes,default=1)<=B,vals,norms,primes))
    fields=[ex.Ideals(f,sorted(allprimes)) for f in fs]
    cols={}
    masks=[]
    for F,inside,smooth,vals,norms,primes in records:
        mask=0
        for j in range(2):
            for pr in fields[j].parity(vals[j],norms[j],primes):
                mask |= 1<<cols.setdefault((j,pr),len(cols))
        masks.append((F,inside,smooth,mask))
    seen=[mask for _,inside,smooth,mask in masks if inside and smooth]
    groups=[]
    for F,g in groupby(sorted(masks),key=lambda r:r[0]):
        same=list(g)
        esc=[mask for _,inside,smooth,mask in same if smooth and not inside]
        if esc:
            groups.append((F,len(esc),ex.rank(seen+esc)-ex.rank(seen)))
        seen.extend(r[3] for r in same)
    return dict(bits=bits,seed=seed,n=n,m=3,E=E,groups=groups,
                stronger=sum(g[2] for g in groups))


if __name__=="__main__":
    for args in [(60,41,8,2),(60,19,7,1),(72,19,8,2)]:
        print(json.dumps(compare(*args)),flush=True)
