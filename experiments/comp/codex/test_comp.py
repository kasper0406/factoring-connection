"""Run with the repository Python -B, from any directory."""
import csv
import ctypes
import itertools
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from flint import acb_poly, arb, ctx, fmpz_mat, fmpz_poly

HERE=Path(__file__).resolve().parent
DATA=Path('/Volumes/git/factoring-connection/experiments/dd/codex')
SCALE=2**48
DEC=10**9
ctx.prec=384
LIB=ctypes.CDLL(str(HERE/'rank.so'))
small=ctypes.c_int
word=ctypes.c_int64
LIB.ranks.argtypes=[small,small,*[ctypes.POINTER(word)]*2,*[ctypes.POINTER(small)]*3]
LIB.ranks.restype=None

def rank(rows):
    if not rows:
        return 0
    expanded=[]
    for row in rows:
        expanded.extend((row,[v for a,b in zip(row[::2],row[1::2]) for v in (2*b,a)]))
    return fmpz_mat(expanded).rank()//2

def endpoint(v,upper=False):
    return int((v.upper().ceil() if upper else v.lower().floor()).unique_fmpz())

def gridpoint(v):
    x=int((v.mid()+arb(1)/2).floor().unique_fmpz())
    assert abs(v-x)<1
    return x

def geometry(cell,records):
    m=cell['m']
    rr=list(map(fmpz_poly,cell['R']))
    h=arb(0)
    loss_sum=arb(0)
    high_loss=0
    sites=[]
    for side in 'PQ':
        X=fmpz_poly(cell[side])
        gs=[fmpz_poly([1]),fmpz_poly([0,1]),*[X//r for r in rr]]
        pol=acb_poly(X.coeffs())
        roots=(X*X-2).complex_roots()
        assert len(roots)==2*cell['n'] and all(k==1 for _,k in roots)
        for j,(z,_) in enumerate(roots):
            r=records[side,j]
            for part,key in [(z.real,'floor_re_1e12'),(z.imag,'floor_im_1e12')]:
                assert endpoint(part*10**12)==int(r[key])
            sign=int(r['sqrt_sign'])
            s=sign*arb(2).sqrt()
            assert pol(z).overlaps(s) and sign*pol(z).real>0
            row=[acb_poly(g.coeffs())(z) for g in gs]
            M=arb(1)
            for v in row:
                M=M.max(abs(v))
            h+=M.log()
            lo,hi=[int(r['loss_'+k+'_1e9']) for k in ('lo','hi')]
            ell=arb(lo+hi,hi-lo)/(2*DEC)
            loss_sum+=ell
            cut=arb(m)*arb(m).log()/4
            assert ell>=cut or ell<cut
            high_loss+=int(ell>=cut)
            if z.imag>=0:
                coords=[gridpoint(component*SCALE) for v in row for u in (v/M,s*v/M) for component in (u.real,u.imag)]
                assert max(map(abs,coords))<2**49
                sites.append((coords,ell))
    return sites,h,h-loss_sum,high_loss

def intersect(m,sites,w):
    data=[c for row,_ in sites for c in row]
    bd=[]
    for _,ell in sites:
        radius=(w-ell).exp()*SCALE
        error=2*m*arb(2).sqrt()
        l=min(2**56,endpoint(radius-error))
        h=min(2**56,endpoint(radius+error,True))
        assert l>0
        bd.extend((l,h))
    mat=small*(2*m*m)
    lo,hi,out=mat(),mat(),(small*2)()
    LIB.ranks(m,len(sites),(word*len(data))(*data),(word*len(bd))(*bd),lo,hi,out)
    certs=[]
    for a,r in [(lo,out[0]),(hi,out[1])]:
        rows=[list(a[j*2*m:(j+1)*2*m]) for j in range(r)]
        assert rank(rows)==r
        if a is lo:
            for row in rows:
                for (c,_),bound in zip(sites,bd[::2]):
                    dist=sum(sum(x*y for x,y in zip(c[k::2],row))**2 for k in (0,1))
                    assert dist<=bound**2
        def has(j):
            e=[0]*(2*m)
            e[2*j]=1
            return rank(rows+[e])==r
        certs.append((r,has(0),has(1)))
    return certs

def decide(a,b):
    assert a>b or a<b or a==b
    return int(a>=b)

def test(job):
    c,N,records=job
    m,n=c['m'],c['n']
    sites,H,Q,count=geometry(c,records)
    x=arb(N).log()
    X=x/(2*((2*n+m-1)//m)+1)
    D=4*n
    P,B=[fmpz_poly(c[s]) for s in 'PQ']
    norms={str(i):abs(int((P(i)**2-2)*(B(i)**2-2))) for i in [-1,0,1] if math.gcd(c['b']-i,N)==1}
    best=min(norms.values())
    lognorm=arb(best).log()
    Ylog=(x*x.log()).sqrt()
    out={k:c[k] for k in ('bits','m','p','n')}
    out.update(high=count,D=D,prem_height=decide(4*X,H),prem_count=int(4*count>=D),
               Q_at_quarter=decide(Q,X/4),Q_nonpos=decide(arb(0),Q),
               affine_min=best,affine_norms=norms,affine_below_X=1-decide(lognorm,X),affine_capped=decide(Ylog,lognorm))
    for label,val in [('H',H),('Q',Q),('H_over_X',H/X),('Q_over_X',Q/X),('X',X)]:
        out[label]=[endpoint(val*DEC),endpoint(val*DEC,True)]
    for label,w in [('half',arb(m)/2),('one',arb(m)),('two',arb(2*m))]:
        out[label]=intersect(m,sites,w)
    return out

def main():
    cells=[json.loads(s) for s in (DATA/'cells.jsonl').read_text().splitlines()]
    fx=json.loads((DATA/'fixtures.json').read_text())
    grouped={}
    for r in csv.DictReader((DATA/'embeddings.csv').open()):
        k=tuple(int(r[c]) for c in ('bits','m','p'))
        grouped.setdefault(k,{})[r['side'],int(r['root'])]=r
    assert len(cells)==len(grouped)==195 and sum(map(len,grouped.values()))==9096
    jobs=[(c,fx[str(c['bits'])]['N'],grouped[c['bits'],c['m'],c['p']]) for c in cells]
    with ProcessPoolExecutor(4) as pool, (HERE/'cells_comp.jsonl').open('w') as f:
        for r in pool.map(test,jobs):
            f.write(json.dumps(r)+'\n')
            f.flush()

if __name__=='__main__':
    main()
