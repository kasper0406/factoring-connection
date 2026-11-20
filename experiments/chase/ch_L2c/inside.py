"""Test the post-sweep binary-coordinate hypothesis."""
from sweep import ex, itertools, json, math, orbit, pair_units, primerange


def audit(record):
    n,m=record["n"],record["m"]
    ex.n,ex.d=n,2*n
    N=int(record["N"])
    _,b,P,Q=ex.select(N)
    fs=[X*X-2 for X in (P,Q)]
    gs=[ex.expand(X,f,ex.generators(X,"divided",m)) for X,f in zip((P,Q),fs)]
    residues=[int(a(b))%N for a in gs[0]]
    primes=list(primerange(record["B"]+1))
    primorial=math.prod(primes)
    fields=[ex.Ideals(f,primes) for f in fs]
    units=pair_units()
    points={orbit(z,units,1) for z in itertools.product((-1,0,1),repeat=2*m)
            if next((v for v in z if v),0)>0 and
            math.gcd(sum(a*v for a,v in zip(z,residues)),N)==1}
    binary=[]
    other=[]
    outside=capped=0
    cols={}
    for z in points:
        vals=[ex.element(g,z) for g in gs]
        norms=[abs(int(f.resultant(a))) for f,a in zip(fs,vals)]
        F=math.prod(norms)
        if not 0<F<=int(record["cap"]):
            continue
        capped+=1
        outside+=any(z[4:])
        if ex.smooth(F,primorial):
            mask=0
            for j in range(2):
                for pr in fields[j].parity(vals[j],norms[j],primes):
                    mask |= 1<<cols.setdefault((j,pr),len(cols))
            (other if any(z[4:]) else binary).append(mask)
    assert capped==record["capped"] and ex.rank(binary+other)==record["covered_rank"]+record["fresh"]
    return dict(bits=record["bits"],seed=record["seed"],n=n,m=m,capped=capped,
                outside=outside,binary_rank=ex.rank(binary),
                fresh=ex.rank(binary+other)-ex.rank(binary))


if __name__=="__main__":
    for seed in (19,41,83):
        with open(f"seed{seed}.jsonl") as records:
            for line in records:
                r=json.loads(line)
                if r["bits"]>=84 and r["E"]==1:
                    print(json.dumps(audit(r)),flush=True)
