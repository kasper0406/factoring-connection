import json
from collections import Counter
from math import comb, gcd, isqrt
from pathlib import Path
from random import Random
from flint import fmpz_mod_poly_ctx
from gmpy2 import is_prime
from samples import prime
from cubic_carry import carries

def half(x,d):
    return (x%d) > d//2

def selector(N,a,d):
    # a is the lower of two candidate lifts. Test a and b=N/a mod d together.
    inv = pow(a,-1,d)
    c = N*inv*inv%d
    assert c%8 == 3
    end = d//8
    assert pow(c,end,d) == 1+d//2 and end%2 == 0
    def phase(j):
        return half(pow(c,j,d),d) ^ (j%2)
    lo,hi = 0,end
    while hi-lo>1:
        mid = (lo+hi)//2
        if phase(mid) == phase(lo):
            lo = mid
        else:
            hi = mid
    h = pow(c,lo,d)*inv%d
    # Both transformations preserve equality of the two upper-half bits.
    h = min(h,d-h)
    h = min(h,d//2-h)
    assert 0<h<d//4
    direction = half(h*a,d)
    assert direction == half(h*(N*inv%d),d)
    return h,direction

def factor_reduction(N,zero):
    d,a = 8,1
    trace = []
    while 2*d <= isqrt(N//2):
        d *= 2
        h,direction = selector(N,a,d)
        answer = zero(N,h,d)
        if bool(answer) != direction:
            a += d//2
        trace.append((d,h,int(answer),a))
    low = isqrt(N//2)
    high = isqrt(2*N)
    start = low+(a-low)%d
    for s in range(start,high+1,d):
        if s>1 and N%s==0:
            return s,trace
    raise AssertionError('Candidate absent')

def oracle_for(p,q):
    def zero(N,h,d):
        m = h*N//d
        vals = [carries(m,s) for s in (p,q)]
        assert all(v == int(half(h*t,d)) for v,t in zip(vals,[q,p]))
        assert vals[0] == vals[1], ('Unexpected proper binomial gcd',N,p,q,h,d)
        return bool(vals[0])
    return zero

def public_fold(N,h,d,L):
    m=h*N//d
    ring=fmpz_mod_poly_ctx(N)
    X=ring.gen()
    # Reduce the central coefficient into its residue class modulo L.
    val=(1+X).pow_mod(2*m,X**L-1)[m%L]
    return gcd(int(val),N)

def large(n):
    rng = Random(20261126*10000+n*10+7)
    tally = Counter()
    examples = {}
    with Path(f'lift{n}.jsonl').open('w') as out:
        for a in (1,9):
            for _ in range(128):
                while True:
                    p=prime(rng,n//2,a,16)
                    q=prime(rng,n//2,3,8)
                    N=p*q
                    if N.bit_length()==n: break
                oracle=oracle_for(p,q)
                s,trace = factor_reduction(N,oracle)
                assert s == p and all(s%d==b for d,h,ans,b in trace)
                h=1 if N%16==3 else 3
                assert oracle(N,h,16) == (a==9)
                got={}
                for L in (8,16,32):
                    g=public_fold(N,h,16,L)
                    state='unit' if g==1 else ('zero' if g==N else 'proper')
                    tally[f'fold/{L}/{a}/{state}']+=1
                    if (g==N) != (a==9):
                        examples.setdefault(f'{L}/{a}',[N,p,q,h,g])
                    got[L]=g
                tally[f'queries/{len(trace)}']+=1
                out.write(json.dumps(dict(N=N,p=p,q=q,trace=trace,fold_gcd=got),separators=(',',':'))+'\n')
    return dict(counts=dict(tally),fold_counterexamples=examples)

def small():
    # Direct verification of the T1 coefficient and all admissible lifts,
    # including lifts above the generic stopping rule when d < min(p,q).
    tally=Counter()
    for p in range(17,720,8):
        if not is_prime(p): continue
        for q in range(11,720,8):
            if not is_prime(q) or not max(p,q)<2*min(p,q): continue
            N=p*q
            if N>40000: continue
            d=16
            while d<min(p,q):
                a=p%(d//2)
                h,direction=selector(N,a,d)
                m=h*N//d
                g=gcd(comb(2*m,m),N)
                assert g in (1,N) and ((g==N)^direction) == (p%d!=a)
                tally['direct']+=1
                d*=2
            h=1 if N%16==3 else 3
            m=h*N//16
            g=gcd(comb(2*m,m),N)
            assert g==(N if p%16==9 else 1)
            tally['T1_direct']+=1
    return dict(tally)

if __name__=='__main__':
    stats={n:large(n) for n in (40,60,80)}
    stats['small']=small()
    Path('lift_counts.json').write_text(json.dumps(stats,indent=2)+'\n')
