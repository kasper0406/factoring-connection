"""Independent/post-registration checks for the half square and denominators."""
from collections import Counter
from fractions import Fraction
from math import gcd
import json
from pathlib import Path
from experiment import lattices, lattice_sum, local, sp

def check_parity():
    counts = Counter()
    for line in Path('records.jsonl').read_text().splitlines():
        r = json.loads(line)
        if r['D'] != -23 or r['split'] != 'hold' or r['jac'] != -1:
            continue
        p,q,n = r['p'],r['q'],r['n']
        if len(local(p,-23)) == 3:
            p,q = q,p
        terms = lattices(p,q,-23)
        def count(k):
            return lattice_sum(terms,k,False)[0]%2
        lo,hi = 1,(n+1)//2
        sign = int(sp.jacobi_symbol(6,q))
        assert count(hi) == (sign==1)
        counts[sign] += 1
        if sign == -1:
            while hi-lo>1:
                mid = (hi+lo)//2
                if count(mid):
                    lo = mid
                else:
                    hi = mid
            assert gcd(hi-1,n) == p
    return counts

def check_denominators():
    count = 0
    for k in range(1,221):
        for a in range(1,k+1):
            if gcd(a,k)>1:
                continue
            numer = sum((2*i-k)*(2*((a*i)%k)-k) for i in range(1,k))
            s = Fraction(numer,4*k*k)
            if k%2:
                den2 = 2 if k%4==3 else 1
                assert (12*k*s-k-1+2*int(sp.jacobi_symbol(a,k)))%8==0
            elif k%4==2:
                den2 = 2 if k%8==6 else 1
            else:
                den2 = 2*(k & -k)
            assert s.denominator & -s.denominator == den2
            count += 1
    return count

def check_weighted_strip():
    from math import isqrt
    nt = 0
    for line in Path('records.jsonl').read_text().splitlines():
        r = json.loads(line)
        if r['split'] != 'hold':
            continue
        p,q,n,D = (r[key] for key in ('p','q','n','D'))
        P,Q = (2*int(len(local(v,D))==3) for v in (p,q))
        terms = lattices(p,q,D)
        k = isqrt(n)
        def delta(x):
            return sum(c*lattice_sum(terms,x+j*2*n)[1] for c,j in ((1,2),(-2,1),(1,0)))
        H = P*Q*k + (1-P)*Q*(k//p) + (1-Q)*P*(k//q)
        assert delta(k+1)-delta(1)==6*n*H
        if D == -4:
            assert r['short']==1+4*(P==Q==2)
        if P==Q==2:
            assert r['full']==1+4*(n%4==3) and r['floor']==2*(n%4==3)
        nt += 1
    return nt

def check_folded_factor():
    counts = Counter()
    for line in Path('records.jsonl').read_text().splitlines():
        r = json.loads(line)
        if r['split'] != 'hold':
            continue
        n,D = r['n'],r['D']
        terms = lattices(r['p'],r['q'],D)
        def prefix(k):
            small = lattice_sum(terms,k+1,False)[0]
            comp = lattice_sum(terms,n-k,False)[0]
            return (r['full']-1+small-comp)%8
        factor = None
        for j in range(1,n.bit_length()):
            hi = n>>j
            if prefix(hi):
                lo = 0
                while hi-lo>1:
                    mid = (hi+lo)//2
                    if prefix(mid):
                        hi = mid
                    else:
                        lo = mid
                factor = gcd(hi,n)
                assert factor in (r['p'],r['q'])
                break
        active = any(len(local(v,D))==3 for v in (r['p'],r['q']))
        assert (factor is not None)==active
        counts[D,active] += 1
    return counts

if __name__ == '__main__':
    print(check_parity(), check_denominators(), check_weighted_strip(), check_folded_factor())
