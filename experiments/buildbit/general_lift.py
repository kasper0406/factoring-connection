import json
from math import comb
from pathlib import Path
from gmpy2 import is_prime
from cubic_carry import carries
from lift import selector

def factor_with_zero(N, zero):
    """Use only the Boolean test and ordinary trial division by lifted candidates."""
    a,d = 1,8
    while d < 2*N:
        half=d
        d*=2
        b=N*pow(a,-1,half)%half
        for v in (a,a+half,b,b+half):
            if 1<v<N and N%v==0:
                return v if v%8==1 else N//v
        # Now both primes exceed d. The larger one always has a zero upper
        # carry when both low digits lie in the lower half.
        h,direction=selector(N,a,d)
        if bool(zero(N,h,d)) != direction:
            a+=half
    raise AssertionError('No factor found')

def validate():
    cases = queries = direct = proper = 0
    unbalanced = 0
    for r in range(17,400,8):
        if not is_prime(r): continue
        for q in range(3,400,8):
            if not is_prime(q): continue
            N=r*q
            def zero(N,h,d):
                nonlocal queries,direct,proper
                queries+=1
                m=h*N//d
                vals=[bool(carries(m,s)) for s in (r,q)]
                proper+=int(vals[0]!=vals[1])
                if N<=12000:
                    direct+=1
                    assert (comb(2*m,m)%N==0) == all(vals)
                return all(vals)
            assert factor_with_zero(N,zero)==r
            cases+=1
            unbalanced+=int(max(r,q)>=2*min(r,q))
    return dict(cases=cases,unbalanced=unbalanced,queries=queries,exact=direct,
                oracle_proper_gcds_ignored=proper)

if __name__=='__main__':
    Path('general_counts.json').write_text(json.dumps(validate(),indent=2)+'\n')
