"""Held-out schedule. Set GRDEG_REL_DEGREE in a fresh process per size."""
import json
import sys
from itertools import combinations_with_replacement
from math import prod
from flint import fmpz_mat
from sympy import primerange, integer_nthroot
from experiment import (n,d,T,one,here,coeffs,fixture,select,measure,Ideals,
                        coordinate_columns)
from sparse import choose


def deficient(N,b,P,Q,m):
    A = [one,T]+[T**i-b**i for i in range(2,m)]
    B = [one,T]+[one*0]*(m-2)
    f,g = P*P-2,Q*Q-2
    cols = [(one,one)]
    # Keep a Q(sqrt(2))-basis after every multiplication for the span test.
    # Compute the lattice index from those selected original monomials.
    from experiment import pivots
    cap = 2*((2*n+m-1)//m)+1
    for t in range(1,cap+1):
        pairs = [(x*u%f,y*v%g) for x,y in cols for u,v in zip(A,B)]
        mats = []
        for x,y in pairs:
            for v,w in ((x,y),(x*P%f,y*Q%g)):
                mats.append(coordinate_columns(v,d)+coordinate_columns(w,d))
        ix = pivots(mats)
        assert all(i%2==0 and j==i+1 for i,j in zip(ix[::2],ix[1::2]))
        cols = [pairs[i//2] for i in ix[::2]]
        M = fmpz_mat([mats[i] for i in ix])
        rright = fmpz_mat([mats[i][d:] for i in ix]).rank()
        if M.rank()-rright == d:
            snf = M.snf()
            idx = prod(abs(int(snf[j,j])) for j in range(M.rank()))
            assert idx%N == 0 and rright<d
            return dict(m=m,first_degree=t,right_rank=rright//2,selected_index=str(idx),cap=cap)
    raise AssertionError("failed to contain first field")


def run(bits):
    degree,m = {60:(7,3),72:(8,3),84:(9,4),96:(10,4)}[bits]
    assert degree==n
    with (here/f'diagonal_{bits}.jsonl').open('w') as out:
        for seed in (2,3):
            N = fixture(bits,seed)
            a,b,P,Q = select(N)
            fourth,exact = integer_nthroot(bits*((bits-1).bit_length())**3,4)
            bound = 2**(int(fourth)+int(not exact))
            primes = list(primerange(bound+1))
            fields = [Ideals(X*X-2,primes) for X in (P,Q)]
            pdata = primes,prod(primes)
            A = next(A for j,A in choose(P,Q) if j==m)
            for mode in ('power','divided','sparse'):
                rec = measure(N,b,P,Q,m,mode,pdata,fields,(A,A) if mode=='sparse' else None)
                rec.update(B=bound,seed=seed)
                if mode=='sparse': rec['generators']=list(map(coeffs,A))
                if bits==60 and seed==2 and mode=='power':
                    rec['deficient_checks'] = [deficient(N,b,P,Q,r) for r in (4,n)]
                out.write(json.dumps(rec)+'\n')
                out.flush()
                print(bits,seed,m,mode,rec['aggregate']['new'],rec['all']['new'],flush=True)

if __name__=='__main__': run(int(sys.argv[1]))
