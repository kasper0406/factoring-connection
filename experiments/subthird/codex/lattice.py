"""Preregistered lattice rescue. Operational selection uses a public square root or an N-only fixture."""
from itertools import combinations, product
from math import prod
import json
import random

from flint import arb, acb, fmpz_mat, ctx
from sympy import nextprime, sqrt_mod, primerange, integer_nthroot
from experiment import (n, d, T, one, coeffs, element, expand, select, fixture, measure,
                        here, coordinate_columns, pivots, Ideals, pari)
from src.rank3 import forms


def split_fixture(bits,seed):
    rng = random.Random(133109*bits+seed)
    scale = 1 << (bits//2)
    def split_prime():
        while True:
            p = int(nextprime(rng.randrange(73*scale//100, scale-100)))
            if p%8 in (1,7): return p
    p,q = split_prime(),split_prime()
    while q==p: q=split_prime()
    v,w = (int(sqrt_mod(2,j)) for j in (p,q))
    a = v+p*((w-v)*pow(p,-1,q)%q)
    N = p*q
    assert N.bit_length()==bits and (a*a-2)%N==0
    return N,min(a,N-a)


def order_basis(P):
    return expand(P,P*P-2,[T**j for j in range(n)])


def grid(P,basis):
    ctx.prec = 384
    zs = forms.acb_roots(coeffs(P*P-2))
    # Each embedding is included. Duplicated conjugate pairs give the usual complex weight.
    def rounded(x):
        q = (x*256+arb(1)/2).floor().unique_fmpz()
        if q is None: raise ArithmeticError("ambiguous embedding grid")
        return int(q)
    rows = []
    for u in basis:
        v = [forms._ev(coeffs(u),z) for z in zs]
        rows.append([rounded(c) for w in v for c in (w.real,w.imag)])
    return rows


def square_dimension(P,A):
    f = P*P-2
    cols = []
    for i,u in enumerate(A):
        for v in A[i:]:
            w = u*v%f
            cols.extend(coordinate_columns(a,d) for a in (w,w*P%f))
    return len(pivots(cols))//2


def reduced_pool(N,b,P,Q):
    bases = [order_basis(X) for X in (P,Q)]
    grids = [grid(X,U) for X,U in zip((P,Q),bases)]
    ids = [[int(i==j) for j in range(2*d)] for i in range(2*d)]
    lattice = [row[:] for row in ids]
    lattice[0][0] = N
    metric = []
    for i,(U,Vgrid) in enumerate(zip(bases,grids)):
        for j,u in enumerate(U):
            jj = i*d+j
            if jj:
                residue = int(u(b)) * (1 if i==0 else -1)
                lattice[jj][0] = -((residue+N//2)%N-N//2)
            metric.append((([0]*(2*d) if i else [])+Vgrid[j]+([0]*(2*d) if not i else []))+ids[jj])
    L = fmpz_mat(lattice)
    M = L*fmpz_mat(metric)
    gram = M*M.transpose()
    H = [[int(gram[i,j]) for j in range(2*d)] for i in range(2*d)]
    # PARI works with columns, the change of basis has determinant +/-1.
    C = pari.matrix(2*d,2*d,[v for row in H for v in row]).qflllgram(1)
    change = fmpz_mat([[int(C[j,i]) for j in range(2*d)] for i in range(2*d)])
    assert abs(change.det()) == 1
    reduced = (change*M).tolist()
    rows = sorted([tuple(map(int,r)) for r in reduced],key=lambda r:sum(x*x for x in r))
    def parse(row):
        c = row[-2*d:]
        u,v = element(bases[0],c[:d]),element(bases[1],c[d:])
        assert (u(b)-v(b))%N==0
        return u,v
    extra = [tuple(a+s*b for a,b in zip(u,v)) for u,v in combinations(rows[:16],2) for s in (1,-1)]
    extra.sort(key=lambda r:sum(x*x for x in r))
    return [parse(r) for r in rows+extra]


def tuples(N,b,P,Q):
    A,B = [one,T],[one,T]
    pool = reduced_pool(N,b,P,Q)
    for m in (3,4):
        for i,(u,v) in enumerate(pool):
            if u==v or not u or not v: continue
            if square_dimension(P,A+[u])<=2*m-1: continue
            if square_dimension(Q,B+[v])<=2*m-1: continue
            A,B = A+[u],B+[v]
            yield m,(A,B),i
            break
        else: return


def run():
    with (here/'lll_results.jsonl').open('w') as outf:
        for kind,bits,seed in product(('special','advice'),(60,72,84,96),range(2)):
            if kind == 'special':
                N,root = fixture(bits,seed),None
            else:
                N,root = split_fixture(bits,seed)
            a,b,P,Q = select(N,root)
            fourth,exact = integer_nthroot(bits*((bits-1).bit_length())**3,4)
            smooth_bound = 2**(int(fourth)+int(not exact))
            primes = list(primerange(smooth_bound+1))
            fields = [Ideals(X*X-2,primes) for X in (P,Q)]
            pdata = primes,prod(primes)
            def save(res):
                res.update(B=smooth_bound,seed=seed,kind=kind,root=str(a))
                outf.write(json.dumps(res)+'\n')
                outf.flush()
            if kind=='advice':
                for m in (3,4):
                    save(measure(N,b,P,Q,m,'power',pdata,fields))
            reached = set()
            for m,(A,B),idx in tuples(N,b,P,Q):
                rec = measure(N,b,P,Q,m,'lattice',pdata,fields,(A,B))
                rec.update(row_index=idx,generators=[list(map(coeffs,U)) for U in (A,B)])
                save(rec)
                reached.add(m)
                print(kind,bits,seed,m,idx,rec['aggregate']['new'],rec['all']['new'],flush=True)
            for m in {3,4}-reached:
                save(dict(N=str(N),bits=bits,m=m,mode='lattice',unreached=True))
                print(kind,bits,seed,m,'unreached',flush=True)

if __name__ == '__main__':
    run()
