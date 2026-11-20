"""Exact divided-difference diagnostic. Run with the repository's Python and PYTHONPATH."""
from collections import Counter
from functools import reduce
from itertools import product, combinations_with_replacement
from math import gcd, isqrt, prod, log
import json
import os
import random
from pathlib import Path

from cypari2 import Pari
from flint import arb, ctx, fmpq_mat, fmpz_mat, fmpz_poly as Poly
from sympy import integer_nthroot, nextprime, isprime, sqrt_mod, primerange
from src.grdeg.ledger import aggregate_alpha, r4_alpha
from src.rank3 import forms
from src.rankdegen.fh import rows_at

here = Path(__file__).resolve().parent
pari = Pari()
pari.allocatemem(128_000_000, silent=True)
ctx.prec = 192
n = int(os.environ.get("GRDEG_REL_DEGREE", "8"))
d = 2*n
T = Poly([0, 1])
one = Poly([1])


def fixture(bits, seed):
    """Diagnostic only. Generate known primes, without factoring any integer."""
    rng = random.Random(10201*bits + seed)
    scale = 1 << (bits//2)
    while True:
        p = int(nextprime(rng.randrange(7*scale//10, 95*scale//100)))
        if p % 8 not in (1, 7):
            continue
        root = int(sqrt_mod(2, p))
        for a in (root, p-root, root+p, 2*p-root):
            N = a*a - 2
            if N.bit_length() != bits or a % 2 == 0:
                continue
            q = N//p
            if q != p and max(p,q) <= 2*min(p,q) and isprime(q):
                return N


def coeffs(f):
    return [int(c) for c in f]


def irreducible(f):
    factors = f.factor()[1]
    return len(factors) == 1 and factors[0][1] == 1


def select(N, root=None):
    a = isqrt(N+2) if root is None else root
    assert (a*a-2)%N == 0
    base = int(integer_nthroot(a,n)[0])
    for offset in range(base-1):
        b = base-offset
        rem = a-b**n
        dig = []
        for _ in range(n-1):
            digit = (rem+b//2) % b-b//2
            dig.append(digit)
            rem = (rem-digit)//b
        P = Poly(dig+[rem,1])
        assert P(b) == a
        if not irreducible(P*P-2):
            continue
        for j in (1,-1,2,-2,3,-3,4,-4):
            Q = P+j*(T-b)*(T*T+1)*(T-1)
            G = Q*Q-2
            if irreducible(G):
                assert Q(b) == a
                return a,b,P,Q
    raise ValueError("no polynomial pair found")


def generators(P, mode, m):
    if mode == "power":
        return [T**j for j in range(m)]
    return [one,T,P//(T*T+1),P//(T-1)][:m]


def coordinate_columns(f, h):
    return [int(f[i]) for i in range(h)]


def pivots(columns):
    M = fmpq_mat([list(row) for row in zip(*columns)])
    R,rank = M.rref()
    return [next(j for j in range(R.ncols()) if R[i,j]) for i in range(rank)]


def product_check(P,Q,A,B,N):
    f,g = P*P-2,Q*Q-2
    m = len(A)
    history, individual = [], []
    for t in range(1,2*((2*n+m-1)//m+1)+2):
        cols = []
        for mon in combinations_with_replacement(range(m),t):
            v = reduce(lambda u,i: (u*A[i])%f,mon,one)
            w = reduce(lambda u,i: (u*B[i])%g,mon,one)
            for u,z in ((v,w),((v*P)%f,(w*Q)%g)):
                cols.append(coordinate_columns(u,d)+coordinate_columns(z,d))
        ii = pivots(cols)
        assert len(ii)%2 == 0
        history.append(len(ii)//2)
        if t == 2:
            individual = [len(pivots([col[i:i+d] for col in cols]))//2 for i in (0,d)]
        if len(ii) == 2*d:
            det = abs(int(fmpz_mat([cols[i] for i in ii]).det()))
            assert det % N == 0
            return dict(degrees=history, square_dimensions=individual, order_index=str(det))
    raise AssertionError("joint powers failed to span")


def direction(z):
    """Canonical projective coordinates via multiplication by the conjugate of the pivot."""
    a,b = next((a,b) for a,b in zip(z[::2],z[1::2]) if a or b)
    key = []
    for c,e in zip(z[::2],z[1::2]):
        key += [c*a-2*e*b,e*a-c*b]
    common = reduce(gcd,key)
    if a*a-2*b*b < 0: common = -common
    return tuple(j//common for j in key)


def expand(P,f,A):
    return [v for a in A for v in (a,(a*P)%f)]


def element(basis,z):
    return sum((v*s for v,s in zip(basis,z) if s),Poly([]))


def smooth(v, primorial):
    while v != 1:
        part = gcd(v,primorial)
        if part == 1: return False
        v //= part
    return True


class Ideals:
    def __init__(self, P, primes):
        self.poly = pari.Polrev(coeffs(P))
        # Supply every small discriminant prime: the basis is maximal at every tested prime.
        disc = abs(int(P.discriminant()))
        bad = [p for p in primes if disc%p == 0]
        self.field = None
        self.bad = bad or [2]
        self.decomps = {}
        self.results = {}

    def parity(self, a, v, primes):
        key = tuple(coeffs(a))
        if key in self.results: return self.results[key]
        if self.field is None:
            self.field = pari.nfinit([self.poly,self.bad])
        basis = self.field.nfalgtobasis(pari.Polrev(list(key)))
        out = []
        for p in primes:
            if v == 1: break
            if v%p: continue
            exp = 0
            while v%p == 0:
                v //= p
                exp += 1
            if p not in self.decomps:
                self.decomps[p] = self.field.idealprimedec(p)
            total = 0
            for i,pr in enumerate(self.decomps[p]):
                w = int(self.field.idealval(basis,pr))
                total += w*int(pr.pr_get_f())
                if w%2: out.append((p,i))
            assert total == exp
        assert v == 1
        self.results[key] = out
        return out


def rank(rows):
    basis = {}
    for v in rows:
        while v:
            j = v.bit_length()-1
            if j not in basis:
                basis[j] = v
                break
            v ^= basis[j]
    return len(basis)


def kranks(points):
    if not points: return 0
    mat = []
    for z in points:
        mat.append(list(z))
        mat.append([u for a,b in zip(z[::2],z[1::2]) for u in (2*b,a)])
    return fmpz_mat(mat).rank()//2


def measure(N,b,P,Q,m,mode,prime_data,fields,chosen=None):
    primes,primorial = prime_data
    f,g = P*P-2,Q*Q-2
    A,B = chosen if chosen is not None else (generators(P,mode,m),generators(Q,mode,m))
    assert all((a(b)-z(b))%N == 0 for a,z in zip(A,B))
    U,V = expand(P,f,A),expand(Q,g,B)
    relU,relV = list(map(coeffs,A)),list(map(coeffs,B))
    rows,_ = rows_at(coeffs(f),relU,coeffs(g),relV,192)
    h = arb(0)
    for row in rows:
        mx = abs(row[0])
        for val in row[1:]: mx = mx.max(abs(val))
        h += mx.log()/arb(2).log()
    index = product_check(P,Q,A,B,N)
    low = 3**(2*d)*2**((N.bit_length()+d-1)//d)
    high = 3**(2*d)*2**(((2*m-1)*N.bit_length()+2*d-2)//(2*d-1))
    dirs = {}
    allnorms = []
    nonunits = Counter()
    values_mod = [int(u(b))%N for u in U]
    for z in forms.box_points(2*m,1):
        if not any(z): continue
        div = gcd(sum(c*v for c,v in zip(z,values_mod)),N)
        if div != 1:
            nonunits[div] += 1
            continue
        a,w = element(U,z),element(V,z)
        v1,v2 = abs(int(f.resultant(a))),abs(int(g.resultant(w)))
        assert v1 and v2
        if not allnorms:
            raw = forms.F_value(coeffs(f),list(map(coeffs,U)),coeffs(g),list(map(coeffs,V)),z)
            assert tuple(map(abs,raw)) == (v1,v2)
        v = v1*v2
        allnorms.append(v)
        k = direction(z)
        if k not in dirs or v<dirs[k][0]:
            dirs[k] = (v,z,v1,v2)
    colids, hits = {}, []
    for v,z,v1,v2 in dirs.values():
        if not smooth(v,primorial): continue
        mask = 0
        for i,(W,vside) in enumerate(((U,v1),(V,v2))):
            for prime in fields[i].parity(element(W,z),vside,primes):
                lab = (i,prime)
                mask |= 1<<colids.setdefault(lab,len(colids))
        mixed = all(a or c for a,c in zip(z[::2],z[1::2]))
        hits.append((v,z,mask,mixed))
    deleted = [mask for v,z,mask,mix in hits if not mix]
    rdel = rank(deleted)
    def metrics(cap):
        sel = [row for row in hits if cap is None or row[0]<=cap]
        mix = [row for row in sel if row[3]]
        masks = [row[2] for row in sel]
        dim = rank(masks)
        inc = rank(deleted+masks)-rdel
        return dict(points=sum(1 for v in allnorms if cap is None or v<=cap),
                    directions=sum(1 for v,*_ in dirs.values() if cap is None or v<=cap),
                    smooth_directions=len(sel),rank=dim,new=inc,
                    new_within_cap=dim-rank(row[2] for row in sel if not row[3]),
                    mixed_k_span=kranks([row[1] for row in mix]),
                    eta=log(1+inc)/(2*m*log(3)))
    return dict(N=str(N),bits=N.bit_length(),n=n,m=m,mode=mode,b=b,P=coeffs(P),Q=coeffs(Q),
                B=max(primes),log2_H=str(h),span=index,nonunits=dict(nonunits),
                X0=str(low),X1=str(high),aggregate=metrics(low),gp=metrics(high),all=metrics(None),
                median_product=str(sorted(allnorms)[len(allnorms)//2]))


def run():
    assert str(aggregate_alpha()["alpha"])=="1/4" and str(r4_alpha()["alpha"])=="1/3"
    with (here/"results.jsonl").open("w") as out:
        for bits,seed in product((60,72,84,96),range(2)):
            N = fixture(bits,seed)
            a,b,P,Q = select(N)
            fourth,exact = integer_nthroot(bits*((bits-1).bit_length())**3,4)
            smooth_bound = 2**(int(fourth)+int(not exact))
            primes = list(primerange(smooth_bound+1))
            factors = [Ideals(X*X-2,primes) for X in (P,Q)]
            pdata = (primes,prod(primes))
            for m,mode in [(2,"power")]+list(product((3,4),("divided","power"))):
                rec = measure(N,b,P,Q,m,mode,pdata,factors)
                rec.update(B=smooth_bound,seed=seed)
                out.write(json.dumps(rec)+"\n")
                out.flush()
                print(bits,seed,m,mode,rec["aggregate"]["new"],rec["all"]["new"],flush=True)

if __name__ == "__main__":
    run()
