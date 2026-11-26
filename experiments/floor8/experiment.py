"""Exact fixed-form semiprime experiment. Run with the repo venv and -B."""
from collections import Counter, defaultdict
from itertools import product, combinations
from functools import lru_cache
from math import gcd, isqrt
from pathlib import Path
from random import Random
import json
import sys

import sympy as sp
from cypari2 import Pari

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from src.classE.fd import floor_sums

gp = Pari()
cubic = gp("x^3-x+1")
FORMS = {-23: (1, 6), -4: (0, 1)}


def prime(rng, bits):
    while True:
        p = rng.getrandbits(bits) | (1 << (bits-1)) | 1
        if sp.isprime(p):
            return p


@lru_cache(maxsize=None)
def local(p, D):
    b, c = FORMS[D]
    root = sp.sqrt_mod(D, p)
    if root is None:
        return [(1, p, 1, 0)]
    # y = slope*x mod p; subtract the origin counted twice.
    inv = pow(2*c, -1, p)
    return [(1, 1, p, ((-b + r)*inv) % p) for r in (root, -root)] + [(-1, p, 1, 0)]


def lattices(p, q, D):
    result = []
    for (sg, d, m, a), (sh, e, n, b) in product(local(p, D), local(q, D)):
        slope = a + m * ((b-a)*pow(m, -1, n) % n)
        result.append((sg*sh, d*e, m*n, slope))
    return result


def lattice_sum(terms, X, weighted=True):
    S = W = 0
    for sign, d, m, slope in terms:
        Y = (X+d-1)//d
        # Count y below Y in the residue class slope*j, via two floors.
        a = (-slope) % m
        hi = floor_sums(Y, a, Y-1, m)
        lo = floor_sums(Y, a, m-1, m)
        count = hi[0]-lo[0]+Y
        S += sign*count
        if weighted:
            moment = hi[1]-lo[1]+Y*(Y-1)//2
            odd = Y//2
            odd_count = (floor_sums(odd, 2*a, a+Y-1, m)[0]
                         - floor_sums(odd, 2*a, a+m-1, m)[0] + odd)
            W += sign*((d*moment-odd_count)//2)
    return S, W


def representation(p, q):
    types = [len(gp.factormod(cubic, r)[0]) for r in (p, q)]
    if types == [1, 1]:
        return 4
    return 8 if types == [3, 3] else 0


def small_checks():
    for p, q in combinations([5, 7, 11, 13, 17, 19], 2):
        n = p*q
        for D, (b, c) in FORMS.items():
            terms = lattices(p, q, D)
            for X in [1, 3, isqrt(n)+1, n, n+isqrt(n)+1]:
                xs = [x for x in range(X) for y in range(X) if (x*x+b*x*y+c*y*y)%n == 0]
                assert lattice_sum(terms, X) == (len(xs), sum(x//2 for x in xs))
        lim = isqrt(4*n//23)
        R = 0
        for t in range(-lim, lim+1):
            value = 4*n-23*t*t
            z = isqrt(value)
            if z*z == value and (z-t)%2 == 0:
                R += 1 + (z != 0)
        assert representation(p, q) == R


def draw(p, q, D, R):
    n = p*q
    k = isqrt(n)
    terms = lattices(p, q, D)
    signs = [1 if len(local(v,D))==3 else -1 for v in (p,q)]
    jac = int(sp.jacobi_symbol(D,n))
    assert jac == signs[0]*signs[1]
    r, w = lattice_sum(terms, n)
    s, z = lattice_sum(terms, k+1)
    P, Q = 1+signs[0], 1+signs[1]
    assert r == (1+P*(p-1))*(1+Q*(q-1))
    assert r%8 == 1+4*(sum(t==1 and v%4==3 for t,v in zip(signs,(p,q)))%2)
    assert w == (r-1)*(n-1)//4

    # Going from [0,j+1) to [0,n+j+1) adds a period and two strips.
    def prefix(j):
        short = lattice_sum(terms, j+1, False)[0]
        long = lattice_sum(terms, n+j+1, False)[0]
        return (long-short-r-2)%8

    expected = 2*(k*P*Q + (1-P)*Q*(k//p) + (1-Q)*P*(k//q))
    assert prefix(k) == expected%8
    recovered = None
    if P == Q == 2:
        assert prefix(k) == 4
        lo, hi = 0, k
        while hi-lo>1:
            mid = (hi+lo)//2
            if prefix(mid) == 0:
                lo = mid
            else:
                hi = mid
        recovered = gcd(hi, n)
        assert recovered == min(p,q)

    # Classical reciprocity at a=2, since s(n,2)=0.
    K = (n-1)*(n-5)//2
    return dict(D=D,p=p,q=q,n=n,R=R,jac=jac,K=K%8, full=r%8, floor=w%8,
                short=s%8, shortfloor=z%8,factor=recovered)


def generate():
    with open('records.jsonl','w') as out:
        for split, seed, size in [('train',823040,256),('hold',824040,512)]:
            for bits in (40,60,80):
                rng = Random(seed+bits)
                for i in range(size):
                    while True:
                        p,q = prime(rng,bits//2), prime(rng,bits//2)
                        if p != q and (p*q).bit_length()==bits:
                            break
                    R = representation(p,q)
                    for D in FORMS:
                        row = draw(p,q,D,R)
                        row.update(split=split,bits=bits)
                        out.write(json.dumps(row)+'\n')
                print(split,bits,flush=True)


def predict(train, test, target, richer):
    tables = [defaultdict(Counter) for _ in range(3)]
    def keys(r):
        return [(), (r['jac'],), (r['jac'],r['n']%8,r['K'])]
    for row in train:
        for table,key in zip(tables,keys(row)):
            table[key][row[target]] += 1
    hits,total = Counter(),Counter()
    for row in test:
        for i in range(2 if richer else 1,-1,-1):
            counts = tables[i][keys(row)[i]]
            if counts:
                pred = min(counts,key=lambda v:(-counts[v],v))
                break
        true = row[target]
        total[true] += 1
        hits[true] += pred == true
    return sum(hits.values()), {v:[hits[v],total[v]] for v in sorted(total)}


def analyze():
    rows = [json.loads(line) for line in Path('records.jsonl').read_text().splitlines()]
    summary = {}
    for D in FORMS:
        train = [r for r in rows if r['D']==D and r['split']=='train']
        test = [r for r in rows if r['D']==D and r['split']=='hold']
        stats = {}
        for target in ['full','floor','short','shortfloor']:
            direct = [sum(r[t]%8 == r[target] for r in test) for t in ('jac','K')]
            accuracy = [predict(train,test,target,rich)[0] for rich in (False,True)]
            table = {R:[sum(r['R']==R and r[target]==v for r in test) for v in range(8)] for R in [0,4,8]}
            stats[target] = dict(direct=direct,accuracy=accuracy,contingency=table)
        stats['recovered'] = sum(r['factor'] is not None for r in test)
        if D == -23:
            for sample in (train,test):
                for r in sample:
                    r['plant'] = (r['n']%4==3)
            train,test = ([r for r in sample if r['R'] in (4,8)] for sample in (train,test))
            for target in ('R','plant'):
                stats[target] = [predict(train,test,target,rich) for rich in (False,True)]
        summary[D]=stats
    Path('summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))

if __name__ == '__main__':
    small_checks()
    generate()
    analyze()
