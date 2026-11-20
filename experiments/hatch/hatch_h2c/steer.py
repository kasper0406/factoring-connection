"""Registered cubic-base experiment. Run with the repository virtual environment."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
from math import gcd, log
from pathlib import Path
from time import perf_counter
import json
import random
import sys

import numpy as np
from flint import fmpz_poly as Poly
from sympy import integer_nthroot, nextprime

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'experiments/subthird/codex'), str(ROOT)]
from experiment import Ideals, coeffs, irreducible, pari, rank

T = Poly([0,1])
H = T**3 + T**2 - 2*T - 1
ONE = (1,0,0)
ZERO = (0,0,0)


def triple(p):
    return tuple(int(p[i]) for i in range(3))


def poly(a):
    return Poly(list(a))


def dictionary():
    units = [T, T*T-2]
    powers = []
    for u in units:
        # Each is a conjugate of z with inverse u^2+u-2.
        inv = (u*u+u-2) % H
        assert (u*inv) % H == 1
        powers.append([(inv**(-i) if i < 0 else u**i) % H for i in range(-2,3)])
    values = [triple(sign*a*b*(2-T)**s % H) for a,b,s,sign in
              product(*powers, range(2), (1,-1))]
    assert len(set(values)) == len(values) == 100
    assert all(abs(int(H.resultant(poly(a)))) in (1,7) for a in values)
    return values


def fixture(bits, seed):
    rng = random.Random(seed + 10201*bits)
    limit = 1 << (bits//2)
    def prime():
        while True:
            p = int(nextprime(rng.randrange(7*limit//10, 95*limit//100)))
            if p % 7 in (1,6):
                return p
    while True:
        p,q = prime(),prime()
        if p != q and (p*q).bit_length() == bits:
            break
    roots = [min(map(int, pari.polrootsmod(pari.Polrev(coeffs(H)), ell))) for ell in (p,q)]
    a = roots[0] + p*((roots[1]-roots[0])*pow(p,-1,q) % q)
    assert int(H(a)) % (p*q) == 0
    return p*q,a


def fields(a,n):
    for b in range(int(integer_nthroot(a,n)[0]),1,-1):
        rem = a - b**n
        dig = []
        for _ in range(n-1):
            dig.append((rem+b//2)%b-b//2)
            rem = (rem-dig[-1])//b
        p = Poly(dig+[rem,1]); f = H(p)
        if not irreducible(f):
            continue
        for j in (1,-1,2,-2,3,-3,4,-4):
            q = p + j*(T-b); g = H(q)
            if irreducible(g):
                assert p(b) == q(b) == a
                return b,(p,q),(f,g)
    raise ValueError('selector exhausted')


def generators(U, sides, bits, seed):
    arms = defaultdict(set)
    rng = random.Random(600001+bits*100+seed)
    zroots = sorted(np.roots(coeffs(H)[::-1]).real)
    def root_pair(p,z):
        co = np.array(coeffs(p)[::-1], dtype=float); co[-1] -= z
        roots = sorted(np.roots(co), key=lambda r: (round(r.real,8),round(r.imag,8)))
        return roots[0],roots[-1]
    roots = [[root_pair(p,z) for z in zroots] for p in sides]
    ev = [np.array([a+b*z+c*z*z for a,b,c in U]) for z in zroots]
    scores = ties = 0
    per_pivot = []
    for piv in range(3):
        j,k = [i for i in range(3) if i != piv]
        def point(a,b):
            t = [ONE]*3
            t[j],t[k] = a,b
            return tuple(t)
        def loss(base, t):
            mon = (1,t,t**3)
            v,w = mon[j]*ev[base][:,None],mon[k]*ev[base][None,:]
            ratio = np.abs(mon[piv]+v+w)/(abs(mon[piv])+np.abs(v)+np.abs(w))
            return np.log(np.maximum(ratio,1e-15))
        selected = set()
        for i,shift,h in product(range(3),range(2),range(2)):
            i2 = (i+shift)%3
            score = loss(i,roots[0][i][h]) + loss(i2,roots[1][i2][1-h])
            ordered = np.argsort(score,axis=None,kind='stable')
            ties += int(score.flat[ordered[7]] == score.flat[ordered[8]])
            selected.update(map(int,ordered[:8])); scores += score.size
        per_pivot.append(len(selected))
        for arm,indices in [('steer',sorted(selected)),('random',rng.sample(range(len(U)**2),len(selected)))]:
            for idx in indices:
                u,v = divmod(idx,len(U)); arms[arm].add(point(U[u],U[v]))
        arms['proper'].add(point(ZERO,ZERO))
        for a in U:
            arms['proper'].update((point(a,ZERO),point(ZERO,a)))
    factors = {}
    for a0,c0 in product(U,repeat=2):
        a,c = poly(a0),poly(c0)
        cand = (triple(a*c % H),triple((c-a*a) % H),ONE)
        factors[cand] = a0,c0
        arms['product'].add(cand)
        if ZERO in cand: arms['proper'].add(cand)
    return arms,factors,dict(scores=scores,cutoff_ties=ties,steered_per_pivot=per_pivot,
                             raw_product_pairs=len(U)**2,cone_slots=3*12*8)


def run(bits,seed,n,U):
    start = perf_counter()
    N,A = fixture(bits,seed); b,sides,abs_polys = fields(A,n)
    arms,factors,work = generators(U,sides,bits,seed)
    @lru_cache(None)
    def lift(a,side):
        return poly(a)(sides[side])
    def element(v,side):
        return (lift(v[0],side)+lift(v[1],side)*T+lift(v[2],side)*T**3) % abs_polys[side]
    def norm(pair):
        return tuple(abs(int(f.resultant(el))) for f,el in zip(abs_polys,pair))
    def evaluated(pair):
        res = [int(el(b)) % N for el in pair]
        assert res[0] == res[1]
        return gcd(res[0],N)
    def pack(pair): return tuple(tuple(coeffs(p)) for p in pair)
    norm_cache = {}
    @lru_cache(None)
    def factor(value):
        if value == 1: return ()
        fa = pari.factorint(value)
        return tuple(int(p) for p in fa[0])
    needed = [set(),set()]
    def keep(pair,ns):
        key = pack(pair); norm_cache[key] = ns
        for i,v in enumerate(ns): needed[i].update(factor(v))
        return key
    valid,excluded,capped,constituents = {},{}, {},{}
    for v in sorted(set.union(*arms.values())):
        pair = tuple(element(v,i) for i in (0,1))
        g = evaluated(pair)
        if g != 1:
            excluded[v] = g; continue
        ns = norm(pair); F = ns[0]*ns[1]
        assert F > 0
        valid[v] = F
        if F >= N*N: continue
        capped[v] = keep(pair,ns)
        if v in factors:
            a,c = factors[v]; separate = []
            for idx in (0,1):
                ps = tuple((T+lift(a,i) if idx==0 else T*T-lift(a,i)*T+lift(c,i)) % abs_polys[i] for i in (0,1))
                assert evaluated(ps)==1
                separate.append(keep(ps,norm(ps)))
            assert all((poly(separate[0][i])*poly(separate[1][i]))%abs_polys[i]==pair[i] for i in (0,1))
            constituents[v] = separate
    ideals = [Ideals(f,sorted(ps)) for f,ps in zip(abs_polys,needed)]
    cols = {}; masks = {}
    for key,ns in norm_cache.items():
        row = 0
        for i,(a,value) in enumerate(zip(key,ns)):
            for label in ideals[i].parity(poly(a),value,factor(value)):
                idx = cols.setdefault((i,label),len(cols)); row |= 1<<idx
        masks[key] = row
    for v,(s,t) in constituents.items(): assert masks[capped[v]] == masks[s]^masks[t]
    def prime_max(key): return max((p for v in norm_cache[key] for p in factor(v)),default=1)
    records = []
    for v,key in capped.items():
        records.append(dict(v=v,norms=norm_cache[key],mask=hex(masks[key]),groups=sorted(k for k in arms if v in arms[k])))
    summary = {}
    for group,points in list(arms.items())+[('union',set.union(arms['steer'],arms['random'],arms['product']))]:
        full = {v for v in points if ZERO not in v}
        use = points if group=='proper' else full
        tested = use & valid.keys(); hits = tested & capped.keys()
        summary[group] = dict(tuples=len(points),full=len(full),valid=len(tested),cap=len(hits),
                              min_full=min((valid[v] for v in full if v in valid),default=None),
                              nonunits=dict(Counter(excluded[v] for v in use if v in excluded)))
        for name,bound in [('all',None),('smooth97',97)]:
            ok = lambda v: bound is None or prime_max(capped[v])<=bound
            base = [masks[capped[v]] for v in arms['proper'] & capped.keys() if ok(v)]
            good = [v for v in hits if ok(v)]; rows = [masks[capped[v]] for v in good]
            row_info = dict(count=len(good),rank=rank(rows),increment=rank(base+rows)-rank(base))
            if group=='product':
                charged = base+[masks[t] for v in good for t in constituents[v]]
                row_info['after_factors']=rank(charged+rows)-rank(charged)
            summary[group][name] = row_info
    result = dict(bits=bits,n=n,seed=seed,N=N,A=A,b=b,sides=[coeffs(p) for p in sides],work=work | dict(unique_evaluations=len(valid)+len(excluded),exact_norm_pairs=len(valid),capped=len(capped),columns=len(cols)),summary=summary)
    raw = dict(result=result,columns=sorted(((j,i,p,index) for (i,(p,index)),j in cols.items())),rows=records)
    Path(f'rows_{bits}_{seed}.json').write_text(json.dumps(raw,separators=(',',':'))+'\n')
    result['seconds']=round(perf_counter()-start,3)
    return result

if __name__ == '__main__':
    U = dictionary()
    for bits,count,n in [(44,6,4),(60,2,5)]:
        for seed in range(count):
            print(json.dumps(run(bits,seed,n,U),separators=(',',':')),flush=True)
