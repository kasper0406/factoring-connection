"""Exhaust the coefficient grid and audit the primary run with relative resultants."""
from functools import lru_cache
from itertools import product
from math import gcd, prod
from pathlib import Path
import json

from steer import H, ONE, ZERO, Poly, coeffs, dictionary, generators, pari, poly, rank


def audit(result, U):
    bits, seed = result['bits'],result['seed']
    n_mod, b, a = result['N'], result['b'], result['A']
    p, q = map(Poly, result['sides'])
    arms, _, _ = generators(U, (p,q), bits, seed)
    grid = {(u,v,ONE) for u,v in product(U,repeat=2)}
    grid.update((u,ONE,v) for u,v in product(U,repeat=2))
    grid.update((ONE,u,v) for u,v in product(U,repeat=2))
    arms['grid'] = grid
    y, t = pari('y'), pari('x')
    z = pari.Mod(y,y**3+y**2-2*y-1)
    fields = [pari.Polrev(coeffs(f))-z for f in (p,q)]
    maximal = [pari.Polrev(coeffs(H(f))).nfinit() for f in (p,q)]
    @lru_cache(None)
    def coef(entry):
        a,b,c = entry
        return a + b*z + c*z*z
    @lru_cache(None)
    def residue(entry):
        return int(poly(entry)(a)) % n_mod
    full = {}; capped = {}; excluded = 0
    for v in sorted(set.union(*arms.values())):
        common = (residue(v[0])+b*residue(v[1])+b**3*residue(v[2])) % n_mod
        if gcd(common,n_mod) != 1:
            excluded += int(v in grid)
            continue
        el = coef(v[0]) + coef(v[1])*t + coef(v[2])*t**3
        ns = tuple(abs(int(f.polresultant(el).norm())) for f in fields)
        assert min(ns) > 0
        if ZERO not in v: full[v] = prod(ns)
        if prod(ns) < n_mod*n_mod: capped[v] = ns
    records = json.loads(Path(f'rows_{bits}_{seed}.json').read_text())['rows']
    original = set.union(*(pts for name,pts in arms.items() if name!='grid'))
    assert {v:ns for v,ns in capped.items() if v in original} == {
        tuple(map(tuple,item['v'])):tuple(item['norms']) for item in records}
    labels, masks, largest = {},{},{}
    for v,ns in capped.items():
        mask = 0; large = 1
        for side,(f,nf) in enumerate(zip((p,q),maximal)):
            el = poly(v[0])(f) + poly(v[1])(f)*Poly([0,1]) + poly(v[2])(f)*Poly([0,0,0,1])
            ideals = nf.idealfactor(pari.Polrev(coeffs(el)))
            factors = []
            for pr,ev in zip(ideals[0],ideals[1]):
                prime, inertia, valuation = int(pr.pr_get_p()), int(pr.pr_get_f()), int(ev)
                factors.append(prime**(inertia*valuation)); large = max(large,prime)
                if valuation%2:
                    label = side,str(nf.idealhnf(pr))
                    mask |= 1<<labels.setdefault(label,len(labels))
            assert prod(factors)==ns[side]
        masks[v],largest[v] = mask,large
    for name,summary in result['summary'].items():
        if name=='proper': continue
        cand = original if name=='union' else arms[name]
        assert min(full[v] for v in cand if v in full)==summary['min_full']
    metrics = {}
    for name,bound in [('all',None),('smooth97',97)]:
        keep = {v for v in capped if bound is None or largest[v]<=bound}
        base = [masks[v] for v in keep & arms['proper']]
        for group,pts in arms.items():
            use = pts if group=='proper' else {v for v in pts if ZERO not in v}
            rows = [masks[v] for v in keep & use]
            inc = rank(base+rows)-rank(base)
            if group=='grid':
                metrics[name]=dict(hits=len(rows),rank=rank(rows),increment=inc)
            else:
                assert dict(count=len(rows),rank=rank(rows),increment=inc).items() <= result['summary'][group][name].items()
    best = min(grid & full.keys(),key=lambda v:(full[v],v))
    el = coef(best[0])+coef(best[1])*t+coef(best[2])*t**3
    best_ns = [abs(int(f.polresultant(el).norm())) for f in fields]
    return dict(bits=bits,seed=seed,N=n_mod,grid=len(grid),excluded=excluded,best=best,
                best_norms=best_ns,original_capped_checked=len(records),metrics=metrics)


if __name__ == '__main__':
    U = dictionary()
    for line in Path('results.jsonl').read_text().splitlines():
        print(json.dumps(audit(json.loads(line),U),separators=(',',':')),flush=True)
