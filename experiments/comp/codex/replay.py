import itertools
import math
import sys
from pathlib import Path
sys.path.insert(0,'/Volumes/git/factoring-connection/codex/scratch/dd_g')
import run_comp_adversary as g
ex=g.ex
ex.n,ex.d=12,24
N=ex.fixture(112,19)
spec=g.build_M2_mix(N,12,4)
P,Q,b=spec['P'],spec['Q'],spec['b']
polys=[X*X-2 for X in (P,Q)]
gens=[ex.expand(X,f,spec[name]) for X,f,name in zip((P,Q),polys,('Ag','Bg'))]
res=[int(t(b))%N for t in gens[0]]
assert res==[int(t(b))%N for t in gens[1]]
units=g.pair_units()
cap=3**48*2**math.ceil(112/24)
orbits={g.orbit_rep(z,units) for z in itertools.product((-1,0,1),repeat=8)
        if math.gcd(sum(s*t for s,t in zip(z,res)),N)==1}
capped=[]
ones=[]
for z in sorted(orbits):
    norm=math.prod(abs(int(f.resultant(ex.element(basis,z)))) for f,basis in zip(polys,gens))
    if 0<norm<=cap and any(z[2:]):  # 1,T and the two pole quotients are k-independent
        capped.append(z)
        if norm==1:
            ones.append(z)
print(f'M2_mix bits=112 seed=19 n=12 m=4 b={b}')
print(f'nonscalar capped orbits={len(capped)} rank={ex.kranks(capped)}')
print(f'nonscalar norm-one orbits={len(ones)} rank={ex.kranks(ones)}')
for z in ones:
    print(z)
