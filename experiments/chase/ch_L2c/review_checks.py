import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import sys

sys.path.insert(0,"/Volumes/git/factoring-connection/experiments/masullo_ext/mext_c")
from cover import ex, pair_units
from flint import acb, arb, ctx
from sympy import cyclotomic_poly, factorint, kronecker_symbol

ctx.prec=256
N=ex.fixture(72,19)
_,b,P,Q=ex.select(N)
U,V=[ex.expand(X,X*X-2,ex.generators(X,"divided",3)) for X in (P,Q)]
rows,_=ex.rows_at(ex.coeffs(P*P-2),list(map(ex.coeffs,U)),
                  ex.coeffs(Q*Q-2),list(map(ex.coeffs,V)),256)
w=(-2,1,-1,1,-1,0)
reps=set()
for u,v in pair_units():
    z=tuple(s for a,c in zip(w[::2],w[1::2]) for s in (a*u+2*c*v,a*v+c*u))
    if max(map(abs,z))<=2:
        reps.add(z)
vals=lambda z: [abs(sum((t*u for t,u in zip(z,row)),acb(0))) for row in rows]
sgns=[1 if row[1].real>0 else -1 for row in rows]
den=arb(3)**3*arb(3).sqrt()*((1+arb(2).sqrt())*2)**4
pairs=[]
tau_product=arb(1)
for i,j in itertools.combinations(range(32),2):
    if sgns[i]!=sgns[j]:
        continue
    x,y=rows[i][::2],rows[j][::2]
    norm=sum(abs(x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3])**2
             for k in range(3)).sqrt()
    assert norm>0
    two_tau=norm/den
    tau_product*=two_tau/2
    pairs.append((i,j,two_tau))
assert len(pairs)==240
for z in sorted(reps):
    vs=vals(z)
    ratios=[vs[i]*vs[j]/t for i,j,t in pairs]
    assert all((r>1) or (r<1) for r in ratios)
    print("rep",z,"analytic tubes",sum(r<1 for r in ratios),"min product/(2tau)",min(ratios,key=float))
vs=vals(w)
inds=sorted(range(16),key=lambda i:float(vs[i]))[:2]
print("smallest first-side tags and base signs",[(i,sgns[i]) for i in inds])
print("global analytic product cap",tau_product**(arb(1)/15))
nodes=(0,1,-1,2)
valid=[]
for size in range(1,4):
    for subset in itertools.combinations(nodes,size):
        coeffs=[Fraction(0) if r in subset else Fraction(1,math.prod(r-s for s in nodes if s!=r and s not in subset)) for r in nodes]
        if all(c.denominator==1 for c in coeffs):
            valid.append(subset)
print("integral Lagrange subset products",len(valid),valid)
K=4
print("nonzero vectors",(vecs:=[(u,v) for u,v in itertools.product(range(-K,K+1), repeat=2) if (u or v) and (u-5*v)%15==0]))
assert all(math.gcd(u,v)>1 for u,v in vecs)
assert 13*45-2**2==7*83 and (13*25+4*5+45)%15==0 and math.gcd(15,2*13*581)==1
cm=int(cyclotomic_poly(15).as_poly().eval(600))
print("CM modulus",cm,"factors",factorint(cm),"base-root symbols",[(p,kronecker_symbol(2,p)) for p in factorint(cm)])
