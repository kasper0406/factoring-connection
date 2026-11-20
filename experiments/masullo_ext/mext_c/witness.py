"""Inspect the one uncovered smooth class specified in witness_prereg.txt."""
from cover import ex, itertools, json, math, pair_units, prefix, unit_key
from sympy import factorint

N = ex.fixture(72, 19)
_, b, P, Q = ex.select(N)
f, g = P*P-2, Q*Q-2
U, V = [ex.expand(X, X*X-2, ex.generators(X, "divided", 3)) for X in (P,Q)]
residues = [int(u(b)) % N for u in U]
points = [z for z in itertools.product(range(-2,3),repeat=6)
          if next((j for j in z if j),0) > 0 and
          math.gcd(sum(t*c for t,c in zip(z,residues)), N) == 1]
rows, _ = ex.rows_at(ex.coeffs(f), list(map(ex.coeffs,U)),
                     ex.coeffs(g), list(map(ex.coeffs,V)), 192)
spaces = [prefix(row,points) for row in rows]
units = pair_units()
records = []
all_primes = set()
for z in sorted({unit_key(z,units) for z in points}):
    a, w = ex.element(U,z), ex.element(V,z)
    v1, v2 = abs(int(f.resultant(a))), abs(int(g.resultant(w)))
    if not 0 < v1*v2 <= 3**32 * 32:
        continue
    factors = sorted(set(factorint(v1)) | set(factorint(v2)))
    all_primes.update(factors)
    covered = any(ex.kranks(s+[z]) == len(s) for s in spaces)
    records.append((z,a,w,v1,v2,factors,covered,max(factors,default=1)<=8192))
fields = [ex.Ideals(X,sorted(all_primes)) for X in (f,g)]
columns = {}
answers = []
for z,a,w,v1,v2,primes,covered,smooth in records:
    mask = 0
    for side,(alpha,norm) in enumerate(((a,v1),(w,v2))):
        for prime in fields[side].parity(alpha,norm,primes):
            mask |= 1 << columns.setdefault((side,prime),len(columns))
    answers.append((z,v1,v2,covered,smooth,mask))
(z,v1,v2,_,_,mask), = [a for a in answers if not a[3] and a[4]]
others = [row[-1] for row in answers if row[1]*row[2] < v1*v2 or (row[3] and row[4])]
print(json.dumps(dict(z=z,norms=[v1,v2],coordinate_deletion=any(
    a==b==0 for a,b in zip(z[::2],z[1::2])), compared_rows=len(others),
    comparison_rank=ex.rank(others), rank_with_witness=ex.rank(others+[mask]))))
