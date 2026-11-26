import json
from math import comb, gcd, prod
from pathlib import Path
from random import Random
from flint import fmpz_mod_poly_ctx, fmpz_poly
from gmpy2 import jacobi, is_prime

RATIONAL_F = fmpz_poly([-1,-2,1,1])

def lift(poly):
    return fmpz_poly([int(c) for c in poly.coeffs()])

def trace(poly, N):
    return (3*int(poly[0]) - int(poly[1]) + 5*int(poly[2])) % N

def cubic_probe(N, rng):
    ring = fmpz_mod_poly_ctx(N)
    X = ring.gen()
    F = ring([-1,-2,1,1])
    while True:
        a = ring([rng.randrange(N) for _ in range(3)])
        norm = int(RATIONAL_F.resultant(lift(a))) % N
        div = gcd(norm,N)
        if div == 1: break
        if div < N: return dict(splits={'norm':div})
    b = (a.pow_mod(3,F) * pow(norm,-1,N)) % F
    T = trace(b,N)
    S = ((T*T - trace((b*b)%F,N))*pow(2,-1,N))%N
    H = ring([T*T-4*S,-8,-2*T,0,1])
    hi = lift(H)
    delta = int(hi.resultant(hi.derivative()))
    g = gcd(delta,N)
    if g != 1:
        return dict(splits={'discriminant':g} if g < N else {},degenerate=True)
    v = X.pow_mod(N,H)-X
    z = int(hi.resultant(lift(v)))%N
    g = gcd(N,z)
    return dict(splits={'resultant':g} if 1<g<N else {},degenerate=False,
                zero=z==0,sign=int(jacobi(z,N)),H=[int(c) for c in H.coeffs()],
                b=[int(c) for c in b.coeffs()])

def verify_cubic(s, row, inert):
    ring = fmpz_mod_poly_ctx(s)
    X = ring.gen()
    H, b = ring(row['H']), ring(row['b'])
    F = ring([-1,-2,1,1])
    nroots = (X.pow_mod(s,H)-X).gcd(H).degree()
    if inert:
        square = b.pow_mod((s**3-1)//2,F) == ring.one()
        assert square and nroots == 1
    else:
        roots = F.roots()
        assert len(roots) == 3
        square = all(jacobi(int(b(x)),s) == 1 for x,_ in roots)
        assert nroots == (4 if square else 0)
    return square

def cubic(n):
    rng = Random(20261126*10000+n*10+4)
    with Path(f'cubic{n}_out.jsonl').open('w') as out:
        for line in Path(f'cubic{n}.tsv').read_text().splitlines():
            N,p,q,label = map(int,line.split())
            row = cubic_probe(N,rng)
            if 'H' in row:
                left = verify_cubic(p,row,label)
                right = verify_cubic(q,row,label)
                row['oracle_square'] = bool(left and right)
            assert all(v in (p,q) for v in row['splits'].values())
            row.update(N=N,label=label)
            out.write(json.dumps(row,separators=(',',':'))+'\n')

def carries(m,s):
    count = carry = 0
    while m or carry:
        m,digit = divmod(m,s)
        carry = int(2*digit + carry >= s)
        count += carry
    return count

def carry_file(n,topic):
    d = 4 if topic=='cf' else 16
    tally = {}
    for line in Path(f'{topic}{n}.tsv').read_text().splitlines():
        N,p,q,label = map(int,line.split())
        vals = [carries(N//d,s) for s in (p,q)]
        assert vals == [int(q % d > d//2), int(p % d > d//2)]
        state = ['unit','proper','zero'][sum(vals)]
        key = f'{label}/{N%d}/{state}'
        tally[key] = tally.get(key,0)+1
    return tally

def small_checks():
    primes = [p for p in range(5,160) if is_prime(p)]
    total = 0
    for p in primes:
        for q in primes:
            if not p < q < 2*p: continue
            N = p*q
            for d in [4,8,16,32,64]:
                if d >= p: continue
                m = N//d
                predicted = prod(s for s,t in [(p,q),(q,p)] if t%d>d//2)
                assert gcd(comb(2*m,m),N) == predicted
                total += 1
    return total

if __name__ == '__main__':
    for n in (40,60,80): cubic(n)
    stats = {f'{topic}{n}':carry_file(n,topic) for topic,sizes in
             [('cf',(40,52,64,80)),('carry',(40,60,80))] for n in sizes}
    stats['direct_small_checks'] = small_checks()
    Path('carry_counts.json').write_text(json.dumps(stats,indent=2)+'\n')
