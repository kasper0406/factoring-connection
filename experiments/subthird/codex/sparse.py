"""Search common, sparse polynomial recipes by their joint conjugate envelope."""
import json
from itertools import product
from math import prod

from flint import arb
from sympy import integer_nthroot, primerange
from experiment import (T, one, Poly, here, fixture, select, coeffs,
                        Ideals, measure)
from lattice import square_dimension
from src.rankdegen.fh import rows_at


def choose(P,Q):
    f,g = coeffs(P*P-2), coeffs(Q*Q-2)
    basis = [one,T]
    for m in (3,4):
        best = None
        for degree in range(2,m+2):
            for tail in product((-1,0,1),repeat=degree-2):
                u = Poly([0,0,*tail,1])
                candidate = basis+[u]
                if any(square_dimension(X,candidate)<=2*m-1 for X in (P,Q)):
                    continue
                lists = list(map(coeffs,candidate))
                rows,_ = rows_at(f,lists,g,lists,192)
                h = arb(1)
                for row in rows:
                    mx = abs(row[0])
                    for y in row[1:]: mx = mx.max(abs(y))
                    h *= mx
                if best is None or h<best[0]: best = (h,u)
        assert best is not None
        basis = basis+[best[1]]
        yield m,basis


def run():
    with (here/'sparse_results.jsonl').open('w') as out:
        for bits,seed in product((60,72,84,96),range(2)):
            N = fixture(bits,seed)
            a,b,P,Q = select(N)
            fourth,exact = integer_nthroot(bits*((bits-1).bit_length())**3,4)
            bound = 2**(int(fourth)+int(not exact))
            primes = list(primerange(bound+1))
            fields = [Ideals(X*X-2,primes) for X in (P,Q)]
            pdata = primes,prod(primes)
            for m,A in choose(P,Q):
                rec = measure(N,b,P,Q,m,'sparse',pdata,fields,(A,A))
                rec.update(B=bound,seed=seed,generators=list(map(coeffs,A)))
                out.write(json.dumps(rec)+'\n')
                out.flush()
                print(bits,seed,m,rec['aggregate']['new'],rec['all']['new'],flush=True)
if __name__ == '__main__':
    run()
