"""Blind replay using the original published fixture loop and its exact LLL output."""
from collections import Counter
from itertools import combinations
from math import prod, gcd
from pathlib import Path
import inspect
import sys
from sympy import Matrix

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "hatch_h1g"))
import toy_experiments as g


def recover(rows, scales, primes, word, bound):
    basis = [[int(i == j) for i in range(4)] for j in range(4)]
    equations = []
    for row in rows[:4]:
        c = [x // t for x,t in zip(row,scales)]
        def linear(x):
            v = [0]*4
            for j in range(3):
                z = g.z5_mul(c[4+4*j:8+4*j],g.z5_sigma(x,j))
                v = g.z5_add(v,z)
            return v
        a = Matrix([linear(x) for x in basis]).T
        equations.append(a.row_join(Matrix([-z for z in c[:4]])))
    answers, full = set(), set()
    for size in range(1,5):
        for subset in combinations(equations,size):
            red, pivots = Matrix.vstack(*subset).rref()
            if pivots != (0,1,2,3):
                continue
            x = tuple(red[i,4] for i in range(4))
            if any(z.q!=1 or abs(z)>bound for z in x):
                continue
            x = tuple(map(int,x))
            if prod(p for p in primes if all((z-w)%p==0 for z,w in zip(x,word))) >= prod(sorted(primes)[:5]):
                answers.add(x)
                if size==4:
                    full.add(x)
    return answers,full


def replay():
    # Execute the published construction unchanged up to its first LLL call.
    # The truth label is referenced only after recover() returns.
    source = inspect.getsource(g.run_E3).split('        # Check shortest vector')[0]
    source += '''
        found, full = recover(red16, scales16, qs, w_coords, B_c)
        observations.append((tuple(alpha) in found,len(found),tuple(alpha) in full))
'''
    observed = []
    env = {**vars(g),'recover':recover,'observations':observed}
    exec(source,env)
    env['run_E3']()
    print('trials / planted recovered / recovered with all four:',len(observed),
          sum(t[0] for t in observed),sum(t[2] for t in observed))
    print('verified list sizes:', dict(Counter(t[1] for t in observed)))


if __name__=='__main__':
    assert 3*7*11==231 and gcd(230-1,231)==1 and pow(230,2,231)==1
    assert 4**2-15==1
    # Additional explicit two-root example for the covolume claim.
    H,d,e = 100,2**20,3**100
    P = d*e
    w = 1+d*pow(d,-1,e)
    f = lambda y: (y-1)*(y-2)
    assert f(w)%P == 0 and gcd(w-1,P)==d
    assert max(abs(f(y)) for y in range(H+1)) < d
    assert 20**2 < 6*120 and H>2**6 and P>2**120
    print('Verified the PIN, PLANT and two-root Int(Z) counterexamples.')
    replay()
