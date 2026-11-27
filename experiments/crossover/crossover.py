"""Hidden bits of C(N) mod b^k beyond public data, versus precision k (see prereg.md)."""
import sys, json, random, time
from collections import defaultdict
from math import log2
from sympy import primerange, jacobi_symbol as J
import flint
sys.path.insert(0, '/Volumes/git/fc-crossover/experiments/modrec')
from hull_count import short_square

random.seed(20261127)
PR = list(primerange(2**13, 2**14))
MOD = 2**16 * 3**8 * 5**6

def tau_table(L):
    # Delta = q * prod (1-q^n)^24 mod MOD, coefficients up to L
    e = [0]*(L+1)
    k = 0
    while True:
        done = True
        for s in (1, -1):
            g = k*(3*k - s)//2 if k else 0
            if g <= L:
                e[g] = (e[g] + (-1)**k) % MOD; done = False
        if k == 0:
            e[0] = 1
        k += 1
        if done and k > 2: break
    P = flint.nmod_poly(e, MOD)
    def mul(a, b):
        c = (a*b).coeffs()[:L+1]
        return flint.nmod_poly([int(x) for x in c], MOD)
    P2 = mul(P, P); P4 = mul(P2, P2); P8 = mul(P4, P4); P16 = mul(P8, P8); P24 = mul(P8, P16)
    c = [int(x) for x in P24.coeffs()] + [0]*(L+2)
    return {n: c[n-1] for n in range(1, L+1)}   # tau(n) = coeff of q^(n-1) in prod

def principal23(p):
    # p split in Q(sqrt-23) and in Hilbert class field iff x^3-x-1 has 3 roots mod p
    return sum(1 for x in range(p) if (x*x*x - x - 1) % p == 0) == 3 if p < 50 else \
        len(flint.nmod_poly([-1, -1, 0, 1], p).roots()) == 3

def families(p, q, tau):
    N = p*q
    c4 = lambda s: J(-4, s); c8 = lambda s: J(-8, s); c23 = lambda s: J(-23, s)
    r2 = 4*(1+c4(p))*(1+c4(q))
    rq2 = 2*(1+c8(p))*(1+c8(q))
    if c23(p) == 1 and c23(q) == 1:
        a, b = principal23(p), principal23(q)
        r23 = 8 if (a and b) else (4 if (not a and not b) else 0)
    else:
        r23 = 0
    sig = (1+p)*(1+q); sig3 = (1+p**3)*(1+q**3)
    return {'ctrl': N*N + 3*N, 'r2': r2, 'rQ2': rq2, 'R23': r23, 'r4': 8*sig, 'r8': 16*sig3,
            'KH': 2*(N+q), 'tau': (tau[p]*tau[q]) % MOD}

_JAC = {}
def public(N, b, k):
    j = _JAC.get(N)
    if j is None:
        j = _JAC[N] = (J(-1, N), J(2, N), J(-3, N), J(5, N), J(-23, N))
    return (N % b**(k+1),) + j

def collision_entropy(rows, fam, b, k):
    cells = defaultdict(lambda: defaultdict(int))
    for N, C in rows:
        cells[public(N, b, k)][C[fam] % b**k] += 1
    same = tot = 0
    for cnt in cells.values():
        n = sum(cnt.values())
        if n < 2: continue
        tot += n*(n-1); same += sum(v*(v-1) for v in cnt.values())
    return (-log2(same/tot) if same else float('inf')), tot

def main(nsamp=200000, nS=20000):
    t0 = time.time(); tau = tau_table(2**14 + 2)
    rows = []
    for _ in range(nsamp):
        p, q = sorted(random.sample(PR, 2))
        rows.append((p*q, families(p, q, tau)))
    print('families done', round(time.time()-t0,1), flush=True)
    # S(N) from N alone (hull evaluator), smaller sample
    tS = time.time(); rowsS = []
    for N, C in rows[:nS]:
        rowsS.append((N, {'S': short_square(N, 1, 6, mod_bits=40)}))
    tS = (time.time() - tS)/nS
    print('S done', flush=True)
    out = {'n': nsamp, 'nS': nS, 'S_sec_per_N': tS, 'Y': {}}
    plan = {2: range(1, 9), 3: range(1, 6), 5: range(1, 5)}
    for b, ks in plan.items():
        for fam in ['ctrl', 'r2', 'rQ2', 'R23', 'r4', 'r8', 'KH', 'tau']:
            out['Y'][f'{fam}|b{b}'] = [round(collision_entropy(rows, fam, b, k)[0], 3) for k in ks]
        kS = {2: range(1, 6), 3: range(1, 4), 5: range(1, 3)}[b]
        out['Y'][f'S|b{b}'] = [round(collision_entropy(rowsS, 'S', b, k)[0], 3) for k in kS]
    out['sec'] = round(time.time() - t0, 1)
    json.dump(out, open('/Volumes/git/fc-crossover/experiments/crossover/results.json', 'w'), indent=1)
    for key, v in out['Y'].items(): print(f'{key:10s}', v)
    print('S sec/N', round(tS, 4), 'total', out['sec'])

if __name__ == '__main__':
    main(*map(int, sys.argv[1:]))
