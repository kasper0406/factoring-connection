"""E2b: which x_r are determined by singleton 4-ranks, by (r/N), and whether pair multipliers k=rs add bits.
Brute force over all sigma for B=23 (8 primes -> 2^11 assignments)."""
import sys, random, itertools, math, json
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from ledger import r4_symbolic, true_sigma, family
from sympy import randprime, primerange, jacobi_symbol
random.seed(5)
rs_all = list(primerange(3, 24))
det_by_jac = {1: [0, 0], -1: [0, 0]}
res = []
for it in range(60):
    p = randprime(2**19, 2**20); q = randprime(2**19, 2**20)
    if p == q: continue
    N = p * q
    rs = [r for r in rs_all if N % r]
    sig = true_sigma(p, q, rs)
    out = {}
    for label, fam in (('single', family(rs)), ('pairs', family(rs, pairs=True))):
        obs = {f: r4_symbolic(N, f[0], f[1], sig)[0] for f in fam}
        cons = []
        for e, b, mu in itertools.product((1, -1), repeat=3):
            for xs in itertools.product((1, -1), repeat=len(rs)):
                x = dict(zip(rs, xs))
                s = (e, b, mu, x)
                if all(r4_symbolic(N, f[0], f[1], s)[0] == obs[f] for f in fam):
                    cons.append(s)
        H = 3 + len(rs)
        out[label] = H - math.log2(len(cons))
        if label == 'single':
            for r in rs:
                vals = {c[3][r] for c in cons if c[:3] == sig[:3]}
                det_by_jac[int(jacobi_symbol(r, N))][0 if len(vals) == 1 else 1] += 1
    out['H'] = 3 + len(rs)
    res.append(out)
summ = {'trials': len(res), 'H_mean': sum(r['H'] for r in res) / len(res),
        'single_bits_mean': sum(r['single'] for r in res) / len(res),
        'pairs_bits_mean': sum(r['pairs'] for r in res) / len(res),
        'pairs_gain_max': max(r['pairs'] - r['single'] for r in res),
        'x_r_det_vs_undet_given_true_B0_by_(r/N)': det_by_jac}
print(json.dumps(summ))
json.dump(summ, open('/Volumes/git/fc-redei/experiments/redei/e2b_pairs.json', 'w'), indent=1)
