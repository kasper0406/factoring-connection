# exploratory (disclosed): matched controls for E3 (same p, q uniform prime in q_range(p))
import sys, json, random, math, numpy as np, gmpy2
sys.path.insert(0, '/Volumes/git/factoring-connection')
import src.inversesieve.experiments as E, src.inversesieve.analyze as A
rs = [json.loads(l) for l in open('/Volumes/git/factoring-connection/experiments/raw/inversesieve_e3.jsonl')]
unc = [r for r in rs if r['kind'] == 'uncovered']
cfgs = {(c['n'], c['shape']): E._with_h0(c) for c in E.configs()}
rng = random.Random(77); ctl = []
for r in unc:
    c = cfgs[(r['n'], r['shape'])]; lo, hi = E.q_range(r['p'], r['n'])
    while True:
        q = rng.randint(lo, hi)
        if gmpy2.is_prime(q): break
    f = E.pair_features(r['p'], q, r['n'], c, c['h0']); ctl.append(f)
jump = lambda r: math.log(r['cfB2'] / r['cfB'])
print('n', len(unc))
print('Lstar KS', A.ks_2samp([math.log(r['Lstar']) for r in unc], [math.log(r['Lstar']) for r in ctl]))
print('cfjump KS', A.ks_2samp([jump(r) for r in unc], [jump(r) for r in ctl]))
print('dq KS', A.ks_2samp([r['dq'] for r in unc], [r['dq'] for r in ctl]))
print('lehman', sum(r['lehman_in'] for r in unc)/len(unc), sum(r['lehman_in'] for r in ctl)/len(ctl))
print('median logL* unc/ctl', np.median([math.log(r['Lstar']) for r in unc]), np.median([math.log(r['Lstar']) for r in ctl]))
print('mean cfjump unc/ctl', np.mean([jump(r) for r in unc]), np.mean([jump(r) for r in ctl]))
