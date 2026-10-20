# exploratory replication on fresh seeds (disclosed): uncovered vs MATCHED control (same p)
import sys, json, random, math, numpy as np, gmpy2, collections
sys.path.insert(0, '/Volumes/git/factoring-connection')
import src.inversesieve.experiments as E, src.inversesieve.analyze as A
from src.inversesieve.design import weights
from src.inversesieve.usets import uncovered_partners
from multiprocessing import Pool
def dead(p, q, H, Y):
    return sum(1 for (a, B) in weights(H, 'signed') if 0 < abs(B*p - a*q) < Y)
def job(args):
    c, i = args
    rng = random.Random(f"inversesieve|replic|{c['n']}|{c['shape']}|{i}")
    plo, phi = E.p_range(c['n']); p = E.random_prime(rng, plo, phi); lo, hi = E.q_range(p, c['n'])
    V = weights(c['H'], 'signed')
    Q = [int(x) for x in uncovered_partners(p, c['primes'], V, lo, hi) if gmpy2.is_prime(int(x))]
    rng.shuffle(Q); out = []
    for q in Q[:5]:
        while True:
            q2 = rng.randint(lo, hi)
            if gmpy2.is_prime(q2): break
        fu = E.pair_features(p, q, c['n'], c, c['h0']); fc = E.pair_features(p, q2, c['n'], c, c['h0'])
        out.append((c['n'], c['shape'], fu['lehman_in'], fc['lehman_in'], math.log(fu['cfB2']/fu['cfB']), math.log(fc['cfB2']/fc['cfB']),
                    dead(p, q, c['H'], c['Y']), dead(p, q2, c['H'], c['Y'])))
    return out
if __name__ == '__main__':
    cf = [E._with_h0(c) for c in E.configs()]
    with Pool(9) as pool:
        res = [x for r in pool.map(job, [(c, i) for c in cf for i in range(80)]) for x in r]
    json.dump(res, open('/Volumes/git/factoring-connection/experiments/raw/inversesieve_x_replicate.json', 'w'))
    for shp in ('M', 'S', None):
        rr = [r for r in res if shp is None or r[1] == shp]
        a = sum(r[2] for r in rr); b = sum(r[3] for r in rr); m = len(rr)
        print(shp, m, 'lehman RR', [round(x,3) for x in A.rr_ci(a, m, b, m)], 'dlogjump', round(np.mean([r[4]-r[5] for r in rr]),3),
              '+-', round(np.std([r[4]-r[5] for r in rr])/math.sqrt(m),3), 'dead unc/ctl', round(np.mean([r[6] for r in rr]),3), round(np.mean([r[7] for r in rr]),3))
    for n in (24,28,32,36,40):
        for shp in ('M','S'):
            rr=[r for r in res if r[0]==n and r[1]==shp]; m=len(rr)
            print(n, shp, m, round(sum(r[2] for r in rr)/m,3), round(sum(r[3] for r in rr)/m,3), round(np.mean([r[6] for r in rr]),2), round(np.mean([r[7] for r in rr]),2))
