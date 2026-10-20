# exploratory: co-coverage of a partner q by two band primes, by weight relation
import sys, random, numpy as np, collections
sys.path.insert(0, '/Volumes/git/factoring-connection')
import src.inversesieve.experiments as E
from src.inversesieve.design import weights
from math import gcd
n, shape = int(sys.argv[1]), sys.argv[2]
c = E.config(n, shape); pr = c['primes']; V = weights(c['H'], 'signed')
rng = random.Random(5); plo, phi = E.p_range(n)
obs = collections.Counter(); exp = collections.Counter()
for t in range(6):
    p = E.random_prime(rng, plo, phi); lo, hi = E.q_range(p, n)
    q = np.arange(lo, hi + 1, dtype=np.int64); W = q.size
    sub = rng.sample(pr, 60)
    lab = {}
    for l in sub:
        # label each q by the weight index covering it at l (-1 none); at most one weight per prime
        m = -np.ones(l, np.int64)
        for i, (A, B) in enumerate(V):
            m[(B * p * pow(A % l, -1, l)) % l] = i
        lab[l] = m[q % l]
    for a in range(len(sub)):
        for b in range(a + 1, len(sub)):
            la, lb = lab[sub[a]], lab[sub[b]]
            both = (la >= 0) & (lb >= 0)
            same = both & (la == lb)
            obs['both'] += int(both.sum()); obs['same'] += int(same.sum())
            na, nb = int((la >= 0).sum()), int((lb >= 0).sum())
            exp['both'] += na * nb / W
            # expected same-weight co-coverage under independence: sum_v n_a(v) n_b(v) / W
            ca = np.bincount(la[la >= 0], minlength=len(V)); cb = np.bincount(lb[lb >= 0], minlength=len(V))
            exp['same'] += float((ca * cb).sum()) / W
print(n, shape, 'both obs/exp', round(obs['both'] / exp['both'], 4), 'same obs', obs['same'], 'same exp', round(exp['same'], 1),
      'diff-weight obs/exp', round((obs['both'] - obs['same']) / (exp['both'] - exp['same']), 4))
