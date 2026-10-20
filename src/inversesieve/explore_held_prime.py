import sys, random, numpy as np, gmpy2
sys.path.insert(0, '/Volumes/git/factoring-connection')
import src.inversesieve.experiments as E
from src.inversesieve.design import weights
from src.inversesieve.usets import uncovered_partners
n, shape = int(sys.argv[1]), sys.argv[2]
c = E.config(n, shape); pr = c['primes']; V = weights(c['H'], 'signed'); h1, h2 = pr[0::2], pr[1::2]
d1 = np.prod([1 - len(V)/l for l in h1]); d2 = np.prod([1 - len(V)/l for l in h2])
rng = random.Random(9); plo, phi = E.p_range(n); t = [0]*6
for _ in range(int(sys.argv[3])):
    p = E.random_prime(rng, plo, phi); lo, hi = E.q_range(p, n)
    q = np.arange(lo, hi+1, dtype=np.int64)
    isp = np.array([bool(gmpy2.is_prime(int(x))) for x in q])
    Q1 = uncovered_partners(p, h1, V, lo, hi); Q = uncovered_partners(p, pr, V, lo, hi)
    s1 = set(Q1.tolist()); s = set(Q.tolist())
    P0 = q[isp]; P1 = [x for x in P0.tolist() if x in s1]; P2 = [x for x in P1 if x in s]
    t[0]+=len(P0); t[1]+=len(P1); t[2]+=len(P2); t[3]+=q.size; t[4]+=Q1.size; t[5]+=Q.size
print(n, shape, 'primes: half1/model', round(t[1]/(t[0]*d1),3), 'held', round(t[2]/(t[1]*d2),3), 'full/model', round(t[2]/(t[0]*d1*d2),3),
      '| integers: half1/model', round(t[4]/(t[3]*d1),3), 'held', round(t[5]/(t[4]*d2),3))
