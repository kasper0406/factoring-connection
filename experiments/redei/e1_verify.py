"""E1: Redei 4-rank formula vs PARI on random (N, k)."""
import random, json, sys
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import *
from sympy import randprime, primerange
random.seed(1)
small = list(primerange(3, 60))
res = {'n': 0, 'fail': 0, 'by_t': {}}
for it in range(600):
    p = randprime(10**2, 10**3); q = randprime(10**3, 10**4)
    if p == q: continue
    k = 1
    for r in random.sample(small, random.randint(0, 2)):
        if r not in (p, q): k *= r
    if random.random() < 0.5: k *= 2
    m = k * p * q
    D = fund_disc(m)
    a = r4_redei(D); b = rank_2k(two_sylow(D), 2)
    res['n'] += 1
    t = len(prime_disc_factors(D))
    res['by_t'].setdefault(t, [0, 0]); res['by_t'][t][0] += 1
    if a != b:
        res['fail'] += 1; res['by_t'][t][1] += 1
        if res['fail'] < 5: print('FAIL', p, q, k, D, a, b)
print(json.dumps(res))
