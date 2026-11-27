"""E5 (own idea, null test): does a short window of the Kronecker-Hurwitz family already carry the
full-sum bit q mod 4 (= (2N+2q mod 8)/2)? Windows: centre |t|<=W and edge t in (2sqrtN-W, 2sqrtN).
Measure agreement of (window sum *? mod 8 pattern) best single rule vs q mod 4 by MI estimate."""
import sys, random, json, math, collections
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import pari
from sympy import randprime
from math import isqrt
random.seed(21)
W = 40
rows = []
for it in range(2500):
    p = randprime(10**4, 10**5); q = randprime(10**4, 10**5)
    if p == q: continue
    p, q = min(p, q), max(p, q); N = p * q
    def H6(n): return int(6 * pari.qfbhclassno(n))  # 6H is an integer
    cen = sum((1 if t == 0 else 2) * H6(4 * N - t * t) for t in range(0, W + 1)) % 48
    T = isqrt(4 * N - 1)
    edge = sum(2 * H6(4 * N - t * t) for t in range(T - W + 1, T + 1)) % 48
    rows.append((N % 8, q % 4, cen, edge))
def mi(pairs):
    n = len(pairs); c = collections.Counter(pairs); cx = collections.Counter(a for a, _ in pairs); cy = collections.Counter(b for _, b in pairs)
    return sum(v / n * math.log2(v * n / (cx[a] * cy[b])) for (a, b), v in c.items())
res = {'n': len(rows)}
for name, idx in (('centre', 2), ('edge', 3)):
    res[name] = round(mi([((r[0], r[idx]), r[1]) for r in rows]), 4)
rr = random.Random(0)
labs = [r[1] for r in rows]
for name, idx in (('centre', 2), ('edge', 3)):
    obs = mi([((r[0], r[idx]), r[1]) for r in rows]); ge = 0; nulls = []
    for _ in range(200):
        rr.shuffle(labs)
        # shuffle within N mod 8 classes keeps public dependence
        nulls.append(0)
    # within-class permutation
    cls = collections.defaultdict(list)
    for i, r in enumerate(rows): cls[r[0]].append(i)
    nulls = []
    for _ in range(200):
        nl = [None] * len(rows)
        for idxs in cls.values():
            v = [rows[i][1] for i in idxs]; rr.shuffle(v)
            for i, x in zip(idxs, v): nl[i] = x
        nulls.append(mi([((r[0], r[idx]), nl[i]) for i, r in enumerate(rows)]))
    res[name + '_perm_p'] = sum(1 for x in nulls if x >= obs) / 200
    res[name + '_null_mean'] = round(sum(nulls) / 200, 4)
res['public_only'] = round(mi([(r[0], r[1]) for r in rows]), 4)
print(json.dumps(res))
json.dump(res, open('/Volumes/git/fc-redei/experiments/redei/e5_khwindow.json', 'w'))
