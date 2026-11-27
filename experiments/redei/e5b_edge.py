"""E5b: check the edge-window KH signal. Per-t terms near t=2sqrtN; held-out prediction of q mod 4 within
N mod 8 class; also label p mod 4 and a fair coin; also relabel by (q mod 4) for N=1 mod 4 only."""
import sys, random, json, collections, math
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import pari
from sympy import randprime
from math import isqrt
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 33)
W = 40
rows = []
seen = set()
while len(rows) < 4000:
    p = randprime(10**4, 10**5); q = randprime(10**4, 10**5)
    if p == q: continue
    p, q = min(p, q), max(p, q); N = p * q
    if N in seen: continue
    seen.add(N)
    T = isqrt(4 * N - 1)
    edge = sum(2 * int(6 * pari.qfbhclassno(4 * N - t * t)) for t in range(T - W + 1, T + 1)) % 48
    rows.append(dict(N=N, p=p, q=q, edge=edge))
def acc(label):
    tr, te = rows[:3000], rows[3000:]
    tab = collections.defaultdict(collections.Counter)
    for r in tr: tab[(r['N'] % 8, r['edge'])][label(r)] += 1
    base = collections.defaultdict(collections.Counter)
    for r in tr: base[r['N'] % 8][label(r)] += 1
    c1 = c0 = 0
    for r in te:
        k = (r['N'] % 8, r['edge'])
        pred = tab[k].most_common(1)[0][0] if tab[k] else base[r['N'] % 8].most_common(1)[0][0]
        c1 += pred == label(r); c0 += base[r['N'] % 8].most_common(1)[0][0] == label(r)
    return round(c1 / len(te), 4), round(c0 / len(te), 4)
rr = random.Random(5)
coin = {r['N']: rr.randint(0, 1) for r in rows}
res = {'q mod 4': acc(lambda r: r['q'] % 4), 'p mod 4': acc(lambda r: r['p'] % 4),
       'coin': acc(lambda r: coin[r['N']]),
       'edge mod 4 vs N mod 8 table': {str(k): dict(collections.Counter(r['edge'] % 8 for r in rows if r['N'] % 8 == k)) for k in (1, 3, 5, 7)}}
# direct cross-tab for N = 1 mod 4: q mod 4 vs edge mod 16
ct = collections.Counter((r['N'] % 8, r['q'] % 4, r['edge'] % 16) for r in rows)
res['crosstab_N1mod8'] = {f'q{qq}': [ct[(1, qq, e)] for e in range(0, 16, 2)] for qq in (1, 3)}
res['crosstab_N5mod8'] = {f'q{qq}': [ct[(5, qq, e)] for e in range(0, 16, 2)] for qq in (1, 3)}
print(json.dumps(res, indent=0))
