"""E3b: in hidden cells (bit 4 of p not fixed by public data + p mod 16), which single class-group
feature carries it? Target here: bit4 = [p mod 32 >= 16] for the prime p with p mod 16 = min(p,q mod 16)
(well-defined when p != q mod 16) -- swap-symmetric target; plus target s bit 7."""
import sys, random, json, collections, math
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import two_sylow, fund_disc, rank_2k
from sympy import randprime, jacobi_symbol
random.seed(12)
KS = [1, 2, 3, 6, 5, 10, 7, 14]
rows = []
while len(rows) < 40000:
    p = randprime(2**11, 2**14); q = randprime(2**11, 2**14)
    if p == q: continue
    N = p * q
    if N % 8 != 1: continue
    syl = [two_sylow(fund_disc(k * N)) for k in KS]
    rows.append((p, q, syl))

def mi(pairs):
    n = len(pairs); c = collections.Counter(pairs); cx = collections.Counter(a for a, _ in pairs); cy = collections.Counter(b for _, b in pairs)
    return sum(v / n * math.log2(v * n / (cx[a] * cy[b])) for (a, b), v in c.items())

res = {}
# target: s bit 7 given (N mod 512, s mod 128), only hidden cells
base = collections.defaultdict(list)
for r in rows:
    N = r[0] * r[1]; s = r[0] + r[1]
    base[(N % 512, s % 128)].append(r)
hidden = [r for k, v in base.items() if len({((x[0] + x[1]) >> 7) & 1 for x in v}) == 2 for r in v]
print('hidden rows', len(hidden))
for i, k in enumerate(KS):
    for j in (2, 3, 4):
        # conditional MI: average over base cells, weighted
        tot = 0
        for key, v in base.items():
            labs = {((x[0] + x[1]) >> 7) & 1 for x in v}
            if len(labs) < 2: continue
            pairs = [(rank_2k(x[2][i], j), ((x[0] + x[1]) >> 7) & 1) for x in v]
            tot += len(v) * mi(pairs)
        res[f'k{k}_r{2**j}'] = round(tot / len(hidden), 4)
# null scale: MI of a random fair coin feature
tot = 0
rr = random.Random(1)
for key, v in base.items():
    labs = {((x[0] + x[1]) >> 7) & 1 for x in v}
    if len(labs) < 2: continue
    pairs = [(rr.randint(0, 1), ((x[0] + x[1]) >> 7) & 1) for x in v]
    tot += len(v) * mi(pairs)
res['null_coin'] = round(tot / len(hidden), 4)
print(json.dumps(res, indent=0))
json.dump(res, open('/Volumes/git/fc-redei/experiments/redei/e3b_which.json', 'w'), indent=1)
