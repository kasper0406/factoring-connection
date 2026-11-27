"""E3: do 8-/16-ranks of a small multiplier family carry bit (k-1) of s=p+q mod 2^k (k=7,8)?
Exact PARI 2-Sylow types; held-out lookup accuracy vs within-cell permutation null."""
import sys, random, json, collections
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import two_sylow, fund_disc, rank_2k
from sympy import randprime
random.seed(11)
KS = [1, 2, 3, 6, 5, 10, 7, 14]
rows = []
seen = set()
while len(rows) < 16000:
    p = randprime(2**11, 2**15); q = randprime(2**11, 2**15)
    if p == q or p in (3, 5, 7) or q in (3, 5, 7): continue
    N = p * q
    if N % 8 != 1 or N in seen: continue        # tower setting; distinct N (dedupe fix)
    seen.add(N)
    syl = [two_sylow(fund_disc(k * N)) for k in KS]
    rows.append((N, p + q, syl))
print('rows', len(rows))

def feats(syl, level):
    return tuple(tuple(rank_2k(s, j) for j in range(2, level + 1)) for s in syl)

out = {}
for k in (7, 8):
    for level, name in ((2, 'r4'), (3, 'r4+r8'), (4, 'r4..r16'), (9, 'full')):
        def key(r):
            N, s, syl = r
            return (N % 2**(k + 1), s % 2**(k - 1), feats(syl, level))
        def lab(r): return (r[1] >> (k - 1)) & 1
        def evaluate(rs, labels):
            tr = rs[:12000]; te = rs[12000:]
            tab = collections.defaultdict(collections.Counter)
            for r, l in zip(tr, labels[:12000]): tab[key(r)][l] += 1
            hit = cor = 0
            for r, l in zip(te, labels[12000:]):
                c = tab.get(key(r))
                if not c: continue
                hit += 1
                cor += (c.most_common(1)[0][0] == l)
            return hit, cor / max(hit, 1)
        labels = [lab(r) for r in rows]
        h, acc = evaluate(rows, labels)
        # null: permute labels within base key (N mod, s mod) -> destroys class-group dependence only
        base = collections.defaultdict(list)
        for i, r in enumerate(rows): base[(r[0] % 2**(k + 1), r[1] % 2**(k - 1))].append(i)
        nl = list(labels)
        rr = random.Random(3)
        for idx in base.values():
            vals = [labels[i] for i in idx]; rr.shuffle(vals)
            for i, v in zip(idx, vals): nl[i] = v
        hn, accn = evaluate(rows, nl)
        # fraction of base-key cells where the bit is already public
        pub = sum(1 for idx in base.values() if len({labels[i] for i in idx}) == 1) / len(base)
        out[f'k{k}_{name}'] = {'test_hits': h, 'acc': round(acc, 4), 'null_acc': round(accn, 4),
                               'frac_base_cells_public': round(pub, 3)}
        print(k, name, out[f'k{k}_{name}'], flush=True)
json.dump(out, open('/Volumes/git/fc-redei/experiments/redei/e3_deep.json', 'w'), indent=1)
