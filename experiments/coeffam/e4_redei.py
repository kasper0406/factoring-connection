"""E4 (HARNESS): h(-4N) mod 16 for N = pq = 1 mod 4 vs ((-1/p),(2/p),(2/q),(p/q)); and h(-4Nm) mod 8 vs (m/p)."""
import collections, json
import cypari2
from sympy import primerange, legendre_symbol as L
pari = cypari2.Pari()
prs = list(primerange(3, 400))
tab = collections.defaultdict(collections.Counter)
for i, p in enumerate(prs):
    for q in prs[i + 1:]:
        N = p * q
        if N % 4 != 1: continue
        h = int(pari.qfbclassno(-4 * N))
        key = (N % 8, L(-1 % p, p), L(2, p), L(2, q), L(p, q))
        tab[key][h % 16] += 1
out = {}
for key, cnt in sorted(tab.items()):
    m8 = sorted(set(v % 8 for v in cnt)); out[str(key)] = dict(mod16=dict(cnt), mod8_values=m8)
    print("N%8, (-1/p),(2/p),(2/q),(p/q) =", key, " h mod 8 values:", m8, " (mod 16:", dict(cnt), ")")
json.dump(out, open("/Volumes/git/fc-coeffam/experiments/coeffam/e4_redei.json", "w"), indent=1)
