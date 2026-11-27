"""E5b (HARNESS): P_N(1) = [q^N] Delta/(1-Delta) mod 2 (= Tr H_{-7N^2} mod 2 by E5). Is it public (function of
N mod 2^J) or factor-dependent? Exact series via S_n = sum_{m odd, m^2<=n} (1+S)_{n-m^2}."""
import collections, json
import numpy as np
from sympy import primerange
M = 1 << 16
S = np.zeros(M + 1, np.uint8); one = np.zeros(M + 1, np.uint8); one[0] = 1
sq = [m * m for m in range(1, 400, 2) if m * m <= M]
for n in range(1, M + 1):
    t = 0
    for s in sq:
        if s > n: break
        t ^= one[n - s] ^ S[n - s]
    S[n] = t
primes = list(primerange(3, M // 3 + 1))
pairs = [(p, q) for i, p in enumerate(primes) for q in primes[i + 1:] if p * q <= M and p != 7 and q != 7]
out = {}
for J in (4, 6, 8, 10):
    m = 1 << J
    pub = collections.defaultdict(set); loc = collections.defaultdict(set)
    for p, q in pairs:
        v = int(S[p * q]); pub[(p * q) % m].add(v); loc[tuple(sorted((p % m, q % m)))].add(v)
    out[J] = dict(pub_ambiguous=sum(len(s) > 1 for s in pub.values()), pub_classes=len(pub),
                  loc_ambiguous=sum(len(s) > 1 for s in loc.values()), loc_classes=len(loc))
    print(J, out[J])
print("fraction P_N(1)=1 on semiprimes:", sum(int(S[p * q]) for p, q in pairs) / len(pairs), len(pairs))
json.dump(out, open("/Volumes/git/fc-coeffam/experiments/coeffam/e5b_PN1.json", "w"), indent=1)
