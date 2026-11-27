"""E1: level-1 mod-2 family Delta^k at semiprimes. HARNESS (uses factorisations to analyse dependence).
Series over F_2 stored as Python-int bitsets: bit n = coefficient of q^n."""
import json, sys, math, itertools, collections
from sympy import primerange

M = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 17
K = int(sys.argv[2]) if len(sys.argv) > 2 else 1 << 12
MASK = (1 << (M + 1)) - 1

def sparse_sq(shift):          # Delta^{2^shift} = sum_{m odd} q^{2^shift m^2}
    out, m = [], 1
    while (m * m << shift) <= M:
        out.append(m * m << shift); m += 2
    return out

SP = {}
def mul_sparse(A, i):
    if i not in SP: SP[i] = sparse_sq(i)
    r = 0
    for s in SP[i]:
        r ^= (A << s)
    return r & MASK

def h(k):
    n3 = n5 = 0
    for i in range(1, k.bit_length()):
        if k >> i & 1:
            j = (i + 1) // 2
            if i % 2: n3 += 1 << (j - 1)
            else: n5 += 1 << (j - 1)
    return n3, n5

# Delta^k for odd k < K via Delta^k = Delta^{k - 2^top} * Delta^{2^top}
D = {1: mul_sparse(1, 0)}       # Delta itself
for k in range(3, K, 2):
    top = k.bit_length() - 1
    D[k] = mul_sparse(D[k - (1 << top)], top)

primes = [p for p in primerange(3, M // 3 + 1)]
semis = []
for i, p in enumerate(primes):
    for q in primes[i + 1:]:
        if p * q > M: break
        semis.append((p, q))
print("semiprimes", len(semis), file=sys.stderr)

def det_by(keyf, vals):
    """is vals (dict (p,q)->bit) a function of keyf(p,q)? returns (#conflicting keys, #keys)"""
    seen = {}; bad = set()
    for pq, v in vals.items():
        kk = keyf(*pq)
        if kk in seen and seen[kk] != v: bad.add(kk)
        seen.setdefault(kk, v)
    return len(bad), len(seen)

res = []
for k in range(1, K, 2):
    A = D[k]
    vals = {(p, q): (A >> (p * q)) & 1 for (p, q) in semis}
    ones = sum(vals.values())
    n3, n5 = h(k)
    row = dict(k=k, n3=n3, n5=n5, h=n3 + n5, ones=ones)
    if ones:
        for j in (3, 4, 5, 6, 7, 8):
            m = 1 << j
            row[f"pub{j}"] = det_by(lambda p, q: (p * q) % m, vals)[0]
            row[f"loc{j}"] = det_by(lambda p, q: tuple(sorted((p % m, q % m))), vals)[0]
    res.append(row)
json.dump(dict(M=M, K=K, nsemi=len(semis), rows=res), open(sys.argv[3] if len(sys.argv) > 3 else
          "/Volumes/git/fc-coeffam/experiments/coeffam/e1_rows.json", "w"))
for r in res:
    if r["ones"]:
        print(r)
