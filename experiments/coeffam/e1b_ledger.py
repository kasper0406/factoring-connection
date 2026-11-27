"""E1b (HARNESS): matrix V[k, (p,q)] = a(pq, Delta^k) mod 2 for odd k < K, semiprimes pq <= M.
(1) per k: is it a function of (p mod 2^j, q mod 2^j)? of that plus non-abelian 2-adic Frobenius features?
(2) ledger: entropy of (a(pq,Delta^k))_{h(k)<=H} inside a public class (N mod 2^J), vs H."""
import sys, json, math, collections
import numpy as np
from sympy import primerange

M = int(sys.argv[1]); K = int(sys.argv[2]); OUT = sys.argv[3]
MASK = (1 << (M + 1)) - 1
def sparse_sq(s):
    out, m = [], 1
    while (m * m << s) <= M: out.append(m * m << s); m += 2
    return out
SP = {}
def mul_sparse(A, i):
    if i not in SP: SP[i] = sparse_sq(i)
    r = 0
    for s in SP[i]: r ^= (A << s)
    return r & MASK
def h(k):
    n3 = n5 = 0
    for i in range(1, k.bit_length()):
        if k >> i & 1:
            j = (i + 1) // 2
            if i % 2: n3 += 1 << (j - 1)
            else: n5 += 1 << (j - 1)
    return n3, n5
def tobits(A):
    b = A.to_bytes((M + 8) // 8 + 1, "little")
    return np.unpackbits(np.frombuffer(b, dtype=np.uint8), bitorder="little")[: M + 1]

primes = list(primerange(3, M // 3 + 1))
P, Q = [], []
for i, p in enumerate(primes):
    for q in primes[i + 1:]:
        if p * q > M: break
        P.append(p); Q.append(q)
P = np.array(P, dtype=np.int64); Q = np.array(Q, dtype=np.int64); Nn = P * Q
ks = list(range(1, K, 2))
D = {1: mul_sparse(1, 0)}
V = np.zeros((len(ks), len(Nn)), dtype=np.uint8)
for idx, k in enumerate(ks):
    if k > 1:
        top = k.bit_length() - 1
        D[k] = mul_sparse(D[k - (1 << top)], top)
    V[idx] = tobits(D[k])[Nn]
print("semiprimes", len(Nn), file=sys.stderr)

# non-abelian 2-adic features of a prime r: root counts mod r of polys whose splitting fields are 2-groups
# unramified outside 2 (D4 / Q16-type), beyond Q(zeta_{2^j}).
POLYS = {"x4-2": [1, 0, 0, 0, -2], "x4+2": [1, 0, 0, 0, 2], "x4-2x2-1": [1, 0, -2, 0, -1],
         "x4+2x2-1": [1, 0, 2, 0, -1], "x4-4x2+2": [1, 0, -4, 0, 2], "x8-2": [1, 0, 0, 0, 0, 0, 0, 0, -2]}
def nroots(c, r):
    return sum(1 for x in range(r) if sum(ci * pow(x, len(c) - 1 - e, r) for e, ci in enumerate(c)) % r == 0) if r < 3000 else None
import cypari2
pari = cypari2.Pari()
PF = [pari("Pol(%s)" % str(c)) for c in POLYS.values()]
def nroots_fast(i, r):
    return int(pari.matsize(pari.polrootsmod(PF[i], r))[1]) if False else len(pari.polrootsmod(PF[i], r))
uniq = sorted(set(P.tolist()) | set(Q.tolist()))
feat = {r: tuple(nroots_fast(i, r) for i in range(len(PF))) for r in uniq}

def det_conflicts(keys, vals):
    seen = {}; bad = set()
    for kk, v in zip(keys, vals):
        if kk in seen and seen[kk] != v: bad.add(kk)
        else: seen.setdefault(kk, v)
    return len(bad), len(seen)

res = []
Pl, Ql = P.tolist(), Q.tolist()
for idx, k in enumerate(ks):
    vals = V[idx].tolist()
    n3, n5 = h(k); row = dict(k=k, n3=n3, n5=n5, ones=int(V[idx].sum()))
    if row["ones"] and n3 + n5 <= 6:
        for j in (4, 6, 8, 10):
            m = 1 << j
            row[f"loc{j}"] = det_conflicts([tuple(sorted((p % m, q % m))) for p, q in zip(Pl, Ql)], vals)
            row[f"locF{j}"] = det_conflicts([tuple(sorted(((p % m,) + feat[p], (q % m,) + feat[q]))) for p, q in zip(Pl, Ql)], vals)
    res.append(row)

# ledger: H -> entropy of the vector over {k: h(k) <= H} within public class N mod 2^J (J=10), sample-weighted
def entropy(cnt):
    t = sum(cnt.values()); return -sum(c / t * math.log2(c / t) for c in cnt.values())
J = 8
pub = (Nn % (1 << J)).tolist()
ledger = []
for H in range(1, 31):
    rows = [i for i, k in enumerate(ks) if sum(h(k)) <= H]
    if max(ks[i] for i in rows) >= K - 1 and H > 1: pass
    sub = V[rows]
    vecs = [sub[:, c].tobytes() for c in range(sub.shape[1])]
    byc = collections.defaultdict(collections.Counter)
    for c, v in zip(pub, vecs): byc[c][v] += 1
    # weighted conditional entropy H(vec | N mod 2^J), only classes with >= 200 samples
    tot = sum(sum(cn.values()) for cn in byc.values() if sum(cn.values()) >= 200)
    Hc = sum(sum(cn.values()) / tot * entropy(cn) for cn in byc.values() if sum(cn.values()) >= 200)
    complete = all(k < K for k in [kk for kk in range(1, 4 * K, 2) if sum(h(kk)) <= H])
    big = max(byc.values(), key=lambda cn: sum(cn.values()))
    ledger.append(dict(H=H, nk=len(rows), dimA=(H + 1) * (H + 2) // 2 - 1, condH_bits=round(Hc, 3), complete=complete,
                       bigclass_n=sum(big.values()), bigclass_distinct=len(big)))
    print(ledger[-1], file=sys.stderr)
json.dump(dict(M=M, K=K, nsemi=len(Nn), polys=list(POLYS), rows=res, ledger=ledger), open(OUT, "w"), indent=1)
