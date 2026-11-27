"""E1f (HARNESS): is the family's entropy about the factors or about public N? Miller-Madow conditional entropy of
v_H(N) = (a(N,Delta^k))_{h(k)<=H} given (i) public N mod 2^J, (ii) local unordered (p mod 2^J, q mod 2^J)."""
import sys, math, collections, json
import numpy as np
from sympy import primerange
M = 1 << 20; K = 512; MASK = (1 << (M + 1)) - 1
def sp(s):
    o, m = [], 1
    while (m * m << s) <= M: o.append(m * m << s); m += 2
    return o
def mul(A, i):
    r = 0
    for s in sp(i): r ^= A << s
    return r & MASK
def h(k):
    n3 = n5 = 0
    for i in range(1, k.bit_length()):
        if k >> i & 1:
            j = (i + 1) // 2
            if i % 2: n3 += 1 << (j - 1)
            else: n5 += 1 << (j - 1)
    return n3 + n5
primes = list(primerange(3, M // 3 + 1)); P, Q = [], []
for i, p in enumerate(primes):
    for q in primes[i + 1:]:
        if p * q > M: break
        P.append(p); Q.append(q)
P = np.array(P); Q = np.array(Q); Nn = P * Q
ks = list(range(1, K, 2)); D = {1: mul(1, 0)}; V = np.zeros((len(ks), len(Nn)), np.uint8)
for t, k in enumerate(ks):
    if k > 1: D[k] = mul(D[k - (1 << (k.bit_length() - 1))], k.bit_length() - 1)
    b = D[k].to_bytes((M + 8) // 8 + 1, "little")
    V[t] = np.unpackbits(np.frombuffer(b, np.uint8), bitorder="little")[: M + 1][Nn]
def condH(keys, vecs):
    g = collections.defaultdict(collections.Counter)
    for kk, v in zip(keys, vecs): g[kk][v] += 1
    tot = len(vecs); Hs = 0.0
    for cn in g.values():
        n = sum(cn.values()); Hc = -sum(c / n * math.log2(c / n) for c in cn.values()) + (len(cn) - 1) / (2 * n * math.log(2))
        Hs += n / tot * Hc
    return round(Hs, 3), round(tot / len(g), 1)
out = []
for H in range(2, 9):
    rows = [t for t, k in enumerate(ks) if h(k) <= H]
    vecs = [V[rows][:, c].tobytes() for c in range(len(Nn))]
    row = {"H": H}
    for J in (4, 8, 12):
        m = 1 << J
        row[f"pub{J}"] = condH((Nn % m).tolist(), vecs)
        row[f"loc{J}"] = condH([tuple(sorted(x)) for x in zip((P % m).tolist(), (Q % m).tolist())], vecs)
    print(row); out.append(row)
json.dump(out, open("/Volumes/git/fc-coeffam/experiments/coeffam/e1f_public_vs_local.json", "w"), indent=1)
