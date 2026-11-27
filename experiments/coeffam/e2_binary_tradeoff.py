"""E2 (factor-blind evaluator, harness only samples N): binary members Delta^{1+2^m}: a = parity #{x,y odd>0: x^2+2^m y^2=N}.
Measure Pr[a=1] and entropy vs enumeration cost sqrt(N/2^m)/2 at 40 bits."""
import random, math, json
from math import isqrt
import numpy as np
from sympy import randprime
random.seed(7)
bits = 40
Ns = []
while len(Ns) < 3000:
    p = randprime(2 ** 19, 2 ** 20); q = randprime(2 ** 19, 2 ** 20)
    if p != q: Ns.append(p * q)
def a_m(N, m):
    ymax = isqrt(N >> m)
    y = np.arange(1, ymax + 1, 2, dtype=np.int64)
    r = N - (y * y << m)
    s = np.floor(np.sqrt(r.astype(np.float64))).astype(np.int64)
    ok = ((s * s == r) | ((s + 1) * (s + 1) == r) | ((s - 1) * (s - 1) == r))
    x = np.where((s + 1) * (s + 1) == r, s + 1, np.where((s - 1) * (s - 1) == r, s - 1, s))
    return int(np.sum(ok & (x % 2 == 1)) & 1), len(y)
rows = []
for m in range(3, 39, 1):
    sample = Ns if m >= 14 else Ns[:600]
    vals = []; cost = 0
    for N in sample:
        v, c = a_m(N, m); vals.append(v); cost += c
    pi = sum(vals) / len(vals)
    H = 0 if pi in (0, 1) else -(pi * math.log2(pi) + (1 - pi) * math.log2(1 - pi))
    rows.append(dict(m=m, n=len(sample), pr1=round(pi, 5), entropy_bits=round(H, 5), mean_cost=cost / len(sample),
                     cost_over_sqrtN=round(cost / len(sample) / 2 ** 20, 6)))
    print(rows[-1])
json.dump(rows, open("/Volumes/git/fc-coeffam/experiments/coeffam/e2_binary_tradeoff.json", "w"), indent=1)
