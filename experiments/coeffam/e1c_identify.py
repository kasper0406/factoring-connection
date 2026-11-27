"""E1c (HARNESS): which single 2-adic feature explains a(pq, Delta^k) for h(k)=3, and the exact rule."""
import sys, itertools, collections
import numpy as np, cypari2
from sympy import primerange
pari = cypari2.Pari()
M = 1 << 17
MASK = (1 << (M + 1)) - 1
def ser(k):
    A = 1
    for i in range(k.bit_length()):
        if k >> i & 1:
            r, m = 0, 1
            while (m * m << i) <= M: r ^= A << (m * m << i); m += 2
            A = r & MASK
    return A
primes = list(primerange(3, M // 3 + 1))
pairs = [(p, q) for i, p in enumerate(primes) for q in primes[i + 1:] if p * q <= M]
POLYS = {"x4-2": "x^4-2", "x4+2": "x^4+2", "x4-2x2-1": "x^4-2*x^2-1", "x4+2x2-1": "x^4+2*x^2-1"}
F = {nm: {r: len(pari.polrootsmod(pari(s), r)) for r in primes} for nm, s in POLYS.items()}
for k in (11, 13, 67):
    A = ser(k)
    vals = [(A >> (p * q)) & 1 for p, q in pairs]
    for nm in POLYS:
        tab = collections.defaultdict(set)
        for (p, q), v in zip(pairs, vals):
            tab[tuple(sorted(((p % 16, F[nm][p]), (q % 16, F[nm][q]))))].add(v)
        bad = sum(1 for s in tab.values() if len(s) > 1)
        print(k, nm, "conflicts", bad, "keys", len(tab))
        if bad == 0:
            # print rule restricted to N = 1 mod 16 stratum, where it is hidden
            for key, s in sorted(tab.items()):
                if (key[0][0] * key[1][0]) % 16 in (1, 9): print("   ", key, s)
            break

print("=== k=11 rule, grouped by public N mod 16; feature = #roots of x^4-2 mod r")
A = ser(11); tab = collections.defaultdict(set)
for (p, q) in pairs:
    tab[tuple(sorted(((p % 16, F["x4-2"][p]), (q % 16, F["x4-2"][q]))))].add((A >> (p * q)) & 1)
byN = collections.defaultdict(list)
for key, s in tab.items(): byN[(key[0][0] * key[1][0]) % 16].append((key, s))
for c, lst in sorted(byN.items()):
    vs = set().union(*[s for _, s in lst])
    if len(vs) > 1:
        print("N=%d mod 16:" % c, sorted(lst))
