"""E1d (HARNESS): exact rules for h(k)<=2 members as functions of (p mod 8, q mod 8); P1a/P1b check."""
import collections
from sympy import primerange
M = 1 << 17; MASK = (1 << (M + 1)) - 1
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
for k in (1, 3, 5, 7, 9, 17):
    A = ser(k); tab = collections.defaultdict(set)
    for p, q in pairs: tab[tuple(sorted((p % 8, q % 8)))].add((A >> (p * q)) & 1)
    print("k=%d" % k, {key: sorted(v) for key, v in sorted(tab.items()) if v != {0}} or "all zero",
          "| ambiguous keys:", [key for key, v in tab.items() if len(v) > 1])
