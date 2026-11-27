"""HARNESS: validate factor-blind evaluators against rules derived at N <= 2^17 (E1c/E1d), at 36-40 bits; time them."""
import sys, time, random, ast, json, collections
sys.path.insert(0, "/Volumes/git/fc-coeffam/experiments/coeffam")
import coeffam_eval as ev
import cypari2
pari = cypari2.Pari()
from sympy import randprime, primerange
# derive the k=11 rule table from small data
M = 1 << 17; MASK = (1 << (M + 1)) - 1
def ser(k):
    A = 1
    for i in range(k.bit_length()):
        if k >> i & 1:
            r, m = 0, 1
            while (m * m << i) <= M: r ^= A << (m * m << i); m += 2
            A = r & MASK
    return A
def feat(r): return (r % 16, len(pari.polrootsmod(pari("x^4-2"), r)))
A11 = ser(11); primes = list(primerange(3, M // 3 + 1)); tab = collections.defaultdict(set)
for i, p in enumerate(primes):
    for q in primes[i + 1:]:
        if p * q > M: break
        tab[tuple(sorted((feat(p), feat(q))))].add((A11 >> (p * q)) & 1)
assert all(len(v) == 1 for v in tab.values())
rule11 = {k: v.pop() for k, v in tab.items()}
rule9 = lambda p, q: int(p % 8 == 3 and q % 8 == 3)
random.seed(20261127)
out = []
for bits in (36, 38, 40):
    agree9 = agree11 = n = 0; t9 = t11 = 0.0; hits = collections.Counter()
    for _ in range(12):
        while True:
            p = randprime(2 ** (bits // 2 - 1), 2 ** (bits // 2)); q = randprime(2 ** (bits // 2 - 1), 2 ** (bits // 2))
            if p != q and (p * q) % 16 in (3, 9, 11): break   # strata where k=9 / k=11 are non-constant
        N = p * q
        s = time.time(); v9 = ev.eval_d9(N); t9 += time.time() - s
        s = time.time(); v11 = ev.eval_d11(N); t11 += time.time() - s
        key = tuple(sorted((feat(p), feat(q))))
        agree9 += v9 == rule9(min(p, q), max(p, q)) if True else 0
        agree11 += (key in rule11 and v11 == rule11[key]); n += 1; hits[(v9, v11)] += 1
    row = dict(bits=bits, n=n, agree_rule9=agree9, agree_rule11=agree11, sec_per_eval_d9=t9 / n,
               sec_per_eval_d11=t11 / n, value_counts={str(k): v for k, v in hits.items()})
    print(row); out.append(row)
json.dump(out, open("/Volumes/git/fc-coeffam/experiments/coeffam/e1e_evaluator.json", "w"), indent=1)
