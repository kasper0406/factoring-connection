"""E4d: Kronecker-Hurwitz family. sum_{t^2<4N} H(4N-t^2) = 2N + 2q for N=pq (p<q), so the 2-parts of
Hurwitz class numbers over the family {4N-t^2} give q mod 2^(k-1) -> factor (Coppersmith at k ~ log2 N/4).
Check the identity (PARI qfbhclassno, exact rationals), and how many family members are nonzero mod 2^k
(terms with many prime factors vanish mod 2^k: 2-rank >= k)."""
import sys, json, random
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from redei import pari
from sympy import randprime
from fractions import Fraction
random.seed(4)
out = []
for it in range(40):
    p = randprime(200, 3000); q = randprime(3000, 20000)
    N = p * q
    S = Fraction(0); nz = {3: 0, 5: 0, 8: 0}; T = 0
    t = 0
    while t * t < 4 * N:
        h = pari.qfbhclassno(4 * N - t * t)
        hf = Fraction(int(pari.numerator(h)), int(pari.denominator(h)))
        w = 1 if t == 0 else 2
        S += w * hf; T += w
        for k in nz:
            if (hf * 6) % (2**k) != 0 if hf.denominator == 1 else True:
                nz[k] += w
        t += 1
    out.append({'N': N, 'ok': S == 2 * N + 2 * q, 'family': T, 'nonzero_mod_2^k(6H)': nz})
print(all(o['ok'] for o in out), out[:3])
json.dump(out, open('/Volumes/git/fc-redei/experiments/redei/e4_kh.json', 'w'), indent=0)
