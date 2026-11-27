"""E6: how many family bits pin p? Legendre fingerprint ((r/p))_{r<=B} (with (-1/p),(2/p)) on ALL primes in
[2^b, 2^(b+1)]: smallest B making it injective, vs log2(#primes)."""
import json, math
from sympy import primerange, legendre_symbol
out = {}
for b in (14, 16, 18, 20):
    ps = list(primerange(2**b, 2**(b + 1)))
    rs = [-1, 2] + list(primerange(3, 400))
    keys = {p: () for p in ps}
    for i, r in enumerate(rs):
        for p in ps:
            v = (p % 4 == 1) if r == -1 else (legendre_symbol(r % p, p) == 1)
            keys[p] = keys[p] + (v,)
        if len(set(keys.values())) == len(ps):
            out[b] = {'nprimes': len(ps), 'log2_nprimes': round(math.log2(len(ps)), 2), 'chars_needed': i + 1, 'B': r}
            break
    print(b, out[b], flush=True)
json.dump(out, open('/Volumes/git/fc-redei/experiments/redei/e6_fingerprint.json', 'w'))
