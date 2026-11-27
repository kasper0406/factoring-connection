"""E2a: half-period parabola floor-sum parity = class number mod 4 = Mordell's factorial sign.
F(p) = sum_{y=1}^{(p-1)/2} floor(2y^2/p);  claim (proved in notes): h(-p) = 3 - 2F(p) mod 4,
for primes p = 3 mod 4, p > 3; Mordell: ((p-1)/2)! = (-1)^((h+1)/2) mod p."""
import numpy as np, json
from math import factorial
from cypari2 import Pari
from sympy import primerange
gp = Pari()
ok = okm = n = nm = 0
for p in primerange(7, 200000):
    if p % 4 != 3: continue
    y = np.arange(1, (p-1)//2 + 1, dtype=np.int64)
    F = int(((2*y*y)//p).sum())
    h = int(gp.qfbclassno(-p))
    n += 1; ok += (h % 4 == (3 - 2*F) % 4)
    if p < 20000:
        s = factorial((p-1)//2) % p
        sign = 1 if s == 1 else -1
        nm += 1; okm += (sign == (-1)**((h+1)//2))
res = dict(primes=n, identity_ok=ok, mordell_checked=nm, mordell_ok=okm)
print(res)
json.dump(res, open('classno_parity.json', 'w'))
