"""First collisions between positive squarefree representatives in the 131 wheel."""
import math
import numpy as np
from benchmark import ODD, signature

pbase = ODD[:31]
limit = 1 << 22
ok = np.ones(limit+1, dtype=bool)
ok[0] = False
for p in pbase:
    ok[p::p] = False
for d in range(2, math.isqrt(limit)+1):
    ok[d*d::d*d] = False
ns = np.flatnonzero(ok)
codes = np.zeros(len(ns), dtype=np.uint32)
for i,p in enumerate(pbase):
    qr = np.ones(p, dtype=np.uint32)
    qr[np.arange(1,p)**2 % p] = 0
    codes |= qr[ns % p] << i
for mode in ('all units', 'odd units', 'odd units with mod 8'):
    seen = {}
    for n,s in zip(ns.tolist(),codes.tolist()):
        if mode != 'all units' and not n % 2:
            continue
        key = (s, n%8 if mode.endswith('8') else 0)
        if key in seen:
            u = seen[key]
            assert signature(u,pbase) == signature(n,pbase)
            print(mode, u, n, s)
            break
        seen[key] = n
