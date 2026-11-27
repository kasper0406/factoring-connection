"""E2c: a large prime multiplier a (~2^40) plus pairs with small r determines (a/p) up to the swap:
so the 4-rank family is a full Legendre oracle (./p) = QR_N oracle. Validated vs factored Redei."""
import sys, random, itertools, json
sys.path.insert(0, '/Volumes/git/fc-redei/experiments/redei')
from ledger import r4_symbolic, true_sigma, family
from redei import r4_redei
from sympy import randprime, primerange
random.seed(9)
ok = 0; tot = 0
for it in range(40):
    p = randprime(2**39, 2**40); q = randprime(2**39, 2**40); a = randprime(2**40, 2**41)
    N = p * q
    rs = [3, 5, 7, 11, a]
    sig = true_sigma(p, q, rs)
    fam = family(rs, pairs=True)
    obs = {}
    for f in fam:
        v, D = r4_symbolic(N, f[0], f[1], sig); obs[f] = v
        if it < 5: assert v == r4_redei(D, None) if abs(D) < 10**40 else True
    cons = []
    for e, b, mu in itertools.product((1, -1), repeat=3):
        for xs in itertools.product((1, -1), repeat=len(rs)):
            s = (e, b, mu, dict(zip(rs, xs)))
            if all(r4_symbolic(N, f[0], f[1], s)[0] == obs[f] for f in fam):
                cons.append(s)
    # up to swap: 2 assignments expected
    tot += 1; ok += (len(cons) == 2)
print(json.dumps({'trials': tot, 'exactly_sigma_and_swap': ok}))
