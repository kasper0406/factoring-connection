"""E3 (HARNESS): exponent-N Hasse-Witt HW_N = [x^{N-1}] f^{(N-1)/2} mod N, f = x^3+ax+b.
Check HW_N = H_p * C_q (mod p), H_p = [x^{p-1}] f^{(p-1)/2} mod p (Hasse invariant = a_p mod p),
C_q = [x^{q-1}] f^{(q-1)/2} in Z reduced mod p. Also Legendre-family truncation F_{N-1}(lam) = sum_{k<N} C(2k,k)^2 (lam/16)^k."""
import json, random
import flint
from sympy import primerange
def coef(n_mod, e, j, a, b):
    f = flint.nmod_poly([b % n_mod, a % n_mod, 0, 1], n_mod)
    return int((f ** e)[j])
random.seed(1)
a, b = -1, 1
prs = list(primerange(11, 140))
ok = tot = splits = 0; bad = []
for i, p in enumerate(prs):
    for q in prs[i + 1:]:
        N = p * q; M = (N - 1) // 2
        hw = coef(N, M, N - 1, a, b)
        Hp = coef(p, (p - 1) // 2, p - 1, a, b); Cq = coef(p, (q - 1) // 2, q - 1, a, b)
        Hq = coef(q, (q - 1) // 2, q - 1, a, b); Cp = coef(q, (p - 1) // 2, p - 1, a, b)
        good = (hw - Hp * Cq) % p == 0 and (hw - Hq * Cp) % q == 0
        ok += good; tot += 1
        if not good and len(bad) < 5: bad.append((p, q))
        from math import gcd
        g = gcd(hw, N); splits += 1 < g < N
# Legendre family truncation mod N with lam = 3 (no 16-inverse issues for odd N)
def F(n_mod, T, lam):
    s, c, inv16 = 0, 1, pow(16, -1, n_mod)
    t = 1
    from math import comb
    for k in range(T):
        s = (s + comb(2 * k, k) ** 2 % n_mod * t) % n_mod; t = t * lam * inv16 % n_mod
    return s
lok = ltot = 0
for i, p in enumerate(prs[:12]):
    for q in prs[i + 1:14]:
        N = p * q
        lhs = F(N, N, 3) % p
        rhs = F(p, p, 3) * F(p, q, 3) % p     # predicted: H_p(lam) * (integer truncation at q, reduced mod p)
        lok += lhs == rhs; ltot += 1
res = dict(curve=[a, b], pairs=tot, product_formula_holds=ok, counterexamples=bad, gcd_splits=splits,
           legendre_pairs=ltot, legendre_product_formula_holds=lok)
print(res)
json.dump(res, open("/Volumes/git/fc-coeffam/experiments/coeffam/e3_hw.json", "w"), indent=1)
