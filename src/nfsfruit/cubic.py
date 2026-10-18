"""OPERATIONAL: arithmetic in Z[alpha] = Z[x]/(f), f monic of small degree d, and an exact
square-root routine (inert-prime Hensel lifting + exact verification).

Elements are lists of d ints (coefficients of 1, alpha, ..., alpha^{d-1})."""
from __future__ import annotations
from typing import List, Optional
import gmpy2
from flint import nmod_poly

class Order:
    def __init__(self, f: List[int]):
        """f = [f0, f1, ..., f_{d-1}, 1] (monic, low -> high)."""
        assert f[-1] == 1
        self.f = [int(c) for c in f]
        self.d = len(f) - 1
        self.R = 1 + max(abs(c) for c in self.f[:-1])     # Cauchy root bound

    def mul(self, a, b, mod=None):
        d = self.d
        prod = [0] * (2 * d - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    prod[i + j] += x * y
        f = self.f
        for k in range(2 * d - 2, d - 1, -1):
            c = prod[k]
            if c:
                prod[k] = 0
                for t in range(d):
                    prod[k - d + t] -= c * f[t]
        out = prod[:d]
        if mod is not None:
            out = [x % mod for x in out]
        return out

    def linear(self, a, b):
        """a - b*alpha"""
        e = [0] * self.d
        e[0] = a
        e[1] = -b
        return e

    def deriv_at_alpha(self):
        """f'(alpha) as an element."""
        d = self.d
        e = [0] * d
        for i in range(1, d + 1):
            c = i * self.f[i]
            if i - 1 < d:
                e[i - 1] += c
        return e

    def prod(self, elems):
        elems = list(elems)
        if not elems:
            one = [0] * self.d
            one[0] = 1
            return one
        while len(elems) > 1:
            nxt = []
            for i in range(0, len(elems) - 1, 2):
                nxt.append(self.mul(elems[i], elems[i + 1]))
            if len(elems) % 2:
                nxt.append(elems[-1])
            elems = nxt
        return elems[0]

    def eval_mod(self, e, m, n):
        acc = 0
        for c in reversed(e):
            acc = (acc * m + c) % n
        return acc

    def powmod(self, a, k, mod):
        res = [0] * self.d
        res[0] = 1
        base = [x % mod for x in a]
        while k:
            if k & 1:
                res = self.mul(res, base, mod)
            base = self.mul(base, base, mod)
            k >>= 1
        return res

    def inert_primes(self, start: int, count: int):
        """Primes l = 3 (mod 4), l >= start, with f irreducible mod l."""
        out = []
        l = int(gmpy2.next_prime(start))
        while len(out) < count:
            if l % 4 == 3:
                fac = nmod_poly(self.f, l).factor()[1]
                if len(fac) == 1 and fac[0][1] == 1 and fac[0][0].degree() == self.d:
                    out.append(l)
            l = int(gmpy2.next_prime(l))
        return out

    def sqrt(self, g, inert_l: int) -> Optional[List[int]]:
        """Exact square root of g in Z[alpha] if one exists, else None.
        Requires d odd or generally l^d = 3 (mod 4); inert_l from inert_primes."""
        l = inert_l
        q = l ** self.d
        assert q % 4 == 3
        gl = [x % l for x in g]
        if not any(gl):
            raise ValueError("l divides g; choose another inert prime")
        s = self.powmod(gl, (q + 1) // 4, l)
        if self.mul(s, s, l) != gl:
            return None                       # non-square mod an inert prime => non-square
        w = self.powmod(s, q - 2, l)          # w = s^{-1} = g^{-1/2} mod l
        # Exact coefficient bound for a root v in Z[alpha] (notes/nfsfruit_results.md, Lemma S):
        # |v_k| <= d (d-1)! R^{(d-1)^2} (d |g|_inf R^{d-1})^{1/2}, using |det V| = |disc f|^{1/2} >= 1.
        d = self.d
        gbits = max(abs(x) for x in g).bit_length()
        cap_bits = gbits // 2 + ((d - 1) ** 2 + d) * self.R.bit_length() + 2 * d + 8
        k = 1
        mod = l
        inv2 = None
        while True:
            k2 = 2 * k
            mod2 = l ** k2
            inv2 = (mod2 + 1) // 2
            gw2 = self.mul(g, self.mul(w, w, mod2), mod2)
            t = [(-x) % mod2 for x in gw2]
            t[0] = (t[0] + 1) % mod2          # 1 - g w^2
            corr = self.mul(w, t, mod2)
            w = [(w[i] + corr[i] * inv2) % mod2 for i in range(self.d)]
            k, mod = k2, mod2
            if mod.bit_length() > cap_bits + 2:
                v = self.mul(g, w, mod)
                v = [x - mod if x > mod // 2 else x for x in v]
                if self.mul(v, v) == list(g):
                    return v
                return None
