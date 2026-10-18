"""Cantor's algorithm for genus-2 Jacobians y^2 = F(X), F monic quintic, over Z/N.

Polynomials are lists of ints mod N, lowest degree first, with no trailing zeros
(zero polynomial = []). Every time the algorithm needs the leading coefficient of a
polynomial, ``_norm`` strips coefficients that are 0 mod N and then checks that the top
coefficient c is a unit. If 1 < gcd(c, N) < N it raises ``NonUnit(d)`` with d a proper
divisor of N (rule 3: every nonunit is output). When N is prime nothing is ever raised,
so the same code is the field arithmetic over F_N.

Reduction lemma (see notes/genus2rh_results.md, Lemma 4): if a computation finishes
without NonUnit, then its reduction mod any prime r | N equals the same computation
over F_r, and all polynomials have the same degree mod N and mod r.
"""
from math import gcd

__all__ = ["NonUnit", "Jac", "IDENTITY"]


class NonUnit(Exception):
    def __init__(self, d):
        super().__init__(d)
        self.d = d


IDENTITY = ((1,), ())


class Jac:
    """Arithmetic in J(y^2 = F) over Z/N; F monic of degree 5 given low-first."""

    def __init__(self, N, F):
        self.N = N
        F = [c % N for c in F]
        assert len(F) == 6 and F[5] == 1
        self.F = F
        self.ops = {"compose": 0, "reduce_step": 0, "inv": 0}

    # --- polynomial helpers -------------------------------------------------
    def _norm(self, a):
        N = self.N
        a = [c % N for c in a]
        while a and a[-1] == 0:
            a.pop()
        if a:
            g = gcd(a[-1], N)
            if g != 1:
                raise NonUnit(g)  # 1 < g < N because a[-1] != 0 mod N
        return a

    def _inv(self, c):
        self.ops["inv"] += 1
        return pow(c, -1, self.N)  # caller guarantees c is a unit

    def add(self, a, b):
        n = max(len(a), len(b))
        return self._norm([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                           for i in range(n)])

    def sub(self, a, b):
        n = max(len(a), len(b))
        return self._norm([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
                           for i in range(n)])

    def mul(self, a, b):
        if not a or not b:
            return []
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] += x * y
        return self._norm(r)

    def divmod(self, a, b):
        """Division by b (b already normalised, so its lead is a unit)."""
        N = self.N
        a = list(a)
        db = len(b) - 1
        assert db >= 0
        inv = self._inv(b[-1])
        q = [0] * max(len(a) - db, 0)
        for k in range(len(a) - 1 - db, -1, -1):
            c = (a[k + db] * inv) % N
            q[k] = c
            if c:
                for j in range(db + 1):
                    a[k + j] = (a[k + j] - c * b[j]) % N
        return self._norm(q), self._norm(a[:db] if db > 0 else [])

    def monic(self, a):
        inv = self._inv(a[-1])
        return self._norm([(c * inv) for c in a])

    def xgcd(self, a, b):
        """Monic g and s, t with s a + t b = g (Euclid; every lead tested)."""
        r0, r1 = self._norm(a), self._norm(b)
        s0, s1, t0, t1 = [1], [], [], [1]
        while r1:
            q, r2 = self.divmod(r0, r1)
            r0, r1 = r1, r2
            s0, s1 = s1, self.sub(s0, self.mul(q, s1))
            t0, t1 = t1, self.sub(t0, self.mul(q, t1))
        if not r0:
            return [], s0, t0
        inv = self._inv(r0[-1])
        return (self._norm([c * inv for c in r0]),
                self._norm([c * inv for c in s0]), self._norm([c * inv for c in t0]))

    # --- group law -------------------------------------------------------------
    def compose(self, D1, D2):
        self.ops["compose"] += 1
        u1, v1 = self._norm(D1[0]), self._norm(D1[1])
        u2, v2 = self._norm(D2[0]), self._norm(D2[1])
        F = self.F
        d1, e1, e2 = self.xgcd(u1, u2)
        d, c1, c2 = self.xgcd(d1, self.add(v1, v2))
        s1, s2, s3 = self.mul(c1, e1), self.mul(c1, e2), c2
        u, _ = self.divmod(self.mul(u1, u2), self.mul(d, d))
        num = self.add(self.add(self.mul(self.mul(s1, u1), v2),
                                self.mul(self.mul(s2, u2), v1)),
                       self.mul(s3, self.add(self.mul(v1, v2), F)))
        v, _ = self.divmod(num, d)
        _, v = self.divmod(v, u)
        while len(u) - 1 > 2:
            self.ops["reduce_step"] += 1
            up, _ = self.divmod(self.sub(F, self.mul(v, v)), u)
            up = self.monic(up)
            _, v = self.divmod(self.sub([], v), up)
            u = up
        u = self.monic(u)
        _, v = self.divmod(v, u)
        return (tuple(u), tuple(v))

    def neg(self, D):
        return (D[0], tuple(self.sub([], list(D[1]))))

    def smul(self, n, D):
        """[n]D by left-to-right double-and-add (n >= 0)."""
        R = IDENTITY
        for bit in bin(n)[2:] if n > 0 else "":
            R = self.compose(R, R)
            if bit == "1":
                R = self.compose(R, D)
        return R

    @staticmethod
    def is_identity(D):
        return tuple(D[0]) == (1,)
