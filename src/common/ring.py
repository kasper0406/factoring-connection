"""Exact arithmetic over R_N = Z/NZ and polynomial quotient rings, with op counters.

Shared infrastructure (owner: infraA).  Operational-safe: nothing in this module
uses or computes a factorisation of N.  The only gcds computed are the
gcd-checked inversions, which raise `NonUnit(divisor)` instead of failing.

Rings
-----
* ``Zmod(N, m)``                         R_N, elements are Python ints in [0, N).
* ``PolyQuotient(base, f, m)``           base[X]/(f), f monic (list low->high).
      - ``truncated(base, K)``           base[X]/(X^K)   (truncated power series)
      - ``binomial_quotient(base, r, c)``base[X]/(X^r - c)
      - ``ext_ring(base, g)``            base[t]/(g(t))  (extension-ring helper)
  Towers are allowed: ``PolyQuotient(ext_ring(Zmod(N), [-c, 0, 1]), ...)``.
* ``FlintQuotient(N, f, m)``             R_N[X]/(f) using python-flint
  (``nmod_poly`` if N < 2^63 else ``fmpz_mod_poly``); same element interface
  for the ops it supports (add/sub/mul/pow/coeffs/is_zero), agrees with
  ``PolyQuotient(Zmod(N), f)`` (see tests/test_common_ring.py).

Element representation: tuples of base elements of length deg f (low->high) for
``PolyQuotient``; flint polynomials for ``FlintQuotient`` (use ``.coeffs(a)``).

Counters (``Metrics`` keys)
  rn_mul, rn_add          -- multiplications / additions in R_N (pure Python rings)
  gcd, inv                -- gcd calls and successful inversions
  q_mul, q_sqr            -- multiplications / squarings in a quotient ring
  q_mul_dimsum            -- sum of quotient dimensions over q_mul+q_sqr
  max dims                -- maxima 'quot_dim', 'ext_degree', 'max_operand_bits'
For ``FlintQuotient`` the R_N-level work is not counted inside flint; use
``q_mul``/``q_mul_dimsum`` (each is one product of two degree<d polys + one
reduction; schoolbook-equivalent cost d^2 R_N-mults).

Flint audit: see ``flint_audit()`` and src/common/README.md.  Multiplication,
mul_low, mulmod by a monic modulus never need inverses or a factorisation;
we verified (a) results agree with the pure-Python reference on random inputs,
(b) contexts for a 4096-bit RSA-shaped modulus are created and used in
milliseconds (any factorisation would be infeasible), and (c) ``ctx.is_prime()``
is only evaluated lazily when called (we never call it in operational code).
"""
from __future__ import annotations

from math import gcd
from typing import Iterable, List, Optional, Sequence, Tuple

from .metrics import Metrics, NullMetrics


class NonUnit(ArithmeticError):
    """Raised when an inversion modulo N hits a nonunit.

    ``divisor`` = gcd(value, N) (> 1).  If 1 < divisor < N it is a proper factor
    of N; divisor == N means the value was 0 mod N.
    """

    def __init__(self, divisor: int, modulus: int, value: Optional[int] = None):
        self.divisor = divisor
        self.modulus = modulus
        self.value = value
        super().__init__(f"nonunit: gcd({value}, {modulus}) = {divisor}")

    @property
    def is_proper(self) -> bool:
        return 1 < self.divisor < self.modulus


def inv_mod(a: int, N: int, m: Optional[Metrics] = None) -> int:
    """gcd-checked inverse of a mod N; raises NonUnit(gcd(a,N)) if not a unit."""
    a %= N
    if m is not None:
        m.inc("gcd")
    g = gcd(a, N)
    if g != 1:
        raise NonUnit(g, N, a)
    if m is not None:
        m.inc("inv")
    return pow(a, -1, N)


# ---------------------------------------------------------------------------
# R_N
# ---------------------------------------------------------------------------
class Zmod:
    """The ring R_N = Z/NZ.  Elements are ints in [0, N)."""

    def __init__(self, N: int, m: Optional[Metrics] = None):
        if N < 2:
            raise ValueError("N must be >= 2")
        self.N = N
        self.m = m if m is not None else NullMetrics()
        self.m.hi("max_operand_bits", N.bit_length())

    # structure
    zero = 0
    one = 1

    def dim(self) -> int:
        return 1

    def total_dim(self) -> int:
        return 1

    def from_int(self, a: int) -> int:
        return a % self.N

    def add(self, a: int, b: int) -> int:
        self.m.inc("rn_add")
        return (a + b) % self.N

    def sub(self, a: int, b: int) -> int:
        self.m.inc("rn_add")
        return (a - b) % self.N

    def neg(self, a: int) -> int:
        return (-a) % self.N

    def mul(self, a: int, b: int) -> int:
        self.m.inc("rn_mul")
        return a * b % self.N

    def inv(self, a: int) -> int:
        return inv_mod(a, self.N, self.m)

    def pow(self, a: int, e: int) -> int:
        if e < 0:
            return self.pow(self.inv(a), -e)
        r = 1
        a %= self.N
        while e:
            if e & 1:
                r = self.mul(r, a)
            e >>= 1
            if e:
                a = self.mul(a, a)
        return r

    def is_zero(self, a: int) -> bool:
        return a % self.N == 0

    def eq(self, a: int, b: int) -> bool:
        return (a - b) % self.N == 0

    def __repr__(self) -> str:
        return f"Zmod({self.N})"


# ---------------------------------------------------------------------------
# base[X]/(f)
# ---------------------------------------------------------------------------
class PolyQuotient:
    """base[X]/(f) for monic f, pure Python reference implementation.

    f is given as a sequence of base elements, low degree first, with
    f[-1] == base.one.  Elements are tuples of length d = deg f.
    Reduction uses only the nonzero lower terms of f (fast for X^r - c, X^K).
    """

    def __init__(self, base, f: Sequence, m: Optional[Metrics] = None, name: str = "X"):
        f = list(f)
        if len(f) < 2:
            raise ValueError("modulus must have degree >= 1")
        if not base.eq(f[-1], base.one):
            raise ValueError("modulus must be monic")
        self.base = base
        self.f = f
        self.d = len(f) - 1
        self.name = name
        self.m = m if m is not None else base.m
        self._low = [(j, f[j]) for j in range(self.d) if not base.is_zero(f[j])]
        self._is_zmod = isinstance(base, Zmod)
        self.zero = tuple(base.zero for _ in range(self.d))
        self.one = self._embed(base.one)
        self.m.hi("quot_dim", self.d)
        self.m.hi("ext_degree", self.total_dim())

    # structure ----------------------------------------------------------
    def dim(self) -> int:
        return self.d

    def total_dim(self) -> int:
        return self.d * self.base.total_dim()

    def _embed(self, b) -> tuple:
        return (b,) + tuple(self.base.zero for _ in range(self.d - 1))

    def from_base(self, b) -> tuple:
        return self._embed(b)

    def from_int(self, a: int) -> tuple:
        return self._embed(self.base.from_int(a))

    def gen(self) -> tuple:
        """The class of X."""
        if self.d == 1:
            return self.reduce([self.base.zero, self.base.one])
        z = [self.base.zero] * self.d
        z[1] = self.base.one
        return tuple(z)

    def from_coeffs(self, coeffs: Iterable) -> tuple:
        return self.reduce(list(coeffs))

    def coeffs(self, a: tuple) -> List:
        return list(a)

    def is_zero(self, a: tuple) -> bool:
        return all(self.base.is_zero(x) for x in a)

    def eq(self, a: tuple, b: tuple) -> bool:
        return all(self.base.eq(x, y) for x, y in zip(a, b))

    # arithmetic ---------------------------------------------------------
    def add(self, a: tuple, b: tuple) -> tuple:
        B = self.base
        return tuple(B.add(x, y) for x, y in zip(a, b))

    def sub(self, a: tuple, b: tuple) -> tuple:
        B = self.base
        return tuple(B.sub(x, y) for x, y in zip(a, b))

    def neg(self, a: tuple) -> tuple:
        return tuple(self.base.neg(x) for x in a)

    def scale(self, b, a: tuple) -> tuple:
        return tuple(self.base.mul(b, x) for x in a)

    def reduce(self, prod: List) -> tuple:
        """Reduce an arbitrary-length coefficient list modulo f."""
        B, d = self.base, self.d
        prod = list(prod)
        if self._is_zmod:
            N = B.N
            nm = 0
            for i in range(len(prod) - 1, d - 1, -1):
                c = prod[i] % N
                if c:
                    s = i - d
                    for j, fj in self._low:
                        prod[s + j] -= c * fj
                        nm += 1
            self.m.inc("rn_mul", nm)
            out = [x % N for x in prod[:d]]
        else:
            for i in range(len(prod) - 1, d - 1, -1):
                c = prod[i]
                if not B.is_zero(c):
                    s = i - d
                    for j, fj in self._low:
                        prod[s + j] = B.sub(prod[s + j], B.mul(c, fj))
            out = prod[:d]
        out += [B.zero] * (d - len(out))
        return tuple(out)

    def mul(self, a: tuple, b: tuple, _kind: str = "q_mul") -> tuple:
        B, d = self.base, self.d
        self.m.inc(_kind)
        self.m.inc("q_mul_dimsum", d)
        n = 2 * d - 1
        if self._is_zmod:
            prod = [0] * n
            nz_b = [(j, y) for j, y in enumerate(b) if y]
            cnt = 0
            for i, x in enumerate(a):
                if x:
                    for j, y in nz_b:
                        prod[i + j] += x * y
                    cnt += len(nz_b)
            self.m.inc("rn_mul", cnt)
            return self.reduce(prod)
        prod = [B.zero] * n
        nz_b = [(j, y) for j, y in enumerate(b) if not B.is_zero(y)]
        for i, x in enumerate(a):
            if B.is_zero(x):
                continue
            for j, y in nz_b:
                prod[i + j] = B.add(prod[i + j], B.mul(x, y))
        return self.reduce(prod)

    def sqr(self, a: tuple) -> tuple:
        return self.mul(a, a, _kind="q_sqr")

    def pow(self, a: tuple, e: int) -> tuple:
        if e < 0:
            raise ValueError("negative exponent")
        r = self.one
        first = True
        for bit in bin(e)[2:]:
            if not first:
                r = self.sqr(r)
            if bit == "1":
                r = a if first else self.mul(r, a)
                first = False
        return r if not first else self.one

    def __repr__(self) -> str:
        return f"PolyQuotient({self.base!r}, deg {self.d})"


def truncated(base, K: int, m: Optional[Metrics] = None) -> PolyQuotient:
    """base[X]/(X^K): power series truncated mod X^K."""
    return PolyQuotient(base, [base.zero] * K + [base.one], m)


def binomial_quotient(base, r: int, c, m: Optional[Metrics] = None) -> PolyQuotient:
    """base[X]/(X^r - c)."""
    f = [base.zero] * r + [base.one]
    f[0] = base.neg(c) if r > 0 else f[0]
    return PolyQuotient(base, f, m)


def ext_ring(base, g: Sequence, m: Optional[Metrics] = None) -> PolyQuotient:
    """Extension-ring helper base[t]/(g(t)), g monic (low->high)."""
    return PolyQuotient(base, g, m, name="t")


def power_of_linear(ring, e: int, a=None, b=None):
    """(a + b*X)^e in `ring` (defaults a=b=1, i.e. (1+X)^e), by repeated squaring."""
    one = ring.one
    x = ring.gen()
    if a is None and b is None:
        lin = ring.add(one, x)
    else:
        base = ring.base
        a = base.one if a is None else a
        b = base.one if b is None else b
        lin = ring.add(ring.from_base(a), ring.scale(b, x))
    return ring.pow(lin, e)


# ---------------------------------------------------------------------------
# python-flint accelerated R_N[X]/(f)
# ---------------------------------------------------------------------------
_NMOD_LIMIT = 1 << 62


class FlintQuotient:
    """R_N[X]/(f) backed by python-flint.  f: list of ints (low->high), monic.

    Uses nmod_poly when N < 2^62 (word-size modulus), else fmpz_mod_poly.
    Only multiplication / mulmod / pow by repeated squaring are used; these
    never invert anything other than the (unit) leading coefficient 1 of f.
    """

    def __init__(self, N: int, f: Sequence[int], m: Optional[Metrics] = None):
        import flint

        f = [int(x) % N for x in f]
        if f[-1] != 1 % N or len(f) < 2:
            raise ValueError("modulus must be monic of degree >= 1")
        self.N = N
        self.d = len(f) - 1
        self.m = m if m is not None else NullMetrics()
        self._small = N < _NMOD_LIMIT
        if self._small:
            self._mk = lambda cs: flint.nmod_poly([int(c) % N for c in cs], N)
        else:
            ctx = flint.fmpz_mod_poly_ctx(N)
            self._mk = lambda cs: ctx([int(c) % N for c in cs])
        self.fpoly = self._mk(f)
        self._trunc_K = self.d if all(c == 0 for c in f[:-1]) else None
        # binomial modulus X^d - c: reduce by folding (no division needed)
        self._binom_c = (-f[0]) % N if all(c == 0 for c in f[1:-1]) else None
        self.zero = self._mk([0])
        self.one = self._mk([1])
        self.m.hi("quot_dim", self.d)
        self.m.hi("max_operand_bits", N.bit_length())

    def dim(self) -> int:
        return self.d

    def total_dim(self) -> int:
        return self.d

    def gen(self):
        return self.reduce([0, 1])

    def from_int(self, a: int):
        return self._mk([a])

    def from_base(self, a: int):
        return self._mk([a])

    def _red(self, p):
        d = self.d
        if p.degree() < d:
            return p
        if self._trunc_K is not None:
            return p.truncate(d) if hasattr(p, "truncate") else p % self.fpoly
        if self._binom_c is not None:
            c = self._binom_c
            while p.degree() >= d:
                hi = p.right_shift(d)
                lo = p.truncate(d) if hasattr(p, "truncate") else p - hi.left_shift(d)
                p = lo + hi * c
            return p
        return p % self.fpoly

    def reduce(self, coeffs: Sequence[int]):
        return self._red(self._mk(list(coeffs)))

    def from_coeffs(self, coeffs):
        return self.reduce(coeffs)

    def coeffs(self, a) -> List[int]:
        cs = [int(c) for c in a.coeffs()]
        return cs + [0] * (self.d - len(cs))

    def is_zero(self, a) -> bool:
        return a.is_zero()

    def eq(self, a, b) -> bool:
        return (a - b).is_zero()

    def add(self, a, b):
        return a + b

    def sub(self, a, b):
        return a - b

    def neg(self, a):
        return -a

    def scale(self, b: int, a):
        return a * (int(b) % self.N)

    def mul(self, a, b, _kind: str = "q_mul"):
        self.m.inc(_kind)
        self.m.inc("q_mul_dimsum", self.d)
        if self._trunc_K is not None:
            return a.mul_low(b, self._trunc_K)
        return self._red(a * b)

    def sqr(self, a):
        return self.mul(a, a, _kind="q_sqr")

    def pow(self, a, e: int):
        r = self.one
        first = True
        for bit in bin(e)[2:]:
            if not first:
                r = self.sqr(r)
            if bit == "1":
                r = a if first else self.mul(r, a)
                first = False
        return r if not first else self.one

    def __repr__(self) -> str:
        return f"FlintQuotient(N={self.N}, deg {self.d})"


def flint_binomial_power(N: int, r: int, c: int, e: int, m: Optional[Metrics] = None) -> List[int]:
    """Coefficients of (1+X)^e mod (N, X^r - c), via FlintQuotient."""
    f = [(-c) % N] + [0] * (r - 1) + [1]
    Q = FlintQuotient(N, f, m)
    return Q.coeffs(power_of_linear(Q, e))


def flint_audit(bits: int = 4096, seed: int = 0) -> dict:
    """Empirical audit that flint does not factor the modulus.

    Builds an RSA-shaped modulus with two random (bits/2)-bit primes (harness
    only: the primes are generated here and discarded), creates contexts and
    performs mulmod / powering.  Returns timings.  A factorisation of a
    4096-bit RSA modulus is infeasible, so millisecond timings show none occurs.
    """
    import random
    import time

    import gmpy2

    rng = random.Random(f"flint_audit:{bits}:{seed}")

    def rp(b):
        while True:
            x = rng.getrandbits(b) | (1 << (b - 1)) | 1
            if gmpy2.is_prime(x):
                return x

    N = rp(bits // 2) * rp(bits // 2)
    t0 = time.perf_counter()
    Q = FlintQuotient(N, [N - 3] + [0] * 63 + [1])
    t1 = time.perf_counter()
    x = power_of_linear(Q, N)
    t2 = time.perf_counter()
    P = PolyQuotient(Zmod(N), [N - 3] + [0] * 63 + [1])
    ref = power_of_linear(P, 1 << 20)
    agree = Q.coeffs(power_of_linear(Q, 1 << 20)) == list(ref)
    return {"bits": N.bit_length(), "ctx_s": t1 - t0, "pow_s": t2 - t1, "agree": agree,
            "nonzero": not Q.is_zero(x)}
