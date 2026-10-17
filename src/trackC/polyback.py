"""Instrumented dense polynomial arithmetic over R_N = Z/NZ (track C, owner factC).

Operational-safe: no factor information is used.  The only divisions are by
*monic* polynomials (via Newton inversion of the reversed divisor, whose
constant term is 1), so no inverse of a possibly-nonunit scalar is ever taken
inside polynomial code.  Scalar inverses go through ``src.common.ring.inv_mod``
(gcd-checked, raises ``NonUnit``).

Flint audit (factC, 2026-10-17, python-flint 0.9.0), complementing
src/common/README.md:
  * ``nmod_poly`` divmod by a polynomial whose leading coefficient is a nonunit
    mod a composite N ABORTS THE PROCESS (C-level ``FLINT exception (Impossible
    inverse)``; not catchable from Python).  ``fmpz_mod_poly`` raises ValueError
    instead and does not report the gcd.  Hence we never call flint division;
    we implement monic division ourselves on top of multiplication.
  * multiplication / ``mul_low`` / add / sub agree with the pure-Python
    Karatsuba backend (tests/test_trackC_baselines.py).

Backends
--------
``PolyRing(N, backend="flint")``: python-flint ``nmod_poly`` (N < 2^62) or
``fmpz_mod_poly``.  Counts polynomial-level operations only.
``PolyRing(N, backend="py")``: pure Python lists, Karatsuba multiplication with
an exact count of R_N multiplications (``rn_mul``) and additions (``rn_add``).

Counters (Metrics keys)
  pmul            number of polynomial products (incl. mullow)
  pmul_len        sum over products of (len a + len b)
  pmul_model      sum of M(n) = n * ceil(log2 n), n = len a + len b  (FFT model)
  rn_mul, rn_add  exact scalar counts (py backend; also scalar ops done outside
                  polynomial products in both backends, via ``smul``)
  maxima: max_degree, stored_coeffs (largest single polynomial length)
"""
from __future__ import annotations

from typing import List, Optional, Sequence

from src.common.metrics import Metrics, NullMetrics

_NMOD_LIMIT = 1 << 62
_KARA_CUTOFF = 24


def _clog2(n: int) -> int:
    return max(1, (n - 1).bit_length())


class PolyRing:
    """Dense univariate polynomials over Z/NZ with op counters."""

    def __init__(self, N: int, backend: str = "flint", m: Optional[Metrics] = None):
        if N < 2:
            raise ValueError("N >= 2")
        self.N = N
        self.backend = backend
        self.m = m if m is not None else NullMetrics()
        self.m.hi("max_operand_bits", N.bit_length())
        if backend == "flint":
            import flint

            if N < _NMOD_LIMIT:
                self._mk = lambda cs: flint.nmod_poly([int(c) % N for c in cs], N)
            else:
                ctx = flint.fmpz_mod_poly_ctx(N)
                self._mk = lambda cs: ctx([int(c) % N for c in cs])
        elif backend != "py":
            raise ValueError("backend must be 'flint' or 'py'")

    # ------------------------------------------------------------------ basics
    def make(self, coeffs: Sequence[int]):
        if self.backend == "py":
            out = [int(c) % self.N for c in coeffs]
            while out and out[-1] == 0:
                out.pop()
            return out
        return self._mk(list(coeffs))

    def coeffs(self, a, length: Optional[int] = None) -> List[int]:
        cs = list(a) if self.backend == "py" else [int(c) for c in a.coeffs()]
        if length is not None:
            cs = cs[:length] + [0] * (length - len(cs))
        return cs

    def length(self, a) -> int:
        return len(a) if self.backend == "py" else a.length()

    def _note(self, a):
        L = self.length(a)
        self.m.hi("stored_coeffs", L)
        self.m.hi("max_degree", L - 1)
        return a

    def add(self, a, b):
        if self.backend == "py":
            N = self.N
            if len(a) < len(b):
                a, b = b, a
            out = list(a)
            for i, y in enumerate(b):
                out[i] = (out[i] + y) % N
            self.m.inc("rn_add", len(b))
            while out and out[-1] == 0:
                out.pop()
            return out
        return a + b

    def sub(self, a, b):
        if self.backend == "py":
            return self.add(a, [(-y) % self.N for y in b])
        return a - b

    def scale(self, c: int, a):
        c %= self.N
        self.m.inc("rn_mul", self.length(a))
        if self.backend == "py":
            return self.make([c * x for x in a])
        return a * c

    def truncate(self, a, n: int):
        if self.backend == "py":
            out = a[:n]
            while out and out[-1] == 0:
                out.pop()
            return out
        return a.truncate(n) if hasattr(a, "truncate") else self.make(self.coeffs(a)[:n])

    def reverse(self, a, length: int):
        """x^(length-1) a(1/x) for a of length <= length."""
        cs = self.coeffs(a, length)
        return self.make(cs[::-1])

    # -------------------------------------------------------------- products
    def _count_mul(self, la: int, lb: int) -> None:
        n = la + lb
        self.m.inc("pmul")
        self.m.inc("pmul_len", n)
        self.m.inc("pmul_model", n * _clog2(n))

    def mul(self, a, b):
        la, lb = self.length(a), self.length(b)
        if la == 0 or lb == 0:
            return self.make([])
        self._count_mul(la, lb)
        if self.backend == "py":
            return self._note(self.make(self._kara(list(a), list(b))))
        return self._note(a * b)

    def mullow(self, a, b, n: int):
        la, lb = self.length(a), self.length(b)
        if la == 0 or lb == 0 or n <= 0:
            return self.make([])
        la, lb = min(la, n), min(lb, n)
        self._count_mul(la, lb)
        if self.backend == "py":
            return self._note(self.make(self._kara(list(a[:n]), list(b[:n]))[:n]))
        return self._note(a.mul_low(b, n))

    # exact-count Karatsuba over Z/NZ (py backend)
    def _school(self, a: List[int], b: List[int]) -> List[int]:
        N = self.N
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    out[i + j] += x * y
        self.m.inc("rn_mul", len(a) * len(b))
        self.m.inc("rn_add", max(0, len(a) * len(b) - len(out)))
        return [c % N for c in out]

    def _kara(self, a: List[int], b: List[int]) -> List[int]:
        if not a or not b:
            return []
        if min(len(a), len(b)) <= _KARA_CUTOFF:
            return self._school(a, b)
        N = self.N
        n = max(len(a), len(b))
        h = n // 2
        a0, a1 = a[:h], a[h:]
        b0, b1 = b[:h], b[h:]
        z0 = self._kara(a0, b0)
        z2 = self._kara(a1, b1) if a1 and b1 else []
        sa = [((a0[i] if i < len(a0) else 0) + (a1[i] if i < len(a1) else 0)) % N
              for i in range(max(len(a0), len(a1)))]
        sb = [((b0[i] if i < len(b0) else 0) + (b1[i] if i < len(b1) else 0)) % N
              for i in range(max(len(b0), len(b1)))]
        self.m.inc("rn_add", len(sa) + len(sb))
        z1 = self._kara(sa, sb)
        for i, c in enumerate(z0):
            z1[i] -= c
        for i, c in enumerate(z2):
            z1[i] -= c
        out = [0] * (len(a) + len(b) - 1)
        for i, c in enumerate(z0):
            out[i] += c
        for i, c in enumerate(z1):
            if i + h < len(out):
                out[i + h] += c
        for i, c in enumerate(z2):
            out[i + 2 * h] += c
        self.m.inc("rn_add", len(z0) + 2 * len(z1) + len(z2))
        return [c % N for c in out]

    # ---------------------------------------------------------- monic division
    def inv_series(self, f, n: int):
        """g with f*g = 1 mod x^n; requires f(0) == 1 (no scalar inverse taken)."""
        if self.coeffs(f, 1)[0] != 1 % self.N:
            raise ValueError("inv_series needs constant term 1 (monic reversal)")
        g = self.make([1])
        prec = 1
        while prec < n:
            prec = min(2 * prec, n)
            e = self.mullow(self.truncate(f, prec), g, prec)       # f g
            e = self.sub(self.make([1]), e)                        # 1 - f g
            g = self.add(g, self.mullow(g, e, prec))               # g + g(1-fg)
        return self.truncate(g, n)

    def divrem_monic(self, a, b, b_rev_inv=None):
        """(q, r) with a = q b + r, deg r < deg b; b monic.

        ``b_rev_inv``: optional precomputed inverse series of reverse(b) with
        precision >= deg a - deg b + 1.
        """
        la, lb = self.length(a), self.length(b)
        if lb == 0:
            raise ZeroDivisionError
        if self.coeffs(b)[lb - 1] != 1 % self.N:
            raise ValueError("divisor must be monic")
        if la < lb:
            return self.make([]), a
        k = la - lb + 1
        if b_rev_inv is None:
            b_rev_inv = self.inv_series(self.reverse(b, lb), k)
        ar = self.reverse(a, la)
        qr = self.mullow(self.truncate(ar, k), b_rev_inv, k)
        q = self.reverse(qr, k)
        r = self.truncate(self.sub(a, self.mul(q, b)), lb - 1)
        return q, r

    def evaluate(self, a, x: int) -> int:
        """Horner evaluation (counted as scalar ops)."""
        cs = self.coeffs(a)
        N = self.N
        r = 0
        for c in reversed(cs):
            r = (r * x + c) % N
        self.m.inc("rn_mul", len(cs))
        self.m.inc("rn_add", len(cs))
        return r


# ---------------------------------------------------------------------------
# subproduct / remainder trees
# ---------------------------------------------------------------------------
def product_tree_linear(R: PolyRing, roots_neg: Sequence[int]):
    """Subproduct tree of the monic linear factors (x + c), c in roots_neg.

    Returns list of levels; levels[0] = leaves, levels[-1] = [root].
    """
    level = [R.make([c, 1]) for c in roots_neg]
    if not level:
        return [[R.make([1])]]
    levels = [level]
    while len(level) > 1:
        nxt = []
        for i in range(0, len(level) - 1, 2):
            nxt.append(R.mul(level[i], level[i + 1]))
        if len(level) % 2:
            nxt.append(level[-1])
        levels.append(nxt)
        level = nxt
    return levels


def multipoint_eval(R: PolyRing, f, points: Sequence[int]) -> List[int]:
    """Values f(x) at the given points via a monic remainder tree."""
    if not points:
        return []
    tree = product_tree_linear(R, [(-x) % R.N for x in points])
    # remainder tree, top-down.  Node i at level depth+1 has children 2i, 2i+1
    # at level depth (only 2i if it was an odd pass-through).
    root = tree[-1][0]
    rems = [R.divrem_monic(f, root)[1] if R.length(f) >= R.length(root) else f]
    for depth in range(len(tree) - 2, -1, -1):
        level = tree[depth]
        new = []
        for i, r in enumerate(rems):
            kids = [2 * i] if 2 * i + 1 >= len(level) else [2 * i, 2 * i + 1]
            for k in kids:
                node = level[k]
                if R.length(r) >= R.length(node):
                    new.append(R.divrem_monic(r, node)[1])
                else:
                    new.append(r)
        rems = new
    return [R.coeffs(r, 1)[0] for r in rems]
