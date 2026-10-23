"""PRIM programs: exponent gates x -> x^{E(N)} that may be retained, nested and summed.  E1 of `prim`.

A node is an expression over Z (rank-1 algebra) built from integer constants, ring operations, pure powers
x^M (M a fixed integer) and exponent gates x^{aN+b} (affine menu) or x^{floor((N-c)/2^j)} (floor menu).  Gates may
take any expression as base (nesting); gate outputs may be reused (retained) and added (summed).

Three evaluations of the same expression:
  * `ev(N, m)`      : the value mod m on input N (m = N is OPERATIONAL; m = p, q prime is the diagonal value).
  * `vec(p, r)`     : DIAGNOSTIC.  Vectorised value mod a prime p when every affine gate exponent is a*r + b with
                      r in Z/(p-1) (the "uniform-r" evaluation; at the true partner q, r = N mod (p-1) = q mod (p-1)).
  * `NestModel`     : the fresh-random model for nodes A * u^{E} - B: u replaced by a uniform element of F_p^*.

`LucasNode` is the one rank-2 node (A = Z[x]/(x^2 - 3x + 1)); it supports `ev` only.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ----------------------------------------------------------------------------- vectorised modular arithmetic
def vpow(base: np.ndarray, e: np.ndarray, p: int) -> np.ndarray:
    """base^e mod p elementwise, base in [0, p), e >= 0 (int64); p < 2^31.  0^e is reported as 0 (gate exponents
    are >= 1 on the support, so a reduced exponent 0 with base 0 still means 0)."""
    base = np.asarray(base, dtype=np.int64) % p
    e = np.asarray(e, dtype=np.int64)
    base, e = np.broadcast_arrays(base, e)
    res = np.ones(base.shape, dtype=np.int64)
    b = base.copy()
    emax = int(e.max()) if e.size else 0
    bit = 0
    while (emax >> bit) > 0:
        mask = ((e >> bit) & 1).astype(bool)
        res = np.where(mask, (res * b) % p, res)
        b = (b * b) % p
        bit += 1
    return np.where(base == 0, 0, res)


# ----------------------------------------------------------------------------- expression DSL
class Ex:
    def ev(self, N: int, m: int) -> int:
        raise NotImplementedError

    def vec(self, p: int, r: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def bases(self) -> set:
        """pure integer constants used as gate/power bases (for the SMALL flag)."""
        return set()

    def affine(self) -> bool:
        return True

    # sugar
    def __add__(self, o): return Add(self, lift(o))
    def __radd__(self, o): return Add(lift(o), self)
    def __sub__(self, o): return Sub(self, lift(o))
    def __rsub__(self, o): return Sub(lift(o), self)
    def __mul__(self, o): return Mul(self, lift(o))
    def __rmul__(self, o): return Mul(lift(o), self)


def lift(x) -> Ex:
    return x if isinstance(x, Ex) else C(int(x))


@dataclass(frozen=True, eq=False)
class C(Ex):
    c: int

    def ev(self, N, m):
        return self.c % m

    def vec(self, p, r):
        return np.full(r.shape, self.c % p, dtype=np.int64)

    def __repr__(self):
        return str(self.c)


@dataclass(frozen=True, eq=False)
class Add(Ex):
    x: Ex
    y: Ex

    def ev(self, N, m): return (self.x.ev(N, m) + self.y.ev(N, m)) % m
    def vec(self, p, r): return (self.x.vec(p, r) + self.y.vec(p, r)) % p
    def bases(self): return self.x.bases() | self.y.bases()
    def affine(self): return self.x.affine() and self.y.affine()
    def __repr__(self): return f"({self.x} + {self.y})"


@dataclass(frozen=True, eq=False)
class Sub(Ex):
    x: Ex
    y: Ex

    def ev(self, N, m): return (self.x.ev(N, m) - self.y.ev(N, m)) % m
    def vec(self, p, r): return (self.x.vec(p, r) - self.y.vec(p, r)) % p
    def bases(self): return self.x.bases() | self.y.bases()
    def affine(self): return self.x.affine() and self.y.affine()
    def __repr__(self): return f"({self.x} - {self.y})"


@dataclass(frozen=True, eq=False)
class Mul(Ex):
    x: Ex
    y: Ex

    def ev(self, N, m): return (self.x.ev(N, m) * self.y.ev(N, m)) % m
    def vec(self, p, r): return (self.x.vec(p, r) * self.y.vec(p, r)) % p
    def bases(self): return self.x.bases() | self.y.bases()
    def affine(self): return self.x.affine() and self.y.affine()
    def __repr__(self): return f"{self.x}*{self.y}"


def _base_consts(x: Ex) -> set:
    return {x.c} if isinstance(x, C) else x.bases()


@dataclass(frozen=True, eq=False)
class PPow(Ex):
    """pure power x^M, M a fixed positive integer."""
    x: Ex
    M: int

    def ev(self, N, m): return pow(self.x.ev(N, m), self.M, m)

    def vec(self, p, r):
        e = self.M % (p - 1)
        if e == 0:
            e = p - 1
        return vpow(self.x.vec(p, r), np.full(r.shape, e, dtype=np.int64), p)

    def bases(self): return _base_consts(self.x)
    def affine(self): return self.x.affine()
    def __repr__(self): return f"{self.x}^M"


@dataclass(frozen=True, eq=False)
class G(Ex):
    """affine exponent gate x^{aN + b} (a >= 1, aN + b >= 1 on the support)."""
    x: Ex
    a: int = 1
    b: int = 0

    def E(self, N): return self.a * N + self.b

    def ev(self, N, m): return pow(self.x.ev(N, m), self.E(N), m)

    def vec(self, p, r):
        e = (self.a * r + self.b) % (p - 1)
        return vpow(self.x.vec(p, r), e, p)

    def bases(self): return _base_consts(self.x)
    def affine(self): return self.x.affine()

    def __repr__(self):
        s = "N" if self.a == 1 else f"{self.a}N"
        if self.b:
            s += f"{self.b:+d}"
        return f"{self.x}^({s})"


@dataclass(frozen=True, eq=False)
class FG(Ex):
    """floor exponent gate x^{floor((N - c) / 2^j)}."""
    x: Ex
    c: int
    j: int

    def E(self, N): return (N - self.c) >> self.j

    def ev(self, N, m): return pow(self.x.ev(N, m), self.E(N), m)
    def vec(self, p, r): raise NotImplementedError("floor gates have no uniform-r form mod p-1")
    def bases(self): return _base_consts(self.x)
    def affine(self): return False
    def __repr__(self): return f"{self.x}^((N-{self.c})/2^{self.j})"


# ----------------------------------------------------------------------------- the fresh-random model
@dataclass(frozen=True, eq=False)
class NestModel:
    """node = A * u^{E(N)} - B.  Model: u uniform in F_x^*.  Probability of a zero at prime x is
    [B/A in (F_x^*)^d] * d/(x-1), d = gcd(e, x-1), e = E(N) mod (x-1)  (A != 0; if A == 0 the node is -B)."""
    A: Ex
    u: Ex
    gate: Ex          # G or FG or PPow with base u (only its exponent is used)
    B: Ex

    def exponent_mod(self, N: Optional[int], x: int, r: Optional[np.ndarray] = None):
        g = self.gate
        if isinstance(g, PPow):
            e = g.M % (x - 1)
            return e if r is None else np.full(r.shape, e, dtype=np.int64)
        if r is None:
            return g.E(N) % (x - 1)
        return (g.a * r + g.b) % (x - 1)

    def prob_scalar(self, N: int, x: int) -> float:
        n = x - 1
        A = self.A.ev(N, x)
        B = self.B.ev(N, x)
        if A == 0:
            return float(B == 0)
        t = B * pow(A, -1, x) % x
        if t == 0:
            return 0.0
        e = self.exponent_mod(N, x)
        d = gcd(e, n)
        return d / n if pow(t, n // d, x) == 1 else 0.0

    def prob_vec(self, p: int, r: np.ndarray) -> np.ndarray:
        n = p - 1
        A = self.A.vec(p, r)
        B = self.B.vec(p, r)
        e = self.exponent_mod(None, p, r)
        d = np.gcd(e, n)
        Ainv = vpow(A, np.full(r.shape, p - 2, dtype=np.int64), p)
        t = (B * Ainv) % p
        ok = vpow(t, n // d, p) == 1
        out = np.where(ok & (t != 0), d / n, 0.0)
        return np.where(A == 0, (B == 0).astype(float), out)


# ----------------------------------------------------------------------------- rank-2 Lucas node
class LucasNode:
    """element test of (alpha^{N+1} + 2)^{N-1} - 1 in Z/m[x]/(x^2 - 3x + 1), alpha = x.  Elements (u, v) = u + v x."""
    name = "L1"

    @staticmethod
    def mul(e, f, m):
        u1, v1 = e
        u2, v2 = f
        vv = v1 * v2
        # x^2 = 3x - 1
        return ((u1 * u2 - vv) % m, (u1 * v2 + u2 * v1 + 3 * vv) % m)

    @classmethod
    def pw(cls, e, k, m):
        r = (1 % m, 0)
        while k:
            if k & 1:
                r = cls.mul(r, e, m)
            e = cls.mul(e, e, m)
            k >>= 1
        return r

    @classmethod
    def ev(cls, N, m):
        a = cls.pw((0, 1), N + 1, m)
        a = ((a[0] + 2) % m, a[1])
        z = cls.pw(a, N - 1, m)
        return ((z[0] - 1) % m, z[1])

    @staticmethod
    def lucas_V(N, m):
        """V_N = alpha^N + alpha^{-N} = trace(alpha^N) = 2u + 3v."""
        u, v = LucasNode.pw((0, 1), N, m)
        return (2 * u + 3 * v) % m


# ----------------------------------------------------------------------------- the registered node menu
TWO, THREE = C(2), C(3)
M64 = 1
for _k in range(2, 65):
    M64 = M64 * _k // gcd(M64, _k)


def _nest(A, u, gate_ab, B, floor=None):
    if floor is not None:
        gate = FG(u, *floor)
    else:
        gate = G(u, *gate_ab)
    node = Sub(Mul(lift(A), gate), lift(B)) if not (isinstance(A, int) and A == 1) else Sub(gate, lift(B))
    return node, NestModel(lift(A), u, gate, lift(B))


def nested_nodes() -> Dict[str, Tuple[Ex, NestModel]]:
    out = {}
    for w in (1, 2, 3):
        out[f"D2a_{w}"] = _nest(1, G(TWO) + 1, (1, 0), w)
    out["D2b"] = _nest(1, G(TWO, 1, -1) + 1, (1, -1), 1)
    out["D2c"] = _nest(1, G(THREE, 1, 1) + 2, (1, -1), 4)
    out["D2d"] = _nest(G(TWO), G(TWO) + 3, (1, 0), 5)
    out["D2e"] = _nest(1, G(TWO) + G(THREE), (1, 1), 1)
    inner = G(G(TWO) + 1) + 1
    out["D3"] = _nest(1, inner, (1, 0), 1)
    out["D2f"] = _nest(1, FG(THREE, 1, 1) + 2, None, 1, floor=(-1, 1))   # (3^{(N-1)/2}+2)^{(N+1)/2} - 1
    return out


def summed_nodes() -> Dict[str, Tuple[Ex, int]]:
    """name -> (node, nominal sparsity t on the group)."""
    s = {}
    s["S1"] = (G(TWO) + G(THREE) - G(C(5)), 3)
    s["S2"] = (G(TWO, 1, -1) + G(THREE, 1, -1) - 2, 3)
    s["S3"] = (G(TWO) + G(THREE, 1, 1) + G(C(5), 1, -1) + 7, 4)
    s["S4"] = (G(TWO) - 2 * G(THREE) + 3 * G(C(5)) - 4 * G(C(7)) + 5 * G(C(11)) - 6 * G(C(13)) + 7, 7)
    s["S5"] = ((G(TWO) + 1) * (G(THREE) + 1) - 2, 4)
    return s


def control_nodes() -> Dict[str, Ex]:
    c = {"C_F": G(TWO, 1, -1) - 1}
    for k in range(1, 17):
        c[f"C_A_{k}"] = G(TWO) - (2 ** k)
    c["C_P"] = PPow(G(TWO) + 1, M64) - 1
    c["C_PM1"] = PPow(TWO, M64) - 1
    return c


def c_p_model() -> NestModel:
    return NestModel(C(1), G(TWO) + 1, PPow(G(TWO) + 1, M64), C(1))


# Jacobi nodes: name -> callable N -> c mod N
def jacobi_values(b: int):
    def J1(N): return (pow(2, N, N) + 1) % N
    def J2(N):
        u = (pow(2, N, N) + 1) % N
        return (pow(u, N, N) + 3) % N
    def J3(N): return (pow(2, (N - 1) // 2, N) + 2) % N
    def J4(N): return (pow(3, N, N) + pow(2, N + 1, N) + 1) % N
    def J5(N): return (pow(2, N - 1, N) + 3) % N
    def J6(N): return (LucasNode.lucas_V(N, N) + 1) % N
    out = {"J1": J1, "J2": J2, "J3": J3, "J4": J4, "J5": J5, "J6": J6}
    for k in range(b, b + 4):
        out[f"J7_{k}"] = (lambda kk: (lambda N: (pow(2, 1 << kk, N) + 1) % N))(k)
    out["J8"] = lambda N: (pow(2, M64, N) - 3) % N
    return out
