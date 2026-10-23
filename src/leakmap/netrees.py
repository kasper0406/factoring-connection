"""One-layer N-exponent trees (class G of factorbits §3.3) and their instrumentation.  E1 of leakmap.

Model (notes/leakmap_results.md §2).  A finite decision tree whose node v computes an element of a fixed free
Z-algebra A (rank 1: Z; rank 2: Z[x]/(x^2 - b x - c)) reduced mod N, and branches on whether the element is 0 in
A/NA (all coordinates; the "element test" of the CRT criterion).  The instrumented run reports a split when the gcd
of the coordinates with N is proper.  Node elements are
  PureNode:  a value given by an N-independent straight-line program (pure, Theorem PT),
  GateNode:  pre * (y^{E(N)} - w), with y, w, pre pure and E an N-dependent exponent function (one exponent gate).
Pure values are callables m -> element mod m; they never see N except as the modulus.

OPERATIONAL: `run_on_N` (input N only).
DIAGNOSTIC (reads a prime factor x of N): `diag_transcript`, `formal_transcript`, `node_components`,
`deviation_class`, `order_mod_prime`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from math import gcd
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import sympy

Elem = Tuple[int, ...]


# ----------------------------------------------------------------------------- algebras
@dataclass(frozen=True)
class Alg:
    """rank 1 (Z) if b is None, else Z[x]/(x^2 - b x - c); elements are tuples (u,) or (u, v) = u + v x."""
    b: Optional[int] = None
    c: Optional[int] = None

    @property
    def rank(self) -> int:
        return 1 if self.b is None else 2

    def one(self) -> Elem:
        return (1,) if self.rank == 1 else (1, 0)

    def red(self, e: Elem, m: int) -> Elem:
        return tuple(x % m for x in e)

    def sub(self, e: Elem, f: Elem, m: int) -> Elem:
        return tuple((x - y) % m for x, y in zip(e, f))

    def mul(self, e: Elem, f: Elem, m: int) -> Elem:
        if self.rank == 1:
            return (e[0] * f[0] % m,)
        u1, v1 = e
        u2, v2 = f
        vv = v1 * v2
        return ((u1 * u2 + self.c * vv) % m, (u1 * v2 + u2 * v1 + self.b * vv) % m)

    def pow(self, e: Elem, k: int, m: int) -> Elem:
        if k < 0:
            raise ValueError("negative exponent")
        if self.rank == 1:
            return (pow(e[0], k, m),)
        r = self.one()
        r = self.red(r, m)
        base = self.red(e, m)
        while k:
            if k & 1:
                r = self.mul(r, base, m)
            base = self.mul(base, base, m)
            k >>= 1
        return r

    def is_zero(self, e: Elem, m: int) -> bool:
        return all(x % m == 0 for x in e)

    @property
    def disc(self) -> int:
        return self.b * self.b + 4 * self.c


Z1 = Alg()


# ----------------------------------------------------------------------------- exponent functions
@dataclass(frozen=True)
class Expo:
    """An N-dependent exponent E(N) >= 0.  `affine` = (alpha, beta) when E(N) = alpha*N + beta, else None."""
    name: str
    f: Callable[[int], int] = field(compare=False, hash=False)
    affine: Optional[Tuple[int, int]] = None

    def __call__(self, N: int) -> int:
        return self.f(N)


def expo_affine(alpha: int, beta: int) -> Expo:
    return Expo(f"{alpha}N{beta:+d}", lambda N, a=alpha, b=beta: a * N + b, (alpha, beta))


def expo_floor_dyadic(j: int) -> Expo:
    return Expo(f"floor((N-1)/2^{j})", lambda N, j=j: (N - 1) >> j, None)


def expo_half_np1() -> Expo:
    return Expo("(N+1)/2", lambda N: (N + 1) >> 1, None)


# ----------------------------------------------------------------------------- pure values
@dataclass(frozen=True)
class Pure:
    """An N-independent pure value: name + callable m -> element mod m."""
    name: str
    f: Callable[[int], Elem] = field(compare=False, hash=False)

    def __call__(self, m: int) -> Elem:
        return self.f(m)


def pure_const(alg: Alg, *coords: int) -> Pure:
    e = tuple(coords) if alg.rank == len(coords) else tuple(coords) + (0,) * (alg.rank - len(coords))
    return Pure(f"const{e}", lambda m, e=e: tuple(x % m for x in e))


def pure_pow(alg: Alg, base: Elem, M: int, minus_one: bool = True) -> Pure:
    def f(m, base=base, M=M):
        r = alg.pow(base, M, m)
        if minus_one:
            r = alg.sub(r, alg.one(), m)
        return r
    return Pure(f"{base}^{M}{'-1' if minus_one else ''}", f)


# ----------------------------------------------------------------------------- nodes
@dataclass(frozen=True)
class PureNode:
    alg: Alg
    value: Pure

    def element(self, N: int, m: int, cache: dict) -> Elem:
        return self.value(m)

    def formal_zero(self, x: int) -> bool:
        return self.alg.is_zero(self.value(x), x)


@dataclass(frozen=True)
class GateNode:
    """element = pre * (y^{E(N)} - w); pre is a rank-1 pure scalar (None = 1).

    t_target: integer t with w = y^t as pure values (so the node is 'y^{E(N)-t} = 1'), or None.
    ambient: 'Fp' (y in F_p^*), 'torus' (norm-one torus of a rank-2 algebra), 'units' (all units of A_p);
    used only to decide UNIVERSAL deviations (exp(G_p) | E(N) - t)."""
    alg: Alg
    y: Pure
    E: Expo
    w: Pure
    pre: Optional[Pure] = None
    t_target: Optional[int] = None
    ambient: str = "Fp"
    label: str = ""

    def gate(self, N: int, m: int, cache: dict) -> Elem:
        key = (self.alg, self.y.name, self.E.name, m)
        v = cache.get(key)
        if v is None:
            v = self.alg.pow(self.y(m), self.E(N), m)
            cache[key] = v
        return v

    def element(self, N: int, m: int, cache: dict) -> Elem:
        z = self.alg.sub(self.gate(N, m, cache), self.w(m), m)
        if self.pre is not None:
            s = self.pre(m)[0]
            z = tuple(s * x % m for x in z)
        return z

    def formal_zero(self, x: int) -> bool:
        """Formal reference: the gate output replaced by a formal variable Z (y is a unit at x).
        pre*(Z - w) is the zero polynomial iff pre == 0 mod x."""
        yx = self.y(x)
        if self.alg.rank == 1 and yx[0] % x == 0:
            raise ValueError("gate base not a unit at x")
        if self.pre is None:
            return False
        return self.pre(x)[0] % x == 0


# ----------------------------------------------------------------------------- trees
def leaf(label) -> tuple:
    return ("leaf", label)


def node(nd, if_nonzero, if_zero) -> tuple:
    return ("node", nd, if_nonzero, if_zero)


def nonadaptive(nodes: Sequence) -> tuple:
    """Lazy non-adaptive tree: every node is evaluated in order; the leaf label is the answer tuple."""
    return ("seq", tuple(nodes))


def _step(t, a):
    """Child of tree t after answer a (1 = zero branch).  Supports lazy ("seq", nodes, i) trees."""
    if t[0] == "seq":
        nodes, i = t[1], (t[2] if len(t) > 2 else 0)
        return ("seq", nodes, i + 1)
    return t[3] if a else t[2]


def _current(t):
    """("node", nd) or ("leaf", label) view of t."""
    if t[0] == "seq":
        nodes, i = t[1], (t[2] if len(t) > 2 else 0)
        return ("node", nodes[i]) if i < len(nodes) else ("leaf", None)
    if t[0] == "node":
        return ("node", t[1])
    return ("leaf", t[1])


def run_on_N(tree: tuple, N: int):
    """OPERATIONAL.  Instrumented run mod N.  Returns ("split", divisor, depth) or ("leaf", transcript)."""
    t, tr, cache = tree, (), {}
    while True:
        cur = _current(t)
        if cur[0] == "leaf":
            return ("leaf", tr)
        nd = cur[1]
        e = nd.element(N, N, cache)
        g = N
        for c in e:
            g = gcd(g, c)
        if 1 < g < N:
            return ("split", g, len(tr))
        a = 1 if g == N else 0
        tr += (a,)
        t = _step(t, a)


def diag_transcript(tree: tuple, N: int, x: int) -> tuple:
    """DIAGNOSTIC.  Diagonal transcript of the prime x | N: answers [element == 0 mod x] with N's exponents."""
    t, tr, cache = tree, (), {}
    while True:
        cur = _current(t)
        if cur[0] == "leaf":
            return tr
        nd = cur[1]
        a = 1 if nd.alg.is_zero(nd.element(N, x, cache), x) else 0
        tr += (a,)
        t = _step(t, a)


def formal_transcript(tree: tuple, x: int) -> tuple:
    """DIAGNOSTIC.  Reference class kappa(x): p-only transcript with gate outputs formal."""
    t, tr = tree, ()
    while True:
        cur = _current(t)
        if cur[0] == "leaf":
            return tr
        a = 1 if cur[1].formal_zero(x) else 0
        tr += (a,)
        t = _step(t, a)


def path_nodes(tree: tuple, transcript: tuple) -> List:
    t, out = tree, []
    for a in transcript:
        out.append(_current(t)[1])
        t = _step(t, a)
    return out


def leaf_of(tree: tuple, transcript: tuple):
    t = tree
    for a in transcript:
        t = _step(t, a)
    cur = _current(t)
    return transcript if t[0] == "seq" else cur[1]


# ----------------------------------------------------------------------------- orders (DIAGNOSTIC)
@lru_cache(maxsize=None)
def _factor(n: int) -> Tuple[Tuple[int, int], ...]:
    return tuple(sorted(sympy.factorint(n).items()))


def order_generic(mulpow: Callable[[int], bool], group_order: int) -> int:
    """Order of an element given a predicate is_identity(element^k) and a multiple of its order."""
    o = group_order
    for l, e in _factor(group_order):
        for _ in range(e):
            if o % l == 0 and mulpow(o // l):
                o //= l
            else:
                break
    return o


def order_mod_prime(a: int, x: int) -> int:
    a %= x
    if a == 0:
        raise ValueError("zero has no order")
    return order_generic(lambda k: pow(a, k, x) == 1, x - 1)


def legendre(a: int, x: int) -> int:
    a %= x
    if a == 0:
        return 0
    return 1 if pow(a, (x - 1) // 2, x) == 1 else -1


def sqrt_mod_prime(a: int, x: int) -> int:
    r = sympy.sqrt_mod(a % x, x)
    if r is None:
        raise ValueError("non-residue")
    return int(r)


@dataclass
class Component:
    order: int          # order of the gate base in this component
    ambient_exp: int    # exponent of the ambient group (for UNIVERSAL)


def node_components(nd: GateNode, x: int) -> List[Component]:
    """DIAGNOSTIC.  Split the gate base of nd into the components of A/xA; return orders and ambient exponents."""
    alg = nd.alg
    y = nd.y(x)
    if alg.rank == 1:
        return [Component(order_mod_prime(y[0], x), x - 1)]
    d = alg.disc % x
    s = legendre(d, x)
    if s == 0:
        raise ValueError("ramified prime")
    if s == 1:
        sq = sqrt_mod_prime(d, x)
        inv2 = pow(2, -1, x)
        roots = [((alg.b + sq) * inv2) % x, ((alg.b - sq) * inv2) % x]
        comps = []
        for r in roots:
            val = (y[0] + y[1] * r) % x
            o = order_mod_prime(val, x)
            amb = x - 1
            comps.append(Component(o, amb))
        return comps
    # inert: F_{x^2}
    G = x * x - 1
    o = order_generic(lambda k: alg.pow(y, k, x) == alg.red(alg.one(), x), G)
    amb = (x + 1) if nd.ambient == "torus" else G
    return [Component(o, amb)]


def deviation_class(nd: GateNode, N: int, x: int, L: int) -> str:
    """DIAGNOSTIC.  For a node that deviates at x (gate hit its target although the formal test is nonzero):
    SMALL if some component order <= L, UNIVERSAL if exp(ambient) | E(N) - t (every base fooled),
    else SUBGROUP (cross gcd coincidence)."""
    comps = node_components(nd, x)
    if min(c.order for c in comps) <= L:
        return "SMALL"
    if nd.t_target is not None:
        En = nd.E(N)
        if all((En - nd.t_target) % c.ambient_exp == 0 for c in comps):
            return "UNIVERSAL"
    return "SUBGROUP"


# ----------------------------------------------------------------------------- the E1 tree menu
def _lcm_upto(B: int) -> int:
    M = 1
    for k in range(2, B + 1):
        M = M * k // gcd(M, k)
    return M


M64 = _lcm_upto(64)


def tree_fermat() -> tuple:
    n2 = GateNode(Z1, pure_const(Z1, 2), expo_affine(1, -1), pure_const(Z1, 1), t_target=0, label="2^(N-1)=1")
    n3 = GateNode(Z1, pure_const(Z1, 3), expo_affine(1, -1), pure_const(Z1, 1), t_target=0, label="3^(N-1)=1")
    return node(n2, leaf("composite"), node(n3, leaf("composite"), leaf("probable prime")))


def tree_mr2(jmax: int = 3) -> tuple:
    nodes = []
    for j in range(jmax, -1, -1):
        E = expo_floor_dyadic(j)
        nodes.append(GateNode(Z1, pure_const(Z1, 2), E, pure_const(Z1, 1), t_target=0, label=f"b{j}=1"))
        nodes.append(GateNode(Z1, pure_const(Z1, 2), E, pure_const(Z1, -1), t_target=None, label=f"b{j}=-1"))
    return nonadaptive(nodes)


LUCAS_ALG = Alg(3, -1)      # x^2 = 3x - 1: alpha = x has norm 1, trace 3; disc 5
QF_ALG = Alg(1, 3)          # x^2 = x + 3: disc 13; sigma(x) = 1 - x


def tree_lucas() -> tuple:
    a = pure_const(LUCAS_ALG, 0, 1)
    one = pure_const(LUCAS_ALG, 1, 0)
    mone = pure_const(LUCAS_ALG, -1, 0)
    n1 = GateNode(LUCAS_ALG, a, expo_affine(1, 1), one, t_target=0, ambient="torus", label="a^(N+1)=1")
    n2 = GateNode(LUCAS_ALG, a, expo_half_np1(), one, t_target=0, ambient="torus", label="a^((N+1)/2)=1")
    n3 = GateNode(LUCAS_ALG, a, expo_half_np1(), mone, t_target=None, ambient="torus", label="a^((N+1)/2)=-1")
    return nonadaptive([n1, n2, n3])


def tree_qfrob() -> tuple:
    a = pure_const(QF_ALG, 0, 1)
    sig = pure_const(QF_ALG, 1, -1)
    n1 = GateNode(QF_ALG, a, expo_affine(1, 0), sig, t_target=None, ambient="units", label="a^N=sigma(a)")
    return nonadaptive([n1])


def tree_arch(K: int = 32) -> tuple:
    nodes = [GateNode(Z1, pure_const(Z1, 2), expo_affine(1, 0), Pure(f"2^{k}", lambda m, k=k: (pow(2, k, m),)),
                      t_target=k, label=f"2^N=2^{k}") for k in range(1, K + 1)]
    return nonadaptive(nodes)


def tree_repair() -> tuple:
    pre = Pure("3^M64-1", lambda m: ((pow(3, M64, m) - 1) % m,))
    n1 = GateNode(Z1, pure_const(Z1, 2), expo_affine(1, -1), pure_const(Z1, 1), pre=pre, t_target=0,
                  label="(2^(N-1)-1)(3^M-1)")
    return nonadaptive([n1])


def tree_pm1() -> tuple:
    return nonadaptive([PureNode(Z1, pure_pow(Z1, (2,), M64))])


TREES: Dict[str, Callable[[], tuple]] = {
    "FERMAT": tree_fermat, "MR2": tree_mr2, "LUCAS": tree_lucas, "QFROB": tree_qfrob,
    "ARCH": tree_arch, "REPAIR": tree_repair, "PM1": tree_pm1,
}
