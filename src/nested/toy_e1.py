"""E1 (notes/nested_prereg.md): toy endogeneity test for multi-level special-q descent in ONE field.

DIAGNOSTIC code.  It never factors N and never uses a factor of N: N enters only through make_f (Lemma GB on N*w).
Smoothness is decided exactly (batch smooth parts over a primorial, then trial division of known-smooth values).

Targets:
  rational  ('r', q, lam) : lattice {(a, b) : a = lam*b mod q}, lam = g0/g1 mod q   (q | G(a, b))
  algebraic ('a', q, r)   : lattice {(a, b) : a = r*b mod q}, r a simple root of f mod q (q | F(a, b))
D1(target): the box {|a| <= V, 0 < b <= V, gcd(a, b) = 1, G != 0, F != 0} of the lattice, V = isqrt(q*S*S),
contains a point whose target-side value divided by q once and whose other value are both B-smooth.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

from src.thmA import arith
from src.thmA import numfield as nf
from src.thmA.smooth import primorial, smooth_parts


# ------------------------------------------------------------------ field
def sample_field(N: int, rng: random.Random, d: int = 3, H: int = 4, wlo: int = 2 ** 5, whi: int = 2 ** 6
                 ) -> Dict:
    n = N.bit_length()
    M = 1 << ((n + 3) // 4)
    while True:
        while True:
            g0, g1 = rng.randrange(M // 2, M + 1), rng.randrange(M // 2, M + 1)
            if gcd(g0, g1) == 1:
                break
        while True:
            w = rng.randrange(wlo + 1, whi + 1)
            if arith.is_prime(w):
                break
        c = [rng.randrange(-H, H + 1) for _ in range(d)]
        f = nf.make_f(N, w, d, g0, g1, c)
        e = 0
        for x in f:
            e = gcd(e, x)
        if e != 1 or f[d] == 0:
            continue
        if not any(nf.irreducible_mod(f, r) for r in arith.primes_upto(200) if f[d] % r):
            continue
        D = nf.disc(f)
        if D == 0:
            continue
        return {"g0": g0, "g1": g1, "w": w, "c": c, "f": f, "d": d, "disc": D, "M": M}


def G_val(A: Dict, a: int, b: int) -> int:
    return A["g1"] * a - A["g0"] * b


def F_val(A: Dict, a: int, b: int) -> int:
    return nf.hom_eval(A["f"], a, b)


# ------------------------------------------------------------------ targets
def rational_target(A: Dict, q: int) -> Optional[Tuple]:
    if A["g1"] % q == 0:
        return None
    return ("r", q, A["g0"] * pow(A["g1"], -1, q) % q)


def algebraic_targets(A: Dict, q: int) -> List[Tuple]:
    f = A["f"]
    if (f[A["d"]] * A["disc"]) % q == 0:
        return []
    return [("a", q, r) for r in nf.roots_mod(f, q)]


def box_points(q: int, lam: int, V: int) -> List[Tuple[int, int]]:
    """All (a, b) with a = lam*b mod q, |a| <= V, 0 < b <= V, gcd(a, b) = 1."""
    e1, e2 = arith.gauss_reduce((q, 0), (lam, 1))
    det = abs(e1[0] * e2[1] - e1[1] * e2[0])
    assert det == q
    # coordinates (x, y) of (a, b) = x e1 + y e2:  |x| <= (|e2_b| V + |e2_a| V)/det, same for y with e1
    X = (abs(e2[0]) + abs(e2[1])) * V // det + 1
    Yc = (abs(e1[0]) + abs(e1[1])) * V // det + 1
    out = []
    for x in range(-X, X + 1):
        for y in range(-Yc, Yc + 1):
            a = x * e1[0] + y * e2[0]
            b = x * e1[1] + y * e2[1]
            if 0 < b <= V and abs(a) <= V and gcd(a, b) == 1:
                out.append((a, b))
    out.sort(key=lambda P: (P[1], P[0]))
    return out


class Ctx:
    def __init__(self, B: int, Y: int):
        self.B, self.Y = B, Y
        self.primes_B = arith.primes_upto(B)
        self.primes_Y = arith.primes_upto(Y)
        self.PB = primorial(self.primes_B)
        self.PY = primorial(self.primes_Y)


def d1(A: Dict, tgt: Tuple, S: int, ctx: Ctx, cache: Dict) -> bool:
    key = (tgt, S)
    if key in cache:
        return cache[key]
    kind, q, lam = tgt
    V = isqrt(q * S * S)
    pts = box_points(q, lam, V)
    xs, ys = [], []
    for (a, b) in pts:
        G, F = G_val(A, a, b), F_val(A, a, b)
        if G == 0 or F == 0:
            continue
        if kind == "r":
            assert G % q == 0
            xs.append(abs(G) // q); ys.append(abs(F))
        else:
            assert F % q == 0
            xs.append(abs(G)); ys.append(abs(F) // q)
    ok = False
    if xs:
        sx = smooth_parts(xs, ctx.PB)
        cand = [i for i in range(len(xs)) if sx[i] == xs[i]]
        if cand:
            sy = smooth_parts([ys[i] for i in cand], ctx.PB)
            ok = any(sy[j] == ys[i] for j, i in enumerate(cand))
    cache[key] = ok
    return ok


def _factor_known_smooth(x: int, primes: Sequence[int]) -> Dict[int, int]:
    x = abs(x)
    fac = {}
    for p in primes:
        if x == 1:
            break
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            fac[p] = e
    assert x == 1
    return fac


def parent_node(A: Dict, ell: int, S1: int, ctx: Ctx) -> Optional[Dict]:
    """First box point of the rational target ell whose cofactors are Y-smooth and which has >= 1 child."""
    tgt = rational_target(A, ell)
    if tgt is None:
        return None
    V = isqrt(ell * S1 * S1)
    pts = box_points(ell, tgt[2], V)
    vals = []
    for (a, b) in pts:
        G, F = G_val(A, a, b), F_val(A, a, b)
        if G == 0 or F == 0:
            continue
        vals.append((a, b, abs(G) // ell, abs(F)))
    if not vals:
        return None
    sx = smooth_parts([v[2] for v in vals], ctx.PY)
    sy = smooth_parts([v[3] for v in vals], ctx.PY)
    fdd = A["f"][A["d"]] * A["disc"]
    for i, (a, b, gx, fy) in enumerate(vals):
        if sx[i] != gx or sy[i] != fy:
            continue
        fg = _factor_known_smooth(gx, ctx.primes_Y)
        ff = _factor_known_smooth(fy, ctx.primes_Y)
        rat = [p for p in fg if p > ctx.B]
        alg, degen = [], 0
        for p in ff:
            if p <= ctx.B:
                continue
            if b % p == 0 or fdd % p == 0:
                degen += 1
                continue
            alg.append(p)
        if not rat and not alg:
            continue
        children = []
        for p in rat:
            t = rational_target(A, p)
            if t is None:
                degen += 1
            else:
                children.append(t)
        for p in alg:
            r = a * pow(b, -1, p) % p
            children.append(("a", p, r))
        if not children:
            continue
        return {"ell": ell, "a": a, "b": b, "children": children, "degenerate": degen}
    return None


class Uniform:
    """Uniform targets of a type in a norm window (B, Y], for one field."""

    def __init__(self, A: Dict, ctx: Ctx):
        self.A = A
        self.rat = [p for p in ctx.primes_Y if p > ctx.B and A["g1"] % p]
        self.alg = []
        for p in ctx.primes_Y:
            if p > ctx.B:
                self.alg.extend(algebraic_targets(A, p))
        self.alg.sort(key=lambda t: t[1])

    def draw(self, kind: str, p: int, rng: random.Random) -> Tuple:
        lo, hi = p / 1.25, p * 1.25
        if kind == "r":
            pool = [q for q in self.rat if lo <= q <= hi]
            return rational_target(self.A, rng.choice(pool))
        pool = [t for t in self.alg if lo <= t[1] <= hi]
        return rng.choice(pool)
