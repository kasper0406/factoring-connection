"""multidescent toy model: level-1 -> level-2 children in Theorem-U-type special-q descent.

DIAGNOSTIC / HARNESS code (notes/multidescent_prereg.md).  N enters only through make_f (Lemma GB on N*w,
the shared root g0/g1 mod N); N is never factored and its factors are never used.  Values of the forms G, F at
lattice points are tested for smoothness exactly (Bernstein smooth parts over a primorial, then trial division of
values already known to be smooth).

Helpers copied (not imported) from src/thmA/{arith,numfield,smooth}.py and src/nested/toy_e1.py, per
notes/AGENT_RULES.md rule 7 (stay inside your own files).

Conventions (as thmA / nested):
  G(x, y) = g1 x - g0 y;  F = F_hat_{g, N w} + G * R_c, F(g0, g1) = N w;  f(x) = F(x, 1), f[i] coeff of x^i.
  rational target  ('r', q, lam): lattice {(a, b): a = lam b mod q}, lam = g0/g1 mod q   (q | G(a, b))
  algebraic target ('a', q, r)  : lattice {(a, b): a = r b mod q}, r a simple root of f mod q (q | F(a, b))
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint
import gmpy2


# ------------------------------------------------------------------ arithmetic (copied from src/thmA/arith.py)
def is_prime(n: int) -> bool:
    return n >= 2 and bool(gmpy2.is_prime(n, 30))


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def gauss_reduce(u: Tuple[int, int], v: Tuple[int, int]) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    def n2(x):
        return x[0] * x[0] + x[1] * x[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        nu = n2(u)
        if nu == 0:
            raise ValueError("degenerate basis")
        dot = u[0] * v[0] + u[1] * v[1]
        m = (2 * dot + nu) // (2 * nu)
        v = (v[0] - m * u[0], v[1] - m * u[1])
        if n2(v) >= nu:
            return u, v
        u, v = v, u


# ------------------------------------------------------------------ smoothness (copied from src/thmA/smooth.py)
def primorial(primes: Sequence[int]) -> int:
    level = list(primes) or [1]
    while len(level) > 1:
        level = [level[i] * level[i + 1] if i + 1 < len(level) else level[i] for i in range(0, len(level), 2)]
    return level[0]


def smooth_parts(xs: Sequence[int], P: int) -> List[int]:
    """s_i = the P-smooth part of |x_i| (x_i != 0); x_i is smooth iff s_i == |x_i|."""
    xs = [abs(x) for x in xs]
    if not xs:
        return []
    assert all(x > 0 for x in xs)
    tree = [xs]
    while len(tree[-1]) > 1:
        cur = tree[-1]
        tree.append([cur[i] * cur[i + 1] if i + 1 < len(cur) else cur[i] for i in range(0, len(cur), 2)])
    rems = [P % tree[-1][0]]
    for lvl in range(len(tree) - 2, -1, -1):
        cur = tree[lvl]
        rems = [rems[i // 2] % cur[i] for i in range(len(cur))]
    e = max(1, max(x.bit_length() for x in xs).bit_length())
    out = []
    for x, r in zip(xs, rems):
        y = r
        for _ in range(e):
            y = y * y % x
        out.append(gcd(x, y) if x > 1 else 1)
    return out


def factor_smooth(x: int, primes: Sequence[int]) -> Dict[int, int]:
    x = abs(x)
    fac: Dict[int, int] = {}
    for p in primes:
        if x == 1:
            break
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            fac[p] = e
    assert x == 1, "factor_smooth called on a non-smooth value"
    return fac


# ------------------------------------------------------------------ forms (copied from src/thmA/numfield.py)
def _centered(x: int, m: int) -> int:
    r = x % m
    return r - m if 2 * r > m else r


def gbe(T: int, d: int, g0: int, g1: int) -> List[int]:
    assert gcd(g0, g1) == 1 and g1 >= 2
    a = [0] * (d + 1)
    for i in range(d, 0, -1):
        a[i] = _centered(T * pow(g0, -i, g1), g1)
        T -= a[i] * g0 ** i
        assert T % g1 == 0
        T //= g1
    a[0] = T
    return a


def hom_eval(f: Sequence[int], a: int, b: int) -> int:
    d = len(f) - 1
    return sum(f[i] * a ** i * b ** (d - i) for i in range(d + 1))


def make_f(N: int, w: int, d: int, g0: int, g1: int, c: Sequence[int]) -> List[int]:
    f = gbe(N * w, d, g0, g1)
    for i, ci in enumerate(c):
        j = d - 1 - i
        f[j + 1] += g1 * ci
        f[j] -= g0 * ci
    assert hom_eval(f, g0, g1) == N * w
    return f


def disc(f: Sequence[int]) -> int:
    d = len(f) - 1
    P = flint.fmpz_poly(list(f))
    res = int(P.resultant(P.derivative()))
    sign = -1 if (d * (d - 1) // 2) % 2 else 1
    assert res % f[d] == 0
    return sign * res // f[d]


def roots_mod(f: Sequence[int], r: int) -> List[int]:
    P = flint.nmod_poly([x % r for x in f], r)
    if P.degree() <= 0:
        return []
    return sorted(int(z) for z, m in P.roots() if m == 1)


def irreducible_mod(f: Sequence[int], r: int) -> bool:
    d = len(f) - 1
    if f[d] % r == 0:
        return False
    P = flint.nmod_poly([x % r for x in f], r)
    _, facs = P.factor()
    return len(facs) == 1 and facs[0][1] == 1 and facs[0][0].degree() == d


def poly_mul(x: Sequence[int], y: Sequence[int]) -> List[int]:
    out = [0] * (len(x) + len(y) - 1)
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                out[i + j] += a * b
    return out


# ------------------------------------------------------------------ fields
def _irreducible_certified(f: Sequence[int]) -> bool:
    return any(irreducible_mod(f, r) for r in primes_upto(400) if f[-1] % r)


def finish_field(N: int, g0: int, g1: int, w: int, f: List[int], extra: Dict | None = None) -> Optional[Dict]:
    d = len(f) - 1
    e = 0
    for x in f:
        e = gcd(e, x)
    if e != 1 or f[d] == 0 or not _irreducible_certified(f):
        return None
    D = disc(f)
    if D == 0:
        return None
    A = {"g0": g0, "g1": g1, "w": w, "f": list(f), "d": d, "disc": D, "N": N}
    if extra:
        A.update(extra)
    return A


def sample_field(N: int, rng: random.Random, d: int = 3, Hc: int = 4, M_bits: int | None = None) -> Dict:
    """Design (a''') toy field: g uniform in [M/2, M]^2 coprime, w a prime in (32, 64], c uniform in [-Hc, Hc]^d."""
    n = N.bit_length()
    if M_bits is None:
        M_bits = -(-(n + 6) // (d + 1))
    M = 1 << M_bits
    while True:
        while True:
            g0, g1 = rng.randrange(M // 2, M + 1), rng.randrange(M // 2, M + 1)
            if gcd(g0, g1) == 1:
                break
        while True:
            w = rng.randrange(33, 65)
            if is_prime(w):
                break
        c = [rng.randrange(-Hc, Hc + 1) for _ in range(d)]
        A = finish_field(N, g0, g1, w, make_f(N, w, d, g0, g1, c), {"c": c, "M": M})
        if A is not None:
            return A


def fiber_field(A: Dict, P0: Tuple[int, int], rng: random.Random, Hh: int) -> Optional[Dict]:
    """Constrained resampling (codex nested_2): F_H = F + G * L0 * H, L0 = b0 x - a0 y, H a uniform binary form of
    degree d - 2 with coefficients in [-Hh, Hh].  Keeps F(P0), F(g0, g1) = N w, hence the common root and every child
    label (q, a0/b0 mod q) of P0.  Returns None if the resampled form fails the field checks."""
    a0, b0 = P0
    d = A["d"]
    h = [rng.randrange(-Hh, Hh + 1) for _ in range(d - 1)]       # H(x, 1) = sum h_i x^i, degree <= d - 2
    Gx = [-A["g0"], A["g1"]]                                     # G(x, 1) = g1 x - g0
    L0 = [-a0, b0]                                               # L0(x, 1) = b0 x - a0
    add = poly_mul(poly_mul(Gx, L0), h)
    f = list(A["f"])
    for i, v in enumerate(add):
        f[i] += v
    B = finish_field(A["N"], A["g0"], A["g1"], A["w"], f, {"h": h, "parent_P0": P0, "M": A.get("M")})
    if B is not None:
        assert hom_eval(B["f"], a0, b0) == hom_eval(A["f"], a0, b0)
        assert hom_eval(B["f"], A["g0"], A["g1"]) == A["N"] * A["w"]
    return B


def G_val(A: Dict, a: int, b: int) -> int:
    return A["g1"] * a - A["g0"] * b


def F_val(A: Dict, a: int, b: int) -> int:
    return hom_eval(A["f"], a, b)


# ------------------------------------------------------------------ targets and boxes
def rational_target(A: Dict, q: int) -> Optional[Tuple]:
    if A["g1"] % q == 0:
        return None
    return ("r", q, A["g0"] * pow(A["g1"], -1, q) % q)


def algebraic_targets(A: Dict, q: int) -> List[Tuple]:
    f = A["f"]
    if (f[A["d"]] * A["disc"]) % q == 0:
        return []
    return [("a", q, r) for r in roots_mod(f, q)]


def box_points(q: int, lam: int, S: int) -> List[Tuple[int, int]]:
    """Reduced-basis box (nested §8): P = x e1 + y e2, |x|, |y| <= S, primitive, normalised to b > 0 (or b = 0, a > 0).
    e1, e2 = Gauss-reduced basis of {(a, b) : a = lam b mod q}."""
    e1, e2 = gauss_reduce((q, 0), (lam % q, 1))
    out = set()
    for x in range(-S, S + 1):
        for y in range(-S, S + 1):
            a = x * e1[0] + y * e2[0]
            b = x * e1[1] + y * e2[1]
            if b < 0 or (b == 0 and a < 0):
                a, b = -a, -b
            if (a, b) == (0, 0) or gcd(a, b) != 1:
                continue
            out.add((a, b))
    return sorted(out, key=lambda P: (P[0] * P[0] + P[1] * P[1], P))


class Ctx:
    def __init__(self, B: int, B2: int):
        self.B, self.B2 = B, B2
        self.primes_B = primes_upto(B)
        self.primes_B2 = primes_upto(B2)
        self.PB = primorial(self.primes_B)
        self.PB2 = primorial(self.primes_B2)


def good_points(A: Dict, tgt: Tuple, S: int, ctx: Ctx) -> int:
    """Number of box points of the target whose two cofactors (target side divided by q once) are B-smooth."""
    kind, q, lam = tgt
    xs, ys = [], []
    for (a, b) in box_points(q, lam, S):
        G, F = G_val(A, a, b), F_val(A, a, b)
        if G == 0 or F == 0:
            continue
        if kind == "r":
            assert G % q == 0
            xs.append(abs(G) // q)
            ys.append(abs(F))
        else:
            assert F % q == 0
            xs.append(abs(G))
            ys.append(abs(F) // q)
    if not xs:
        return 0
    sx = smooth_parts(xs, ctx.PB)
    cand = [i for i in range(len(xs)) if sx[i] == xs[i]]
    if not cand:
        return 0
    sy = smooth_parts([ys[i] for i in cand], ctx.PB)
    return sum(1 for j, i in enumerate(cand) if sy[j] == ys[i])


def level1_nodes(A: Dict, ell: int, S1: int, ctx: Ctx, max_nodes: int = 1) -> List[Dict]:
    """Box points of the rational level-1 target ell (in order) whose cofactors G/ell and F are B2-smooth, with their
    children (rational primes and simple degree-one prime ideals of norm in (B, B2]); at most max_nodes of them."""
    tgt = rational_target(A, ell)
    if tgt is None:
        return []
    vals = []
    for (a, b) in box_points(ell, tgt[2], S1):
        G, F = G_val(A, a, b), F_val(A, a, b)
        if G == 0 or F == 0:
            continue
        vals.append((a, b, abs(G) // ell, abs(F)))
    if not vals:
        return []
    sx = smooth_parts([v[2] for v in vals], ctx.PB2)
    cand = [i for i, v in enumerate(vals) if sx[i] == v[2]]
    if not cand:
        return []
    sy = smooth_parts([vals[i][3] for i in cand], ctx.PB2)
    fdd = A["f"][A["d"]] * A["disc"]
    out = []
    for j, i in enumerate(cand):
        a, b, gx, fy = vals[i]
        if sy[j] != fy:
            continue
        fg = factor_smooth(gx, ctx.primes_B2)
        ff = factor_smooth(fy, ctx.primes_B2)
        rat, alg, degen = [], [], 0
        for p in fg:
            if p <= ctx.B:
                continue
            if A["g1"] % p == 0:
                degen += 1
            else:
                rat.append(p)
        for p in ff:
            if p <= ctx.B:
                continue
            if b % p == 0 or fdd % p == 0 or ff[p] > 1:
                degen += 1
                continue
            alg.append((p, a * pow(b, -1, p) % p))
        out.append({"ell": ell, "a": a, "b": b, "rat": sorted(rat), "alg": sorted(alg), "degenerate": degen,
                    "rat_exp": {p: fg[p] for p in rat}})
        if len(out) >= max_nodes:
            break
    return out
