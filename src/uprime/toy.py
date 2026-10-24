"""uprime toy: level-1 acceptance with MIXED smoothness (rational B, algebraic B2 with squarefree (B, B2]-part)
and square child boxes (notes/uprime_prereg.md, E2).  DIAGNOSTIC / HARNESS.

N enters only through make_f (Lemma GB on N*w); N is never factored.  The gcd test q | N for children is a gcd
with N (a factor exit in the real algorithm), not a use of the factors.

Helpers (arithmetic, smooth parts, base expansion, fields) are COPIED verbatim from src/multidescent/toy.py
(itself copied from src/thmA, src/nested), per notes/AGENT_RULES.md rule 7.
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



# ------------------------------------------------------------------ uprime additions
def lattice_box(q: int, lam: int, V: int) -> List[Tuple[int, int]]:
    """Primitive points (mod +-1, normalised b > 0 or b = 0 < a) of {(a, b): a = lam b mod q} in [-V, V]^2."""
    e1, e2 = gauss_reduce((q, 0), (lam % q, 1))
    # every lattice point in the square has |coeffs| <= bound (reduced basis: |x e1 + y e2| >= c max(|x||e1|, |y||e2|))
    n1 = (e1[0] ** 2 + e1[1] ** 2) ** 0.5
    n2 = (e2[0] ** 2 + e2[1] ** 2) ** 0.5
    X = int(2 * 1.5 * V / n1) + 2
    Yb = int(2 * 1.5 * V / n2) + 2
    out = set()
    for x in range(-X, X + 1):
        for y in range(-Yb, Yb + 1):
            a = x * e1[0] + y * e2[0]
            b = x * e1[1] + y * e2[1]
            if abs(a) > V or abs(b) > V:
                continue
            if b < 0 or (b == 0 and a < 0):
                a, b = -a, -b
            if (a, b) == (0, 0) or gcd(a, b) != 1:
                continue
            out.add((a, b))
    return sorted(out)


def lambda1(q: int, lam: int) -> int:
    """Squared length of a shortest nonzero vector of {(a, b): a = lam b mod q}."""
    e1, _ = gauss_reduce((q, 0), (lam % q, 1))
    return e1[0] ** 2 + e1[1] ** 2


def reduced_box(q: int, lam: int, S: int) -> List[Tuple[int, int]]:
    """Reduced-basis box (as multidescent / nested §8), primitive, normalised."""
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
    return sorted(out)


class Ctx:
    def __init__(self, B: int, B2: int):
        self.B, self.B2 = B, B2
        self.primes_B = primes_upto(B)
        self.primes_B2 = primes_upto(B2)
        self.large = [p for p in self.primes_B2 if p > B]
        self.PB = primorial(self.primes_B)
        self.PB2 = primorial(self.primes_B2)


def classify_level1(A: Dict, ell: int, S1: int, ctx: Ctx) -> Dict:
    """Level-1 box of the rational target ell (reduced-basis box radius S1).  Returns counts and the accepted
    points with their algebraic children.
      rat:  G/ell B-smooth;  alg2: F B2-smooth;  both: rat and alg2;
      acc:  both, and the (B, B2]-part of F squarefree, coprime to b * f_d * disc(f) and to N  (Acc_ell).
      removed_sq / removed_disc / removed_N: why a 'both' point is not 'acc' (first reason in this order)."""
    g0, g1 = A["g0"], A["g1"]
    if g1 % ell == 0:
        return {"n": 0}
    lam = g0 * pow(g1, -1, ell) % ell
    pts = []
    for (a, b) in reduced_box(ell, lam, S1):
        G, F = G_val(A, a, b), F_val(A, a, b)
        if G == 0 or F == 0:
            continue
        assert G % ell == 0
        pts.append((a, b, abs(G) // ell, abs(F)))
    res = {"n": len(pts), "rat": 0, "alg2": 0, "both": 0, "acc": 0, "removed_sq": 0, "removed_disc": 0,
           "removed_N": 0, "sqdiv_any": 0, "nodes": []}
    if not pts:
        return res
    sx = smooth_parts([p[2] for p in pts], ctx.PB)
    sy = smooth_parts([p[3] for p in pts], ctx.PB2)
    fdd = A["f"][A["d"]] * A["disc"]
    N = A["N"]
    for i, (a, b, gx, fy) in enumerate(pts):
        r = sx[i] == gx
        al = sy[i] == fy
        res["rat"] += r
        res["alg2"] += al
        # unconditional: does p^2 | F for some p in (B, B2]  (exact, over all box points)
        if al and any(fy % (p * p) == 0 for p in ctx.large if p * p <= fy):
            res["sqdiv_any"] += 1                 # B2-smooth with p^2 | F, p in (B, B2]  (rational side ignored)
        if not (r and al):
            continue
        res["both"] += 1
        fac = factor_smooth(fy, ctx.primes_B2)
        bigp = {p: e for p, e in fac.items() if p > ctx.B}
        if any(e > 1 for e in bigp.values()):
            res["removed_sq"] += 1
            continue
        if any(b % p == 0 or fdd % p == 0 for p in bigp):
            res["removed_disc"] += 1
            continue
        if any(N % p == 0 for p in bigp):
            res["removed_N"] += 1
            continue
        res["acc"] += 1
        res["nodes"].append({"a": a, "b": b, "children": sorted((p, a * pow(b, -1, p) % p) for p in bigp)})
    return res


def child_count(A: Dict, q: int, r: int, S2: int, ctx: Ctx) -> Dict:
    """Square child box: primitive points of {a = r b mod q} in [-V2, V2]^2, V2 = floor(sqrt(q) S2).
    Good(P'): q || F(P'), G(P') B-smooth, F(P')/q B-smooth.  Also returns lambda1^2 (irregularity check)."""
    V2 = int(q ** 0.5 * S2)
    pts = lattice_box(q, r, V2)
    xs, ys = [], []
    for (a, b) in pts:
        G, F = G_val(A, a, b), F_val(A, a, b)
        if G == 0 or F == 0:
            continue
        if F % q != 0:
            continue                              # cannot happen when (q, r) is a root of f
        if (F // q) % q == 0:
            continue                              # q^2 | F(P'): not q || F(P')
        xs.append(abs(G))
        ys.append(abs(F) // q)
    good = 0
    if xs:
        sx = smooth_parts(xs, ctx.PB)
        cand = [i for i in range(len(xs)) if sx[i] == xs[i]]
        if cand:
            sy = smooth_parts([ys[i] for i in cand], ctx.PB)
            good = sum(1 for j, i in enumerate(cand) if sy[j] == ys[i])
    return {"good": good, "npts": len(pts), "lam1sq": lambda1(q, r), "irregular": 4 * S2 * S2 * lambda1(q, r) < q}


def random_label(A: Dict, lo: float, hi: float, rng: random.Random, ctx: Ctx) -> Optional[Tuple[int, int]]:
    """A uniformly random simple degree-one prime ideal (q, r) of A with lo <= q <= hi, q not dividing f_d disc."""
    fdd = A["f"][A["d"]] * A["disc"]
    cands = [p for p in ctx.large if lo <= p <= hi and fdd % p]
    rng.shuffle(cands)
    for p in cands:
        rts = roots_mod(A["f"], p)
        if rts:
            # weight labels uniformly: accept p with probability (#roots)/d, then a uniform root
            if rng.random() < len(rts) / A["d"]:
                return p, rng.choice(rts)
    return None
