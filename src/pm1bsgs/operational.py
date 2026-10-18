"""OPERATIONAL algorithms for the pm1bsgs study.  Input: N, public parameters, seed.

Nothing here uses p, q, phi(N), group orders or any factorisation routine.  The
only "number theory" used is a sieve for the public primes <= B (stage-1
exponent / prime-only stage 2) and gcds with N.  All polynomial arithmetic is
over Z/N via python-flint ``fmpz_mod_poly`` using ONLY multiplication and ``%``
by a monic polynomial (the audited operations, see src/common/README.md); the
derivative is formed by hand from coefficients (multiplication by small ints).

Methods (all share the stage-1 routine and the base selection):
  (a) ``pm1_bsgs``     : stage 1 (B1) + BSGS stage 2 over a window of size X=m*J,
                         hash-free collision test by Pollard-Strassen multipoint
                         evaluation  prod_j F(v_j),  F(T)=prod_i (T-u_i).
                         mode "N"     : u_i=b^(N+i), v_j=b^(j m), j=0..J-1
                                        (the user's  b^N == b^x (mod p),  x=jm-i);
                         mode "order" : u_i=b^i,     v_j=b^(j m), j=1..J
                                        (plain order finding, ord_p(b) | jm-i).
                         ``naive=True`` replaces multipoint evaluation by the
                         plain product of all m*J differences (linear cost).
  (b) ``pm1_prime``    : stage 1 (B1) + classical prime-only stage 2:
                         prod_{B1<r<=B2, r prime} (b^r - 1).
  (c) ``pm1_birthday`` : stage 1 (B1) + birthday-paradox stage 2 of Brent 1986:
                         value-independent walk x_{j+1} = x_j^2 b^{c_j}, c_j in {0,1}
                         (K points), all pairs tested at once via
                         prod_i F'(x_i) = +- prod_{i<j}(x_i-x_j)^2.
  (d) plain BSGS       : ``pm1_bsgs`` with B1 <= 1 (no stage 1), i.e. search a
                         multiple of ord_p(a) (Pollard/Strassen N^(1/4) type).

Cost accounting (Metrics keys):
  rn_mul    scalar multiplications mod N.  A modular power a^e is charged the
            binary-method count  (bitlen(e)-1) + (popcount(e)-1).
  gcd       gcd calls with N.
  poly_mul / poly_rem         number of polynomial products / remainders.
  poly_eq   ring-op equivalents of polynomial work under the explicit model
            M(d) = d*ceil(log2 d)  (d = larger operand degree) for a product and
            2*M(max(a-b+1, b)) for (deg a) mod (deg b); remainder degrees are
            charged at their nominal (generic) value so counts are data-free.
  eq_ops    = rn_mul + gcd + poly_eq   (the single budget unit used throughout).
  backtrack_* extra work after a gcd == N.
Maxima: max_degree, stored_coeffs, max_operand_bits.
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from src.common.metrics import Metrics, NullMetrics

STAGE1_CHUNK = 64      # stage-1 prime-power steps per gcd
PRIME2_CHUNK = 256     # prime-only stage 2: primes per gcd


# --------------------------------------------------------------------------
# public helpers (no factor knowledge)
# --------------------------------------------------------------------------
def primes_upto(n: int) -> List[int]:
    """Public primes <= n (sieve of Eratosthenes)."""
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= n:
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
        i += 1
    return [i for i in range(n + 1) if s[i]]


def pow_cost(e: int) -> int:
    """Binary-method multiplication count for x^e (e >= 0)."""
    if e <= 1:
        return 0
    return (e.bit_length() - 1) + (bin(e).count("1") - 1)


def mcost(d: int) -> int:
    """Model cost of a degree-d polynomial product: d*ceil(log2 d) (>=1)."""
    if d <= 1:
        return 1
    return d * (d - 1).bit_length()


def rem_cost(a: int, b: int) -> int:
    """Model cost of (degree a) mod (monic degree b), a >= b."""
    return 2 * mcost(max(a - b + 1, b))


def stage1_steps(B1: int) -> List[int]:
    """Prime-power steps of stage 1: each prime l <= B1 repeated floor(log_l B1) times."""
    out: List[int] = []
    for l in primes_upto(B1):
        pk = l
        while pk <= B1:
            out.append(l)
            pk *= l
    return out


def stage1_exponent(B1: int) -> int:
    E = 1
    for l in stage1_steps(B1):
        E *= l
    return E


def choose_base(N: int, seed: int) -> int:
    rng = random.Random(f"pm1bsgs-base|{N}|{seed}")
    return rng.randrange(2, N - 1)


def shift_target(N: int, seed: int) -> int:
    """A public pseudo-random exponent T with the bit length of N (control for
    mode "N": the window [T-(J-1)m, T+m-1] replaces [N-(J-1)m, N+m-1])."""
    rng = random.Random(f"pm1bsgs-shift|{N}|{seed}")
    k = N.bit_length()
    return rng.getrandbits(k - 1) | (1 << (k - 1))


def birthday_schedule(N: int, seed: int, K: int) -> List[int]:
    """Public, value-independent bits c_1..c_{K-1} of Brent's (1986, Sec. 6)
    walk  x_1 = b,  x_{j+1} = x_j^2 * b^{c_j}  (exponent e_{j+1} = 2 e_j + c_j)."""
    rng = random.Random(f"pm1bsgs-bday|{N}|{seed}")
    return [rng.getrandbits(1) for _ in range(max(0, K - 1))]


# --------------------------------------------------------------------------
# polynomial arithmetic over Z/N (audited flint operations only)
# --------------------------------------------------------------------------
class PolyZN:
    def __init__(self, N: int, m: Metrics):
        self.N = N
        self.ctx = flint.fmpz_mod_poly_ctx(N)
        self.m = m

    def linear(self, r: int):
        return self.ctx([(-r) % self.N, 1])

    def mul(self, f, g, df: int, dg: int):
        self.m.inc("poly_mul")
        self.m.inc("poly_eq", mcost(max(df, dg)))
        self.m.hi("max_degree", df + dg)
        return f * g

    def rem(self, f, g, a: int, b: int):
        self.m.inc("poly_rem")
        self.m.inc("poly_eq", rem_cost(a, b))
        return f % g

    def const_value(self, f) -> int:
        c = f.coeffs()
        return int(c[0]) if len(c) else 0

    def derivative(self, f, d: int):
        c = [int(x) for x in f.coeffs()]
        self.m.inc("rn_mul", max(0, len(c) - 1))
        return self.ctx([(i * c[i]) % self.N for i in range(1, len(c))]), d - 1


def product_tree(ring: PolyZN, roots: Sequence[int]):
    """Levels of the subproduct tree; level 0 = linear leaves.  Each level is a
    list of (poly, degree).  Odd nodes are carried up unchanged."""
    level = [(ring.linear(r), 1) for r in roots]
    levels = [level]
    while len(level) > 1:
        nxt = []
        for i in range(0, len(level) - 1, 2):
            (f, df), (g, dg) = level[i], level[i + 1]
            nxt.append((ring.mul(f, g, df, dg), df + dg))
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
        levels.append(level)
    ring.m.hi("stored_coeffs", sum(d + 1 for lv in levels for _, d in lv))
    return levels


def remainder_tree(ring: PolyZN, F, degF: int, levels) -> List[int]:
    """Values F(r_i) at the leaves of the subproduct tree `levels`."""
    top = levels[-1][0]
    R, dR = F, degF
    if dR >= top[1]:
        R, dR = ring.rem(R, top[0], dR, top[1]), top[1] - 1
    cur = [(R, dR)]
    for k in range(len(levels) - 2, -1, -1):
        lv = levels[k]
        nxt = []
        for idx, (g, dg) in enumerate(lv):
            parent = idx // 2
            PR, dPR = cur[parent]
            single = (idx == len(lv) - 1) and (len(lv) % 2 == 1)
            if single or dPR < dg:
                nxt.append((PR, dPR))
            else:
                nxt.append((ring.rem(PR, g, dPR, dg), dg - 1))
        cur = nxt
    return [ring.const_value(R) if dR <= 0 else _eval_small(ring, R, levels[0][i][0])
            for i, (R, dR) in enumerate(cur)]


def _eval_small(ring: PolyZN, R, leaf) -> int:  # pragma: no cover - only if deg R >= 1 at a leaf
    return ring.const_value(R % leaf)


# --------------------------------------------------------------------------
# stage 1
# --------------------------------------------------------------------------
def stage1(N: int, a: int, B1: int, m: Metrics) -> Tuple[str, int]:
    """b = a^E, E = prod_{l<=B1} l^floor(log_l B1).  gcd(b-1, N) every
    STAGE1_CHUNK steps; on gcd == N the chunk is replayed step by step.
    Returns ("factor", d) | ("stage1_N", 0) | ("ok", b)."""
    steps = stage1_steps(B1) if B1 > 1 else []
    b = a
    for c0 in range(0, len(steps), STAGE1_CHUNK):
        chunk = steps[c0:c0 + STAGE1_CHUNK]
        b_start = b
        for l in chunk:
            b = pow(b, l, N)
            m.inc("rn_mul", pow_cost(l))
        g = gcd(b - 1, N)
        m.inc("gcd")
        if g == 1:
            continue
        if g < N:
            return "factor", g
        b = b_start                         # replay chunk with per-step gcds
        for l in chunk:
            b = pow(b, l, N)
            m.inc("backtrack_rn_mul", pow_cost(l))
            g = gcd(b - 1, N)
            m.inc("backtrack_gcd")
            if g == 1:
                continue
            if g < N:
                return "factor", g
            return "stage1_N", 0
    return "ok", b


# --------------------------------------------------------------------------
# stage 2 variants
# --------------------------------------------------------------------------
def bsgs_stage2(N: int, b: int, mbaby: int, J: int, mode: str, m: Metrics,
                block: Optional[int] = None, naive: bool = False,
                target: Optional[int] = None) -> Optional[int]:
    """Search t in the window W (|W| = mbaby*J) with ord_p(b) | t, hash-free.
    mode "N": W = [N-(J-1)m, N+m-1] (t = N+i-jm);  mode "order": W = [1, Jm]."""
    assert mode in ("N", "order") and mbaby >= 1 and J >= 1
    block = block or mbaby
    ring = PolyZN(N, m)
    if mode == "N":
        T = N if target is None else target   # target != N only for the shift control
        u = pow(b, T, N)
        m.inc("rn_mul", pow_cost(T))
    else:
        u = 1
    baby = [u]
    for _ in range(mbaby - 1):
        u = u * b % N
        baby.append(u)
    m.inc("rn_mul", mbaby - 1)
    bm = pow(b, mbaby, N)
    m.inc("rn_mul", pow_cost(mbaby))
    if not naive:
        F, degF = product_tree(ring, baby)[-1][0]
    v = 1 if mode == "N" else bm
    js = list(range(0, J)) if mode == "N" else list(range(1, J + 1))
    first = True
    for c0 in range(0, J, block):
        pts = []
        for _ in js[c0:c0 + block]:
            if not first:
                v = v * bm % N
            first = False
            pts.append(v)
        m.inc("rn_mul", len(pts) - (1 if c0 == 0 else 0))
        if naive:
            vals = []
            for w in pts:
                acc = 1
                for x in baby:
                    acc = acc * (w - x) % N
                vals.append(acc)
            m.inc("rn_mul", len(pts) * len(baby))
        else:
            levels = product_tree(ring, pts)
            vals = remainder_tree(ring, F, degF, levels)
        acc = 1
        for val in vals:
            acc = acc * val % N
        m.inc("rn_mul", len(vals))
        g = gcd(acc, N)
        m.inc("gcd")
        if g == 1:
            continue
        if g < N:
            return g
        for w, val in zip(pts, vals):         # backtrack: per point, then per pair
            g = gcd(val, N)
            m.inc("backtrack_gcd")
            if 1 < g < N:
                return g
            if g == N:
                for x in baby:
                    g = gcd(w - x, N)
                    m.inc("backtrack_gcd")
                    if 1 < g < N:
                        return g
    return None


def prime_stage2(N: int, b: int, B1: int, B2: int, m: Metrics) -> Optional[int]:
    """Classical stage 2: prod over primes B1 < r <= B2 of (b^r - 1)."""
    ps = [r for r in primes_upto(B2) if r > B1]
    if not ps:
        return None
    x = pow(b, ps[0], N)
    m.inc("rn_mul", pow_cost(ps[0]))
    table = {0: 1, 1: b}                      # gap 1 only occurs for 2 -> 3
    b2 = b * b % N
    m.inc("rn_mul")
    table[2] = b2
    maxg = 2

    def step(g: int) -> int:
        nonlocal maxg
        while maxg < g:
            table[maxg + 2] = table[maxg] * b2 % N
            m.inc("rn_mul")
            maxg += 2
        return table[g]

    acc = 1
    for c0 in range(0, len(ps), PRIME2_CHUNK):
        chunk = ps[c0:c0 + PRIME2_CHUNK]
        x_start = x
        for k, r in enumerate(chunk):
            acc = acc * (x - 1) % N
            m.inc("rn_mul")
            nxt = chunk[k + 1] if k + 1 < len(chunk) else (ps[c0 + len(chunk)] if c0 + len(chunk) < len(ps) else None)
            if nxt is not None:
                x = x * step(nxt - r) % N
                m.inc("rn_mul")
        g = gcd(acc, N)
        m.inc("gcd")
        if g == 1:
            continue
        if g < N:
            return g
        x = x_start                            # replay chunk per prime
        for k, r in enumerate(chunk):
            g = gcd(x - 1, N)
            m.inc("backtrack_gcd")
            if 1 < g < N:
                return g
            if k + 1 < len(chunk):
                x = x * step(chunk[k + 1] - r) % N
                m.inc("backtrack_rn_mul")
        return None                            # every hit in this chunk was common to p and q
    return None


def birthday_stage2(N: int, b: int, K: int, seed: int, m: Metrics) -> Optional[int]:
    """Birthday stage 2 (Brent 1986): K points of the value-independent walk
    x_{j+1} = x_j^2 b^{c_j} in <b>; gcd(prod_i F'(x_i), N) tests all K(K-1)/2
    pairs at once (F = prod (T - x_i), F'(x_i) = prod_{j != i} (x_i - x_j))."""
    cs = birthday_schedule(N, seed, K)
    x = b
    pts = [x]
    for c in cs:
        x = x * x % N
        if c:
            x = x * b % N
        pts.append(x)
    m.inc("rn_mul", len(cs) + sum(cs))
    if K < 2:
        return None
    ring = PolyZN(N, m)
    levels = product_tree(ring, pts)
    F, degF = levels[-1][0]
    dF, ddF = ring.derivative(F, degF)
    vals = remainder_tree(ring, dF, ddF, levels)
    acc = 1
    for val in vals:
        acc = acc * val % N
    m.inc("rn_mul", len(vals))
    g = gcd(acc, N)
    m.inc("gcd")
    if g == 1:
        return None
    if g < N:
        return g
    for i, val in enumerate(vals):
        g = gcd(val, N)
        m.inc("backtrack_gcd")
        if 1 < g < N:
            return g
        if g == N:
            for j in range(len(pts)):
                if j != i:
                    g = gcd(pts[i] - pts[j], N)
                    m.inc("backtrack_gcd")
                    if 1 < g < N:
                        return g
    return None


# --------------------------------------------------------------------------
# top-level operational entry points
# --------------------------------------------------------------------------
def _finish(m: Metrics, **res) -> Dict:
    c = m.counts
    m.set("eq_ops", c.get("rn_mul", 0) + c.get("gcd", 0) + c.get("poly_eq", 0))
    m.set("eq_ops_with_backtrack", m.fields["eq_ops"] + c.get("backtrack_rn_mul", 0) + c.get("backtrack_gcd", 0))
    m.hi("max_operand_bits", m.fields.get("N_bits", 0))
    out = {"factor": None, "stage": None, "cls": None}
    out.update(res)
    return out


def _start(N: int, seed: int, B1: int, m: Metrics, base: Optional[int] = None):
    m.set("N_bits", N.bit_length())
    a = choose_base(N, seed) if base is None else base
    g = gcd(a, N)
    m.inc("gcd")
    if g > 1:
        return ("done", _finish(m, factor=g, stage="base", cls="ACCIDENTAL"))
    st, val = stage1(N, a, B1, m)
    if st == "factor":
        return ("done", _finish(m, factor=val, stage="stage1", cls="STRUCTURAL"))
    if st == "stage1_N":
        return ("done", _finish(m, factor=None, stage="stage1_N"))
    return ("ok", val)


def pm1_bsgs(N: int, seed: int, B1: int, mbaby: int, J: int, mode: str = "N",
             m: Optional[Metrics] = None, naive: bool = False, block: Optional[int] = None,
             base: Optional[int] = None, target: Optional[int] = None) -> Dict:
    """Method (a) (B1 > 1) or (d) (B1 <= 1).  `base` overrides the seeded base;
    `target` replaces N in mode "N" (shift control)."""
    m = m if m is not None else NullMetrics()
    st, val = _start(N, seed, B1, m, base)
    if st == "done":
        return val
    d = bsgs_stage2(N, val, mbaby, J, mode, m, block=block, naive=naive, target=target)
    return _finish(m, factor=d, stage="stage2" if d else None, cls="STRUCTURAL" if d else None)


def pm1_prime(N: int, seed: int, B1: int, B2: int, m: Optional[Metrics] = None,
              base: Optional[int] = None) -> Dict:
    """Method (b)."""
    m = m if m is not None else NullMetrics()
    st, val = _start(N, seed, B1, m, base)
    if st == "done":
        return val
    d = prime_stage2(N, val, B1, B2, m)
    return _finish(m, factor=d, stage="stage2" if d else None, cls="STRUCTURAL" if d else None)


def pm1_birthday(N: int, seed: int, B1: int, K: int, m: Optional[Metrics] = None) -> Dict:
    """Method (c)."""
    m = m if m is not None else NullMetrics()
    st, val = _start(N, seed, B1, m)
    if st == "done":
        return val
    d = birthday_stage2(N, val, K, seed, m)
    return _finish(m, factor=d, stage="stage2" if d else None, cls="STRUCTURAL" if d else None)
