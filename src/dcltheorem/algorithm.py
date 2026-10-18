"""OPERATIONAL: DCL-det, the deterministic divisor-cover Lehman algorithm.

Receives only N and a public Design (a function of the bit length n; see
params.py).  Never uses p, q, phi(N), group orders or any factoring /
primality routine.  Base beta = 2 (fixed, no randomness).  Every inversion
mod N is gcd-checked, every nontrivial gcd is returned.

Algorithm (notes/dcltheorem_results.md, section 2), for each block k = 0..K-1:
  1. L = prod_{l in P_k} l^2;  u_l = 2^{L/l^2} for all l by a balanced
     all-but-one exponent tree;  gamma = 2^L.
  2. For each l in P_k and (A, B) in V: a = A N B^{-1} mod l; if a is a
     nonzero square mod l (Euler), r = a^{(l+1)/4} (l == 3 mod 4) and the
     roots are {r, l - r}.  For each root: s = N r^{-1} mod l^2,
     c = (B r + A s) mod l^2, whole-range window [g_lo, g_hi] from the exact
     integer bounds Lo_v <= B p + A N/p <= Hi_v on the balanced range, and
     giants  G_{tau,j} = u_l^{B + A N - c} gamma^{-(g_lo + j m)},
     0 <= j < ceil(W/m).
  3. Babies gamma^i, 0 <= i < m;  f(X) = prod (X - gamma^i);  f(G) for every
     giant by a remainder tree;  d = gcd(f(G), N);  output d if 1 < d < N.
Output None (FAIL) if no block produced a factor.

Cost counters (Metrics): mulN (multiplications mod N, square-and-multiply
counted exactly), mul_small (mults mod l or l^2), bigmul_bits (sum of operand
bit lengths of integer products in the exponent product tree), poly_cost /
poly_ops (stage 3; see stage3.py), gcdN, invN, isqrt_ops.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint
import gmpy2

from src.common.metrics import Metrics
from .params import Design
from .stage3 import eval_baby_poly_at_giants


class Found(Exception):
    def __init__(self, d: int, how: str):
        self.d, self.how = d, how


def _gcd_check(x: int, N: int, met, how: str) -> int:
    met.inc("gcdN")
    g = gcd(x % N, N)
    if 1 < g < N:
        raise Found(g, how)
    return g


def powN(x: int, e: int, N: int, met) -> int:
    """x^e mod N, e >= 0; charges bits(e)-1 squarings + popcount(e)-1 mults."""
    if e <= 0:
        return 1 % N if e == 0 else None  # negative exponents are not used
    met.inc("mulN", e.bit_length() + bin(e).count("1") - 2)
    return pow(x, e, N)


def pow_small(x: int, e: int, M: int, met) -> int:
    met.inc("mul_small", max(0, e.bit_length() + bin(e).count("1") - 2))
    return pow(x, e, M)


# ---------------------------------------------------------------------------
# balanced all-but-one exponent tree
# ---------------------------------------------------------------------------
def int_product_tree(xs: Sequence[int], met) -> List[List[int]]:
    # gmpy2/GMP multiplication (FFT for large operands), so the tree costs
    # O(M(log L) log #leaves) bit operations; CPython ints would be Karatsuba.
    levels = [[gmpy2.mpz(x) for x in xs]]
    while len(levels[-1]) > 1:
        cur, nxt = levels[-1], []
        for i in range(0, len(cur) - 1, 2):
            met.inc("bigmul_bits", cur[i].bit_length() + cur[i + 1].bit_length())
            nxt.append(cur[i] * cur[i + 1])
        if len(cur) % 2:
            nxt.append(cur[-1])
        levels.append(nxt)
    return levels


def all_but_one(base: int, exps: Sequence[int], N: int, met) -> Tuple[List[int], int]:
    """([base^{E/e_j} for j], E) with E = prod exps, via the balanced product
    tree: a node holding base^{E/prod(node)} passes base^{E/prod(node)*prod(sib)}
    to each child.  Per level the exponent bit lengths sum to <= log2 E + #nodes."""
    levels = int_product_tree(exps, met)
    E = int(levels[-1][0])
    vals = [base % N]
    for lev in range(len(levels) - 2, -1, -1):
        cur = levels[lev]
        parent_count = len(levels[lev + 1])
        new = [None] * len(cur)
        for pi in range(parent_count):
            li, ri = 2 * pi, 2 * pi + 1
            if ri < len(cur):
                new[li] = powN(vals[pi], int(cur[ri]), N, met)
                new[ri] = powN(vals[pi], int(cur[li]), N, met)
            else:
                new[li] = vals[pi]          # odd node carried up unchanged
        vals = new
    return vals, E


# ---------------------------------------------------------------------------
# windows (exact integers)
# ---------------------------------------------------------------------------
def balanced_range(N: int) -> Tuple[int, int]:
    """Integers p with p^2 < N <= 2 p^2 lie in [a0, b0]."""
    a0 = isqrt(N // 2)
    while 2 * a0 * a0 < N:
        a0 += 1
    b0 = isqrt(N)
    return a0, b0


def f_bounds(N: int, A: int, B: int, a0: int, b0: int, met) -> Tuple[int, int]:
    """Lo <= B x + A N / x <= Hi for all real x in [a0, b0]:  Lo = floor(2 sqrt(ABN))
    (AM-GM), Hi = max over the endpoints of ceil(F) (convexity)."""
    met.inc("isqrt_ops")
    lo = isqrt(4 * A * B * N)
    hi = max(B * x + -(-(A * N) // x) for x in (a0, b0))
    return lo, hi


# ---------------------------------------------------------------------------
# one block
# ---------------------------------------------------------------------------
_NONRES: Dict[int, int] = {}


def sqrt_prime(a: int, l: int, met) -> int:
    """A square root of the nonzero square a mod the public prime l.
    l == 3 (mod 4): a^{(l+1)/4}.  Otherwise Tonelli--Shanks with the least
    quadratic nonresidue z (found by Euler's criterion; the search is charged
    once per l -- its total over a band is O(#band) on average by Erdos 1961 /
    Pollack 2012; only the secondary design mod4=False uses this branch)."""
    if l % 4 == 3:
        return pow_small(a, (l + 1) // 4, l, met)
    if l not in _NONRES:
        z = 2
        while pow_small(z, (l - 1) // 2, l, met) != l - 1:
            z += 1
        _NONRES[l] = z
        met.inc("nonres_search", z - 1)
    z = _NONRES[l]
    Q, S = l - 1, 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    c, x, t, R = pow_small(z, Q, l, met), pow_small(a, (Q + 1) // 2, l, met), pow_small(a, Q, l, met), S
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % l
            i += 1
            met.inc("mul_small")
        b = pow_small(c, 1 << (R - i - 1), l, met)
        x, c, t, R = x * b % l, b * b % l, t * b * b % l, i
        met.inc("mul_small", 4)
    return x


def build_targets(N: int, design: Design, k: int, met) -> List[Dict]:
    band = design.bands[k]
    _NONRES.clear()                       # nonresidue searches are charged per run
    a0, b0 = balanced_range(N)
    bounds = {v: f_bounds(N, v[0], v[1], a0, b0, met) for v in design.V}
    out = []
    for l in band:
        if N % l == 0:
            raise Found(l, "band_prime")
        l2 = l * l
        Nl = N % l
        for (A, B) in design.V:
            met.inc("vl_pairs")
            if (A * B) % l == 0:          # finite-size guard (never triggers when Y > H)
                continue
            a = A * Nl * pow(B, -1, l) % l
            met.inc("mul_small", 3)
            if pow_small(a, (l - 1) // 2, l, met) != 1:
                continue
            r = sqrt_prime(a, l, met)
            lo, hi = bounds[(A, B)]
            for root in sorted({r, l - r}):
                s = N * pow(root, -1, l2) % l2
                c = (B * root + A * s) % l2
                met.inc("mul_small", 4)
                g_lo = -((-(lo - c)) // l2)
                g_hi = (hi - c) // l2
                out.append({"A": A, "B": B, "l": l, "r": root, "c": c,
                            "g_lo": g_lo, "W": max(1, g_hi - g_lo + 1)})
    return out


def choose_m(targets: List[Dict], mbar: int) -> int:
    """Public choice of the baby length: minimise m + sum ceil(W/m) over a
    geometric grid in [1, max W] (never above max W; max W <= mbar)."""
    if not targets:
        return 1
    Wmax = max(t["W"] for t in targets)
    best, m = None, 1
    while True:
        mm = min(m, Wmax)
        c = mm + sum(-(-t["W"] // mm) for t in targets)
        if best is None or c < best[0]:
            best = (c, mm)
        if mm == Wmax:
            break
        m = (119 * m) // 100 + 1            # integer-only grid (no float overflow)
    return best[1]


def run_block(N: int, design: Design, k: int, met, record: bool = False,
              m_override: Optional[int] = None) -> Dict:
    band = design.bands[k]
    if not band:
        return {"factor": None, "status": "empty_band"}
    met.start(f"block{k}")
    targets = build_targets(N, design, k, met)
    m = m_override or choose_m(targets, design.mbar[k])
    sq = [l * l for l in band]
    u_list, L = all_but_one(2, sq, N, met)
    u = dict(zip(band, u_list))
    gamma = powN(u_list[0], sq[0], N, met)              # 2^L
    _gcd_check(gamma, N, met, "gamma")
    met.inc("invN")
    ginv = pow(gamma, -1, N)
    gm_inv = powN(ginv, m, N, met)
    Nmod = N
    giants, meta = [], []
    # per-l precomputation w_l = u_l^N (charged once per band prime used), then
    # G = u^{B + A N - c} gamma^{-g_lo} = w^A * u^{B - c + l^2} * gamma^{-(g_lo + 1)}
    # (u^{l^2} = gamma), all exponents nonnegative and O(log(H l^2)) bits.
    w = {}
    for t in targets:
        if t["l"] not in w:
            w[t["l"]] = powN(u[t["l"]], N, N, met)
    for ti, t in enumerate(targets):
        A, B, l, c, g_lo, W = t["A"], t["B"], t["l"], t["c"], t["g_lo"], t["W"]
        G = (powN(w[l], A, N, met) * powN(u[l], B - c + l * l, N, met) % N
             * powN(ginv, g_lo + 1, N, met) % N)
        met.inc("mulN", 2)
        for j in range(-(-W // m)):
            giants.append(G)
            meta.append((ti, j))
            G = G * gm_inv % N
            met.inc("mulN")
    met.inc("giants", len(giants))
    met.inc("babies", m)
    met.inc("targets", len(targets))
    vals = eval_baby_poly_at_giants(N, gamma, m, giants, met)
    res = {"factor": None, "status": "none", "m": m, "n_targets": len(targets),
           "n_giants": len(giants), "L_bits": L.bit_length()}
    gs = []
    for h, v in enumerate(vals):
        met.inc("gcdN")
        g = gcd(v, N)
        if record:
            gs.append(g)
        if 1 < g < N and res["factor"] is None:
            ti, j = meta[h]
            res.update({"factor": g, "status": "collision", "hit": {"target": ti, "j": j}})
            if not record:
                break
    if record:
        res["targets"] = targets
        res["giant_meta"] = meta
        res["giant_gcds"] = gs
        res["gamma"] = gamma
        res["u"] = u
    met.stop(f"block{k}")
    return res


def dcl_det(N: int, design: Design, met: Optional[Metrics] = None, early_exit: bool = True,
            record: bool = False, blocks: Optional[Sequence[int]] = None) -> Dict:
    """Deterministic DCL.  Returns {"factor": d or None, "block": k, "blocks": [...]}.
    With early_exit, stops after the first block that yields a factor; the
    per-block metrics are recorded either way (prefix b{k}_)."""
    met = met or Metrics()
    out = {"factor": None, "status": "fail", "block": None, "blocks": []}
    try:
        _gcd_check(2, N, met, "even")
        if N % 2 == 0:
            return {"factor": 2, "status": "even", "block": None, "blocks": []}
    except Found as e:
        return {"factor": e.d, "status": "gcd_" + e.how, "block": None, "blocks": []}
    for k in (blocks if blocks is not None else range(design.K)):
        bm = Metrics()
        try:
            r = run_block(N, design, k, bm, record=record)
        except Found as e:
            r = {"factor": e.d, "status": "gcd_" + e.how}
        r["k"] = k
        r["met"] = bm.to_record()
        met.merge(bm)
        out["blocks"].append(r)
        if r["factor"] and out["factor"] is None:
            out.update({"factor": r["factor"], "status": r["status"], "block": k})
            if early_exit:
                break
    return out


def metric_E(counts: Dict[str, int], n: int) -> int:
    """Pre-registered bit-operation proxy in n-bit-modular-multiplication units:
    E = mulN + mul_small + poly_cost + ceil(log2 n) * gcdN + ceil(bigmul_bits / n)
        + ceil(log2 n) * invN."""
    lg = max(1, (n - 1).bit_length())
    return (counts.get("mulN", 0) + counts.get("mul_small", 0) + counts.get("poly_cost", 0)
            + lg * (counts.get("gcdN", 0) + counts.get("invN", 0))
            + -(-counts.get("bigmul_bits", 0) // n))
