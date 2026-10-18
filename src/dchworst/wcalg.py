"""OPERATIONAL: DCL-wc, a worst-case variant of DCL-det (notes/dchworst_results.md, section 3).

Sees only N and public parameters (weights, prime blocks, T, B0).  Never uses p, q,
phi(N), group orders or any factoring/primality routine on N.  Every output is
validated (1 < d < N, d | N).  Arithmetic is exact.

Differences from DCL-det (src/dcltheorem/algorithm.py):
  (0) square / even / small-base checks;
  (1) base selection (Lemma O): for b = 2, 3, ..., B0 find the least r <= T with
      gcd(b^r - 1, N) > 1 by baby-step giant-step polynomial evaluation.  A proper
      gcd is output; if no such r exists, b has order > T modulo every prime factor;
      if gcd = N then ord_p(b) = ord_q(b) = r.  If every b <= B0 is of the last
      kind, Lambda = lcm(r_b) divides p - 1 for every prime p | N and the
      candidates 1 + j Lambda are tried;
  (2) per block, a gamma-order pre-step (Lemma G) on the baby list;
  (3) every giant whose gcd(f(G), N) equals N is EXPLORED down the baby
      subproduct tree (Lemma X); at an exact N-collision leaf the candidate
      S = B p + A q = c + l^2 g is tested by solving x^2 - S x + A B N = 0;
  (4) signed weights (A < 0) and the diagonal weight (1, 1) are allowed.
No clean band, no order assumption mod q, no FULL exclusion is needed.

CORRECTNESS PROTOTYPE (reviews, notes section 8): u_l and b^{E/r^e} use the all-but-one
exponent tree (src.dcltheorem.algorithm.all_but_one), but Lemma G scans r linearly and the
exploration evaluates node polynomials one point at a time (O(m) per N-giant instead of the
batched O~(m + #N-giants) per level of Theorem R (iii)).  Outputs are unaffected.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from src.common.metrics import Metrics
from src.dcltheorem.algorithm import all_but_one


class Found(Exception):
    def __init__(self, d: int, how: str):
        self.d, self.how = d, how


def _chk(x: int, N: int, met, how: str) -> int:
    met.inc("gcdN")
    g = gcd(x % N, N)
    if 1 < g < N:
        raise Found(g, how)
    return g


def _pw(x: int, e: int, N: int, met) -> int:
    """x^e mod N for any integer e (x a unit when e < 0); charged square-and-multiply."""
    a = abs(e)
    if a:
        met.inc("mulN", a.bit_length() + bin(a).count("1") - 2)
    if e < 0:
        met.inc("invN")
        x = pow(x, -1, N)
    return pow(x, a, N)


# ---------------------------------------------------------------------------
# polynomial helpers over Z/N (monic only)
# ---------------------------------------------------------------------------
def _pc(d: int) -> int:
    d = max(0, d)
    return (d + 1) * max(1, (d + 1).bit_length())


def subproduct_tree(ctx, pts: List[int], met) -> List[List]:
    lev = [[ctx([(-x) % ctx.modulus(), 1]) for x in pts]]
    while len(lev[-1]) > 1:
        cur, nxt = lev[-1], []
        for i in range(0, len(cur) - 1, 2):
            met.inc("poly_cost", _pc(cur[i].degree() + cur[i + 1].degree()))
            nxt.append(cur[i] * cur[i + 1])
        if len(cur) % 2:
            nxt.append(cur[-1])
        lev.append(nxt)
    return lev


def multipoint(ctx, f, pts: List[int], met) -> List[int]:
    if not pts:
        return []
    lev = subproduct_tree(ctx, pts, met)
    rems = [f % lev[-1][0]]
    met.inc("poly_cost", _pc(f.degree()))
    for k in range(len(lev) - 2, -1, -1):
        new = []
        for i, node in enumerate(lev[k]):
            r = rems[i // 2]
            met.inc("poly_cost", _pc(r.degree()))
            new.append(r % node if r.degree() >= node.degree() else r)
        rems = new
    out = []
    for r in rems:
        c = r.coeffs()
        out.append(int(c[0]) if c else 0)
    return out


def _children(lev, k: int, i: int):
    """Children (level k-1 indices) of node i at level k of a subproduct tree."""
    below = len(lev[k - 1])
    return [c for c in (2 * i, 2 * i + 1) if c < below]


# ---------------------------------------------------------------------------
# Lemma O: least r <= T with gcd(b^r - 1, N) > 1 (BSGS by polynomial evaluation)
# ---------------------------------------------------------------------------
def first_r(b: int, T: int, N: int, met) -> Tuple[str, int]:
    """Scans 1 <= r <= T_scan = s*ceil(T/s), s = isqrt(T) + 1 (T <= T_scan = T + O(sqrt T)).
    Returns ('large', 0)  if gcd(b^r - 1, N) = 1 for all 1 <= r <= T_scan,
            ('equal', r) if the least such r (possibly in (T, T_scan]) has gcd = N.
    Raises Found on a proper divisor.  Cost O~(sqrt T) ring operations."""
    s = isqrt(T) + 1
    nj = -(-T // s)
    ctx = flint.fmpz_mod_poly_ctx(N)
    binv = _pw(b, -1, N, met)
    roots, x = [], 1
    for _ in range(s):                       # b^{-i}, i = 1..s
        x = x * binv % N
        met.inc("mulN")
        roots.append(x)
    f = subproduct_tree(ctx, roots, met)[-1][0]
    bs = _pw(b, s, N, met)
    X, pts = 1, []
    for _ in range(nj):                      # X_j = b^{j s}
        pts.append(X)
        X = X * bs % N
        met.inc("mulN")
    vals = multipoint(ctx, f, pts, met)
    for j, v in enumerate(vals):
        met.inc("gcdN")
        if gcd(v, N) > 1:
            y = pts[j]
            for i in range(1, s + 1):        # b^{js+i}, in increasing r
                y = y * b % N
                met.inc("mulN")
                g = _chk(y - 1, N, met, "first_r")
                if g == N:
                    return "equal", j * s + i
            raise AssertionError("block gcd > 1 but no element gcd > 1")
    return "large", 0


def select_base(N: int, T: int, B0: int, met) -> Tuple[int, Dict]:
    """Lemma O.  Returns a base b with ord_p(b) > T for every prime p | N, or raises
    Found.  If every b in [2, B0] has equal small order, tries 1 + j*Lambda."""
    info = {"equal": {}}
    lam = 1
    for b in range(2, B0 + 1):
        _chk(b, N, met, "small_base")
        kind, r = first_r(b, T, N, met)
        if kind == "large":
            info["base"] = b
            return b, info
        info["equal"][b] = r
        lam = lam * r // gcd(lam, r)
    info["Lambda"] = lam
    j = 1
    while 1 + j * lam <= isqrt(N):
        _chk(1 + j * lam, N, met, "lambda_enum")
        met.inc("lambda_cands")
        j += 1
    # Unreachable for N = pq (Lemma O: Lambda | p - 1, so p = 1 + j*Lambda was tried).
    raise Found(-1, "all_bases_equal_FAIL")


# ---------------------------------------------------------------------------
# targets of one block (signed weights allowed)
# ---------------------------------------------------------------------------
def weights(H: int, signed: bool, diag: bool) -> List[Tuple[int, int]]:
    out = []
    for B in range(1, H + 1):
        for a in range(1, H + 1):
            if gcd(a, B) != 1:
                continue
            if a != B or diag:
                out.append((a, B))
            if signed:
                out.append((-a, B))
    return out


def balanced_range(N: int) -> Tuple[int, int]:
    a0 = isqrt(N // 2)
    while 2 * a0 * a0 < N:
        a0 += 1
    return a0, isqrt(N)


def f_bounds(N: int, A: int, B: int, a0: int, b0: int) -> Tuple[int, int]:
    """Lo <= B x + A N / x <= Hi for x in [a0, b0]."""
    if A > 0:   # convex, AM-GM lower bound
        return isqrt(4 * A * B * N), max(B * x + -(-(A * N) // x) for x in (a0, b0))
    # A < 0: increasing in x
    return B * a0 + (A * N) // a0, B * b0 + -(-(A * N) // b0)


def _sqrt3mod4(a: int, l: int, met) -> int:
    met.inc("mul_small", 2 * l.bit_length())
    return pow(a, (l + 1) // 4, l)


def build_targets(N: int, band: Sequence[int], W: Sequence[Tuple[int, int]], met) -> List[Dict]:
    a0, b0 = balanced_range(N)
    bounds = {v: f_bounds(N, v[0], v[1], a0, b0) for v in W}
    out = []
    for l in band:
        if N % l == 0:
            if 1 < l < N:
                raise Found(l, "band_prime")
            continue
        l2 = l * l
        for (A, B) in W:
            if (A * B) % l == 0:
                continue
            a = A * N * pow(B, -1, l) % l
            met.inc("mul_small", 3)
            if a == 0 or pow(a, (l - 1) // 2, l) != 1:
                continue
            r = _sqrt3mod4(a, l, met)
            lo, hi = bounds[(A, B)]
            for root in sorted({r, l - r}):
                s = N * pow(root, -1, l2) % l2
                c = (B * root + A * s) % l2
                g_lo = -((-(lo - c)) // l2)
                g_hi = (hi - c) // l2
                out.append({"A": A, "B": B, "l": l, "c": c, "g_lo": g_lo,
                            "W": max(1, g_hi - g_lo + 1)})
    return out


def choose_m(targets: List[Dict]) -> int:
    if not targets:
        return 1
    Wmax = max(t["W"] for t in targets)
    best, m = None, 1
    while True:
        mm = min(m, Wmax)
        cst = mm + sum(-(-t["W"] // mm) for t in targets)
        if best is None or cst < best[0]:
            best = (cst, mm)
        if mm == Wmax:
            break
        m = (119 * m) // 100 + 1
    return best[1]


def _candidate(N: int, A: int, B: int, S: int, met) -> Optional[int]:
    """Solve x^2 - S x + A B N = 0 (roots B p, A q when S = B p + A q)."""
    met.inc("isqrt_ops")
    disc = S * S - 4 * A * B * N
    if disc < 0:
        return None
    D = isqrt(disc)
    if D * D != disc:
        return None
    for x in ((S + D) // 2, (S - D) // 2):
        d = gcd(x, N)
        if 1 < d < N:
            return d
    return None


# ---------------------------------------------------------------------------
# Lemma G (gamma pre-step) and the equal-order split
# ---------------------------------------------------------------------------
def _small_factor(o: int) -> Dict[int, int]:
    """Trial division of the PUBLIC integer o (o < m)."""
    f, d = {}, 2
    while d * d <= o:
        while o % d == 0:
            f[d] = f.get(d, 0) + 1
            o //= d
        d += 1
    if o > 1:
        f[o] = f.get(o, 0) + 1
    return f


def equal_order_split(N: int, b: int, band: Sequence[int], o: int, T: int, met) -> None:
    """gamma = b^L has equal order o mod p and q.  E = L*o is a FACTORED multiple of
    ord_N(b).  Either a valuation of ord_p(b), ord_q(b) differs (proper gcd), or
    ord_p(b) = ord_q(b) = D | gcd(p-1, q-1) with D > T (Lemma O) and 1 + jD is tried."""
    fac = {l: 2 for l in band}
    for r, e in _small_factor(o).items():
        fac[r] = fac.get(r, 0) + e
    E = 1
    for r, e in fac.items():
        E *= r ** e
    D = 1
    rs = list(fac)
    ys, _ = all_but_one(b, [r ** fac[r] for r in rs], N, met)   # b^{E / r^e} for all r
    for r, y in zip(rs, ys):
        e = fac[r]
        v = None
        for j in range(e + 1):
            if j:
                y = _pw(y, r, N, met)
            g = _chk(y - 1, N, met, "order_split")
            if g == N and v is None:
                v = j
        D *= r ** v
    met.set("equal_order_D", D)
    j = 1
    while 1 + j * D <= isqrt(N):
        _chk(1 + j * D, N, met, "D_enum")
        j += 1
    # Unreachable for N = pq (Lemma G: D | p - 1, so p = 1 + j*D was tried).
    raise Found(-1, "equal_order_FAIL")


# ---------------------------------------------------------------------------
# one block
# ---------------------------------------------------------------------------
def run_block(N: int, b: int, band: Sequence[int], W, T: int, met, log: Dict) -> None:
    targets = build_targets(N, band, W, met)
    if not targets:
        return
    m = choose_m(targets)
    u_list, L = all_but_one(b, [l * l for l in band], N, met)  # u_l = b^{L/l^2}
    u = dict(zip(band, u_list))
    gamma = _pw(u[band[0]], band[0] ** 2, N, met)
    _chk(gamma, N, met, "gamma")
    # babies and Lemma G: least r in [1, m) with gcd(gamma^r - 1, N) > 1
    babies, x = [], 1
    for _ in range(m):
        babies.append(x)
        x = x * gamma % N
        met.inc("mulN")
    o = None
    for r in range(1, m):        # (a product tree + binary descent in the theorem)
        met.inc("gcdN")
        g = gcd(babies[r] - 1, N)
        if g > 1:
            if g < N:
                raise Found(g, "gamma_order")
            o = r
            break
    if o is not None:
        log["equal_order"] = o
        equal_order_split(N, b, band, o, T, met)
    ginv = _pw(gamma, -1, N, met)
    gm = _pw(ginv, m, N, met)
    giants, meta = [], []
    for ti, t in enumerate(targets):
        A, B, l, c, g_lo, Wt = t["A"], t["B"], t["l"], t["c"], t["g_lo"], t["W"]
        G = _pw(u[l], B + A * N - c, N, met) * _pw(gamma, -g_lo, N, met) % N
        for j in range(-(-Wt // m)):
            giants.append(G)
            meta.append((ti, j))
            G = G * gm % N
            met.inc("mulN")
    met.inc("giants", len(giants))
    met.inc("babies", m)
    met.inc("targets", len(targets))
    ctx = flint.fmpz_mod_poly_ctx(N)
    lev = subproduct_tree(ctx, babies, met)
    vals = multipoint(ctx, lev[-1][0], giants, met)
    Ng = []
    for h, v in enumerate(vals):
        g = _chk(v, N, met, "collision")
        if g == N:
            Ng.append(h)
    log["N_giants"] = log.get("N_giants", 0) + len(Ng)
    # Lemma X: explore N-giants down the tree; at most one N-node per level per
    # giant when ord_p(gamma) >= m (Lemma G); leaves give exact candidates.
    top = len(lev) - 1
    frontier = [(h, 0) for h in Ng]          # (giant, node index at level `k`)
    k = top
    while frontier and k > 0:
        nxt = []
        for (h, i) in frontier:
            for cidx in _children(lev, k, i):
                node = lev[k - 1][cidx]
                met.inc("explore_evals")
                met.inc("poly_cost", node.degree() + 1)
                g = _chk(int(node(giants[h])), N, met, "explore")
                if g == N:
                    nxt.append((h, cidx))
        met.hi("explore_frontier", len(nxt))
        frontier = nxt
        k -= 1
    for (h, i) in frontier:                  # leaves: exact N-collisions at baby i
        ti, j = meta[h]
        t = targets[ti]
        gval = t["g_lo"] + j * m + i
        met.inc("leaf_candidates")
        d = _candidate(N, t["A"], t["B"], t["c"] + t["l"] ** 2 * gval, met)
        if d:
            raise Found(d, "exact_leaf")


def dcl_wc(N: int, blocks: Sequence[Sequence[int]], H: int, T: int, B0: int,
           signed: bool = True, diag: bool = True, met: Optional[Metrics] = None,
           force_base: Optional[int] = None) -> Dict:
    """force_base (TEST MODE): skip Lemma O and use this base; correctness on covered
    pairs then still holds except in the equal-order branch, whose enumeration
    1 + jD is complete but not cost-bounded."""
    met = met or Metrics()
    log: Dict = {}
    if N < 4:
        return {"factor": None, "status": "fail_small", "log": log}
    try:
        if N % 2 == 0:
            raise Found(2, "even")
        s = isqrt(N)
        if s * s == N:
            raise Found(s, "square")
        if force_base is None:
            b, info = select_base(N, T, B0, met)
        else:
            b = force_base
            _chk(b, N, met, "small_base")
        log["base"] = b
        W = weights(H, signed, diag)
        for k, band in enumerate(blocks):
            if band:
                run_block(N, b, band, W, T, met, log)
    except Found as e:
        d = e.d
        ok = isinstance(d, int) and 1 < d < N and N % d == 0
        return {"factor": d if ok else None, "status": e.how if ok else "fail_" + e.how,
                "log": log}
    return {"factor": None, "status": "fail", "log": log}
