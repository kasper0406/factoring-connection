"""Track C baselines for the factorial / interval unit predicate (owner factC).

OPERATIONAL code: every function receives only N, the public interval [A, B)
and configuration.  No factor information is used.  Every scalar inversion is
gcd-checked (``src.common.ring.inv_mod``); a nonunit is returned as output.

Predicates
  U(K)    := [gcd((K-1)!, N) = 1]                 (prefix; = interval [1, K))
  I(A,B)  := [gcd(prod_{A<=j<B} j, N) > 1]        (interval contains a multiple
                                                   of some prime factor of N)
Methods (all compute the product value mod N as well, for cross-checking):
  naive : running product, gcd once per block of `block` factors (batched)
  ps    : Pollard-Strassen: product tree for F(x)=prod_{i<b}(x+A+i), monic
          remainder tree at x = 0, b, 2b, ...  (O(M(sqrt L) log L))
  bgs   : Bostan-Gaudry-Schost style: values of g_d(x)=prod_{j<d}(A+v x+j) at
          x=0..d maintained by Lagrange *shifting* of evaluations (BGS Thm 5),
          doubling d -> 2d needs 3 shifts, each one middle-product (O(M(d))).
          Total O(M(sqrt L)).  BGS needs O(sqrt L) auxiliary integers to be
          units; they are gcd-checked and a hit is reported as ``aux_hit``.

Result dict keys: method, N, A, B, value (prod mod N; None if aborted by an
aux hit), nonunit (bool, = I(A,B)), divisor (gcd found, or None), witness
(an integer whose gcd with N is the divisor), aux_hit (bool), counters.
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import NonUnit, inv_mod

from .polyback import PolyRing, multipoint_eval, product_tree_linear


# ---------------------------------------------------------------------------
# scalar helpers
# ---------------------------------------------------------------------------
def scalar_product_tree(vals: Sequence[int], N: int, m: Metrics) -> List[List[int]]:
    level = [v % N for v in vals]
    levels = [level]
    while len(level) > 1:
        nxt = [level[i] * level[i + 1] % N for i in range(0, len(level) - 1, 2)]
        m.inc("rn_mul", len(level) // 2)
        if len(level) % 2:
            nxt.append(level[-1])
        levels.append(nxt)
        level = nxt
    return levels


def find_nonunit(vals: Sequence[int], N: int, m: Metrics) -> Optional[int]:
    """BGS Lemma 12: index of some vals[i] with gcd(vals[i], N) > 1, or None.

    One product tree (len(vals) multiplications) and O(log len) gcds.
    """
    if not vals:
        return None
    levels = scalar_product_tree(vals, N, m)
    m.inc("gcd")
    if gcd(levels[-1][0], N) == 1:
        return None
    idx = 0
    for depth in range(len(levels) - 2, -1, -1):
        level = levels[depth]
        left = 2 * idx
        if left + 1 >= len(level):
            idx = left
            continue
        m.inc("gcd")
        idx = left if gcd(level[left], N) > 1 else left + 1
    return idx


def batch_inverse(vals: Sequence[int], N: int, m: Metrics) -> List[int]:
    """Montgomery batch inversion with a single gcd-checked inversion.

    Raises NonUnit(divisor, N, value) for a located nonunit element.
    """
    n = len(vals)
    if n == 0:
        return []
    pref = [0] * n
    acc = 1
    for i, v in enumerate(vals):
        acc = acc * v % N
        pref[i] = acc
    m.inc("rn_mul", n)
    try:
        inv = inv_mod(acc, N, m)
    except NonUnit:
        i = find_nonunit(vals, N, m)
        g = gcd(vals[i], N)
        raise NonUnit(g, N, vals[i])
    out = [0] * n
    for i in range(n - 1, 0, -1):
        out[i] = inv * pref[i - 1] % N
        inv = inv * vals[i] % N
    out[0] = inv
    m.inc("rn_mul", 2 * (n - 1))
    return out


def _result(method, N, A, B, value, divisor, witness, aux_hit, m: Metrics, **extra) -> Dict:
    rec = {"method": method, "N": N, "A": A, "B": B, "value": value,
           "nonunit": (divisor is not None and not aux_hit) or (value is not None and gcd(value, N) > 1),
           "divisor": divisor, "witness": witness, "aux_hit": aux_hit}
    rec.update(extra)
    rec.update(m.to_record())
    return rec


def _localize_in_range(N: int, lo: int, hi: int, m: Metrics) -> Tuple[int, int]:
    """Find j in [lo, hi) with gcd(j, N) > 1 (caller guarantees existence)."""
    vals = list(range(lo, hi))
    i = find_nonunit(vals, N, m)
    assert i is not None
    return vals[i], gcd(vals[i], N)


# ---------------------------------------------------------------------------
# (a) naive running product with batched gcds
# ---------------------------------------------------------------------------
def naive_interval(N: int, A: int, B: int, m: Optional[Metrics] = None, block: int = 512) -> Dict:
    m = m if m is not None else Metrics()
    value = 1
    divisor = witness = None
    j = A
    while j < B:
        hi = min(B, j + block)
        bp = 1
        for t in range(j, hi):
            bp = bp * t % N
        m.inc("rn_mul", hi - j)
        m.inc("gcd")
        if divisor is None and gcd(bp, N) > 1:
            witness, divisor = _localize_in_range(N, j, hi, m)
        value = value * bp % N
        m.inc("rn_mul")
        j = hi
    return _result("naive", N, A, B, value, divisor, witness, False, m)


# ---------------------------------------------------------------------------
# (b) Pollard-Strassen
# ---------------------------------------------------------------------------
def ps_interval(N: int, A: int, B: int, m: Optional[Metrics] = None, backend: str = "flint",
                block_len: Optional[int] = None) -> Dict:
    m = m if m is not None else Metrics()
    L = B - A
    if L <= 16:
        return dict(naive_interval(N, A, B, m), method="ps")
    R = PolyRing(N, backend, m)
    b = block_len or isqrt(L)
    nb = L // b
    # F(x) = prod_{i<b} (x + A + i)
    m.start("build")
    F = product_tree_linear(R, [A + i for i in range(b)])[-1][0]
    m.stop("build")
    m.start("eval")
    vals = multipoint_eval(R, F, [j * b for j in range(nb)])
    m.stop("eval")
    tail = list(range(A + nb * b, B))
    allv = vals + tail
    value = 1
    for v in allv:
        value = value * v % N
    m.inc("rn_mul", len(allv))
    m.inc("gcd")
    divisor = witness = None
    if gcd(value, N) > 1:
        i = find_nonunit(allv, N, m)
        if i < nb:
            witness, divisor = _localize_in_range(N, A + i * b, A + (i + 1) * b, m)
        else:
            witness = allv[i]
            divisor = gcd(witness, N)
    return _result("ps", N, A, B, value, divisor, witness, False, m, block=b, nblocks=nb,
                   backend=backend)


# ---------------------------------------------------------------------------
# (c) BGS shifting of evaluations
# ---------------------------------------------------------------------------
class _BGS:
    def __init__(self, N: int, R: PolyRing, m: Metrics, fact_inv: List[int], small_inv: List[int]):
        self.N, self.R, self.m = N, R, m
        self.fact_inv = fact_inv      # 1/i! for i <= dmax
        self.small_inv = small_inv    # 1/i for 1 <= i <= 2*dmax+1

    def shift(self, vals: List[int], a: int, a_int_witness) -> List[int]:
        """Values P(0..d) of deg<=d poly P -> values P(a..a+d) (BGS Thm 5).

        Needs a-d, ..., a+d to be units (checked; NonUnit raised, with
        ``a_int_witness(t)`` giving the integer whose unit status equals that
        of a+t).  1..d must be units (precomputed by caller).
        """
        N, R, m = self.N, self.R, self.m
        d = len(vals) - 1
        # Ptilde_i = P(i) / (i! (d-i)! (-1)^(d-i))
        fi = self.fact_inv
        pt = []
        for i, v in enumerate(vals):
            c = v * fi[i] % N * fi[d - i] % N
            pt.append(c if (d - i) % 2 == 0 else (-c) % N)
        m.inc("rn_mul", 2 * (d + 1))
        # S_j = 1/(a - d + j), j = 0..2d
        S_arg = [(a - d + j) % N for j in range(2 * d + 1)]
        try:
            S = batch_inverse(S_arg, N, m)
        except NonUnit as e:
            t = S_arg.index(e.value % N) - d
            raise NonUnit(e.divisor, N, a_int_witness(t))
        conv = R.mul(R.make(pt), R.make(S))
        cc = R.coeffs(conv, 3 * d + 1)
        Q = cc[d:2 * d + 1]
        # Delta_k = prod_{j=0}^{d} (a + k - j)
        delta = 1
        for j in range(d + 1):
            delta = delta * ((a - j) % N) % N
        out = []
        for k in range(d + 1):
            out.append(Q[k] * delta % N)
            if k < d:
                # Delta_{k+1} = Delta_k * (a+k+1) / (a+k-d);  1/(a+k-d) = S_k
                delta = delta * ((a + k + 1) % N) % N * S[k] % N
        m.inc("rn_mul", (d + 1) + (d + 1) + 2 * d)
        return out


def bgs_interval(N: int, A: int, B: int, m: Optional[Metrics] = None, backend: str = "flint") -> Dict:
    m = m if m is not None else Metrics()
    L = B - A
    if L <= 64:
        return dict(naive_interval(N, A, B, m), method="bgs")
    R = PolyRing(N, backend, m)
    v = isqrt(L)
    # --- auxiliary unit conditions on small integers 1..2v+3 and v --------
    small = list(range(1, 2 * v + 4))
    i = find_nonunit(small, N, m)
    if i is not None:
        w = small[i]
        return _result("bgs", N, A, B, None, gcd(w, N), w, True, m)
    small_inv = [0] + batch_inverse(small, N, m)
    fact_inv = [1] * (v + 2)
    f = 1
    for k in range(1, v + 2):
        f = f * small_inv[k] % N
        fact_inv[k] = f
    m.inc("rn_mul", v + 1)
    vinv = small_inv[v]
    eng = _BGS(N, R, m, fact_inv, small_inv)

    bits = bin(v)[2:]
    d = 1
    vals = [A % N, (A + v) % N]                     # g_1(x) = A + v x at x = 0, 1
    try:
        for bit in bits[1:]:
            # ---- doubling d -> 2d:  g_2d(x) = g_d(x) g_d(x + d/v)
            ext = eng.shift(vals, d + 1, lambda t, d=d: d + 1 + t)
            a1 = d * vinv % N
            s1 = eng.shift(vals, a1, lambda t, d=d: d + v * t)
            a2 = (a1 + d + 1) % N
            s2 = eng.shift(vals, a2, lambda t, d=d: d + v * (d + 1 + t))
            G = vals + ext
            S = s1 + s2
            vals = [G[x] * S[x] % N for x in range(2 * d + 1)]
            m.inc("rn_mul", 2 * d + 1)
            d *= 2
            if bit == "1":
                # ---- increment d -> d+1
                vals = [vals[x] * ((A + v * x + d) % N) % N for x in range(d + 1)]
                last = 1
                base = A + v * (d + 1)
                for j in range(d + 1):
                    last = last * ((base + j) % N) % N
                vals.append(last)
                m.inc("rn_mul", 2 * (d + 1))
                d += 1
    except NonUnit as e:
        w = e.value
        return _result("bgs", N, A, B, None, gcd(w, N), w, True, m)
    assert d == v
    blocks = vals[:v]                               # block x covers [A+vx, A+v(x+1))
    tail = list(range(A + v * v, B))
    allv = blocks + tail
    value = 1
    for t in allv:
        value = value * t % N
    m.inc("rn_mul", len(allv))
    m.inc("gcd")
    divisor = witness = None
    if gcd(value, N) > 1:
        i = find_nonunit(allv, N, m)
        if i < v:
            witness, divisor = _localize_in_range(N, A + i * v, A + (i + 1) * v, m)
        else:
            witness = allv[i]
            divisor = gcd(witness, N)
    return _result("bgs", N, A, B, value, divisor, witness, False, m, block=v, backend=backend)


METHODS = {"naive": naive_interval, "ps": ps_interval, "bgs": bgs_interval}


def interval_predicate(N: int, A: int, B: int, method: str = "bgs", m: Optional[Metrics] = None,
                       backend: str = "flint") -> Dict:
    """I(A,B) with fall-back: if BGS hits an auxiliary nonunit, finish with PS.

    The aux hit is still reported (it is a factor discovery, classified by the
    harness; the integers involved are < 3(B-A) + 4).
    """
    m = m if m is not None else Metrics()
    if method == "naive":
        return naive_interval(N, A, B, m)
    fn = METHODS[method]
    r = fn(N, A, B, m, backend=backend)
    if r["aux_hit"]:
        aux = (r["divisor"], r["witness"])
        r = ps_interval(N, A, B, m, backend=backend)
        r["method"] = method + "+ps_fallback"
        r["aux_hit"] = True
        r["aux_divisor"], r["aux_witness"] = aux
    return r


def U_predicate(N: int, K: int, method: str = "bgs", m: Optional[Metrics] = None,
                backend: str = "flint") -> Dict:
    """U(K) = [gcd((K-1)!, N) = 1]  (True = unit)."""
    r = interval_predicate(N, 1, K, method, m, backend)
    r["U"] = not r["nonunit"]
    return r


def least_factor_by_search(N: int, method: str = "bgs", m: Optional[Metrics] = None,
                           backend: str = "flint") -> Dict:
    """Binary search on U(K) for K in [2, isqrt(N)+1]; returns least prime factor
    (or None if N is prime/1).  Every call's counters accumulate in m.
    The search uses only the predicate's boolean output (no gcd shortcut)."""
    m = m if m is not None else Metrics()
    lo, hi = 2, isqrt(N) + 1           # U(lo) true; find largest K with U(K) true
    calls = 0
    if U_predicate(N, hi + 1, method, m, backend)["U"]:
        calls += 1
        return {"p": None, "calls": calls, **m.to_record()}
    calls += 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        calls += 1
        if U_predicate(N, mid, method, m, backend)["U"]:
            lo = mid
        else:
            hi = mid - 1
    return {"p": lo, "calls": calls, **m.to_record()}
