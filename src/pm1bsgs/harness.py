"""DIAGNOSTIC / HARNESS code for pm1bsgs.  KNOWS THE FACTORS.  Never imported by
operational.py or costmodel.py.

Provides, for an instance N = p q with known p, q:
  * factorisation of p-1, q-1 (sympy, harness only), B-smooth parts S_B(p-1);
  * exact element orders o_p = ord_p(a), and d_p = ord_p(a^E);
  * exact outcome predictors for the operational methods (a)-(d), mirroring the
    gcd / backtracking logic of operational.py exactly (validated in tests and
    experiments/raw/pm1bsgs_validate.jsonl);
  * exact failure probabilities over a uniformly random base.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import numpy as np
from sympy import factorint

from src.common.instances import is_prime
from .operational import birthday_schedule, choose_base


# ------------------------------------------------------------------ basics
def e_l(l: int, B1: int) -> int:
    """floor(log_l B1) computed exactly (0 if l > B1)."""
    e, pk = 0, l
    while pk <= B1:
        e += 1
        pk *= l
    return e


@dataclass
class PrimeInfo:
    r: int
    fac: Dict[int, int]       # factorisation of r-1

    @classmethod
    def of(cls, r: int) -> "PrimeInfo":
        return cls(r, dict(factorint(r - 1)))

    def smooth_part(self, B: int) -> int:
        """S_B = gcd(r-1, E_B) = prod_{l<=B} l^min(v_l(r-1), floor(log_l B))."""
        S = 1
        for l, v in self.fac.items():
            S *= l ** min(v, e_l(l, B))
        return S

    def order(self, a: int) -> int:
        a %= self.r
        if a == 0:
            return 0
        o = self.r - 1
        for l in self.fac:
            while o % l == 0 and pow(a, o // l, self.r) == 1:
                o //= l
        return o


def factor_of_order(o: int, fac: Dict[int, int]) -> Dict[int, int]:
    out = {}
    for l in fac:
        v = 0
        while o % l == 0:
            o //= l
            v += 1
        if v:
            out[l] = v
    assert o == 1
    return out


def stage1_key(o: int, ofac: Dict[int, int], B1: int):
    """Index (as a sortable tuple) of the stage-1 step after which o | (partial
    exponent), or None if o does not divide E_B1."""
    if B1 <= 1:
        return None
    for l, v in ofac.items():
        if l > B1 or v > e_l(l, B1):
            return None
    if o == 1:
        return (2, 1)
    lstar = max(ofac)
    return (lstar, ofac[lstar])


def residual_order(o: int, ofac: Dict[int, int], B1: int) -> int:
    """d = ord(a^E) = o / gcd(o, E_B1)."""
    if B1 <= 1:
        return o
    d = o
    for l, v in ofac.items():
        d //= l ** min(v, e_l(l, B1))
    return d


@dataclass
class BaseState:
    """Per-(instance, base) truth."""
    a: int
    op: int
    oq: int
    opf: Dict[int, int]
    oqf: Dict[int, int]


def base_state(N: int, P: PrimeInfo, Q: PrimeInfo, seed: int) -> BaseState:
    a = choose_base(N, seed)
    op, oq = P.order(a), Q.order(a)
    return BaseState(a, op, oq, factor_of_order(op, P.fac) if op else {},
                     factor_of_order(oq, Q.fac) if oq else {})


def stage1_outcome(bs: BaseState, B1: int):
    """("base",) | ("stage1",) | ("stage1_N",) | ("ok", d_p, d_q)."""
    if bs.op == 0 or bs.oq == 0:
        return ("base",)
    kp = stage1_key(bs.op, bs.opf, B1)
    kq = stage1_key(bs.oq, bs.oqf, B1)
    if kp is None and kq is None:
        return ("ok", residual_order(bs.op, bs.opf, B1), residual_order(bs.oq, bs.oqf, B1))
    if kp is not None and kq is not None and kp == kq:
        return ("stage1_N",)
    return ("stage1",)


# ------------------------------------------------------------------ stage 2 predictors
def _first_multiples(d: int, lo: int, hi: int, k: int = 2) -> List[int]:
    """Up to k smallest multiples of d in [lo, hi]."""
    t = -(-lo // d) * d
    out = []
    while t <= hi and len(out) < k:
        out.append(t)
        t += d
    return out


def window(mbaby: int, J: int, mode: str, N: int) -> Tuple[int, int]:
    if mode == "N":
        return N - (J - 1) * mbaby, N + mbaby - 1
    return 1, J * mbaby


def interval_pure_hit(dp: int, dq: int, lo: int, hi: int) -> bool:
    """Exists t in [lo,hi] divisible by exactly one of dp, dq?  (Checking the
    first two multiples of each suffices; proof in notes/pm1bsgs_analysis.md.)"""
    for t in _first_multiples(dp, lo, hi) + _first_multiples(dq, lo, hi):
        if (t % dp == 0) != (t % dq == 0):
            return True
    return False


def predict_bsgs(dp: int, dq: int, mbaby: int, J: int, mode: str, N: int) -> bool:
    lo, hi = window(mbaby, J, mode, N)
    return interval_pure_hit(dp, dq, lo, hi)


def predict_prime(dp: int, dq: int, B1: int, B2: int) -> bool:
    hp = B1 < dp <= B2 and is_prime(dp)
    hq = B1 < dq <= B2 and is_prime(dq)
    return (hp or hq) and dp != dq


def birthday_first_pure(N: int, seed: int, dp: int, dq: int, Kmax: int) -> Optional[int]:
    """Smallest K <= Kmax such that the first K walk points contain a pair that
    collides mod exactly one of p, q (as seen through dp, dq); None if none."""
    cs = birthday_schedule(N, seed, Kmax)
    K = 1 << 10
    while True:
        K = min(K, Kmax)
        res = _bday_first(cs[:K - 1], dp, dq)
        if res is not None or K >= Kmax:
            return res
        K *= 4


def _keys(cs, d) -> np.ndarray:
    """Exponents e_1 = 1, e_{j+1} = 2 e_j + c_j, reduced mod d."""
    e = 1 % d
    out = [e]
    app = out.append
    for c in cs:
        e = (2 * e + c) % d
        app(e)
    return np.array(out, dtype=np.int64)


def _first_split(ka: np.ndarray, kb: np.ndarray) -> Optional[int]:
    """Smallest index j such that some i < j has ka[i]==ka[j] and kb[i]!=kb[j]."""
    n = len(ka)
    idx = np.arange(n)
    order = np.lexsort((idx, ka))
    ka_s, kb_s, ix_s = ka[order], kb[order], idx[order]
    newgrp = np.ones(n, dtype=bool)
    newgrp[1:] = ka_s[1:] != ka_s[:-1]
    grp_start = np.maximum.accumulate(np.where(newgrp, np.arange(n), 0))
    # within a group (sorted by index) the first split is the first element whose
    # kb differs from the group's first element's kb
    differs = kb_s != kb_s[grp_start]
    if not differs.any():
        return None
    return int(ix_s[differs].min())


def _bday_first(cs, dp, dq) -> Optional[int]:
    kp = _keys(cs, dp)
    kq = _keys(cs, dq)
    a = _first_split(kp, kq)
    b = _first_split(kq, kp)
    c = [x for x in (a, b) if x is not None]
    return min(c) + 1 if c else None


# ------------------------------------------------------------------ outcome labels
def outcome_from_orders(op: int, opf: Dict[int, int], oq: int, oqf: Dict[int, int], N: int,
                        method: str, params: Dict) -> str:
    """Exact outcome of an operational run whose base has ord_p(a)=op, ord_q(a)=oq.
    Labels: "stage1" (factor in stage 1), "stage1_N" (stage-1 gcd = N, failure),
    "success2" (factor in stage 2), "common2" (stage-2 hits exist but every
    hit is shared by p and q -> only gcd = N, failure), "nohit" (no hit)."""
    bs = BaseState(0, op, oq, opf, oqf)
    B1 = params.get("B1", 1)
    st = stage1_outcome(bs, B1)
    if st[0] != "ok":
        return st[0]
    dp, dq = st[1], st[2]
    if method == "bsgs":
        lo, hi = window(params["m"], params["J"], params["mode"], N)
        if interval_pure_hit(dp, dq, lo, hi):
            return "success2"
        if _first_multiples(dp, lo, hi, 1) or _first_multiples(dq, lo, hi, 1):
            return "common2"
        return "nohit"
    if method == "prime":
        B2 = params["B2"]
        if predict_prime(dp, dq, B1, B2):
            return "success2"
        if dp == dq and B1 < dp <= B2 and is_prime(dp):
            return "common2"
        return "nohit"
    raise ValueError(method)


# ------------------------------------------------------------------ exact probabilities
def divisors_with_phi(fac: Dict[int, int]) -> List[Tuple[int, int, Dict[int, int]]]:
    """All (d, phi(d), factorisation of d) for d | prod l^v."""
    out = [(1, 1, {})]
    for l, v in fac.items():
        new = []
        for d, ph, f in out:
            pk = 1
            for k in range(v + 1):
                f2 = dict(f)
                if k:
                    f2[l] = k
                new.append((d * pk, ph * (pk - pk // l if k else 1), f2))
                pk *= l
        out = new
    return out


def exact_outcome_probs(P: PrimeInfo, Q: PrimeInfo, N: int, method: str, params: Dict,
                        max_pairs: int = 2_000_000) -> Optional[Dict]:
    """Exact outcome distribution over a base a uniform in (Z/N)^*:
    ord_p(a) = o with probability phi(o)/(p-1), independently for q (CRT).
    Returns {label: Fraction} or None if too many divisor pairs."""
    from fractions import Fraction
    Dp, Dq = divisors_with_phi(P.fac), divisors_with_phi(Q.fac)
    if len(Dp) * len(Dq) > max_pairs:
        return None
    acc: Dict[str, int] = {}
    for op, php, opf in Dp:
        for oq, phq, oqf in Dq:
            lab = outcome_from_orders(op, opf, oq, oqf, N, method, params)
            acc[lab] = acc.get(lab, 0) + php * phq
    tot = (P.r - 1) * (Q.r - 1)
    return {k: Fraction(v, tot) for k, v in acc.items()}


def user_first_x_failure_prob(P: PrimeInfo, Q: PrimeInfo, N: int, B1: int):
    """The user's single test: x* = N mod d_p (least x >= 0 with d_p | N - x) and
    gcd(b^N - b^x*, N).  Exact P(q also divides, i.e. d_q | N - x*), for b uniform
    in the subgroup of E-th powers mod p and mod q (independent), CONDITIONED on
    d_p > 1 and d_q > 1 (i.e. stage 1 did not already reach either prime).
    Returns None if R_p = 1 or R_q = 1."""
    from fractions import Fraction
    Rp, Rq = residual_group_fac(P, B1), residual_group_fac(Q, B1)
    Rpv, Rqv = _val(Rp), _val(Rq)
    if Rpv == 1 or Rqv == 1:
        return None
    Dp = [x for x in divisors_with_phi(Rp) if x[0] > 1]
    Dq = [x for x in divisors_with_phi(Rq) if x[0] > 1]
    bad = 0
    for dp, php, _ in Dp:
        xs = N % dp
        for dq, phq, _ in Dq:
            if (N - xs) % dq == 0:
                bad += php * phq
    return Fraction(bad, (Rpv - 1) * (Rqv - 1))


def _val(fac: Dict[int, int]) -> int:
    v = 1
    for l, e in fac.items():
        v *= l ** e
    return v


def residual_group_fac(P: PrimeInfo, B1: int) -> Dict[int, int]:
    """Factorisation of R = (r-1)/S_B1 (b = a^E is uniform on the subgroup of order R)."""
    out = {}
    for l, v in P.fac.items():
        w = v - (min(v, e_l(l, B1)) if B1 > 1 else 0)
        if w:
            out[l] = w
    return out


def prob_equal_orders(P: PrimeInfo, Q: PrimeInfo, B1: int, exclude_one: bool = False):
    """P(ord_p(b) = ord_q(b)) = sum_{d | gcd(R_p,R_q)} phi(d)^2 / (R_p R_q);
    with exclude_one, conditioned on both orders > 1 (None if R_p or R_q = 1)."""
    from fractions import Fraction
    Rp, Rq = residual_group_fac(P, B1), residual_group_fac(Q, B1)
    g = {l: min(v, Rq[l]) for l, v in Rp.items() if l in Rq}
    if not exclude_one:
        return Fraction(sum(ph * ph for _, ph, _ in divisors_with_phi(g)), _val(Rp) * _val(Rq))
    if _val(Rp) == 1 or _val(Rq) == 1:
        return None
    return Fraction(sum(ph * ph for d, ph, _ in divisors_with_phi(g) if d > 1), (_val(Rp) - 1) * (_val(Rq) - 1))
