"""Public parameter schedules for the deterministic DCL algorithm (DCL-det).

MODEL / PUBLIC code: everything here is a function of the bit length n (and
of the public knobs eps / lam_target).  Nothing depends on N's factors, and
-- important for the sieve proof (notes/dcltheorem_results.md, Lemma 9) --
nothing depends on N itself: the prime bands and the weight set are fixed by n.

Design D = (n, H, Y, K blocks):
  V      = {(A, B) : 1 <= A, B <= H, gcd(A, B) = 1, A != B}
  P_k    = {l prime : l == 3 (mod 4), 2^k Y < l <= 2^(k+1) Y},  k = 0..K-1
  mbar_k = an n-only upper bound for every whole-range window width in block k
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil, gcd, isqrt, log
from typing import Dict, List, Optional, Tuple


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def weight_set(H: int) -> List[Tuple[int, int]]:
    """V(H): positive coprime pairs, A != B (the symmetric pair (1,1) is two-sided
    by construction -- see Lemma 4 -- and is excluded)."""
    return [(A, B) for A in range(1, H + 1) for B in range(1, H + 1)
            if A != B and gcd(A, B) == 1]


def iroot2_floor(x: int) -> int:
    return isqrt(x)


def pow2_frac(n: int, num: int, den: int) -> int:
    """floor(2^(n*num/den)) exactly (integer root of a power of two)."""
    e_num = n * num
    # 2^(e_num/den) = (2^e_num)^(1/den)
    x = 1 << e_num
    lo, hi = 1, 1 << (e_num // den + 2)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** den <= x:
            lo = mid
        else:
            hi = mid - 1
    return lo


@dataclass
class Design:
    n: int
    H: int
    Y: int
    K: int = 4
    mod4: bool = True                      # primes l == 3 (mod 4) only
    label: str = ""
    V: List[Tuple[int, int]] = field(default_factory=list)
    bands: List[List[int]] = field(default_factory=list)
    mbar: List[int] = field(default_factory=list)

    def __post_init__(self):
        self.V = weight_set(self.H)
        top = (1 << self.K) * self.Y
        ps = primes_upto(top)
        self.bands = []
        for k in range(self.K):
            lo, hi = (1 << k) * self.Y, (1 << (k + 1)) * self.Y
            band = [l for l in ps if lo < l <= hi and l > 2 and (not self.mod4 or l % 4 == 3)]
            self.bands.append(band)
        self.mbar = [span_bound(self.n, self.H) // (b[0] ** 2) + 2 if b else 0
                     for b in self.bands]

    def lam(self, k: int = 0) -> float:
        """Expected sieve multiplicity of block k: |V| * sum_{l in P_k} 1/l (display)."""
        return len(self.V) * sum(1.0 / l for l in self.bands[k])

    def summary(self) -> Dict:
        return {"n": self.n, "H": self.H, "Y": self.Y, "K": self.K, "mod4": self.mod4,
                "label": self.label, "nV": len(self.V),
                "band_sizes": [len(b) for b in self.bands],
                "band_ranges": [[b[0], b[-1]] if b else None for b in self.bands],
                "lam": [round(self.lam(k), 4) for k in range(self.K)],
                "mbar": self.mbar}


def span_bound(n: int, H: int) -> int:
    """n-only upper bound for max_v (Hi_v - Lo_v) over N < 2^n:
    Hi_v <= (B + sqrt2 A) sqrt N + 2H + 2 and Lo_v >= 0, so
    span <= (1 + sqrt 2) H 2^(n/2) + 2H + 3  <=  (5/2) H 2^(n/2) + 2H + 3."""
    return (5 * H * (pow2_frac(n, 1, 2) + 1)) // 2 + 2 * H + 3


# ---------------------------------------------------------------------------
# Theorem schedule (notes/dcltheorem_results.md, section 2)
# ---------------------------------------------------------------------------
def theorem_design(n: int, eps_num: int, eps_den: int) -> Design:
    """H = ceil(n 2^(eps n/2)), Y = ceil(2^(n/6) H^(-1/3)), K = 4, l == 3 (mod 4)."""
    H = n * pow2_frac(n * eps_num, 1, 2 * eps_den)       # ~ n 2^(eps n / 2)
    Y = max(3, -(-pow2_frac(n, 1, 6) // max(1, round(H ** (1.0 / 3)))))
    # (the cube root of H only enters through Y; any constant-factor choice of
    #  Y is allowed by the proof, so rounding here is harmless)
    return Design(n=n, H=H, Y=Y, K=4, mod4=True, label=f"theorem eps={eps_num}/{eps_den}")


# ---------------------------------------------------------------------------
# Finite-size schedule (pre-registered; notes/dcltheorem_prereg.md)
# ---------------------------------------------------------------------------
def _expected_cost(n: int, H: int, Y: int, V, band: List[int]) -> int:
    """n-only estimate of the block-0 list cost m + sum ceil(W/m) with one
    expected target per (v, l) (each (v, l) has 1 root on average), windows
    W_v(l) = span_v/l^2 + 1 with span_v computed at N_ref = 3 * 2^(n-2)."""
    Nref = 3 << (n - 2)
    P = isqrt(Nref)
    widths = []
    for (A, B) in V:
        hi = max(B * P + A * Nref // P, B * isqrt(Nref // 2) + A * Nref // isqrt(Nref // 2)) + 1
        lo = isqrt(4 * A * B * Nref)
        sp = hi - lo
        for l in band:
            widths.append(sp // (l * l) + 2)
    if not widths:
        return 1 << 62
    best = None
    m = 4
    while m <= 4 * max(widths):
        c = m + sum(-(-w // m) for w in widths)
        if best is None or c < best:
            best = c
        m = int(m * 1.19) + 1
    return best


def schedule_design(n: int, lam_target: float, K: int = 4, mod4: bool = True,
                    label: str = "") -> Design:
    """Smallest H whose block-0 multiplicity |V| sum 1/l reaches lam_target, with
    Y minimising the n-only expected block-0 cost on a 2^(1/4) grid."""
    H = 2
    while True:
        V = weight_set(H)
        y0 = pow2_frac(n, 1, 6) / max(1.0, H ** (1.0 / 3))
        best = None
        for j in range(-12, 13):
            Y = max(H + 1, int(y0 * 2 ** (j / 4)))
            ps = [l for l in primes_upto(2 * Y) if l > Y and (not mod4 or l % 4 == 3)]
            if not ps:
                continue
            c = _expected_cost(n, H, Y, V, ps)
            if best is None or c < best[0]:
                best = (c, Y, ps)
        c, Y, ps = best
        lam = len(V) * sum(1.0 / l for l in ps)
        if lam >= lam_target:
            return Design(n=n, H=H, Y=Y, K=K, mod4=mod4,
                          label=label or f"schedule lam>={lam_target:.3f}")
        H += 1


def lam_full(n: int) -> float:
    """Full-coverage multiplicity target ln(#prime candidates) + 2 = ln(P / ln P) + 2."""
    P = 2.0 ** (n / 2)
    return log(P / log(P)) + 2.0
