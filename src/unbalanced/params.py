"""PUBLIC parameter designs for DCL-wide (ratio-band DCL).

Everything here is a function of the bit length n and public knobs only.
Nothing depends on N or its factors.

An *instance* is a weight (A, B) together with a ratio interval
(rho_lo, rho_hi] (rationals, given as integer pairs).  For a given N the
algorithm computes the window of g = (B p + A q - c)/l^2 over the integer
p-range {p : rho_lo < N/p^2 <= rho_hi}, i.e. over the p's whose cofactor
ratio q/p lies in that interval.  Theorem U (notes/unbalanced_results.md):

  Regime I  (ratio band j, R_j = 2^j <= H^2/C0):
      instances (A, B, 2^j, 2^{j+1}) for coprime A != B with
      A <= a_j = floor(H 2^{-j/2}),  B <= b_j = floor(H 2^{j/2});
      every window span is <= (1+sqrt2) H sqrt N + O(1) -- the SAME bound as D',
      so one prime band set, one gamma and one shared baby list serve all bands.
  Regime II (ratio band j with 2^j > H^2/C0; Theorem U_theta):
      instances (1, B, B - Delta, B + Delta) for B in [R_j, 2 R_j),
      with a band-specific Y_j ~ Delta^{2/3} N^{1/6} R_j^{-5/6}.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

from src.dcltheorem.params import (pow2_frac, primes_upto, schedule_design, icbrt,
                                   span_bound)


@dataclass(frozen=True)
class Inst:
    A: int
    B: int
    rlo_num: int      # rho_lo = rlo_num / rlo_den   (exclusive)
    rlo_den: int
    rhi_num: int      # rho_hi = rhi_num / rhi_den   (inclusive)
    rhi_den: int
    band: int         # ratio band label j (bookkeeping only)


@dataclass
class WideDesign:
    n: int
    Y: int
    K: int
    mod4: bool
    insts: List[Inst]
    label: str = ""
    H: int = 0
    Delta: int = 0
    bands: List[List[int]] = field(default_factory=list)

    def __post_init__(self):
        top = (1 << self.K) * self.Y
        ps = primes_upto(top)
        self.bands = []
        for k in range(self.K):
            lo, hi = (1 << k) * self.Y, (1 << (k + 1)) * self.Y
            self.bands.append([l for l in ps if lo < l <= hi and l > 2
                               and (not self.mod4 or l % 4 == 3)])

    def summary(self) -> Dict:
        jb = {}
        for t in self.insts:
            jb[t.band] = jb.get(t.band, 0) + 1
        return {"n": self.n, "H": self.H, "Delta": self.Delta, "Y": self.Y, "K": self.K,
                "mod4": self.mod4, "label": self.label, "n_insts": len(self.insts),
                "insts_per_band": jb,
                "band_sizes": [len(b) for b in self.bands],
                "band_ranges": [[b[0], b[-1]] if b else None for b in self.bands]}


# ---------------------------------------------------------------------------
# Regime I
# ---------------------------------------------------------------------------
def regime1_weights(H: int, j: int) -> List[Tuple[int, int]]:
    """V_j: coprime A != B, A <= floor(H 2^{-j/2}), B <= floor(H 2^{j/2})."""
    a = isqrt((H * H) >> j) if j >= 0 else H
    b = isqrt((H * H) << j)
    return [(A, B) for A in range(1, a + 1) for B in range(1, b + 1)
            if A != B and gcd(A, B) == 1]


def regime1_insts(H: int, js: List[int]) -> List[Inst]:
    out = []
    for j in js:
        for (A, B) in regime1_weights(H, j):
            out.append(Inst(A, B, 1 << j, 1, 1 << (j + 1), 1, j))
    return out


def regime1_Jmax(H: int, C0: int = 1) -> int:
    """Largest j with 2^j <= H^2 / C0 (a_j >= 1 when C0 = 1)."""
    j = 0
    while (1 << (j + 1)) * C0 <= H * H:
        j += 1
    return j


def theorem_design_U(n: int, eps_num: int, eps_den: int, C0: int = 16) -> WideDesign:
    """Theorem U schedule: H, Y, K = 4 as in D' (dcltheorem.params.theorem_design),
    ratio bands j = 0..J with 2^J <= H^2/C0."""
    H = n * pow2_frac(n * eps_num, 1, 2 * eps_den)
    Y = max(3, -(-pow2_frac(n, 1, 6) // max(1, icbrt(H))))
    J = regime1_Jmax(H, C0)
    return WideDesign(n=n, Y=Y, K=4, mod4=True, insts=regime1_insts(H, list(range(J + 1))),
                      label=f"thmU eps={eps_num}/{eps_den}", H=H)


def wide_from_C1(n: int, js: List[int], lam: float = 3.0, K: int = 4) -> WideDesign:
    """Finite-size regime-I design: H and Y from the pre-registered D' schedule C1
    (dcltheorem.params.schedule_design(n, lam)), plus ratio bands js."""
    D = schedule_design(n, lam, K=K, mod4=True)
    return WideDesign(n=n, Y=D.Y, K=K, mod4=True, insts=regime1_insts(D.H, js),
                      label=f"W1 C1-H js={js[0]}..{js[-1]}", H=D.H)


# ---------------------------------------------------------------------------
# Regime II
# ---------------------------------------------------------------------------
def regime2_insts(j: int, Delta: int) -> List[Inst]:
    """Weights (1, B) for B in [2^j, 2^{j+1}) (each B belongs to exactly one dyadic
    band), each with ratio window (B - Delta, B + Delta].  A pair with ratio R in
    band j has >= Delta + 1 such B with |B - R| <= Delta (2 Delta + 1 away from the
    band edges)."""
    R = 1 << j
    return [Inst(1, B, max(1, B - Delta), 1, B + Delta, 1, j) for B in range(R, 2 * R)]


def regime2_Y(n: int, j: int, Delta: int) -> int:
    """Y_j = Delta^{2/3} 2^{n/6} 2^{-5j/6} (balances span/Y^2 against R_j Y)."""
    # (Delta^2 * 2^{n/2} / 2^{5j/2})^{1/3}, integer-only
    x = Delta * Delta * pow2_frac(n, 1, 2)
    y = pow2_frac(5 * j, 1, 2) if j else 1
    return max(3, icbrt(x // max(1, y)))


def regime2_design(n: int, j: int, Delta: int, K: int = 4, Y: Optional[int] = None,
                   label: str = "") -> WideDesign:
    Yj = Y or regime2_Y(n, j, Delta)
    return WideDesign(n=n, Y=Yj, K=K, mod4=True, insts=regime2_insts(j, Delta),
                      label=label or f"W2 j={j} Delta={Delta}", Delta=Delta)


def regime2_span_bound(n: int, j: int, Delta: int) -> int:
    """sqrt(N) Delta^2 / (2 R_j^{3/2}) + 2 for N < 2^n, R_j = 2^j >= 2 Delta (display)."""
    return (pow2_frac(n, 1, 2) + 1) * Delta * Delta // max(1, 2 * pow2_frac(3 * j, 1, 2)) + 2


# ---------------------------------------------------------------------------
# Finite-size regime-II schedule (pre-registered; notes/unbalanced_prereg.md)
# ---------------------------------------------------------------------------
def _band_primes(Y: int) -> List[int]:
    return [l for l in primes_upto(2 * Y) if l > Y and l % 4 == 3]


def _expected_cost_r2(n: int, j: int, Delta: int, ps: List[int]) -> int:
    """n-only estimate of the block-0 list cost m + sum ceil(W/m), one expected
    target per (instance, l), windows from N_ref = 3 * 2^(n-2)."""
    Nref = 3 << (n - 2)
    widths = []
    for t in regime2_insts(j, Delta):
        x_lo, x_hi = _rr(Nref, t)
        if x_lo > x_hi:
            continue
        lo = isqrt(4 * t.B * Nref)
        hi = max(t.B * x + -(-Nref // x) for x in (x_lo, x_hi))
        for l in ps:
            widths.append((hi - lo) // (l * l) + 2)
    if not widths:
        return 1 << 62
    best, m = None, 4
    while m <= 4 * max(widths):
        c = m + sum(-(-w // m) for w in widths)
        best = c if best is None or c < best else best
        m = (119 * m) // 100 + 1
    return best


def _rr(N: int, t: Inst) -> Tuple[int, int]:
    a, b = t.rhi_num, t.rhi_den
    x_lo = isqrt((N * b) // a)
    while a * x_lo * x_lo < N * b:
        x_lo += 1
    a2, b2 = t.rlo_num, t.rlo_den
    x_hi = isqrt((N * b2) // a2) + 1
    while a2 * x_hi * x_hi >= N * b2:
        x_hi -= 1
    return x_lo, x_hi


def schedule_regime2(n: int, j: int, lam_target: float = 3.0, K: int = 4) -> WideDesign:
    """For each Y on a 2^(1/4) grid around regime2_Y: Delta = least integer with
    lam0 = (2 Delta + 1) * sum_{l in P_0} 1/l >= lam_target (the per-pair
    multiplicity of block 0: #{B : |B - R| <= Delta} ~ 2 Delta + 1); keep the
    (Y, Delta) minimising the n-only expected block-0 cost; require Y > 2 Delta + 1
    (distinct classes, Lemma 9).  Display floats are used only to centre the grid."""
    best = None
    y0 = None
    for jj in range(-10, 11):
        Delta0 = 20
        if y0 is None:
            y0 = regime2_Y(n, j, Delta0)
        Y = max(7, int(y0 * 2 ** (jj / 4)))
        ps = _band_primes(Y)
        if not ps:
            continue
        s = sum(1.0 / l for l in ps)
        Delta = 1
        while (2 * Delta + 1) * s < lam_target:
            Delta += 1
        if Y <= 2 * Delta + 1 or (1 << j) < 2 * Delta:
            continue
        c = _expected_cost_r2(n, j, Delta, ps)
        if best is None or c < best[0]:
            best = (c, Y, Delta)
    if best is None:
        raise ValueError(f"regime II infeasible at n={n}, j={j}")
    _, Y, Delta = best
    return regime2_design(n, j, Delta, K=K, Y=Y, label=f"W2 j={j}")
