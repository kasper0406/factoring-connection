"""HARNESS.  Fourier statistics of Farey-arc sets in Z/q and Erdos-Turan certificates.

Farey-arc set S(H', k0) = {Bbar*k mod q : 1 <= B <= H', |k| <= k0}.  For 4 H'^3 k0 < q the map
k/B -> Bbar*k mod q is a Freiman 2-isomorphism from R = {k/B} subset Q onto S (Lemma F0), so the
additive energy of S equals that of R over Q.

Erdos-Turan / phase-discarding certificate for a weight g <= 1_S:
  sum_l g(a lbar) = (J/q) ghat(0) + (1/q) sum_{h != 0} ghat(h) K(a h),  K(h) = sum_l e(h lbar/q),
so  sum_l g(a lbar) >= J*sum(g)/q - kappa * W(g),  W(g) = (1/q) sum_{h!=0} |ghat(h)|,  kappa = max_{h!=0}|K(h)|.
FFTs are floating point; they are used only to report these statistics.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Dict, List

import numpy as np


def inv_band(q: int, L: List[int]) -> np.ndarray:
    return np.array([pow(l, -1, q) for l in L], dtype=np.int64)


def kloosterman_all(q: int, L: List[int]) -> np.ndarray:
    """|K(h)| for all h mod q (FFT of the indicator of {lbar})."""
    x = np.zeros(q)
    np.add.at(x, inv_band(q, L), 1.0)
    return np.abs(np.fft.fft(x))


def farey_tri(q: int, Hp: int, k0: int, shape: str = "tri") -> np.ndarray:
    """g(x) = max over B <= Hp of max(0, 1 - |k|/(k0+1)) with x = Bbar k (tri), or 1_S (sharp)."""
    g = np.zeros(q)
    k = np.arange(-k0, k0 + 1, dtype=np.int64)
    val = 1.0 - np.abs(k) / (k0 + 1.0) if shape == "tri" else np.ones(len(k))
    for B in range(1, Hp + 1):
        x = (pow(B, -1, q) * k) % q
        np.maximum.at(g, x, val) if len(k) < 4096 else _max_scatter(g, x, val)
    return g


def farey_direct(q: int, Hp: int, k0: int, shape: str = "tri") -> np.ndarray:
    """Same function as farey_tri, computed as g(x) = max_B tri(B x mods q) (O(q) per B)."""
    x = np.arange(q, dtype=np.int64)
    g = np.zeros(q)
    for B in range(1, Hp + 1):
        k = (B * x) % q
        k = np.minimum(k, q - k)
        v = (1.0 - k / (k0 + 1.0)) if shape == "tri" else np.ones(q)
        v[k > k0] = 0.0
        np.maximum(g, v, out=g)
    return g


def farey(q: int, Hp: int, k0: int, shape: str = "tri") -> np.ndarray:
    if 2 * k0 + 1 >= q // 64:
        return farey_direct(q, Hp, min(k0, (q - 1) // 2), shape)
    return farey_tri(q, Hp, k0, shape)


def _max_scatter(g, x, val):
    order = np.argsort(x, kind="stable")
    xs, vs = x[order], val[order]
    # positions in one B are distinct when 2k0+1 <= q (k -> Bbar k is injective)
    np.maximum(g[xs], vs, out=vs)
    g[xs] = vs


def stats(g: np.ndarray) -> Dict[str, float]:
    q = len(g)
    gh = np.fft.fft(g)
    a = np.abs(gh)
    W = (a.sum() - a[0]) / q
    E = float((a ** 4).sum() / q)
    s2 = float((g * g).sum())
    LB = s2 ** 1.5 / E ** 0.5 if E > 0 else 0.0
    return {"mass": float(g.sum()), "W": float(W), "E": E, "sum_g2": s2, "LB": float(LB),
            "supp": int((g > 0).sum())}


def certificate(J: int, kappa: float, st: Dict[str, float], q: int) -> float:
    return J * st["mass"] / q - kappa * st["W"]


# ---------------------------------------------------------------------------
# exact (small) checks of Lemma F0 / F1
# ---------------------------------------------------------------------------
def farey_rationals(Hp: int, k0: int):
    return sorted({Fraction(k, B) for B in range(1, Hp + 1) for k in range(-k0, k0 + 1)})


def energy_exact(xs) -> int:
    """Additive energy #{x1+x2 = x3+x4} of a finite set (Fractions or ints)."""
    from collections import Counter
    c = Counter(a + b for a in xs for b in xs)
    return sum(v * v for v in c.values())


def energy_mod_q(q: int, Hp: int, k0: int) -> int:
    S = sorted({(pow(B, -1, q) * k) % q for B in range(1, Hp + 1) for k in range(-k0, k0 + 1)})
    from collections import Counter
    c = Counter((a + b) % q for a in S for b in S)
    return sum(v * v for v in c.values())


def lcm_sum(H: int) -> Fraction:
    """sum over B1..B4 <= H of min(B)/lcm(B1..B4) (the quantity bounded in Lemma F1)."""
    from math import lcm
    tot = Fraction(0)
    for b1 in range(1, H + 1):
        for b2 in range(1, H + 1):
            for b3 in range(1, H + 1):
                for b4 in range(1, H + 1):
                    tot += Fraction(min(b1, b2, b3, b4), lcm(b1, b2, b3, b4))
    return tot
