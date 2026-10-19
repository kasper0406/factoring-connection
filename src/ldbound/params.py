"""Exponent bookkeeping for design (a'''_k) (notes/ldbound_results.md §1, §5).

Design (a'''_k) = pcest's (a''_k) plus one extra field parameter: a uniform random prime
w in (W, 2W], W = N^(s + 1/(4k) + eta0/2), and F(g0, g1) = N w (Lemma GB applied to N w).
All sizes are exponents of N (eta0 -> 0 slack omitted unless stated):
  s = 1/(d+1);  l = N^(1/(2k));  V = N^(1/(4k));  rational modulus t = |G| = N^(s + 1/(4k));
  |a_0| = O(M w) = N^(2s + 1/(4k));  H = N^h.
Constraints on h:
  (W)   H >= w N^(eta0/2)  (FM window, |F^| << t Y_omega; X_F unchanged): h >= s + 1/(4k)
  (lvl) level 1/2 for FM (Lemma 8.5'):   h >= s - (d-2)/(4k)      [implied by (W)]
  (ker) Lemma Kern H >= V^2:            h >= 1/(2k)              [implied by (W) iff k >= (d+1)/4]
theta''' = (s - 1/(4k)) + (s + h + d/(4k)) = 2s + h + (d-1)/(4k) = 3s + d/(4k)  at h = s + 1/(4k).
Constant (Lenstra-Pomerance shape test): c = (1 + theta)/sqrt 2; infimum over k: (1 + 3/(d+1))/sqrt 2,
the same as for (a''_k) (pcest).
"""
from __future__ import annotations

from fractions import Fraction
from math import sqrt
from typing import Dict

from src.pcest.params import theta as theta_pp  # design (a''_k), for comparison


def h_floor(d: int, k: int) -> Fraction:
    s = Fraction(1, d + 1)
    return max(s + Fraction(1, 4 * k), s - Fraction(d - 2, 4 * k), Fraction(1, 2 * k))


def theta3(d: int, k: int) -> Fraction:
    s = Fraction(1, d + 1)
    return 2 * s + h_floor(d, k) + Fraction(d - 1, 4 * k)


def c_const(th) -> float:
    return (1 + float(th)) / sqrt(2)


def c_ps(th) -> float:
    """Theorem A' (Pollard-Strassen smoothness tests), closed form valid for theta <= 1/2."""
    th = float(th)
    assert th <= 0.5
    return 1 / sqrt(2) + (sqrt(2) - sqrt(2 - 4 * th)) / 2


def min_k(d: int, target: float = 1.0, design: str = "a3", kmax: int = 10 ** 5) -> int:
    for k in range(1, kmax):
        if 4 * k <= d + 1:          # FM needs k > (d+1)/4
            continue
        th = theta3(d, k) if design == "a3" else theta_pp(d, k, "pair")
        if th < 1 and c_const(th) < target:
            return k
    raise ValueError("no k found")


def margins(d: int, k: int) -> Dict[str, Fraction]:
    """Exponent (of N) of every polynomial margin used in the LD* proof; all must be > 0
    (each is then at least the fixed slack and absorbs any L^{O(1)} factor)."""
    s = Fraction(1, d + 1)
    h = h_floor(d, k)
    nu = Fraction(1, 4 * k)
    return {
        # FM needs q = s + nu > 1/(2k): k > (d+1)/4
        "FM_k": s + nu - 2 * nu,
        # Lemma 8.5' level: q < xi/2 with xi = s + h + d nu
        "FM_level": (s + h + d * nu) / 2 - (s + nu),
        # Lemma Kern: H >= V^2
        "Kern": h - 2 * nu,
        # diagonal / window: T = H V^(d-1) >> l q (moduli t Delta << X)
        "diag": h + (d - 1) * nu - 2 * nu,
        # step (C) of (***): Y^2 q^2 << T rho rho, i.e. h + (d-1) nu > 4 nu
        "C_R2": h + (d - 1) * nu - 4 * nu,
        # step (C): pi_J = N^(1/(2k)) beats 1/rho_F = L^{O(1)}
        "C_X": 2 * nu,
        # step (B): X'_+ q^2 << X rho rho: X / X' = l T
        "B_R2": 2 * nu + h + (d - 1) * nu,
        # Brun-Titchmarsh: W / t_max = N^(eta0/2) -- by construction (slack only)
    }


def table(ds=(7, 9, 11, 13, 19), ks=(10, 12, 20, 30, 50, 100, 1000)) -> Dict:
    out = {}
    for d in ds:
        rows = {}
        for k in ks:
            if 4 * k <= d + 1:
                continue
            t3, t2 = theta3(d, k), theta_pp(d, k, "pair")
            rows[k] = dict(theta_a3=float(t3), c_a3=c_const(t3),
                           theta_a2=float(t2), c_a2=c_const(t2))
        out[d] = dict(rows=rows, c_limit=c_const(Fraction(3, d + 1)),
                      min_k_a3=min_k(d) if c_const(Fraction(3, d + 1)) < 1 else None,
                      min_k_a2=min_k(d, design="a2") if c_const(Fraction(3, d + 1)) < 1 else None)
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(table(), indent=1))
