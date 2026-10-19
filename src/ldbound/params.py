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


def h_floor_mu(d: int, k: int, mu: Fraction) -> Fraction:
    """Height floor with M = N^mu, w ~ N^(mu + nu) (codex ldbound_2 variant):
    Kern 2nu; FM level mu - (d-2)nu; window |a_0| ~ N w / M^d << M H:  1 - d mu + nu."""
    nu = Fraction(1, 4 * k)
    return max(2 * nu, mu - (d - 2) * nu, 1 - d * mu + nu)


def theta_mu(d: int, k: int, mu: Fraction) -> Fraction:
    nu = Fraction(1, 4 * k)
    return 2 * mu + h_floor_mu(d, k, mu) + (d - 1) * nu


def best_mu(d: int, k: int):
    """Minimise theta_mu over mu > nu (piecewise linear in mu: check the breakpoints)."""
    nu = Fraction(1, 4 * k)
    cands = {Fraction(1, d + 1), (1 + (d - 1) * nu) / (d + 1), (1 - nu) / d, d * nu}
    cands = [m for m in cands if m > nu]
    m = min(cands, key=lambda m: (theta_mu(d, k, m), m))
    return m, theta_mu(d, k, m)


def theta3opt(d: int, k: int) -> Fraction:
    return best_mu(d, k)[1]


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
        th = theta_pp(d, k, "pair") if design == "a2" else {"a3": theta3, "a3opt": theta3opt}[design](d, k)
        if th < 1 and c_const(th) < target:
            return k
    raise ValueError("no k found")


def margins(d: int, k: int, mu: Fraction = None) -> Dict[str, Fraction]:
    """Exponent (of N) of every polynomial margin used in the LD* proof, for M = N^mu (default
    mu = 1/(d+1)) and h = h_floor_mu; all must be > 0 (each then absorbs any L^{O(1)} factor;
    the eta0-slack is added on top)."""
    nu = Fraction(1, 4 * k)
    mu = Fraction(1, d + 1) if mu is None else mu
    h = h_floor_mu(d, k, mu)
    return {
        "FM_k": mu - nu,                              # FM: q = mu + nu > 2 nu; fibres M/V >> 1
        "FM_level": (mu + h + d * nu) / 2 - (mu + nu) + Fraction(1, 10 ** 9),  # >= 0 at floor
        "Kern": h - 2 * nu + Fraction(1, 10 ** 9),    # >= 0 at floor (eta0 added)
        "window": h - (1 - d * mu + nu) + Fraction(1, 10 ** 9),  # |a_0| V^d << M H V^d
        "diag": h + (d - 1) * nu - 2 * nu,            # T >> l q
        "C_R2": h + (d - 1) * nu - 4 * nu,            # class (C): Y^2 q^2 << T rho rho
        "C_X": 2 * nu,                                # class (C): pi_J >> 1/rho_F
        "B_R2": 2 * nu + h + (d - 1) * nu,            # class (B): X/X' = l T
        "theta_le_1": 1 - theta_mu(d, k, mu) + Fraction(1, 10 ** 9),
    }


def table(ds=(7, 9, 11, 13, 19), ks=(10, 12, 20, 30, 50, 100, 1000)) -> Dict:
    out = {}
    for d in ds:
        rows = {}
        for k in ks:
            if 4 * k <= d + 1:
                continue
            t3, t2 = theta3(d, k), theta_pp(d, k, "pair")
            mu, t3o = best_mu(d, k)
            rows[k] = dict(theta_a3=float(t3), c_a3=c_const(t3), mu_opt=str(mu),
                           theta_a3opt=float(t3o), c_a3opt=c_const(t3o),
                           theta_a2=float(t2), c_a2=c_const(t2))
        out[d] = dict(rows=rows, c_limit=c_const(Fraction(3, d + 1)),
                      min_k_a3opt=min_k(d, design="a3opt") if c_const(Fraction(3, d + 1)) < 1 else None,
                      min_k_a3=min_k(d) if c_const(Fraction(3, d + 1)) < 1 else None,
                      min_k_a2=min_k(d, design="a2") if c_const(Fraction(3, d + 1)) < 1 else None)
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(table(), indent=1))
