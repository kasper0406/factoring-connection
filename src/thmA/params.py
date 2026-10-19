"""Exact exponent bookkeeping for Theorem A in design (a'''_k; mu)  (notes/thmA_proof.md §1, §5).

Independent re-implementation (agent thmA); it does NOT import src/ldbound/params.py, so the two can be
compared in tests.  All quantities are exponents of N, as Fractions.

  s = 1/(d+1), nu = 1/(4k), special-q prime l = N^(2 nu), box V = N^nu * S, S = L^(O(1)),
  M = N^mu, w ~ N^(mu + nu + eta0/2), H = N^h,
  h >= h_*(mu) + eta0,  h_*(mu) = max(2 nu [Kern], mu - (d-2) nu [FM level 1/2], 1 - d mu + nu [window]),
  theta = (mu - nu) [rational cofactor G/l] + (mu + h + d nu) [algebraic value] = 2 mu + h + (d-1) nu.
  Cost exponent (Lenstra-Pomerance shape test, B = L[1/2, a]):
     c(a) = max(2a, a + 1/sqrt2, a + theta/(2a));  min_a c(a) = (1 + theta)/sqrt2 at a = theta/sqrt2
  provided theta <= 1.   c < 1  <=>  (1 + theta)^2 < 2  (exact test on rationals).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple


def h_star(d: int, k: int, mu: Fr) -> Fr:
    nu = Fr(1, 4 * k)
    return max(2 * nu, mu - (d - 2) * nu, 1 - d * mu + nu)


def theta(d: int, k: int, mu: Fr, eta0: Fr = Fr(0)) -> Fr:
    """theta including slack: h = h_*(mu) + eta0."""
    nu = Fr(1, 4 * k)
    return 2 * mu + h_star(d, k, mu) + eta0 + (d - 1) * nu


def admissible(d: int, k: int, mu: Fr, eta0: Fr = Fr(0)) -> Dict[str, bool]:
    """Every side condition used somewhere in the chain (zero-slack versions are checked with >=/<=)."""
    nu = Fr(1, 4 * k)
    h = h_star(d, k, mu) + eta0
    th = theta(d, k, mu, eta0)
    return {
        "d_odd_ge3": d >= 3 and d % 2 == 1,
        "k_gt_(d+1)/4": Fr(k) > Fr(d + 1, 4),          # FM: q0 = mu+nu > 2 nu needs mu > nu; pcest form
        "mu_ge_nu+eta0": mu >= nu + eta0,               # FM step 1 (M/V -> inf), Lemma 8.5' Q/l >= x^eps
        "mu_lt_1": mu < 1,
        "theta_le_1": th <= 1,                          # cost formula, |f_d| < N, 2H+1 < N, w < N
        "mu+h_lt_1": mu + h < 1,                        # ||f|| << M H < N
        "w_lt_N": mu + nu + eta0 / 2 < 1,               # w < N (gcd(w, N) is a proper-factor test)
        "h_gt_(5-d)nu": h > (5 - d) * nu,               # L2_fam class (C) second term
        "kern": h >= 2 * nu,
        "level": h >= mu - (d - 2) * nu,
        "window": h >= 1 - d * mu + nu,
    }


def best_mu(d: int, k: int) -> Tuple[Fr, Fr]:
    """Exact minimiser of theta(mu) (zero slack) over the admissible interval [nu, 1).
    theta is convex piecewise linear in mu, so the minimum is at a breakpoint or an endpoint."""
    nu = Fr(1, 4 * k)
    cands = {nu, d * nu, (1 - nu) / d, (1 + (d - 1) * nu) / (d + 1), Fr(1, d + 1)}
    cands = [m for m in cands if nu <= m < 1]
    best = min(cands, key=lambda m: (theta(d, k, m), m))
    return best, theta(d, k, best)


def closed_form(d: int, k: int) -> Tuple[Fr, Fr]:
    """codex ldbound_3 closed form, re-derived:
    k <  (d^2+1)/4: mu = (1-nu)/d,           theta = 2/d + ((d^2+d-2)/d) nu   (Kern + window bind)
    k >= (d^2+1)/4: mu = (1+(d-1)nu)/(d+1),  theta = 3/(d+1) + ((4d-2)/(d+1)) nu (level + window bind)"""
    nu = Fr(1, 4 * k)
    if Fr(k) < Fr(d * d + 1, 4):
        return (1 - nu) / d, Fr(2, d) + Fr(d * d + d - 2, d) * nu
    return (1 + (d - 1) * nu) / (d + 1), Fr(3, d + 1) + Fr(4 * d - 2, d + 1) * nu


def c_lt_1(th: Fr) -> bool:
    """(1 + theta)/sqrt2 < 1 exactly."""
    return (1 + th) ** 2 < 2


def eta0_threshold_sq(th0: Fr) -> Tuple[Fr, Fr]:
    """c(eta0) = (1 + th0 + eta0)/sqrt2 < 1  <=>  eta0 < sqrt2 - 1 - th0.  Return (th0, (1+th0)^2) for
    the exact statement; numerical value via float in reports only."""
    return th0, (1 + th0) ** 2


def min_k(d: int, kmax: int = 2000) -> int | None:
    for k in range(1, kmax + 1):
        if Fr(k) <= Fr(d + 1, 4):
            continue
        mu, th = best_mu(d, k)
        if th <= 1 and c_lt_1(th) and all(v for kk, v in admissible(d, k, mu).items()
                                          if kk not in ("mu_ge_nu+eta0",)):
            return k
    return None


def cost_exponent(a: Fr, th: Fr) -> float:
    """max(2a, a + 1/sqrt2, a + theta/(2a)) as a float (reporting only)."""
    from math import sqrt
    a_ = float(a)
    return max(2 * a_, a_ + 1 / sqrt(2), a_ + float(th) / (2 * a_))


def table(ds=(7, 9, 11, 13, 19), ks=(30, 50, 100, 1000)) -> List[dict]:
    from math import sqrt
    out = []
    for d in ds:
        row = {"d": d, "min_k": min_k(d)}
        mk = row["min_k"]
        mu, th = best_mu(d, mk)
        row["theta_at_min_k"] = th
        row["c_at_min_k"] = (1 + float(th)) / sqrt(2)
        row["eta0_max_at_min_k"] = sqrt(2) - 1 - float(th)
        for k in ks:
            row[f"c_k{k}"] = (1 + float(best_mu(d, k)[1])) / sqrt(2)
        row["inf"] = (1 + 3 / (d + 1)) / sqrt(2)
        out.append(row)
    return out


if __name__ == "__main__":
    for r in table():
        print(r)
