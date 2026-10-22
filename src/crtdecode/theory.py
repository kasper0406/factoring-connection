"""Cost and regime formulas for crtdecode (MODEL; floats only in model evaluation, never in decoders).

- c_of_a(a): factoring constant when the decoder costs L[1/3, a*kappa] per window (Prop RD's softer target):
  min_kappa max((1 + a) kappa + 2/(9 kappa^2), 2 kappa) = 3 (1/18)^{1/3} (1 + a)^{2/3} for a >= 1/3.
- a_star(): the a at which c_of_a meets MNFS ((92 + 26 sqrt 13)/27)^{1/3}.
- dictionary(bits): the CRT-code parameters of Prop SI at kappa*, nu* (o(1) dropped): code length n = pi(y),
  "dimension" k = ln h / ln y, required agreement t = u0 = x / ln y, Johnson agreement sqrt(k n), the q-ary
  capacity rate at agreement fraction t/n, and the average list size (h + 1) rho(u0) (rho by leading order).
"""
from __future__ import annotations

import math
from typing import Dict

MNFS = ((92 + 26 * math.sqrt(13)) / 27) ** (1.0 / 3.0)
KAPPA_STAR = (2.0 / 9.0) ** (1.0 / 3.0)
NU_STAR = 2.0 / (3.0 * KAPPA_STAR)


def c_of_a(a: float) -> float:
    if a < 1.0 / 3.0:
        raise ValueError("closed form valid for a >= 1/3")
    return 3.0 * (1.0 / 18.0) ** (1.0 / 3.0) * (1.0 + a) ** (2.0 / 3.0)


def c_of_a_grid(a: float, grid: int = 200000) -> float:
    best = float("inf")
    for i in range(1, grid):
        k = 0.2 + 1.8 * i / grid
        best = min(best, max((1 + a) * k + 2.0 / (9.0 * k * k), 2.0 * k))
    return best


def a_star() -> float:
    return (MNFS / (3.0 * (1.0 / 18.0) ** (1.0 / 3.0))) ** 1.5 - 1.0


def F(alpha: float, x: float) -> float:
    return x ** alpha * math.log(x) ** (1.0 - alpha)


def Hq(p: float, q: float) -> float:
    """q-ary entropy in base q."""
    if p <= 0:
        return 0.0
    if p >= 1:
        return math.log(q - 1) / math.log(q)
    return (p * math.log(q - 1) - p * math.log(p) - (1 - p) * math.log(1 - p)) / math.log(q)


def ln_rho(u: float) -> float:
    """de Bruijn's expansion ln rho(u) = -u (ln u + ln ln u - 1 + (ln ln u - 1)/ln u) + O(u (lnln u/ln u)^2)."""
    lu = math.log(u)
    llu = math.log(lu)
    return -u * (lu + llu - 1.0 + (llu - 1.0) / lu)


def dictionary(bits: float, kappa: float = KAPPA_STAR, nu: float = NU_STAR) -> Dict[str, float]:
    x = bits * math.log(2)
    lny = kappa * F(1 / 3, x)
    lnh = nu * F(2 / 3, x)
    y = math.exp(lny)
    n = y / lny                       # pi(y), leading order
    u0 = x / lny
    k = lnh / lny
    t = u0
    alpha = t / n                     # agreement fraction (prime coordinates, unweighted)
    R = k / n                         # rate
    johnson_t = math.sqrt(k * n)
    Rcap = 1.0 - Hq(1.0 - alpha, y)   # q-ary capacity rate at agreement fraction alpha (alphabet ~ y)
    ln_list = lnh + ln_rho(u0)        # ln of average list size (h + 1) rho(u0), de Bruijn expansion
    return {"bits": bits, "x": x, "ln_y": lny, "ln_h": lnh, "n": n, "ln_n": math.log(n), "k": k, "t": t,
            "t_over_k2": t / k ** 2, "ln_n_over_klnk": math.log(n) / (k * math.log(k)),
            "alpha": alpha, "R": R, "alpha_over_R": alpha / R, "johnson_t": johnson_t,
            "t_over_johnson": t / johnson_t, "R_cap": Rcap, "R_over_Rcap": R / Rcap,
            "s_logratio": math.log(alpha) / math.log(R), "ln_avg_list": ln_list,
            "ln_h_over_u0lnu0": lnh / (u0 * math.log(u0))}


if __name__ == "__main__":
    print(f"a* = {a_star():.10f}; c(1) = {c_of_a(1.0):.6f} (6^(1/3) = {6 ** (1/3):.6f}); "
          f"grid check c(a*) = {c_of_a_grid(a_star()):.6f} vs MNFS {MNFS:.6f}")
    for b in (512, 1024, 2048, 4096):
        d = dictionary(b)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items()})
