"""Cost formulas for sinonobliv (MODEL; floats only in model evaluation).

- `si_cost(tau)`: factoring constant with a constructor costing L[1/3, tau] per output:
  min_kappa max(kappa + tau + 2/(9 kappa^2), 2 kappa)  (suffmap §1.4; siconstruct §1).
- `tau_threshold()`: the tau at which si_cost meets MNFS ((92 + 26 sqrt 13)/27)^{1/3} = 1.9019.
- `prop_cd_log_hits`: natural log of the expected number of smooth members per window for an m-parameter,
  degree-d identity family with r forced pieces (Prop CD, heuristic H_piece), at x = ln X.
- `ms_case_bound`: the two-case bound of Prop MS (independent windows, Cohn–Heninger region).
- `toy_path`: siconstruct's toy L-path parameters (y, h) at a bit size.
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

MNFS = ((92 + 26 * math.sqrt(13)) / 27) ** (1.0 / 3.0)
KAPPA_STAR = (2.0 / 9.0) ** (1.0 / 3.0)
NU_STAR = 2.0 / (3.0 * KAPPA_STAR)


def si_cost(tau: float, grid: int = 200000) -> Tuple[float, float]:
    """Return (c, kappa_opt) minimising max(kappa + tau + 2/(9 kappa^2), 2 kappa) over kappa in (0.2, 2)."""
    best = (float("inf"), None)
    for i in range(1, grid):
        k = 0.2 + 1.8 * i / grid
        c = max(k + tau + 2.0 / (9.0 * k * k), 2.0 * k)
        if c < best[0]:
            best = (c, k)
    return best


def tau_threshold(target: float = MNFS) -> float:
    lo, hi = 0.0, 2.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if si_cost(mid, grid=20000)[0] < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def F(alpha: float, x: float) -> float:
    return x ** alpha * math.log(x) ** (1.0 - alpha)


def prop_cd_log_hits(x: float, m: int, d: int, r: int, kappa: float = KAPPA_STAR, nu: float = NU_STAR) -> Dict[str, float]:
    """ln E[#smooth members per window] for a family with m parameters, degree d, r pieces, at ln X = x.

    Hits per window ~ (h+1) X^{m/d - 1} (definite forms); smoothness ~ rho(u) * r^u (equal pieces, H_piece);
    rho via -u ln u (leading order). Returns the three log-contributions and the total.
    """
    lny = kappa * F(1.0 / 3.0, x)
    u = x / lny
    lnh = nu * F(2.0 / 3.0, x)
    placement = lnh + (m / d - 1.0) * x
    smooth = -u * math.log(u)
    gain = u * math.log(r) if r > 1 else 0.0
    return {"placement": placement, "smooth": smooth, "gain": gain, "total": placement + smooth + gain}


def ms_case_bound(t: float, pi: float, m: int, ell: float) -> Dict[str, float]:
    """Prop MS case split, as in the note: for a certified divisor with ln d = ell (sigma_m <= ell <= pi),
    case 1 if ell >= j1 = m t/(m-1) + 1 (geometric decay); otherwise case 2, where
    pi/t <= (ell/t)^{m+1} < (m/(m-1) + 1/t)^{m+1} <= 8 e^2 for 2 <= m <= t."""
    sig = t ** (m / (m + 1.0)) * pi ** (1.0 / (m + 1.0))
    j1 = m * t / (m - 1.0) + 1.0
    return {"sigma": sig, "j1": j1, "case": 1 if ell >= j1 else 2,
            "case2_pi_cap": (m / (m - 1.0) + 1.0 / t) ** (m + 1) * t}


def toy_path(bits: int, kappa: float = KAPPA_STAR, nu: float = 1.376) -> Tuple[int, int]:
    """siconstruct's toy L-path (o(1) dropped): y = round(e^{kappa F_1/3}), h = round(e^{nu F_2/3})."""
    x = bits * math.log(2)
    return int(round(math.exp(kappa * F(1 / 3, x)))), int(round(math.exp(nu * F(2 / 3, x))))


def consecutive_product_log_hits(x: float, r: int, kappa: float = KAPPA_STAR, nu: float = NU_STAR) -> float:
    """ln of expected smooth members per window for t(t+1)...(t+r-1): hits (h+1)/(r X^{1-1/r}), smoothness
    rho(u/r)^r (H_piece)."""
    lny = kappa * F(1.0 / 3.0, x)
    u = x / lny
    lnh = nu * F(2.0 / 3.0, x)
    return lnh - math.log(r) - (1 - 1.0 / r) * x - u * math.log(u / r)


if __name__ == "__main__":
    print("MNFS", MNFS)
    print("si_cost(0)", si_cost(0.0))
    print("tau threshold", tau_threshold())
    for a in (1.0, 1.1, 1.2, 2.0):
        # tau = a*kappa: per-output cost polynomial of degree a in pi(y)
        best = min((max(k * (1 + a) + 2 / (9 * k * k), 2 * k), k) for k in [0.3 + i * 1e-4 for i in range(15000)])
        print("tau = %.1f kappa:" % a, best)
    for bits in (512, 1024, 2048, 4096):
        x = bits * math.log(2)
        print(bits, {r: round(consecutive_product_log_hits(x, r), 1) for r in (2, 3, 4, 8)})
