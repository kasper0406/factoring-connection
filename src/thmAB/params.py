"""thmAB: leading-order exponent bookkeeping for the cell-covariance proof of Theorems AB / AB2
(notes/thmAB_results.md §1.3).  Exact rationals only.

A quantity is represented by the leading term of its logarithm, c * D^e * lam^f, as the pair (e, f); lam = log log N,
D = log N / lam.  "A << B" holds at leading order iff (e_A, f_A) < (e_B, f_B) lexicographically (equal pairs are
reported as borderline).  This checks the note's exponent arithmetic; it is not an independent proof (same status as
src/c12alg/params.py and src/c12rat/params.py).

Design: c12alg's 2L_DW(alpha), alpha in [3/8, 2/5).
"""
from __future__ import annotations

from fractions import Fraction as F
from typing import Dict, List, Tuple

Term = Tuple[F, F]


def theta(al: F) -> F:
    """Child loss exponent: max(3a - 1, 1 - 2a, 2a - 1/2)."""
    return max(3 * al - 1, 1 - 2 * al, 2 * al - F(1, 2))


def theta_par(al: F) -> F:
    """Parent loss exponent: max(1 - 2a, 2 - 5a) (no child-filter spread, obligation O3)."""
    return max(1 - 2 * al, 2 - 5 * al)


def zeta_child(al: F) -> F:
    """S0 = exp(D^zeta): exceptional gcd-class threshold for children."""
    return (theta(al) + al) / 2


def zeta_par(al: F) -> F:
    """S0_par = exp(D^zeta_par)."""
    return (theta_par(al) + 4 * al - 1) / 2


P0_EXP = F(7, 10)          # P0 = ceil(D^{7/10})


def comparisons(al: F) -> List[Tuple[str, Term, Term]]:
    one, zero = F(1), F(0)
    s = 3 - 7 * al                      # parent precision exponent / u2 exponent
    th, thp = theta(al), theta_par(al)
    zc, zp = zeta_child(al), zeta_par(al)
    z0 = 3 * al - F(1, 2)               # eta0 log N = D^z0 lam
    return [
        # --- children (y = B, log B = a D^alpha lam, gamma ~ D^-alpha, u ~ D^alpha)
        ("theta = 2a - 1/2 on [3/8, 2/5)", (th, zero), (2 * al - F(1, 2), zero + F(1, 10**9))),
        ("child losses Xi = e^{D^theta lam^2} << saving e^{-c D^a / lam^2}", (th, F(2)), (al, F(-2))),
        ("kappa0 = T0^3 inside Hildebrand range (h_* >= B^{-5/12})", (th, F(2)), (al, one)),
        ("S0: Rankin tail S0^{-1/log P0} beats Xi", (th, F(2)), (zc, F(-1))),
        ("S0: gamma log S0 -> 0", (zc, zero), (al, zero)),
        ("S0: mid-prime tail exp(2 omega(Delta') P0^{-1/2}) << S0^{1/2}", (al - P0_EXP / 2, zero), (zc, zero)),
        ("omega(child Delta') = o(P0)", (al, zero), (P0_EXP, zero)),
        ("child moduli Delta' <= x^{o(1)} (Fact D level)", (al, one), (2 * al, one)),
        ("X^- = S B^{-C_X} in the relevant range (log X^- >= log S / 2)", (al, one), (2 * al, one)),
        ("TV N^{-eta0/2} << every L[alpha]^C", (al, one), (z0, one)),
        # --- parents (y = B2, log B2 = b2 D^{4a-1} lam, gamma2 ~ D^{1-4a}, u2 ~ D^{3-7a})
        ("parent losses Xi_p << parent saving (O3: separate thresholds)", (thp, F(2)), (s, F(-2))),
        ("kappa0_par inside Hildebrand range at B2", (thp, F(2)), (4 * al - 1, one)),
        ("S0_par: Rankin tail beats Xi_p", (thp, F(2)), (zp, F(-1))),
        ("S0_par: gamma2 log S0_par -> 0", (zp, zero), (4 * al - 1, zero)),
        ("S0_par: mid-prime tail exp(2 omega P0^{-1/2}) << S0_par^{1/2}", (1 - al - P0_EXP / 2, zero), (zp, zero)),
        ("omega(parent Delta') = o(P0)", (1 - al, zero), (P0_EXP, zero)),
        ("parent moduli Delta' <= x_F^{o(1)} (Fact D level)", (1 - al, one), (2 - 3 * al, one)),
        ("u2 -> infinity faster than lam^4", (zero, F(4)), (s, zero)),
        ("Hildebrand's range at B2: B2 >= exp((log log x)^{5/3+eps})", (zero, F(2)), (4 * al - 1, one)),
        ("parent FO tail: B2^{-1/2} beats exp(omega P0^{-0.49})", (1 - al - F(49, 100), zero), (4 * al - 1, one)),
        # --- outputs vs budgets
        ("AB(i): second-order child D^{-2a} << first-order need D^{-a}", (-2 * al, F(20)), (-al, F(-20))),
        ("AB(ii): parent D^{2-8a} << need D^{7a-3}", (2 - 8 * al, F(20)), (7 * al - 3, F(-20))),
        ("AB(ii): parent D^{-2(3-7a)} << need D^{7a-3}", (-2 * s, F(20)), (-s, F(-20))),
        ("parent arch. variance D^{2-8a} <= D^{-2a}", (2 - 8 * al, zero), (-2 * al + F(1, 10**9), zero)),
        ("AB2 budget: A/eta -> 0", (-2 * al + s, F(20)), (zero, zero)),
        ("AB2 budget: sqrt(AB)/eta -> 0", (-al - min(al, s) + s, F(20)), (zero, zero)),
    ]


def acceptability_ok(al: F) -> str:
    """Parent acceptability subtraction (log N)^C / B against rho2^2: ok in the open window; at 3/8 it needs a > 2 c2."""
    s = 3 - 7 * al
    return "ok" if s < al else ("needs a > 2 c2 (O4)" if s == al else "FAIL")


def margins(al: F) -> Dict[str, str]:
    out = {}
    for name, a, b in comparisons(al):
        out[name] = "ok" if a < b else ("borderline" if a == b else "FAIL")
    out["parent acceptability"] = acceptability_ok(al)
    return out


def all_ok(al: F) -> bool:
    m = margins(al)
    return all(v == "ok" for k, v in m.items() if k != "parent acceptability")
