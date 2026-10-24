"""Leading-order exponent bookkeeping for Theorem RAT (notes/c12rat_results.md §1.3).  Exact rationals only.

Every quantity is a power D^e (times powers of lambda = log log N, ignored at leading order).  A margin is a
comparison e_small < e_large of exponents of D; it holds for N >= N_0 iff the strict inequality holds.

Designs:
  'W'  c12alg's window design 2L_DW(alpha), 3/8 < alpha < 2/5 (B2 = L[4 alpha - 1], eta0 log N = D^{3 alpha - 1/2});
  'P'  pc3/uprime design 2L' with DW at alpha = 3/8 (B2 = L[1/2], eta0 = D^{-1/2}).
"""
from __future__ import annotations

from fractions import Fraction as F
from typing import Dict, List, Tuple


def exponents(alpha: F, design: str) -> Dict[str, F]:
    a = F(alpha)
    e = {}
    e["logB"] = a                       # log y = log B = a D^alpha lambda
    e["gamma"] = -a                     # 1 - beta  <<  1/(a D^alpha)   (HT (alpha2))
    e["logX_child"] = 2 * a             # log of child rational values  (log M + log V2)
    e["logX_par"] = 2 * a               # log of parent rational cofactor X_r
    e["u"] = a                          # u = log X / log B  (child and parent rational sides)
    e["logY"] = 1 - a                   # special-q l = L[1 - alpha]; V1^2 ~ Y S1^2
    if design == "W":
        e["logB2"] = 4 * a - 1
        e["eta0logN"] = 3 * a - F(1, 2)
    elif design == "P":
        assert a == F(3, 8)
        e["logB2"] = F(1, 2)
        e["eta0logN"] = F(1, 2)
    else:
        raise ValueError(design)
    e["logK_child"] = e["logB2"]        # child moduli K = q Delta <= 4 Q2 S2^2
    e["logK_par"] = e["logY"]           # parent moduli K <= 2 V1^2
    e["logP0"] = F(0)                   # P0 = D^{7/10}: log P0 ~ lambda (exponent 0)
    e["P0"] = F(7, 10)
    # losses (exponent of log of the loss)
    e["Xi_child_moduli"] = e["gamma"] + e["logK_child"]      # (|P'||P''|)^{1-beta}, |P'| <= V2 ~ B2 S2
    e["Xi_par_moduli"] = e["gamma"] + e["logY"]               # (V1^2)^{1-beta}
    e["Xi_filter"] = e["gamma"] + e["eta0logN"]               # (X / X^-)^{1-beta}, X^- = X N^{-eta0/4}
    e["saving"] = e["u"]                                       # log(1/eps) ~ u / log^2 u
    return e


def margins(alpha: F, design: str) -> List[Tuple[str, F, F, bool]]:
    """(name, small, large, holds): holds iff small < large (strict)."""
    e = exponents(alpha, design)
    out = []

    def m(name, small, large):
        out.append((name, F(small), F(large), F(small) < F(large)))

    # Fact D: every modulus is a small conductor, K <= x^{eta_D}
    m("child moduli below x^eta_D", e["logK_child"], e["logX_child"])
    m("parent moduli below x^eta_D", e["logK_par"], e["logX_par"])
    # P0: omega(K) / P0 -> 0 needs P0 >> log K; P0# << interval lengths needs P0 << log X
    m("P0 beats omega(child K)", e["logK_child"], e["P0"])
    m("P0 beats omega(parent K)", e["logK_par"], e["P0"])
    m("P0# below the cell intervals", e["P0"], e["logX_child"])
    # losses below the saving
    for k in ("Xi_child_moduli", "Xi_par_moduli", "Xi_filter"):
        m(k + " < BDH saving", e[k], e["saving"])
    # Hildebrand's range: relative cell widths >= y^{-5/12}: every loss exponent < log B exponent
    for k in ("Xi_child_moduli", "Xi_par_moduli", "Xi_filter"):
        m(k + " inside Hildebrand range", e[k], e["logB"])
    # second-order output D^{-2 alpha} beats the needs: (R1) D^{-alpha}, (R2) D^{7 alpha - 3}, and eta at 3/8
    m("D^-2a below (R1) budget D^-a", -2 * e["u"], -e["u"])
    m("D^-2a below (R2) budget D^(7a-3)", -2 * e["u"], 7 * F(alpha) - 3)
    m("D^-2a below eta = D^-(3-7a)", -2 * e["u"], -(3 - 7 * F(alpha)))
    # parent box much larger than the child primes (label weights flat in r): V1 >> B2
    m("B2 below V1 (V1 >> q)", e["logB2"], e["logY"])
    return out


def all_hold(alpha: F, design: str) -> bool:
    return all(h for *_, h in margins(alpha, design))
