"""Exponent bookkeeping for design (a''_k) (notes/pcest_results.md §1, §5).

All sizes are exponents of N.  s = 1/(d+1); special-q prime l = N^(1/(2k)); box points
|P| <= V = N^(1/(4k)) L^(o(1)); rational cofactor n = G/l = N^(s - 1/(4k));
algebraic value F = N^(s + h + d/(4k)), H = N^h the height of R.
Constraints on h (eta -> 0 slack omitted):
- level 1/2 for FM:     |G| <= |F|^(1/2):    h >= s - (d - 2)/(4k)
- LV kernel (FM):       H >> V:              h >= 1/(4k)
- pair kernel (Lemma Kern, PC_alg):  H >> V^2:  h >= 1/(2k)
theta = (s - 1/(4k)) + (s + h + d/(4k)) = 2s + h + (d - 1)/(4k).
With the Lenstra-Pomerance shape test, c = (1 + theta)/sqrt(2) (hypD Corollary K).
"""
from __future__ import annotations

from fractions import Fraction
from math import sqrt
from typing import Dict


def theta(d: int, k: int, floor: str = "pair") -> Fraction:
    s = Fraction(1, d + 1)
    fl = {"lv": Fraction(1, 4 * k), "pair": Fraction(1, 2 * k)}[floor]
    h = max(fl, s - Fraction(d - 2, 4 * k))
    return 2 * s + h + Fraction(d - 1, 4 * k)


def c_const(th) -> float:
    return (1 + float(th)) / sqrt(2)


def c_limit(d: int) -> float:
    """Infimum over growing fixed k: theta -> 3/(d+1)."""
    return c_const(Fraction(3, d + 1))


def binding_floor_k(d: int, floor: str = "pair") -> Fraction:
    """k below which the kernel floor (not the level condition) sets h."""
    # floor > s - (d-2)/(4k)  <=>  k < (d+1)(c0 + d - 2)/4 with floor = c0/(4k)
    c0 = {"lv": 1, "pair": 2}[floor]
    return Fraction((d + 1) * (c0 + d - 2), 4)


def table(ds=(7, 9, 11, 19), ks=(5, 10, 20, 50, 100, 1000)) -> Dict:
    out = {}
    for d in ds:
        row = {}
        for k in ks:
            t_lv, t_pair = theta(d, k, "lv"), theta(d, k, "pair")
            row[k] = dict(theta_lv=float(t_lv), theta_pair=float(t_pair),
                          c_lv=c_const(t_lv), c_pair=c_const(t_pair))
        out[d] = dict(rows=row, c_limit=c_limit(d),
                      k_floor_lv=float(binding_floor_k(d, "lv")),
                      k_floor_pair=float(binding_floor_k(d, "pair")))
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(table(), indent=1))
