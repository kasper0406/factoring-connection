"""OPERATIONAL code: divisor-cover Lehman (DCL) search.

Receives only N, public parameters (lam, Hmax, optional m) and a seed.  Never
uses p, q, phi or group orders.  Every inversion mod N is gcd-checked and every
gcd hit is returned.  Collision detection reuses cycloharvey's Harvey Alg. 4.1
implementation (imported, not modified).

One shared baby list across all moduli
--------------------------------------
Target (A, B, M, r) knows beta^{M^2 g} = beta^{B + A N - c}  (mod p), where
g = (B p + A q - c)/M^2 (since p == 1, q == N mod p - 1).  With L = prod_M M^2
(over the distinct moduli used) and gamma = beta^L,

      gamma^g = (beta^{L/M^2})^{B + A N - c}   (mod p),

so every target is a giant step against the SAME baby list gamma^i.  All
beta^{L/M^2} are computed by a product tree in the exponent
(`all_but_one_powers`), O(log L * log #moduli) multiplications.
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Sequence

from src.common.metrics import Metrics
from src.cycloharvey.operational import (FactorFound, _check, find_collisions, inv_mod,
                                         recover_from_u)
from .dcl import best_design, design


def _pow(x: int, e: int, N: int, met) -> int:
    met.inc("rn_pow")
    met.inc("rn_mul_pow", max(0, e.bit_length() + bin(e).count("1") - 2))
    return pow(x, e, N)


def all_but_one_powers(x: int, exps: Sequence[int], N: int, met) -> List[int]:
    """[x^{prod_{i != j} exps[i]} for j] by recursive halving (remainder tree in
    the exponent).  Cost O(log(prod exps) * log len) multiplications."""
    if len(exps) == 1:
        return [x]
    h = len(exps) // 2
    left, right = exps[:h], exps[h:]
    pl, pr = 1, 1
    for e in left:
        pl *= e
    for e in right:
        pr *= e
    return (all_but_one_powers(_pow(x, pr, N, met), left, N, met)
            + all_but_one_powers(_pow(x, pl, N, met), right, N, met))


def dcl_search(N: int, lam: float, seed: int = 0, Hmax: int = 12, m: Optional[int] = None,
               met=None, early_exit: bool = True) -> Dict:
    met = met or Metrics()
    d = design(N, lam, m, Hmax) if m is not None else best_design(N, lam, Hmax)
    m, targets = d["m"], d["targets"]
    rng = random.Random(f"multimod-dcl|{N}|{seed}")
    found = None
    try:
        # trial division by the public moduli (p > M is assumed by the window
        # analysis; a hit here is a factor anyway)
        moduli = sorted({t["M"] for t in targets})
        for M in moduli:
            _check(M, N, met, "modulus")
        beta = rng.randrange(2, N - 1)
        _check(beta, N, met, "base")
        sq = [M * M for M in moduli]
        parts = all_but_one_powers(beta, sq, N, met)          # beta^{L/M^2}
        bm = dict(zip(moduli, parts))
        gamma = _pow(parts[0], sq[0], N, met)                  # beta^L
        _check(gamma - 1, N, met, "gamma_order")
        baby, x = {}, 1
        for i in range(m):
            if i:
                _check(x - 1, N, met, "baby_order")
            if x in baby:
                return {"factor": None, "status": "small_order", "design": _dsum(d),
                        "met": met.to_record()}
            baby[x] = i
            x = x * gamma % N
            met.inc("rn_mul")
        ginv = inv_mod(gamma, N, met)
        gm_inv = _pow(ginv, m, N, met)
        giants: List[int] = []
        meta: List[tuple] = []
        for ti, t in enumerate(targets):
            met.inc("handles")
            A, B, M, c, g_lo, W = t["A"], t["B"], t["M"], t["c"], t["g_lo"], t["W"]
            v = _pow(bm[M], B + A * N - c, N, met) * _pow(ginv, g_lo, N, met) % N
            met.inc("rn_mul")
            for k in range(-(-W // m)):
                if v in baby:                                   # two-sided (exact) match
                    met.inc("exact_match")
                    u = c + M * M * (g_lo + baby[v] + k * m)
                    dd = recover_from_u(N, A, B, u)
                    if dd and early_exit:
                        return {"factor": dd, "status": "exact_match", "design": _dsum(d),
                                "hit": {"baby": baby[v], "target": ti, "k": k},
                                "met": met.to_record()}
                    if dd and found is None:
                        found = dd
                else:
                    giants.append(v)
                    meta.append((ti, k))
                v = v * gm_inv % N
                met.inc("rn_mul")
        met.set("baby", m)
        met.set("giant", len(giants) + met.counts.get("exact_match", 0))
        if found is not None:
            return {"factor": found, "status": "exact_match", "design": _dsum(d),
                    "met": met.to_record()}
        g = find_collisions(N, gamma, m, giants, met)
        hit = met.fields.get("hit")
        info = None
        if g and hit:
            i, h = hit
            ti, k = meta[h]
            info = {"baby": i, "target": ti, "k": k}
        return {"factor": g, "status": "collision" if g else "none", "hit": info,
                "design": _dsum(d), "targets": targets, "met": met.to_record()}
    except FactorFound as e:
        return {"factor": e.d, "status": "gcd_" + e.how, "design": _dsum(d),
                "met": met.to_record()}


def _dsum(d: Dict) -> Dict:
    return {k: d[k] for k in ("lam", "m", "giants", "cost", "n_targets", "n_moduli",
                              "Mmin", "Mmax", "Hmax_used")}
