"""OPERATIONAL: DCL-wide, deterministic ratio-band DCL (Theorem U / U_theta).

Receives only N and a public WideDesign (params.py).  Never uses p, q, phi(N),
group orders, factoring or primality routines.  Base 2.  Identical to
dcltheorem.algorithm.dcl_det (Theorem D') except that each weight instance
carries its own ratio interval (rho_lo, rho_hi]: its window is computed over the
integer p-range {p : rho_lo < N/p^2 <= rho_hi} instead of the balanced range.
With the single instance family (A,B,1,2), (A,B) in V(H), the target list is
identical to D''s (tests/test_unbalanced.py checks this).

Per block k: exponent tree u_l = 2^{L/l^2}, gamma = 2^L; for each band prime l
and instance (A,B,rho): roots r of B r^2 == A N (mod l), c = (B r + A N r^{-1})
mod l^2, window [g_lo, g_hi] from Lo = floor(2 sqrt(ABN)) and
Hi = max_{x in {x_lo, x_hi}} ceil(B x + A N/x); one shared baby list
gamma^i (0 <= i < m) for ALL instances of the block; giants
u_l^{B + A N - c} gamma^{-(g_lo + j m)}; d = gcd(f(G), N).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

from src.common.metrics import Metrics
from src.dcltheorem.algorithm import (Found, _gcd_check, _proper, all_but_one, choose_m,
                                      metric_E, powN, pow_small, sqrt_prime, _NONRES)
from src.dcltheorem.stage3 import eval_baby_poly_at_giants

from .params import Inst, WideDesign


def ratio_range(N: int, t: Inst) -> Tuple[int, int]:
    """Integers p with rho_lo < N/p^2 <= rho_hi lie in [x_lo, x_hi] (x_lo > x_hi: empty).
    rho_hi = a/b:  N/p^2 <= a/b  <=>  a p^2 >= N b  ->  x_lo = ceil(sqrt(N b / a)).
    rho_lo = a'/b': N/p^2 > a'/b' <=> a' p^2 < N b'  ->  x_hi = max x with a' x^2 < N b'."""
    a, b = t.rhi_num, t.rhi_den
    x_lo = isqrt((N * b) // a)
    while a * x_lo * x_lo < N * b:
        x_lo += 1
    while x_lo > 1 and a * (x_lo - 1) ** 2 >= N * b:
        x_lo -= 1
    a2, b2 = t.rlo_num, t.rlo_den
    x_hi = isqrt((N * b2) // a2) + 1
    while a2 * x_hi * x_hi >= N * b2:
        x_hi -= 1
    return max(1, x_lo), x_hi


def window_bounds(N: int, t: Inst, met) -> Optional[Tuple[int, int]]:
    x_lo, x_hi = ratio_range(N, t)
    if x_lo > x_hi:
        return None
    met.inc("isqrt_ops", 3)
    lo = isqrt(4 * t.A * t.B * N)
    hi = max(t.B * x + -(-(t.A * N) // x) for x in (x_lo, x_hi))
    return lo, hi


def build_targets(N: int, design: WideDesign, k: int, met) -> List[Dict]:
    band = design.bands[k]
    _NONRES.clear()
    bounds = {}
    for t in design.insts:
        wb = window_bounds(N, t, met)
        if wb is not None:
            bounds[t] = wb
    out = []
    for l in band:
        if N % l == 0:
            if 1 < l < N:
                raise Found(l, "band_prime")
            continue
        l2 = l * l
        Nl = N % l
        for t, (lo, hi) in bounds.items():
            A, B = t.A, t.B
            met.inc("vl_pairs")
            if (A * B) % l == 0:          # no root (and no coverage: l | Bp - Aq => l | pq)
                continue
            a = A * Nl * pow(B, -1, l) % l
            met.inc("mul_small", 3)
            if pow_small(a, (l - 1) // 2, l, met) != 1:
                continue
            r = sqrt_prime(a, l, met)
            for root in sorted({r, l - r}):
                s = N * pow(root, -1, l2) % l2
                c = (B * root + A * s) % l2
                met.inc("mul_small", 4)
                g_lo = -((-(lo - c)) // l2)
                g_hi = (hi - c) // l2
                if g_hi < g_lo:
                    continue
                out.append({"A": A, "B": B, "j": t.band, "rho": (t.rlo_num, t.rlo_den, t.rhi_num, t.rhi_den),
                            "l": l, "r": root, "c": c, "g_lo": g_lo, "W": g_hi - g_lo + 1})
    return out


def run_block(N: int, design: WideDesign, k: int, met, record: bool = False,
              m_override: Optional[int] = None) -> Dict:
    band = design.bands[k]
    if not band:
        return {"factor": None, "status": "empty_band"}
    targets = build_targets(N, design, k, met)
    if not targets:
        return {"factor": None, "status": "no_targets", "m": 0, "n_targets": 0, "n_giants": 0}
    m = m_override or choose_m(targets, 0)
    sq = [l * l for l in band]
    u_list, L = all_but_one(2, sq, N, met)
    u = dict(zip(band, u_list))
    gamma = powN(u_list[0], sq[0], N, met)
    _gcd_check(gamma, N, met, "gamma")
    met.inc("invN")
    ginv = pow(gamma, -1, N)
    gm_inv = powN(ginv, m, N, met)
    giants, meta = [], []
    w = {}
    for t in targets:
        if t["l"] not in w:
            w[t["l"]] = powN(u[t["l"]], N, N, met)
    for ti, t in enumerate(targets):
        A, B, l, c, g_lo, W = t["A"], t["B"], t["l"], t["c"], t["g_lo"], t["W"]
        # G = u^{B + A N - c} gamma^{-g_lo} = w^A u^{B - c + l^2} gamma^{-(g_lo+1)}
        e1 = B - c + l * l
        G = powN(w[l], A, N, met) * powN(u[l], e1, N, met) % N
        if g_lo + 1 >= 0:
            G = G * powN(ginv, g_lo + 1, N, met) % N
        else:
            G = G * powN(gamma, -(g_lo + 1), N, met) % N
        met.inc("mulN", 2)
        for jj in range(-(-W // m)):
            giants.append(G)
            meta.append((ti, jj))
            G = G * gm_inv % N
            met.inc("mulN")
    met.inc("giants", len(giants))
    met.inc("babies", m)
    met.inc("targets", len(targets))
    vals = eval_baby_poly_at_giants(N, gamma, m, giants, met)
    res = {"factor": None, "status": "none", "m": m, "n_targets": len(targets),
           "n_giants": len(giants), "L_bits": L.bit_length()}
    gs = []
    for h, v in enumerate(vals):
        met.inc("gcdN")
        g = gcd(v, N)
        if record:
            gs.append(g)
        if 1 < g < N and res["factor"] is None:
            ti, jj = meta[h]
            res.update({"factor": g, "status": "collision", "hit": {"target": ti, "j": jj}})
            if not record:
                break
    if record:
        res.update({"targets": targets, "giant_meta": meta, "giant_gcds": gs,
                    "gamma": gamma, "u": u})
    return res


def dclw(N: int, design: WideDesign, met: Optional[Metrics] = None, early_exit: bool = True,
         record: bool = False, blocks: Optional[Sequence[int]] = None) -> Dict:
    """DCL-wide.  Returns {"factor": d or None, "block": k, "blocks": [...]};
    every output passes the proper-divisor check."""
    met = met or Metrics()
    out = {"factor": None, "status": "fail", "block": None, "blocks": []}
    if N < 4:
        return {"factor": None, "status": "fail_small", "block": None, "blocks": []}
    met.inc("gcdN")
    if N % 2 == 0:
        return {"factor": 2, "status": "even", "block": None, "blocks": []}
    for k in (blocks if blocks is not None else range(design.K)):
        bm = Metrics()
        try:
            r = run_block(N, design, k, bm, record=record)
        except Found as e:
            r = {"factor": e.d, "status": "gcd_" + e.how}
        if r.get("factor") is not None and not _proper(r["factor"], N):
            r["improper"] = r["factor"]
            r["factor"] = None
            r["status"] = "rejected_improper"
        r["k"] = k
        r["met"] = bm.to_record()
        met.merge(bm)
        out["blocks"].append(r)
        if r["factor"] and out["factor"] is None:
            out.update({"factor": r["factor"], "status": r["status"], "block": k})
            if early_exit:
                break
    return out


__all__ = ["dclw", "build_targets", "run_block", "ratio_range", "window_bounds", "metric_E"]
