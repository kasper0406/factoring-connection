"""c12alg E2: the thmA toy pipeline run with design-DW fields (composite integer multiplier w).

OPERATIONAL part: factor_dw(N, cfg, seed) sees only N, the public config and a seed (as src/thmA/pipeline.py).
It reuses the thmA pipeline by import; the only replaced step is SAMPLE_FIELD, swapped in for the duration of
the call (the thmA module is restored afterwards).  The DW repairs (notes/c12alg_results.md §4) are:
  * w uniform among ALL integers in (W, 2W] (no primality, no coprimality with g1); gcd(w, N) checked;
  * F = w * Fhat_{g,N} + (g1 x - g0 y) R_c  (multiplier outside the base expansion), F(g0, g1) = w N;
  * degree drop f_d = 0 -> reject;  content e > 1 -> factor exit if 1 < gcd(e, N) < N, otherwise reject;
  * irreducibility certificate: unchanged code; primes r | w cannot certify (f mod r has the root m_r or
    degree < d), which the analysis accounts for (Lemma IR-DW);
  * N | f'(m) -> reject (unchanged).
Diagnostic counters added: Lemma CU (Fhat_N(P) a unit mod |G(P)|) checked on every returned atom.
"""
from __future__ import annotations

import contextlib
import random
from math import gcd
from typing import Dict, Optional

from src.thmA import arith
from src.thmA import numfield as nf
from src.thmA import pipeline as tp

from . import dw


def sample_field_dw(N: int, P: Dict, rng: random.Random, ctr: Dict) -> Optional[Dict]:
    d, M, H, W = P["d"], P["M"], P["H"], P["W"]
    while True:
        g0, g1 = rng.randrange(M // 2, M + 1), rng.randrange(max(2, M // 2), M + 1)
        if gcd(g0, g1) == 1:
            break
    w = rng.randrange(W + 1, 2 * W + 1)                 # uniform INTEGER (design DW)
    c = [rng.randrange(-H, H + 1) for _ in range(d)]
    tp._inc(ctr, "dw_w_composite", 0 if arith.is_prime(w) else 1)
    if not tp._gcd_check(g1, N, "field_g1", ctr) or not tp._gcd_check(w, N, "field_w", ctr):
        tp._inc(ctr, "rej_gcdN"); return None
    f = dw.dw_coeffs(N, w, d, g0, g1, c)
    t = f[d]
    if t == 0:
        tp._inc(ctr, "rej_degree_drop"); return None
    if not tp._gcd_check(t, N, "field_fd", ctr):
        tp._inc(ctr, "rej_gcdN"); return None
    e = 0
    for x in f:
        e = gcd(e, x)
    if e != 1:
        if gcd(e, N) not in (1, N):
            if ctr.get("_continue"):
                ctr.setdefault("accidental_continued", []).append("field_content")
                return None
            raise tp.FactorFound(gcd(e, N), "field_content")
        assert w % e == 0 or gcd(e, N) == N, "content must divide w (Lemma CT)"
        tp._inc(ctr, "rej_content"); return None
    usable = [r for r in P["irr_primes"] if (t * N) % r]
    tp._inc(ctr, "ir_primes_dividing_w", sum(1 for r in usable if w % r == 0))
    if not any(nf.irreducible_mod(f, r) for r in usable):
        tp._inc(ctr, "rej_irreducibility_uncertified"); return None
    # primes r | w with r not dividing g1 never certify (f(m_r) = 0 mod r): record violations of that claim
    for r in usable:
        if w % r == 0 and g1 % r and nf.irreducible_mod(f, r):
            tp._inc(ctr, "IR_DW_CLAIM_VIOLATION")
    g1inv = pow(g1, -1, N)
    m = g0 * g1inv % N
    fp = [(i + 1) * f[i + 1] for i in range(d)]
    fpm = 0
    for cc in reversed(fp):
        fpm = (fpm * m + cc) % N
    if not tp._gcd_check(fpm, N, "field_fprime", ctr):
        tp._inc(ctr, "rej_N_divides_fprime"); return None
    D = t * nf.disc(f)
    g = nf.monic_g(f)
    return {"g0": g0, "g1": g1, "w": w, "c": c, "f": f, "t": t, "D": D, "g": g,
            "beta_img": t * m % N, "g1inv": g1inv, "sq_cache": {}, "fhat": dw.fhat(N, d, g0, g1), "N": N}


_orig_sq = tp.sq_search


def _sq_search_checked(N, A, ell, P, ctx, ctr):
    res = _orig_sq(N, A, ell, P, ctx, ctr)
    if res is not None and "fhat" in A:
        tval = abs(res["G"])
        if gcd(tval, N) == 1:
            u = nf.hom_eval(A["fhat"], res["a"], res["b"])
            tp._inc(ctr, "cu_checked")
            if gcd(u, tval) != 1:
                tp._inc(ctr, "CU_VIOLATION")
            # the forced residue under DW: F(P) = w * Fhat_N(P) mod |G(P)|
            if (res["F"] - A["w"] * u) % tval:
                tp._inc(ctr, "FORCED_RESIDUE_VIOLATION")
    return res


@contextlib.contextmanager
def dw_fields():
    saved = (tp.sample_field, tp.sq_search)
    tp.sample_field, tp.sq_search = sample_field_dw, _sq_search_checked
    try:
        yield
    finally:
        tp.sample_field, tp.sq_search = saved


def factor_dw(N: int, cfg: Dict, seed: int) -> Dict:
    with dw_fields():
        return tp.factor(N, cfg, seed)
