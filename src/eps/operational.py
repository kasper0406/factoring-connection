"""The frozen statistic menu of the eps/QR_2 battery.  OPERATIONAL: `compute_stats(N)` sees only N.

Definitions are in notes/eps_prereg.md section 2; the canonical name list is `MENU` (order frozen).
Conventions: jac(v) = Jacobi symbol (v mod N / N) in {-1, 0, 1} (0 iff gcd(v, N) > 1, recorded as a
nonunit); par(v) = (v mod N) & 1 on the representative in [0, N); all statistics are bounded integers.
"""
from __future__ import annotations

import time
from typing import Dict, List, Tuple

import gmpy2

from .algebra import (Nonunits, bm_mod_N, det_mod_N, jac, mult_matrix, rank_mod_N, section_values,
                      sketch, sketch_G)
from .cf import cf_features
from .forms import form_features
from .mr import mr_features

PHI5 = [1, 1, 1, 1, 1]
PHI12 = [1, 0, -1, 0, 1]
X3M2 = [-2, 0, 0, 1]
FIB = [-1, -1, 1]

# (name, group, role) -- role: 'nc' negative control (floor), 'menu', 'pc_eps', 'pc_qr2'
MENU: List[Tuple[str, str, str]] = (
    [(n, "NC", "nc") for n in (
        "nc_N_mod8", "nc_jac2_N", "nc_N_mod256", "nc_bit8", "nc_bit9", "nc_bit10", "nc_bit11",
        "nc_N_mod3_is1", "nc_jacN_5", "nc_jacN_7", "nc_lead16", "nc_chisum1000")]
    + [(f"cf{T}_{k}", "CF", "menu") for T in (64, 1024) for k in (
        "n_a1", "n_a2", "sumbits", "maxbits", "n_Qeven", "n_sq_even", "n_Q1mod4_odd", "n_2sq_odd")]
    + [(f"{k}_T{T}", "PC", "pc_eps") for k in ("pc1", "pc2", "pc3") for T in (64, 1024)]
    + [(n, "PC", "pc_eps") for n in ("pc4_mr", "pc5_ss")]
    + [(f"{k}_T{T}", "PC", "pc_qr2") for k in ("pcq2", "pcq3") for T in (64, 1024)]
    + [(n, "A", "menu") for n in (
        "A1c2_jac", "A1c2_par", "A1c3_jac", "A1c3_par",
        "A2c3_jac0", "A2c3_jac1", "A2c3_jacsum01", "A2c3_jacdif01", "A2c3_par0", "A2c3_par1",
        "A4c2_jac0", "A4c2_jac1", "A4c2_jac2", "A4c2_jac3", "A4c2_par0", "A4c2_par1", "A4c2_par2", "A4c2_par3",
        "A4c5_jacsum", "A4c5_parsum",
        "A8c3_jacsum", "A8c3_parsum", "A8c3_jac0",
        "A16c5_jacsum", "A16c5_parsum", "A16c5_mod3is0")]
    + [(n, "B", "menu") for n in (
        "Bx3m2_jac0", "Bx3m2_jac1", "Bx3m2_jac2", "Bx3m2_par0", "Bx3m2_mod3is0", "Bx3m2_mod5is0",
        "Bx3m2_jacres", "Bx3m2_rank",
        "Bphi5_jac0", "Bphi5_jacsum", "Bphi5_parsum",
        "Bphi12_jac0", "Bphi12_jacsum",
        "Bfib_jac0", "Bfib_jac1", "Bfib_par0",
        "BG1x3m2_jac0", "BG1x3m2_jac1")]
    + [(n, "BM", "menu") for n in ("bm_Fc", "bm_coef64", "bm_euler", "bm_coef64_jacC")]
    + [(n, "CL", "menu") for n in (
        "cl_cnt", "cl_cnt_odd", "cl_cnt_prime", "cl_cnt_le10", "cl_cnt_le25", "cl_min_a")]
)
MENU_NAMES = [m[0] for m in MENU]
S_BONFERRONI = 100

_ODD_PRIMES_1000 = [x for x in range(3, 1001) if all(x % d for d in range(2, int(x ** 0.5) + 1))]


def _par(v: int, N: int) -> int:
    return (v % N) & 1


def compute_stats(N: int) -> Tuple[Dict[str, int], Dict]:
    """All menu statistics for N (a positive odd non-square integer).  Returns (stats, aux).

    aux: 'nonunits' list of (where, divisor), 'gcd_checks', CF first-hit indices, per-group seconds."""
    nu = Nonunits(N)
    st: Dict[str, int] = {}
    tm: Dict[str, float] = {}
    n = N.bit_length()

    t0 = time.perf_counter()
    # ---- NC: floor-type negative controls
    st["nc_N_mod8"] = N % 8
    st["nc_jac2_N"] = int(gmpy2.jacobi(2, N))
    st["nc_N_mod256"] = N % 256
    for b in (8, 9, 10, 11):
        st[f"nc_bit{b}"] = (N >> b) & 1
    st["nc_N_mod3_is1"] = int(N % 3 == 1)
    st["nc_jacN_5"] = int(gmpy2.jacobi(N, 5))
    st["nc_jacN_7"] = int(gmpy2.jacobi(N, 7))
    st["nc_lead16"] = N >> (n - 16)
    st["nc_chisum1000"] = sum(int(gmpy2.kronecker(-N, l)) for l in _ODD_PRIMES_1000)
    tm["NC"] = time.perf_counter() - t0

    # ---- CF prefix (menu + positive controls)
    t0 = time.perf_counter()
    cf = cf_features(N)
    aux_cf = {k: cf.pop(k) for k in list(cf) if k.startswith("first_") or k in ("cf_nsmooth", "cf_nonunit", "cf_fb_ok")}
    st.update(cf)
    mr = mr_features(N)
    aux_cf["mr_split"] = mr.pop("mr_split")
    st.update(mr)
    tm["CF"] = time.perf_counter() - t0

    # ---- Track A section values
    t0 = time.perf_counter()
    for c in (2, 3):
        (v,) = section_values(N, 1, c)
        st[f"A1c{c}_jac"] = jac(v, N, nu, f"A1c{c}")
        st[f"A1c{c}_par"] = _par(v, N)
    s0, s1 = section_values(N, 2, 3)
    st["A2c3_jac0"] = jac(s0, N, nu, "A2c3_0")
    st["A2c3_jac1"] = jac(s1, N, nu, "A2c3_1")
    st["A2c3_jacsum01"] = jac(s0 + s1, N, nu, "A2c3_s")
    st["A2c3_jacdif01"] = jac(s0 - s1, N, nu, "A2c3_d")
    st["A2c3_par0"] = _par(s0, N)
    st["A2c3_par1"] = _par(s1, N)
    V = section_values(N, 4, 2)
    for r in range(4):
        st[f"A4c2_jac{r}"] = jac(V[r], N, nu, f"A4c2_{r}")
    for r in range(4):
        st[f"A4c2_par{r}"] = _par(V[r], N)
    V = section_values(N, 4, 5)
    st["A4c5_jacsum"] = sum(jac(x, N, nu, "A4c5") for x in V)
    st["A4c5_parsum"] = sum(_par(x, N) for x in V)
    V = section_values(N, 8, 3)
    st["A8c3_jacsum"] = sum(jac(x, N, nu, "A8c3") for x in V)
    st["A8c3_parsum"] = sum(_par(x, N) for x in V)
    st["A8c3_jac0"] = jac(V[0], N, nu, "A8c3_0")
    V = section_values(N, 16, 5)
    st["A16c5_jacsum"] = sum(jac(x, N, nu, "A16c5") for x in V)
    st["A16c5_parsum"] = sum(_par(x, N) for x in V)
    st["A16c5_mod3is0"] = int(V[0] % 3 == 0)
    tm["A"] = time.perf_counter() - t0

    # ---- Track B sketches
    t0 = time.perf_counter()
    S = sketch(N, X3M2)
    for j in range(3):
        st[f"Bx3m2_jac{j}"] = jac(S[j], N, nu, f"Bx3m2_{j}")
    st["Bx3m2_par0"] = _par(S[0], N)
    st["Bx3m2_mod3is0"] = int(S[0] % 3 == 0)
    st["Bx3m2_mod5is0"] = int(S[0] % 5 == 0)
    Mx = mult_matrix(N, X3M2, S)
    res = det_mod_N(Mx, N, nu, "Bx3m2_det")
    st["Bx3m2_jacres"] = 0 if res is None else jac(res, N, nu, "Bx3m2_res")
    st["Bx3m2_rank"] = rank_mod_N(Mx, N, nu, "Bx3m2_rank")
    S = sketch(N, PHI5)
    st["Bphi5_jac0"] = jac(S[0], N, nu, "Bphi5_0")
    st["Bphi5_jacsum"] = sum(jac(x, N, nu, "Bphi5") for x in S)
    st["Bphi5_parsum"] = sum(_par(x, N) for x in S)
    S = sketch(N, PHI12)
    st["Bphi12_jac0"] = jac(S[0], N, nu, "Bphi12_0")
    st["Bphi12_jacsum"] = sum(jac(x, N, nu, "Bphi12") for x in S)
    S = sketch(N, FIB)
    st["Bfib_jac0"] = jac(S[0], N, nu, "Bfib_0")
    st["Bfib_jac1"] = jac(S[1], N, nu, "Bfib_1")
    st["Bfib_par0"] = _par(S[0], N)
    G = sketch_G(N, X3M2, 1)
    st["BG1x3m2_jac0"] = jac(G[0], N, nu, "BG1_0")
    st["BG1x3m2_jac1"] = jac(G[1], N, nu, "BG1_1")
    tm["B"] = time.perf_counter() - t0

    # ---- Berlekamp-Massey lengths over a fixed prefix of 64 terms
    t0 = time.perf_counter()
    seq = [(pow(1 + c, N, N) - 1 - pow(c, N, N)) % N for c in range(2, 66)]
    st["bm_Fc"], _ = bm_mod_N(seq, N, nu, "bm_Fc")
    seq = section_values(N, 64, 3)
    L, C = bm_mod_N(seq, N, nu, "bm_coef64")
    st["bm_coef64"] = L
    st["bm_coef64_jacC"] = jac(C[L], N, nu, "bm_coef64_C") if L > 0 else 0
    e = (N - 1) // 2
    seq = [pow(c, e, N) for c in range(2, 66)]
    st["bm_euler"], _ = bm_mod_N(seq, N, nu, "bm_euler")
    tm["BM"] = time.perf_counter() - t0

    # ---- small reduced forms of Cl(-4N)
    t0 = time.perf_counter()
    st.update(form_features(N))
    tm["CL"] = time.perf_counter() - t0

    assert set(st) == set(MENU_NAMES), set(st) ^ set(MENU_NAMES)
    aux = {"nonunits": nu.events, "gcd_checks": nu.checks, "t": tm}
    aux.update(aux_cf)
    return {k: st[k] for k in MENU_NAMES}, aux
