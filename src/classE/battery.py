"""The frozen statistic menu of the classE battery.  OPERATIONAL: `compute_stats(N)` sees only N.

Definitions: notes/classE_prereg.md section 2.  All statistics are exact bounded integers.

Groups
* NC  (5)  : size / residue controls (two are exactly degenerate under the epsilon floor).
* ARG (8 x 13 = 104): class-E statistics of a/N at 8 structured arguments (structure.arguments).
* FD2 (6)  : two-argument shifted d = 2 Fourier-Dedekind sums N^2 s_1(a, b; N) (poly-time, fd.py).
* FW  (7)  : reduced-form walk f^j, j = 1..32, f the prime form of the least split prime ell >= 3 of D = -4N.
* CF3 (5)  : continued fraction of sqrt(3N), 256 steps (infrastructure path of disc 12N).
* PC  (3)  : one-sided positive controls: pc2_T1024, pc3_T1024 (fire => eps = +1; src/eps/cf.py),
             pcA_fermat_T1024 (Fermat's method succeeds within 1024 steps => p/q >= 7/8 for these sizes).
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Tuple

import gmpy2

from src.eps.cf import cf_features
from src.genusleak.forms import ambiguous_split, compose, is_ambiguous, prime_form, reduce_form, small_primes
from src.newfamily.modsym import modsym_plus
from src.newfamily.recip import dedekind6, jacobi
from src.classE.fd import fd_scaled
from src.classE.structure import arguments, convergents, cf, middle_index, path_gcd_hit, remainders

ARGS = ("isqrt", "isqrt2", "cbrt2", "div3", "fermres", "euler3", "invsqrt", "hash2")
PER_ARG = ("len", "n_a1", "maxbits", "altsign", "altbits", "dmod7", "phimod5", "ms_mod3", "rem1mod4",
           "midbits", "prealt_sign", "dr1_mod7", "dr1_bits")
FD2_PAIRS = (("isqrt", "euler3"), ("cbrt2", "fermres"), ("isqrt", "invsqrt"))
FD2_STATS = ("mod7", "bits")
FW_STATS = ("fw_a1mod4", "fw_a0mod3", "fw_bitsum", "fw_bneg", "fw_cflen", "fw_minbits", "fw_sq")
CF3_STATS = ("cf3_n_a1", "cf3_maxbits", "cf3_sumbits", "cf3_Qsq", "cf3_Q1mod5")
NC_STATS = ("nc_N_mod8", "nc_N_mod256", "nc_bit8", "nc_lead12", "nc_lead16")
PC_STATS = (("pc2_T1024", "pc_eps"), ("pc3_T1024", "pc_eps"), ("pcA_fermat_T1024", "pc_tb"))

# role: "menu" counts toward the registered conclusion; "menuG" = genus-expected (Q) statistics, reported
# separately (a hit there is family Q, not a new family); "nc" = controls; pc_* = positive controls.
MENU: List[Tuple[str, str, str]] = (
    [(n, "NC", "nc") for n in NC_STATS]
    + [(f"{a}_{k}", a.upper(), "menu") for a in ARGS for k in PER_ARG]
    + [(f"fd2_{a}_{b}_{k}", "FD2", "menu") for a, b in FD2_PAIRS for k in FD2_STATS]
    + [(n, "FW", "menuG" if n in ("fw_minbits", "fw_sq") else "menu") for n in FW_STATS]
    + [(n, "CF3", "menuG" if n in ("cf3_n_a1", "cf3_maxbits", "cf3_Qsq") else "menu") for n in CF3_STATS]
    + [(n, "PC", role) for n, role in PC_STATS]
)
MENU_NAMES = [m[0] for m in MENU]


def _sgn(x) -> int:
    return (x > 0) - (x < 0)


def per_arg(a: int, N: int) -> Dict[str, int]:
    qs = cf(a, N)
    H, S = convergents(qs)
    m = middle_index(S, N)
    alt = sum(x if i % 2 == 0 else -x for i, x in enumerate(qs))
    prealt = sum(x if i % 2 == 0 else -x for i, x in enumerate(qs[:m]))
    r = len(qs)
    D = dedekind6(a, N)                         # 6 N s(a, N)
    # Rademacher Phi = -alt + 2 - (-1)^r  (newfamily DS2)
    phi = -alt + 2 - (-1) ** r
    ms = modsym_plus(a, N) * 2
    assert ms.denominator == 1
    rs = remainders(a, N)
    dr = fd_scaled(1, [a, 1], N)               # N^2 s_1(a, 1; N): shifted d = 2 (Dedekind-Rademacher type)
    return {
        "len": r, "n_a1": sum(1 for x in qs if x == 1), "maxbits": max(x.bit_length() for x in qs),
        "altsign": _sgn(alt), "altbits": abs(alt).bit_length(),
        "dmod7": D % 7, "phimod5": phi % 5, "ms_mod3": int(ms) % 3,
        "rem1mod4": sum(1 for x in rs if x % 4 == 1),
        "midbits": qs[m].bit_length() if m < r else 0,
        "prealt_sign": _sgn(prealt),
        "dr1_mod7": dr % 7, "dr1_bits": abs(dr).bit_length(),
    }


def form_walk(N: int, J: int = 32, splits: list = None) -> Dict[str, int]:
    D = -4 * N
    ell = None
    for l in small_primes(2000)[1:]:          # odd primes
        if N % l == 0:
            continue
        if jacobi(D, l) == 1:
            ell = l
            break
    f = prime_form(D, ell)
    g = f
    st = {k: 0 for k in FW_STATS}
    minbits = None
    for _ in range(J):
        a, b, c = g
        st["fw_a1mod4"] += a % 4 == 1
        st["fw_a0mod3"] += a % 3 == 0
        st["fw_bitsum"] += a.bit_length()
        st["fw_bneg"] += b < 0
        st["fw_cflen"] += len(cf(abs(b), 2 * a)) if b % (2 * a) else 0
        mb = a.bit_length()
        minbits = mb if minbits is None else min(minbits, mb)
        st["fw_sq"] += gmpy2.is_square(a)
        if splits is not None and is_ambiguous(g):   # rule 3: an ambiguous form on the walk may split N
            h = ambiguous_split(N, g)
            if h is not None:
                splits.append(("form_walk", h))
        g = reduce_form(compose(g, f))
    st["fw_minbits"] = minbits
    return st


def cf3(N: int, T: int = 256) -> Dict[str, int]:
    M = 3 * N
    a0 = int(gmpy2.isqrt(M))
    P, Q, a = 0, 1, a0
    st = {k: 0 for k in CF3_STATS}
    for i in range(1, T + 1):
        P = a * Q - P
        Q = (M - P * P) // Q
        if Q == 0:
            break
        a = (a0 + P) // Q
        st["cf3_n_a1"] += a == 1
        st["cf3_maxbits"] = max(st["cf3_maxbits"], a.bit_length())
        st["cf3_sumbits"] += a.bit_length()
        st["cf3_Qsq"] += bool(gmpy2.is_square(Q))
        st["cf3_Q1mod5"] += Q % 5 == 1
    return st


def fermat_pc(N: int, T: int = 1024) -> int:
    x = int(gmpy2.isqrt(N))
    if x * x < N:
        x += 1
    for _ in range(T):
        y2 = x * x - N
        if gmpy2.is_square(y2):
            y = int(gmpy2.isqrt(y2))
            return int(1 < x - y < N)
        x += 1
    return 0


def compute_stats(N: int) -> Tuple[Dict[str, int], Dict]:
    out: Dict[str, int] = {}
    nonunits = []
    path_hits = []
    n = N.bit_length()
    out["nc_N_mod8"] = N % 8
    out["nc_N_mod256"] = N % 256
    out["nc_bit8"] = (N >> 8) & 1
    out["nc_lead12"] = N >> (n - 12)
    out["nc_lead16"] = N >> (n - 16)
    args = arguments(N)
    for name in ARGS:
        a = args[name]
        g = gcd(a, N)
        if a <= 0 or g != 1:
            if 1 < g < N:
                nonunits.append((name, g))
            for k in PER_ARG:
                out[f"{name}_{k}"] = 0
            continue
        h = path_gcd_hit(a, N)                 # rule 3: any factor revealed along the Euclid path
        if h is not None:
            path_hits.append((name, h))
        for k, v in per_arg(a, N).items():
            out[f"{name}_{k}"] = v
    for a_name, b_name in FD2_PAIRS:
        a, b = args[a_name], args[b_name]
        if a > 0 and b > 0 and gcd(a * b, N) == 1:
            v = fd_scaled(1, [a, b], N)
            out[f"fd2_{a_name}_{b_name}_mod7"] = v % 7
            out[f"fd2_{a_name}_{b_name}_bits"] = abs(v).bit_length()
        else:
            out[f"fd2_{a_name}_{b_name}_mod7"] = 0
            out[f"fd2_{a_name}_{b_name}_bits"] = 0
    out.update(form_walk(N, splits=path_hits))
    out.update(cf3(N))
    f = cf_features(N, T_max=1024, T_list=(1024,))
    out["pc2_T1024"] = int(f["pc2_T1024"])
    out["pc3_T1024"] = int(f["pc3_T1024"])
    out["pcA_fermat_T1024"] = fermat_pc(N)
    return out, {"nonunits": nonunits, "path_hits": path_hits, "cf_nonunit": f.get("cf_nonunit", 0)}
