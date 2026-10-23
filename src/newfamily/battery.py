"""The frozen statistic menu of the newfamily battery.  OPERATIONAL: `compute_stats(N)` sees only N.

Definitions: notes/newfamily_prereg.md section 2.  All statistics are exact bounded integers.
Arguments a(N) (each a deterministic function of N; gcd(a, N) is checked first, a nonunit is recorded and the
argument's statistics are set to 0):
  sqrt  : a = isqrt(N)
  gold  : a = floor(N (sqrt5 - 1)/2)
  cbrt2 : a = floor(N^{2/3})
  ferm  : a = (isqrt(N) + 1)^2 - N             (the first Fermat residue, ~ 2 sqrt N)
  euler : a = 2^{(N-1)/2} mod N                (Euler-criterion element)
  hash  : a = 1 + (SHA-256(decimal N) mod (N - 1))   (a public pseudo-random argument)
Per argument (14): jac, dsign, dbits, dmod7, dmod5 (Dedekind sum D = 6N s(a, N)), len, n_a1, maxbits, altsign,
altbits, rem1mod4, rem0mod3 (CF of a/N and Euclid remainders), ms_sign, ms_mod3 (plus modular symbol of 11a:
m = 2 [a/N]^+ in Z).
"""
from __future__ import annotations

import hashlib
from math import gcd
from typing import Dict, List, Tuple

import gmpy2

from src.eps.cf import cf_features
from src.newfamily.modsym import modsym_plus
from src.newfamily.recip import cf_stats, dedekind6, jacobi

ARGS = ("sqrt", "gold", "cbrt2", "ferm", "euler", "hash")
PER_ARG = ("jac", "dsign", "dbits", "dmod7", "dmod5", "len", "n_a1", "maxbits", "altsign", "altbits",
           "rem1mod4", "rem0mod3", "ms_sign", "ms_mod3")

MENU: List[Tuple[str, str, str]] = (
    [(n, "NC", "nc") for n in ("nc_N_mod8", "nc_N_mod256", "nc_bit8", "nc_bit9", "nc_jacN_5", "nc_lead16",
                               "nc_cubic7", "nc_ded5_bits", "nc_ded5_sign", "nc_ms_a3")]
    + [(f"{a}_{k}", a.upper(), "menu") for a in ARGS for k in PER_ARG]
    + [(n, "PC", "pc_eps") for n in ("pc2_T1024", "pc3_T1024")]
    + [(n, "PC", "pc_qr2") for n in ("pcq2_T1024", "pcq3_T1024")]
)
MENU_NAMES = [m[0] for m in MENU]
S_BONFERRONI = 100


def _sgn(x) -> int:
    return (x > 0) - (x < 0)


def arguments(N: int) -> Dict[str, int]:
    r = int(gmpy2.isqrt(N))
    g5 = int(gmpy2.isqrt(5 * N * N))
    return {
        "sqrt": r,
        "gold": (g5 - N) // 2,
        "cbrt2": int(gmpy2.iroot(N * N, 3)[0]),
        "ferm": (r + 1) ** 2 - N,
        "euler": pow(2, (N - 1) // 2, N),
        "hash": 1 + int.from_bytes(hashlib.sha256(str(N).encode()).digest(), "big") % (N - 1),
    }


def per_arg(a: int, N: int) -> Dict[str, int]:
    D = dedekind6(a, N)
    st = cf_stats(a, N)
    m = modsym_plus(a, N) * 2
    assert m.denominator == 1
    m = int(m)
    alt = st["altsum"]
    return {
        "jac": jacobi(a, N),
        "dsign": _sgn(D), "dbits": abs(D).bit_length(), "dmod7": D % 7, "dmod5": D % 5,
        "len": st["len"], "n_a1": st["n_a1"], "maxbits": st["maxbits"],
        "altsign": _sgn(alt), "altbits": abs(alt).bit_length(),
        "rem1mod4": st["rem1mod4"], "rem0mod3": st["rem0mod3"],
        "ms_sign": _sgn(m), "ms_mod3": m % 3,
    }


def compute_stats(N: int) -> Tuple[Dict[str, int], Dict]:
    out: Dict[str, int] = {}
    nonunits = []
    n = N.bit_length()
    out["nc_N_mod8"] = N % 8
    out["nc_N_mod256"] = N % 256
    out["nc_bit8"] = (N >> 8) & 1
    out["nc_bit9"] = (N >> 9) & 1
    out["nc_jacN_5"] = jacobi(N, 5)
    out["nc_lead16"] = N >> (n - 16)
    # cubic character of N mod 7 = (1 + 3 omega / (N))_3 by Theorem RA (exponent: index of N in (Z/7)^*/cubes)
    out["nc_cubic7"] = {1: 0, 6: 0, 3: 1, 4: 1, 2: 2, 5: 2}[N % 7] if N % 7 else 0
    D5 = dedekind6(5, N) if N % 5 else 0
    out["nc_ded5_bits"] = abs(D5).bit_length()
    out["nc_ded5_sign"] = _sgn(D5)
    out["nc_ms_a3"] = int(2 * modsym_plus(3, N)) if N % 3 else 0
    for name, a in arguments(N).items():
        g = gcd(a, N)
        if g != 1:
            if 1 < g < N:
                nonunits.append((name, g))
            for k in PER_ARG:
                out[f"{name}_{k}"] = 0
            continue
        for k, v in per_arg(a, N).items():
            out[f"{name}_{k}"] = v
    f = cf_features(N, T_max=1024, T_list=(1024,))
    for k in ("pc2_T1024", "pc3_T1024", "pcq2_T1024", "pcq3_T1024"):
        out[k] = int(f[k])
    return out, {"nonunits": nonunits, "cf_nonunit": f.get("cf_nonunit", 0)}
