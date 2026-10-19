"""Algorithms B and D: residues from the smooth part of N^2 - 1, then Coppersmith.

OPERATIONAL. Inputs: N and public parameters (smoothness bound B, lower size bound
on the target modulus, a size band for p). Never uses the factors; no primality test.

Algorithm B (Theorem S). Promise: an integer g >= Dmin with P+(g) <= B divides
p^2 - 1 and N^2 - 1 (e.g. g | gcd(p^2-1, q^2-1) for N = pq).
  1. Pollard–Strassen (block products over Z/M) finds every prime <= B dividing
     M = N - 1 and M = N + 1, with multiplicity: H = B-smooth part of N^2 - 1.
  2. For every divisor d | H with d >= Dmin (largest first) and every r mod d with
     r^2 == 1 (mod d): Coppersmith with known residue p == r (mod d), subdividing
     the band into rigorous chunks.
Algorithm D (Corollary W). g = gcd(p-1, q-1) arbitrary (not smooth): dovetail
B = 2^t, budget 2^{t/2} chunks per divisor, residue r = 1 (or all r^2 == 1).
"""
from __future__ import annotations

from itertools import product as iproduct
from math import gcd, isqrt
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from src.common.metrics import Metrics, NullMetrics
from src.specialfam.copper import best_params, n_chunks, subdivide_search
from src.specialfam.normbox import small_trial
from src.specialfam.polyprod import block_products


def smooth_part_ps(M: int, B: int, mt: Optional[Metrics] = None) -> Dict[int, int]:
    """{ell: v_ell(M)} for every prime ell <= B dividing M, by Pollard–Strassen.

    Block products prod_{j<m}(2 + m v + j) mod M by one remainder tree (polyprod),
    then in-order block gcds; the first x with gcd(x, M_cur) > 1 is the least prime
    factor of M_cur in range, so no primality test is needed. Cost O~(sqrt B).
    """
    mt = mt if mt is not None else NullMetrics()
    M = abs(M)
    out: Dict[int, int] = {}
    if M <= 1 or B < 2:
        return out
    # powers of 2 directly (cheap, and keeps M odd)
    while M % 2 == 0:
        out[2] = out.get(2, 0) + 1
        M //= 2
    if B < 3 or M == 1:
        return out
    lo = 3
    span = B - lo + 1
    m = isqrt(span) + 1
    V = -(-span // m)
    mt.inc("ps_units", m + V)
    if M == 1:
        return out
    vals = block_products(M, list(range(m)), [[lo + m * v, 1] for v in range(V)], mt)
    Mcur = M
    for v, val in enumerate(vals):
        if Mcur == 1:
            break
        mt.inc("gcd")
        if gcd(val % Mcur, Mcur) == 1:
            continue
        for j in range(m):
            x = lo + m * v + j
            if x > B or Mcur == 1:
                break
            mt.inc("gcd")
            if gcd(x, Mcur) > 1:  # x is the least prime factor of Mcur in [lo, B]
                e = 0
                while Mcur % x == 0:
                    Mcur //= x
                    e += 1
                out[x] = out.get(x, 0) + e
                mt.inc("ps_scan", m)
    return out


def merge(a: Dict[int, int], b: Dict[int, int]) -> Dict[int, int]:
    c = dict(a)
    for k, v in b.items():
        c[k] = c.get(k, 0) + v
    return c


def divisors(fac: Dict[int, int]) -> List[Tuple[int, Dict[int, int]]]:
    items = sorted(fac.items())
    out = [(1, {})]
    for ell, e in items:
        nxt = []
        for d, f in out:
            pe = 1
            for k in range(e + 1):
                g = dict(f)
                if k:
                    g[ell] = k
                nxt.append((d * pe, g))
                pe *= ell
        out = nxt
    return out


def unit_roots(fac: Dict[int, int]) -> List[int]:
    """All r mod d (d = prod ell^e) with r^2 == 1 (mod d), by CRT."""
    mods, opts = [], []
    for ell, e in sorted(fac.items()):
        q = ell ** e
        if ell == 2:
            if e == 1:
                o = [1]
            elif e == 2:
                o = [1, 3]
            else:
                o = [1, q - 1, q // 2 - 1, q // 2 + 1]
        else:
            o = [1, q - 1]
        mods.append(q)
        opts.append(o)
    d = 1
    for q in mods:
        d *= q
    res = []
    for combo in iproduct(*opts):
        r = 0
        for q, a in zip(mods, combo):
            Q = d // q
            r = (r + a * Q * pow(Q, -1, q)) % d
        res.append(r)
    return sorted(set(res)) if d > 1 else [0]


def balanced_band(N: int) -> Tuple[int, int]:
    """Smaller factor of a balanced semiprime (q <= 2p): p in [isqrt(N//2), isqrt(N)]."""
    return isqrt(N // 2), isqrt(N)


def smoothcopp_factor(N: int, B: int, Dmin: int, band: Tuple[int, int],
                      Dmax_lat: int = 9, residues: str = "pm",
                      metrics: Optional[Metrics] = None, max_total_chunks: Optional[int] = None) -> Dict:
    """Algorithm B. residues: 'pm' (r^2 == 1, from N^2-1) or 'one' (r = 1, from N-1)."""
    mt = metrics if metrics is not None else NullMetrics()
    t = small_trial(N, mt)
    if t:
        return {"factor": t, "how": "trial"}
    fac = smooth_part_ps(N - 1, B, mt)
    if residues == "pm":
        fac = merge(fac, smooth_part_ps(N + 1, B, mt))
    P_lo, P_hi = band
    params = best_params(N, max(P_lo, 2), Dmax_lat)
    mt.set("lattice", list(params[:2]))
    Y = params[2]
    # Lemma W: g >= Dmin with P+(g) <= B has a divisor in [Dmin, B*Dmin); only those are tried
    divs = [(d, f) for d, f in divisors(fac) if Dmin <= d < Dmin * max(B, 2)]
    divs.sort(key=lambda t: -t[0])
    total = 0
    pairs = 0
    for d, f in divs:
        roots = [1] if residues == "one" else unit_roots(f)
        for r in roots:
            pairs += 1
            left = None if max_total_chunks is None else max_total_chunks - total
            if left is not None and left <= 0:
                return {"factor": None, "chunks": total, "pairs": pairs, "H": fac, "how": "budget"}
            fct, used = subdivide_search(N, d, r, P_lo, P_hi, P_lo, params, mt, max_chunks=left)
            total += used
            if fct:
                mt.inc("chunks", total)
                return {"factor": fct, "d": d, "r": r, "chunks": total, "pairs": pairs,
                        "n_divisors": len(divs), "H": fac, "how": "coppersmith", "Y": Y}
    mt.inc("chunks", total)
    return {"factor": None, "chunks": total, "pairs": pairs, "n_divisors": len(divs), "H": fac}


def common_factor_factor(N: int, band: Tuple[int, int], residues: str = "one",
                         Dmax_lat: int = 9, tmax: Optional[int] = None,
                         metrics: Optional[Metrics] = None) -> Dict:
    """Algorithm D: dovetail over B = 2^t with chunk budget 2^{ceil(t/2)} per divisor.

    Finds p when g = gcd(p-1, q-1) (residues='one') or g | gcd(p^2-1, q^2-1)
    (residues='pm') is large, with cost N^{max(gamma/2, beta(1-beta)-gamma)+o(1)}.
    """
    mt = metrics if metrics is not None else NullMetrics()
    t0 = small_trial(N, mt)
    if t0:
        return {"factor": t0, "how": "trial"}
    P_lo, P_hi = band
    params = best_params(N, max(P_lo, 2), Dmax_lat)
    Y = params[2]
    tried = set()
    total_chunks = 0
    ps_units = 0
    t = 2
    tmax = tmax if tmax is not None else (N.bit_length() // 2 + 2)
    while t <= tmax:
        B = 1 << t
        budget = 1 << ((t + 1) // 2)
        sub = Metrics()
        fac = smooth_part_ps(N - 1, B, sub)
        if residues == "pm":
            fac = merge(fac, smooth_part_ps(N + 1, B, sub))
        ps_units += sub.get("ps_units")
        mt.merge(sub)
        for d, f in sorted(divisors(fac), key=lambda x: -x[0]):
            if d in tried or d < 2:
                continue
            if n_chunks(d, P_lo, P_hi, Y) > budget:
                continue
            tried.add(d)
            roots = [1] if residues == "one" else unit_roots(f)
            for r in roots:
                fct, used = subdivide_search(N, d, r, P_lo, P_hi, P_lo, params, mt)
                total_chunks += used
                if fct:
                    return {"factor": fct, "d": d, "t": t, "chunks": total_chunks,
                            "ps_units": ps_units, "how": "coppersmith"}
        t += 1
    return {"factor": None, "chunks": total_chunks, "ps_units": ps_units}
