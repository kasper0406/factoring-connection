"""Coverage geometry for DCL-type designs (HARNESS / DIAGNOSTIC: uses p, q).

Notation (notes/dchworst_results.md section 1).
  * A weight is a primitive integer vector (A, B) with B >= 1, A != 0.
    - positive weights: 1 <= A, B <= H, (A, B) != (1, 1)   (Theorem D' set V(H))
    - signed weights:   1 <= |A|, B <= H, (A, B) != (1, 1)  (McKee-Pinch +/- form)
    The height of (A, B) is max(|A|, B).
  * For a prime l not dividing pq put t_l = p * q^{-1} (mod l).  Then
        l | Bp - Aq   <=>   A * B^{-1} == t_l (mod l).
  * mu_l(t) = least height of a weight with A B^{-1} == t (mod l)  (the
    first minimum of the lattice Lambda_l = {(A,B) : l | Bp - Aq} restricted to
    the weight set).  Table: hmin_table(l, Hcap, signed); entries > Hcap are
    stored as Hcap + 1 (censored).
  * For a family L of primes, h_L(p, q) = min_{l in L} mu_l(t_l).
    (p, q) is covered by (L, H)  <=>  h_L(p, q) <= H.

Exact integer arithmetic throughout (numpy int64 with operands < 2^31 * 2^31).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from src.dcltheorem.params import primes_upto


# ---------------------------------------------------------------------------
# weights and per-prime tables
# ---------------------------------------------------------------------------
def weights(H: int, signed: bool) -> List[Tuple[int, int]]:
    out = []
    for B in range(1, H + 1):
        for a in range(1, H + 1):
            if gcd(a, B) != 1:
                continue
            if not (a == 1 and B == 1):
                out.append((a, B))
            if signed:
                out.append((-a, B))
    return out


def _weight_arrays(Hcap: int, signed: bool):
    W = weights(Hcap, signed)
    A = np.array([w[0] for w in W], dtype=np.int64)
    B = np.array([w[1] for w in W], dtype=np.int64)
    h = np.maximum(np.abs(A), B).astype(np.int64)
    return A, B, h


_WCACHE: Dict[Tuple[int, bool], tuple] = {}


def weight_arrays(Hcap: int, signed: bool):
    key = (Hcap, signed)
    if key not in _WCACHE:
        _WCACHE[key] = _weight_arrays(Hcap, signed)
    return _WCACHE[key]


def inverse_table(l: int) -> np.ndarray:
    """inv[x] = x^{-1} mod l for 1 <= x < l (inv[0] = 0); l prime < 2^31."""
    x = np.arange(l, dtype=np.int64)
    r = np.ones(l, dtype=np.int64)
    e = l - 2
    b = x.copy()
    while e:
        if e & 1:
            r = r * b % l
        b = b * b % l
        e >>= 1
    r[0] = 0
    return r


def hmin_table(l: int, Hcap: int, signed: bool) -> np.ndarray:
    """mu_l(t) for every t mod l (censored at Hcap + 1).  Requires Hcap < l."""
    assert Hcap < l
    A, B, h = weight_arrays(Hcap, signed)
    invB = np.array([0] + [pow(b, -1, l) for b in range(1, Hcap + 1)], dtype=np.int64)
    idx = (A % l) * invB[B] % l
    tab = np.full(l, Hcap + 1, dtype=np.int64)
    np.minimum.at(tab, idx, h)
    return tab


def nu_counts(tab: np.ndarray, Hcap: int) -> np.ndarray:
    """nu[H] = number of nonzero classes t with mu_l(t) <= H, H = 0..Hcap."""
    c = np.bincount(tab[1:], minlength=Hcap + 2)[: Hcap + 1]
    return np.cumsum(c)


# ---------------------------------------------------------------------------
# populations
# ---------------------------------------------------------------------------
def all_pairs(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Every (p, q) in S_n = {p < q <= 2p primes, 2^(n-1) <= pq < 2^n} (exhaustive)."""
    lo = 1 << (n - 1)
    hi = 1 << n
    pr = np.array(primes_upto(isqrt(hi) * 2 + 2), dtype=np.int64)
    P_list, Q_list = [], []
    plo = isqrt(lo // 2)
    for p in pr[(pr >= max(3, plo - 2)) & (pr * pr < hi)]:
        p = int(p)
        qmin = max(p + 1, -(-lo // p))
        qmax = min(2 * p, (hi - 1) // p)
        if qmin > qmax:
            continue
        i0 = np.searchsorted(pr, qmin, side="left")
        i1 = np.searchsorted(pr, qmax, side="right")
        if i1 > i0:
            qs = pr[i0:i1]
            P_list.append(np.full(len(qs), p, dtype=np.int64))
            Q_list.append(qs)
    return np.concatenate(P_list), np.concatenate(Q_list)


# ---------------------------------------------------------------------------
# the core sweep
# ---------------------------------------------------------------------------
def t_residues(P: np.ndarray, Q: np.ndarray, l: int, inv: np.ndarray = None) -> np.ndarray:
    if inv is None:
        inv = inverse_table(l)
    return (P % l) * inv[Q % l] % l


def sweep(P: np.ndarray, Q: np.ndarray, primes: Sequence[int], Hcap: int, signed: bool,
          record_every: Iterable[int] = ()) -> Dict:
    """Running min over the family `primes` (in the given order).

    Returns h (per pair, censored at Hcap+1), and for each prefix length J in
    record_every: worst-case h*_J, the histogram of h_J, and the Poisson
    prediction Pr[h_J > H] = prod_{l in prefix} (1 - nu_l(H)/(l-1)).
    Primes dividing p or q for some pair are an error (callers use l < p)."""
    rec = set(record_every)
    run = np.full(len(P), Hcap + 1, dtype=np.int64)
    logsurv = np.zeros(Hcap + 1)          # sum log(1 - nu/(l-1)) per H
    out = {"prefix": []}
    for j, l in enumerate(primes, start=1):
        tab = hmin_table(l, Hcap, signed)
        nu = nu_counts(tab, Hcap)
        with np.errstate(divide="ignore"):
            logsurv += np.log(np.maximum(1e-300, 1.0 - nu / (l - 1)))
        t = t_residues(P, Q, l)
        np.minimum(run, tab[t], out=run)
        if j in rec:
            hist = np.bincount(run, minlength=Hcap + 2)
            out["prefix"].append({
                "J": j, "l_max": int(l), "hstar": int(run.max()),
                "argmax": int(np.argmax(run)),
                "hist": hist.tolist(),
                "poisson_surv": np.exp(logsurv).tolist(),
            })
    out["h"] = run
    return out


def covering_events(p: int, q: int, primes: Sequence[int], H: int, signed: bool):
    """Brute force list of (A, B, l) with l | Bp - Aq (for tests / spot checks)."""
    W = weights(H, signed)
    return [(A, B, l) for l in primes for (A, B) in W if (B * p - A * q) % l == 0]


def band(Y: int, top: int = None, mod4: bool = True) -> List[int]:
    top = 2 * Y if top is None else top
    return [l for l in primes_upto(top) if l > Y and l > 2 and (not mod4 or l % 4 == 3)]
