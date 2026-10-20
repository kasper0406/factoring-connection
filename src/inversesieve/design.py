"""DCL designs and their exact forbidden residue sets.  OPERATIONAL-SAFE: functions here see only N
and public design parameters (no factors).

Notation (notes/inversesieve_results.md section 1).
  * band(Y): primes l = 3 (mod 4) in (Y, 2Y].
  * A weight is a primitive (A, B), B >= 1, A != 0, height max(|A|, B).
      kind 'signed'   : 1 <= |A|, B <= H, gcd = 1, INCLUDING (1,1) and (-1,1)  (Theorem R set; symmetric A -> -A)
      kind 'signed_dw': as 'signed' but WITHOUT (1,1) (the dchworst e1 convention)
      kind 'positive' : 1 <= A, B <= H, A != B, gcd = 1 (Theorem D' set V(H))
  * F_V(l) = {A * B^{-1} mod l : (A, B) in V}            (the Farey slope set)
  * R_l(N) = {x mod l : B x^2 = A N (mod l) for some (A,B) in V}  (classes of candidates covered at l)
    For l not dividing N:  x in R_l(N)  <=>  x != 0 and x^2 N^{-1} in F_V(l).
  * A candidate x (integer) is covered at l iff x mod l in R_l(N).  For the true factor p this is
    l | Bp - Aq (Bp^2 - AN = p (Bp - Aq), l not dividing p).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Sequence, Tuple

import numpy as np

from src.dcltheorem.params import primes_upto


def band(Y: int, top: int = None) -> List[int]:
    top = 2 * Y if top is None else top
    return [l for l in primes_upto(top) if l > Y and l % 4 == 3]


def weights(H: int, kind: str = "signed") -> List[Tuple[int, int]]:
    out = []
    for B in range(1, H + 1):
        for a in range(1, H + 1):
            if gcd(a, B) != 1:
                continue
            if kind == "positive":
                if a != B:
                    out.append((a, B))
                continue
            if kind == "signed":
                out.append((a, B)); out.append((-a, B))
            elif kind == "signed_dw":
                if not (a == 1 and B == 1):
                    out.append((a, B))
                out.append((-a, B))
            else:
                raise ValueError(kind)
    return out


def farey_set(l: int, V: Sequence[Tuple[int, int]]) -> np.ndarray:
    """Sorted array of the distinct residues A * B^{-1} mod l."""
    s = {(A * pow(B, -1, l)) % l for (A, B) in V}
    return np.array(sorted(s), dtype=np.int64)


def sqrt_mod_3mod4(a: int, l: int):
    """Square roots of a mod l (l = 3 mod 4 prime). Returns () if a is a non-residue."""
    a %= l
    if a == 0:
        return (0,)
    if pow(a, (l - 1) // 2, l) != 1:
        return ()
    r = pow(a, (l + 1) // 4, l)
    return (r, l - r) if r != l - r else (r,)


def forbidden_set(N: int, l: int, V: Sequence[Tuple[int, int]]) -> np.ndarray:
    """R_l(N): residues x mod l with B x^2 = A N (mod l) for some weight.  Requires l prime, l = 3 mod 4,
    l not dividing N, l not dividing any B.  Exact; uses only N."""
    if N % l == 0:
        raise ValueError("l | N")
    out = set()
    for (A, B) in V:
        a = (A % l) * N % l * pow(B, -1, l) % l
        for r in sqrt_mod_3mod4(a, l):
            out.add(r)
    return np.array(sorted(out), dtype=np.int64)


def forbidden_mask(N: int, l: int, V) -> np.ndarray:
    m = np.zeros(l, dtype=bool)
    m[forbidden_set(N, l, V)] = True
    return m


def candidate_interval(N: int) -> Tuple[int, int]:
    """[a0, b0] = [ceil(sqrt(N/2)), floor(sqrt(N))]: contains p for every (p,q) with p < q <= 2p."""
    a0 = isqrt(N // 2)
    while 2 * a0 * a0 < N:
        a0 += 1
    return a0, isqrt(N)


def local_density(N: int, primes: Sequence[int], V) -> float:
    """Independent-model density prod_l (1 - |R_l(N)|/l)  (float only for reporting)."""
    d = 1.0
    for l in primes:
        d *= 1.0 - len(forbidden_set(N, l, V)) / l
    return d


def lam(primes: Sequence[int], nV: int) -> float:
    """lambda = |V| * sum 1/l (expected covering events of a random pair; report only)."""
    return nV * sum(1.0 / l for l in primes)
