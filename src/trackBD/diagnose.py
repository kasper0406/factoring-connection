"""Track B/D DIAGNOSTIC harness (agent sketchBD).  KNOWS THE FACTORS.

Never imported by operational code (src/trackBD/sketch.py, linalg.py).
Uses flint factorisation of polynomials over F_p for *known* primes p only.

Main entry points
-----------------
* prime_data(N, factors)                 per-prime (p, k=v_p(N), m=N/p^k)
* components(p, f)                       irreducible factorisation of f mod p
* frob_check(p, g, m)                    Frobenius-coincidence test for alpha = X mod g:
                                         exists i<deg g with alpha^m = alpha^(p^i)
                                         and (1+alpha)^m = (1+alpha)^(p^i)
* component_vanishing(N, factors, f, S)  which components of S vanish mod each p,
                                         each labelled FROB / SPORADIC / TRIVIAL / RAMIFIED
* classify_coef_hit, classify_euclid_hit event classifiers (see notes, Section B2)

Criterion being tested (Theorem B1 in notes/trackBD_analysis.md):
  over F_p, F_N(X) = F_m(X)^(p^k)   with  F_m = (1+X)^m - 1 - X^m,  N = p^k m, p !| m,
so for a root alpha of f mod p:  S(alpha) = 0  <=>  (1+alpha)^m = 1 + alpha^m.
"""
from __future__ import annotations

from math import comb, gcd
from typing import Dict, List, Optional, Sequence, Tuple

import flint


def prime_data(N: int, factors: Sequence[int]) -> List[Tuple[int, int, int]]:
    out = []
    for p in sorted(set(factors)):
        k = 0
        M = N
        while M % p == 0:
            M //= p
            k += 1
        out.append((p, k, M))
    return out


def poly_disc_int(f: Sequence[int]) -> int:
    """Discriminant of an integer polynomial (via flint fmpz_poly)."""
    P = flint.fmpz_poly(list(f))
    try:
        return int(P.discriminant())
    except Exception:  # pragma: no cover
        import sympy
        X = sympy.Symbol("X")
        return int(sympy.discriminant(sum(int(c) * X ** i for i, c in enumerate(f)), X))


def public_integers(f: Sequence[int]) -> Dict[str, int]:
    """Integers computable from the public f whose gcd with N is 'trial division'."""
    f0 = int(f[0])
    fm1 = sum(int(c) * (-1) ** i for i, c in enumerate(f))
    # Res(f, X^2+X+1): roots of Phi_3 are universal roots of F_m (Cauchy-Mirimanoff)
    res3 = int(flint.fmpz_poly(list(f)).resultant(flint.fmpz_poly([1, 1, 1])))
    return {"f0": f0, "fm1": fm1, "disc": poly_disc_int(f), "res_phi3": res3}


def td_divides(p: int, pub: Dict[str, int]) -> List[str]:
    return [k for k, v in pub.items() if v != 0 and v % p == 0]


def components(p: int, f: Sequence[int]):
    """[(g, e)] irreducible monic factors of f mod p (harness: p known)."""
    P = flint.nmod_poly([int(c) % p for c in f], p)
    _, facs = P.factor()
    return [(g, e) for g, e in facs]


def is_universal_root_factor(p: int, g) -> bool:
    """g | X (X+1) (X^2+X+1) mod p: roots of F_m for every m prime to 6 (Cauchy-Mirimanoff
    universal factors; alpha = 0, -1 or a primitive cube root of unity)."""
    U = flint.nmod_poly([0, 1], p) * flint.nmod_poly([1, 1], p) * flint.nmod_poly([1, 1, 1], p)
    return (U % g).is_zero()


def frob_check(p: int, g, m: int) -> Dict:
    """Frobenius coincidence for alpha = X mod g in F_p[X]/(g)."""
    d = g.degree()
    x = flint.nmod_poly([0, 1], p) % g
    one = flint.nmod_poly([1], p)
    y = (x + one) % g
    am = x.pow_mod(m, g) if d > 0 else x
    bm = y.pow_mod(m, g)
    ai, bi = x, y
    for i in range(d):
        if am == ai and bm == bi:
            return {"frob": True, "i": i, "d": d}
        ai = ai.pow_mod(p, g)
        bi = bi.pow_mod(p, g)
    return {"frob": False, "i": None, "d": d}


def H_vanishes(p: int, g, m: int) -> bool:
    """(1+alpha)^m == 1 + alpha^m in F_p[X]/(g)."""
    x = flint.nmod_poly([0, 1], p) % g
    one = flint.nmod_poly([1], p)
    lhs = (x + one).pow_mod(m, g)
    rhs = (one + x.pow_mod(m, g)) % g
    return lhs == rhs


def component_vanishing(N: int, factors: Sequence[int], f: Sequence[int], S: Sequence[int],
                        pub: Optional[Dict[str, int]] = None) -> Dict[int, Dict]:
    """For each prime p | N: full vanishing, root-level vanishing per component, labels.

    Label of a vanishing component g (S == 0 mod (p, g^e) / or g | S for root level):
      TRIVIAL   alpha = 0 or alpha = -1 (S(0)=0, S(-1)=0 identically for odd N)
      RAMIFIED  e > 1 (p | disc f)
      PRIMEPOW  m == 1 (N is a power of p: F_N == 0 mod p identically)
      FROB      Frobenius coincidence holds for alpha
      SPORADIC  H_m(alpha) = 0 without Frobenius coincidence
    """
    out = {}
    pub = pub if pub is not None else public_integers(f)
    for p, k, m in prime_data(N, factors):
        Sp = flint.nmod_poly([int(c) % p for c in S], p)
        comps = []
        for g, e in components(p, f):
            root_v = (Sp % g).is_zero()
            full_v = (Sp % (g ** e)).is_zero() if e > 1 else root_v
            lab = None
            if root_v:
                if m == 1:
                    lab = "PRIMEPOW"
                elif is_universal_root_factor(p, g):  # X, X+1, X^2+X+1 factors
                    lab = "TRIVIAL"
                elif e > 1:
                    lab = "RAMIFIED"
                else:
                    fc = frob_check(p, g, m)
                    lab = "FROB" if fc["frob"] else "SPORADIC"
            comps.append({"deg": g.degree(), "e": e, "root_v": root_v, "full_v": full_v,
                          "label": lab})
        out[p] = {"k": k, "m": m, "full": Sp.is_zero(),
                  "deg_gcd": sum(c["deg"] for c in comps if c["root_v"]),
                  "comps": comps, "td": td_divides(p, pub)}
    return out


def _support_count(p: int, k: int, m: int, r: int, c: int, j: int, N: int, cap: int = 10 ** 6) -> Optional[int]:
    """# support terms of F_N mod p folding to position j mod (X^r - c).

    Over F_p, F_N = F_m(X^(p^k)); support = {p^k t : 1<=t<=m-1, binom(m,t) != 0 mod p}.
    Returns None if m is too large to enumerate."""
    if m > cap:
        return None
    P = p ** k
    cnt = 0
    for t in range(1, m):
        if (P * t) % r == j and comb(m, t) % p:
            cnt += 1
    return cnt


def classify_prime_coef(p: int, info: Dict, j: int, f_kind: Tuple, N: int) -> str:
    """Why coefficient j of S is 0 mod p.  f_kind = ('binom', r, c) or ('gen',)."""
    if info["m"] == 1:
        return "PRIMEPOW"
    if info["td"]:
        return "TD"
    if info["full"]:
        labs = {c["label"] for c in info["comps"]}
        if labs <= {"FROB"}:
            return "FULL_FROB"
        if "SPORADIC" in labs:
            return "FULL_SPORADIC"
        return "FULL_" + "+".join(sorted(x for x in labs if x))
    if f_kind[0] == "binom":
        _, r, c = f_kind
        sc = _support_count(p, info["k"], info["m"], r, c, j, N)
        if sc == 0:
            return "SUPPORT"
    return "ACC"


def classify_euclid(p: int, info: Dict, all_info: Dict[int, Dict]) -> str:
    """Why Euclid over Z/N hit a nonunit divisible by p."""
    if info["m"] == 1:
        return "PRIMEPOW"
    if info["td"]:
        return "TD"
    others = [v["deg_gcd"] for q, v in all_info.items() if q != p]
    if any(o != info["deg_gcd"] for o in others):
        labs = {c["label"] for c in info["comps"] if c["root_v"]}
        if "SPORADIC" in labs:
            return "SPLIT_SPORADIC"
        if labs <= {"FROB"} and labs:
            return "SPLIT_FROB"
        if not labs:
            return "SPLIT_OTHER"  # p-side gcd smaller than some other prime's
        return "SPLIT_" + "+".join(sorted(x for x in labs if x))
    return "ACC"


def scalar_roots_census(p: int, m: int) -> Dict[str, int]:
    """Exact enumeration over c in F_p of (1+c)^m == 1 + c^m (F-type) and of
    (1+c)^(m-1) == c^(m-1)   [G_1-type, exponent N-1 = p m - 1 = m - 1 mod (p-1)]

    Classifies F-roots: trivial {0,-1}, FROB (c^(m-1) = 1 = (1+c)^(m-1)), SPORADIC.
    Cost O(p) modexps: use only for p <~ 2^16."""
    F_triv = F_frob = F_spor = 0
    G_roots = 0
    e = m % (p - 1)
    for c in range(p):
        lhs = pow(1 + c, m, p)
        rhs = (1 + pow(c, m, p)) % p
        if lhs == rhs:
            if c in (0, p - 1) or (c * c + c + 1) % p == 0:  # universal points incl. omega
                F_triv += 1
            elif pow(c, m - 1, p) == 1 and pow(1 + c, m - 1, p) == 1:
                F_frob += 1
            else:
                F_spor += 1
        if c not in (0, p - 1):
            # G_1(c) = (1+c)^(N-1) - c^(N-1); mod p exponent N-1 == m-1 (mod p-1) when p !| c(1+c)
            if pow(1 + c, m - 1, p) == pow(c, m - 1, p):
                G_roots += 1
    return {"F_triv": F_triv, "F_frob": F_frob, "F_spor": F_spor, "G1": G_roots,
            "g": gcd(p - 1, m - 1), "e": e}
