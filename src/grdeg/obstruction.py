"""The subfield-linking obstruction (SL) for frame-less type-D families (grdeg E2).

A frame-less type-D linked module is, after Theorem-S reduction, defined over a common subfield k = Q(mu) with
[k:Q] = e >= 2 (a Weil restriction). To give it genuine linking maps phi: O -> Z/N, phi': O' -> Z/N that AGREE on the
k-structure with small-height generators on BOTH sides, one must fix a ring homomorphism O_k -> Z/N, i.e. a root of
the minimal polynomial mu_k of a generator of k modulo N. This module demonstrates:

 - link_requirement(N, mu_k): operational code (sees only N and the public field k) can realise the linking only by
   solving  mu_k(s) == 0 (mod N)  -- a root-finding-mod-N call, which it must treat as an external charged routine.
 - For e = 2 (mu_k = x^2 - D0), root-finding is square-root mod N. CORRECTION (Opus critic + codex): producing ONE
   root of a FIXED D0 is NOT known factoring-equivalent -- Rabin 1979 needs RANDOM squares. rabin_split only shows the
   weaker fact that two CRT-INEQUIVALENT roots give a proper factor by gcd. This whole module is a DEMOTED remark: the
   operative obstruction is the relative-resultant height penalty (relres.py, Proposition RR), not construction
   hardness.

OPERATIONAL/constructive functions receive only (N, public field data, seed): link_requirement, naive_link_attempt.
HARNESS/diagnostic functions (named *_harness) are factor-aware and used only to exhibit the Rabin behaviour and to
classify splits STRUCTURAL vs ACCIDENTAL; they are never called by operational code.
"""
from __future__ import annotations

import math
from fractions import Fraction as Fr
from typing import Dict, List, Optional, Tuple

from sympy import isprime


# ----------------------------------------------------------------- operational (sees only N, public k)

def eigenvector_lemma(C_matrix: List[List[int]], v: List[int], mu_k: List[int], N: int) -> Dict[str, object]:
    """Non-circular form of SL0 (Opus critic's rewrite). theta_k = a generator of k acts on the module M=span(u_j);
    C is the integer matrix with (C v)_i = phi(theta_k * u_i) in linking coordinates v_j = phi(u_j). Let phi be a ring
    map with some v_{j0} a unit mod N. Multiplicativity phi(theta_k u_i) = phi(theta_k) phi(u_i) forces C v = s v
    (mod N) with s = phi(theta_k), and the Cayley-Hamilton relation mu_k(C) = 0 gives mu_k(s) v = 0, hence
    mu_k(s) = 0 (mod N) at the unit coordinate. This verifies, for given C, v, that C v = s*v (mod N) and
    mu_k(s) = 0 (mod N): a genuine root of mu_k mod N is forced, with NO assumed O_k-structure on Z/N (not circular)."""
    rr = len(v)
    Cv = [sum(C_matrix[i][j] * v[j] for j in range(rr)) % N for i in range(rr)]
    # find a unit coordinate j0
    j0 = next((j for j in range(rr) if math.gcd(v[j], N) == 1), None)
    if j0 is None:
        return {"ok": False, "reason": "no unit coordinate (phi(u_j0) must be a unit mod N)"}
    s = (Cv[j0] * pow(v[j0], -1, N)) % N
    is_eigen = all((Cv[i] - s * v[i]) % N == 0 for i in range(rr))
    mu_s = 0
    for c in reversed(mu_k):
        mu_s = (mu_s * s + c) % N
    return {"ok": True, "s": s, "C v == s v (mod N)": is_eigen, "mu_k(s) mod N": mu_s, "root_forced": (mu_s == 0)}


def link_requirement(N: int, mu_k: List[int]) -> Dict[str, object]:
    """Given only N and the minimal polynomial mu_k (low->high) of a generator of the common subfield k,
    state what the linking of a type-D family over k requires. e = deg mu_k.

    Returns the required mod-N congruence and the hardness classification. Does NOT attempt to solve it
    (operational code may not factor); it names the charged sub-routine needed.
    """
    e = len(mu_k) - 1
    assert e >= 2, "type-D obstruction is for a nontrivial subfield e>=2"
    hardness = {
        2: "square-root mod N (CONJECTURED hard for a single fixed D0; Rabin 1979 only shows that TWO "
           "CRT-inequivalent roots split N -- not that producing one fixed root is factoring-equivalent)",
    }.get(e, f"root of degree-{e} polynomial mod N (conjectured hard; RSA-type, not known >= factoring)")
    return {
        "subfield_degree_e": e,
        "required_congruence": f"find s in Z/N with mu_k(s) == 0 (mod N)",
        "reason": "phi(mu_generator) must equal phi'(mu_generator) with small-height generators on both sides, "
                  "so an O_k -> Z/N ring map (a root of mu_k mod N) is needed",
        "charged_subroutine": "root_mod_N(mu_k, N)",
        "hardness": hardness,
    }


def naive_link_attempt(N: int, mu_k: List[int], search_bound: int = 10 ** 6) -> Dict[str, object]:
    """Operational attempt to find a root of mu_k mod N WITHOUT factoring: brute search up to search_bound.
    This is the honest cost an operational algorithm pays if it has no sqrt-mod-N oracle. Returns the found root
    (if any) and the number of trials. For cryptographic N this never succeeds below Theta(sqrt smallest prime)."""
    def ev(s):
        acc = 0
        for c in reversed(mu_k):
            acc = (acc * s + c) % N
        return acc % N
    trials = 0
    for s in range(min(search_bound, N)):
        trials += 1
        if ev(s) == 0:
            return {"found": True, "root": s, "trials": trials}
    return {"found": False, "root": None, "trials": trials}


# ----------------------------------------------------------------- harness (factor-aware; exhibits Rabin)

def _sqrt_mod_p(a: int, p: int) -> Optional[int]:
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def rabin_split_harness(D0: int, p: int, q: int) -> Dict[str, object]:
    """Exhibit Rabin's reduction: pick two CRT-inequivalent square roots of D0 mod N=pq; their difference's gcd with
    N is a proper factor (STRUCTURAL). The 'diagonal' pair gives a trivial gcd (no split). Classifies each gcd hit."""
    N = p * q
    rp = _sqrt_mod_p(D0, p)
    rq = _sqrt_mod_p(D0, q)
    if rp is None or rq is None:
        return {"splittable": False, "reason": "D0 not a QR mod p or q (k inert at a prime)"}

    def crt(xp, xq):
        # x == xp (p), x == xq (q)
        return (xp * q * pow(q, -1, p) + xq * p * pow(p, -1, q)) % N

    s = crt(rp, rq)         # a root
    s_mixed = crt(rp, (-rq) % q)   # CRT-inequivalent root (flip one prime)
    s_diag = (N - s) % N    # = -s: CRT-equivalent sign flip on BOTH primes
    results = []
    for label, other in (("mixed root (CRT-inequivalent)", s_mixed), ("-s (diagonal)", s_diag)):
        gpos = math.gcd((s - other) % N, N)
        gneg = math.gcd((s + other) % N, N)
        g = next((x for x in (gpos, gneg) if 1 < x < N), 1)
        if g in (p, q):
            cls = "STRUCTURAL (explained: Rabin root difference)"
        elif g == 1 or g == N:
            cls = "no split"
        else:
            cls = "ACCIDENTAL"
        results.append({"pair": label, "gcd": g, "class": cls})
    return {"splittable": True, "N": N, "root_s": s, "results": results}


def null_gcd_rate(p: int) -> float:
    """Null expectation for a random gcd hit on N=pq: ~1/p (random integer shares factor p w.p. 1/p)."""
    return 1.0 / p


if __name__ == "__main__":
    # public view: a type-D family over k=Q(sqrt2) needs a root of x^2-2 mod N
    N = 1000003 * 1000033
    print("link requirement (operational, k=Q(sqrt2)):")
    for key, v in link_requirement(N, [-2, 0, 1]).items():
        print(f"  {key}: {v}")
    print("naive brute root search (no factoring), bound 1e5:",
          naive_link_attempt(N, [-2, 0, 1], search_bound=10 ** 5))
    print("\nRabin split (harness, k=Q(sqrt2), D0=2; primes with 2 a QR, i.e. == +-1 mod 8):")
    from sympy import nextprime
    def pick(start):
        p = nextprime(start)
        while p % 8 not in (1, 7):
            p = nextprime(p)
        return p
    p, q = pick(1000000), pick(2000000)
    out = rabin_split_harness(2, p, q)
    print(f"  p={p}, q={q}, N={p*q}")
    print("  null random-gcd rate ~1/p =", null_gcd_rate(p))
    for r in out["results"]:
        print("  ", r)
