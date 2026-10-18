"""OPERATIONAL primality tests with leak instrumentation (owner: primality).

Every function receives only N, public parameters and a base / seed.  It never
uses p, q, phi(N), lambda(N) or any group order, and never calls a factoring
or factor-revealing primality routine.  (``gmpy2.is_square`` / ``isqrt`` are used
for the perfect-square checks that Baillie-PSW and Selfridge's method need.)

Each test returns a dict with the verdict and a list ``leaks`` of every *proper*
divisor 1 < d < N obtained as gcd(N, x) for an intermediate value x that the
test computes (or x +- 1 for the values it compares with +-1).  Leaks are classified:

* ``structural``: the compared quantity (a ring element, possibly a polynomial or
  a quadratic-ring element) vanishes modulo d as a whole.  This is the CRT event
  "vanishes mod one prime but not the other" that the analysis predicts.
* ``accidental``: a single coordinate of a multi-coordinate value is divisible
  by d although the element is not (probability ~ 1/p per coordinate), or a
  trial-division / nonunit-base hit.  Reported separately (AGENT_RULES rule 3).
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Optional, Tuple

import gmpy2

from src.common.metrics import Metrics, NullMetrics

_NULL = NullMetrics()


def v2(x: int) -> int:
    x = abs(x)
    if x == 0:
        raise ValueError("v2(0)")
    return (x & -x).bit_length() - 1


def split2(x: int) -> Tuple[int, int]:
    """x = 2^t * m with m odd."""
    t = v2(x)
    return t, x >> t


def _proper(d: int, N: int) -> bool:
    return 1 < d < N


def _leak(leaks: List[dict], N: int, x: int, kind: str, where, cls: str = "structural",
          metrics=_NULL) -> None:
    metrics.inc("gcd")
    d = gcd(x % N, N)
    if _proper(d, N):
        leaks.append({"divisor": d, "kind": kind, "where": where, "class": cls})


# ---------------------------------------------------------------- Fermat ----
def fermat(N: int, a: int, chain: bool = False, metrics=_NULL) -> dict:
    """Fermat test a^(N-1) == 1.  Observable: gcd(a^(N-1) - 1, N).

    With ``chain=True`` also gcd(a^e - 1, N) for every prefix exponent e of the
    left-to-right binary chain for N-1 (the intermediates a naive implementation holds).
    """
    leaks: List[dict] = []
    _leak(leaks, N, a, "nonunit_base", None, "accidental", metrics)
    e = N - 1
    if chain:
        x, pref = 1, 0
        for bit in bin(e)[2:]:
            x = x * x % N
            pref = 2 * pref
            metrics.inc("rn_mul")
            if bit == "1":
                x = x * a % N
                pref += 1
                metrics.inc("rn_mul")
            _leak(leaks, N, x - 1, "chain_minus1", pref, metrics=metrics)
        y = x
    else:
        y = pow(a, e, N)
        metrics.inc("rn_mul", 2 * e.bit_length())
        _leak(leaks, N, y - 1, "minus1", e, metrics=metrics)
    return {"verdict": y == 1, "residue": y, "leaks": leaks}


# ---------------------------------------------------------- Miller-Rabin ----
def miller_rabin(N: int, a: int, metrics=_NULL) -> dict:
    """Strong (Miller-Rabin) test to base a.  N-1 = 2^t m.

    Transcript b_i = a^(m 2^i), i = 0..t.  Observables: gcd(b_i - 1, N) and
    gcd(b_i + 1, N) for all i (this includes the classical 'nontrivial square
    root of 1' event).  ``pattern`` is the +-1 pattern mod N:
      'S0'  b_0 = 1;  'L<j>' b_(j-1) = -1 (j >= 1, strong liar at level j);
      'NONE' none of these (composite verdict, the overwhelmingly common case).
    """
    leaks: List[dict] = []
    _leak(leaks, N, a, "nonunit_base", None, "accidental", metrics)
    t, m = split2(N - 1)
    b = pow(a, m, N)
    metrics.inc("rn_mul", 2 * m.bit_length())
    seq = [b]
    for _ in range(t):
        b = b * b % N
        metrics.inc("rn_mul")
        seq.append(b)
    for i, x in enumerate(seq):
        _leak(leaks, N, x - 1, "minus1", i, metrics=metrics)
        _leak(leaks, N, x + 1, "plus1", i, metrics=metrics)
    if seq[0] == 1:
        pattern = "S0"
    else:
        pattern = "NONE"
        for i in range(t):
            if seq[i] == N - 1:
                pattern = f"L{i + 1}"
                break
    # classical event: first i with b_i == 1 has b_(i-1) != +-1
    sqrt_split = None
    for i in range(1, t + 1):
        if seq[i] == 1:
            if seq[i - 1] not in (1, N - 1):
                sqrt_split = gcd(seq[i - 1] - 1, N)
            break
    return {"verdict": pattern != "NONE", "t": t, "m": m, "seq": seq, "pattern": pattern,
            "sqrt_split": sqrt_split, "leaks": leaks}


# ----------------------------------------------------- Solovay-Strassen ----
def solovay_strassen(N: int, a: int, metrics=_NULL) -> dict:
    """Euler test a^((N-1)/2) == (a/N).  Observables gcd(e -+ 1, N)."""
    leaks: List[dict] = []
    _leak(leaks, N, a, "nonunit_base", None, "accidental", metrics)
    J = int(gmpy2.jacobi(a, N))
    e = pow(a, (N - 1) // 2, N)
    metrics.inc("rn_mul", 2 * N.bit_length())
    _leak(leaks, N, e - 1, "minus1", None, metrics=metrics)
    _leak(leaks, N, e + 1, "plus1", None, metrics=metrics)
    metrics.inc("gcd")
    d_res = gcd((e - J) % N, N)  # the 'residual' observable gcd(b - J, N)
    if e == 1:
        cls = "+1"
    elif e == N - 1:
        cls = "-1"
    else:
        cls = "other"
    return {"verdict": J != 0 and e == J % N, "jacobi": J, "euler": e, "euler_class": cls,
            "residual_split": d_res if _proper(d_res, N) else None, "leaks": leaks}


# ---------------------------------------------------- quadratic rings ------
class QuadRing:
    """Z/N[x]/(x^2 - P x + Q).  Elements are pairs (u0, u1) = u0 + u1 x."""

    def __init__(self, N: int, P: int, Q: int, metrics=_NULL):
        self.N, self.P, self.Q = N, P % N, Q % N
        self.m = metrics

    def mul(self, A, B):
        N, P, Q = self.N, self.P, self.Q
        a0, a1 = A
        b0, b1 = B
        self.m.inc("q_mul")
        c2 = a1 * b1
        return ((a0 * b0 - Q * c2) % N, (a0 * b1 + a1 * b0 + P * c2) % N)

    def pow(self, A, e: int):
        R = (1, 0)
        B = A
        while e:
            if e & 1:
                R = self.mul(R, B)
            e >>= 1
            if e:
                B = self.mul(B, B)
        return R


def lucas_UV(N: int, P: int, Q: int, k: int, metrics=_NULL) -> Tuple[int, int]:
    """(U_k, V_k) mod N via x^k = U_k x - Q U_(k-1) in Z/N[x]/(x^2 - P x + Q)."""
    R = QuadRing(N, P, Q, metrics)
    u0, u1 = R.pow((0, 1), k)
    return u1 % N, (2 * u0 + P * u1) % N


def selfridge(N: int, metrics=_NULL) -> dict:
    """Selfridge method A: first D in 5, -7, 9, -11, ... with (D/N) = -1; P = 1, Q = (1-D)/4.

    Returns {'D','P','Q'} or {'nonunit': d} if some gcd(D, N) is proper, or
    {'square': True} if N is a perfect square (then no such D exists).
    """
    if gmpy2.is_square(N):
        return {"square": True}
    D = 5
    while True:
        metrics.inc("jacobi")
        j = int(gmpy2.jacobi(D % N, N))
        if j == 0:
            d = gcd(abs(D), N)
            if _proper(d, N):
                return {"nonunit": d, "D": D}
        if j == -1:
            return {"D": D, "P": 1, "Q": (1 - D) // 4}
        D = -D - 2 if D > 0 else -D + 2


def strong_lucas(N: int, P: int, Q: int, metrics=_NULL) -> dict:
    """Strong Lucas test with D = P^2 - 4Q, (D/N) = -1 assumed.  N+1 = 2^t' m'.

    Observables: gcd(U_(m' 2^i), N), gcd(V_(m' 2^i), N) for i = 0..t'
    (U_(m' 2^t') = U_(N+1) is the plain Lucas test), plus gcd(Q D, N).
    U_k = 0 mod p  <=>  gamma^k = 1 and V_k = 0 mod p <=> gamma^k = -1, gamma = alpha/beta.
    Pattern: 'S0' U_m' = 0; 'L<j>' V_(m' 2^(j-1)) = 0; 'NONE'.
    """
    leaks: List[dict] = []
    D = P * P - 4 * Q
    _leak(leaks, N, Q * D, "nonunit_QD", None, "accidental", metrics)
    t, m = split2(N + 1)
    R = QuadRing(N, P, Q, metrics)
    X = R.pow((0, 1), m)  # x^m'
    Us, Vs = [], []
    for i in range(t + 1):
        u0, u1 = X
        Us.append(u1 % N)
        Vs.append((2 * u0 + P * u1) % N)
        if i < t:
            X = R.mul(X, X)
    for i in range(t + 1):
        _leak(leaks, N, Us[i], "U", i, metrics=metrics)
        _leak(leaks, N, Vs[i], "V", i, metrics=metrics)
    if Us[0] == 0:
        pattern = "S0"
    else:
        pattern = "NONE"
        for i in range(t):
            if Vs[i] == 0:
                pattern = f"L{i + 1}"
                break
    return {"verdict": pattern != "NONE", "lucas_verdict": Us[t] == 0, "t": t, "m": m,
            "U": Us, "V": Vs, "pattern": pattern, "leaks": leaks}


def bpsw(N: int, small_primes=(3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37), metrics=_NULL) -> dict:
    """Baillie-PSW: trial division by tiny primes, MR base 2, strong Lucas (Selfridge)."""
    leaks: List[dict] = []
    if N % 2 == 0:
        return {"verdict": N == 2, "leaks": []}
    for sp in small_primes:
        if N % sp == 0 and N != sp:
            leaks.append({"divisor": sp, "kind": "trial_division", "where": sp, "class": "accidental"})
            return {"verdict": False, "leaks": leaks, "stage": "trial"}
    mr = miller_rabin(N, 2, metrics)
    leaks += [dict(l, stage="MR2") for l in mr["leaks"]]
    sel = selfridge(N, metrics)
    if sel.get("square"):
        return {"verdict": False, "leaks": leaks, "stage": "square", "mr": mr["pattern"]}
    if "nonunit" in sel:
        leaks.append({"divisor": sel["nonunit"], "kind": "selfridge_jacobi0", "where": sel["D"],
                      "class": "accidental"})
        return {"verdict": False, "leaks": leaks, "stage": "selfridge", "mr": mr["pattern"]}
    lu = strong_lucas(N, sel["P"], sel["Q"], metrics)
    leaks += [dict(l, stage="SLucas") for l in lu["leaks"]]
    return {"verdict": mr["verdict"] and lu["verdict"], "leaks": leaks, "D": sel["D"],
            "mr": mr["pattern"], "lucas": lu["pattern"]}


# ------------------------------------------------- quadratic Frobenius ----
def frobenius_quadratic(N: int, Delta: int, alpha: Tuple[int, int], metrics=_NULL) -> dict:
    """Frobenius condition alpha^N == sigma(alpha) in R = Z/N[w]/(w^2 - Delta).

    For (Delta/N) = -1 this is Grantham's quadratic Frobenius condition (for
    alpha = (P + w)/2 a root of x^2 - P x + Q it is x^(N+1) == Q mod (N, f)).
    Observable z = alpha^N - sigma(alpha) = z0 + z1 w:
      structural leak  gcd(gcd(z0, z1), N) proper (z vanishes mod d);
      accidental leak  gcd(z0, N) or gcd(z1, N) proper otherwise.
    """
    leaks: List[dict] = []
    a0, a1 = alpha
    _leak(leaks, N, a0 * a0 - Delta * a1 * a1, "nonunit_norm", None, "accidental", metrics)
    R = QuadRing(N, 0, -Delta, metrics)  # w^2 = Delta
    z = R.pow((a0 % N, a1 % N), N)
    z0, z1 = (z[0] - a0) % N, (z[1] + a1) % N
    metrics.inc("gcd")
    d = gcd(gcd(z0, z1), N)
    if _proper(d, N):
        leaks.append({"divisor": d, "kind": "frob", "where": None, "class": "structural"})
    else:
        for c, x in (("z0", z0), ("z1", z1)):
            _leak(leaks, N, x, c, None, "accidental", metrics)
    return {"verdict": z0 == 0 and z1 == 0, "z": (z0, z1), "leaks": leaks}


# ------------------------------------------------------------------ AKS ----
def aks_parameters(N: int) -> Tuple[int, int]:
    """Smallest r with ord_r(N) > (log2 N)^2 (exact integer comparison), and
    l = floor(sqrt(phi(r)) * log2 N) (computed with a safe integer upper bound)."""
    import math
    L2 = math.log2(N)
    bound = int(L2 * L2)  # ord_r(N) > log^2 N  <=>  ord_r(N) >= floor(log^2)+1
    r = 2
    while True:
        if gcd(r, N) == 1:
            o, x = 1, N % r
            ok = True
            y = x
            while y != 1:
                y = y * x % r
                o += 1
                if o > bound:
                    break
            if o > bound:
                break
        r += 1
    phi_r = sum(1 for k in range(1, r) if gcd(k, r) == 1)
    ell = int(math.isqrt(phi_r) * L2) + 1
    return r, ell


def aks(N: int, r: Optional[int] = None, a_values=None, metrics=_NULL) -> dict:
    """AKS steps 3 and 5 instrumented; step 1 (perfect power) skipped for semiprimes.

    Step 3: gcd(a, N) for a <= r  (trial division; leaks are 'accidental', small factor).
    Step 5: Delta_a = (X+a)^N - X^N - a mod (N, X^r - 1), a in ``a_values``.
      structural leak: gcd(content(Delta_a), N) proper (Delta_a vanishes mod d);
      coefficient leak ('coefficient'): some single coefficient gcd proper otherwise.  NOT
      assumed to be a 1/p event: when q < r the mod-p residue is a sparse Lucas/Kummer
      pattern (see notes/primality_analysis.md, AKS section); compare with the null r*l/p.
    Uses python-flint nmod_poly (only multiplication and reduction by the monic X^r - 1).
    """
    import flint
    if r is None:
        r, ell = aks_parameters(N)
    else:
        ell = None
    if a_values is None:
        a_values = range(1, (ell or 1) + 1)
    leaks: List[dict] = []
    for a in range(2, min(r, N - 1) + 1):
        metrics.inc("gcd")
        d = gcd(a, N)
        if _proper(d, N):
            leaks.append({"divisor": d, "kind": "aks_trial", "where": a, "class": "accidental"})
            break
    if N <= r:
        return {"verdict": True, "r": r, "leaks": leaks}
    mod = flint.nmod_poly([N - 1] + [0] * (r - 1) + [1], N)  # X^r - 1
    Xn = [0] * r
    Xn[N % r] = 1
    XN = flint.nmod_poly(Xn, N)
    verdict = True
    n_coeff_hits = 0
    for a in a_values:
        base = flint.nmod_poly([a % N, 1], N)
        P = base.pow_mod(N, mod)
        metrics.inc("q_mul", 2 * N.bit_length())
        Dl = P - XN - a
        cs = [int(c) for c in Dl.coeffs()]
        if any(cs):
            verdict = False
        cont = 0
        for c in cs:
            cont = gcd(cont, c)
        metrics.inc("gcd")
        d = gcd(cont, N)
        if cs and _proper(d, N):
            leaks.append({"divisor": d, "kind": "aks_poly", "where": a, "class": "structural"})
            continue
        for i, c in enumerate(cs):
            metrics.inc("gcd")
            dd = gcd(c, N)
            if _proper(dd, N):
                n_coeff_hits += 1
                if n_coeff_hits <= 5:
                    leaks.append({"divisor": dd, "kind": "aks_coeff", "where": (a, i),
                                  "class": "coefficient"})
    return {"verdict": verdict, "r": r, "ell": ell, "leaks": leaks, "n_coeff_hits": n_coeff_hits}
