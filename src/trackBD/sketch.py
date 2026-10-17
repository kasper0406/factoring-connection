"""Track B/D operational sketches and readouts (agent sketchBD).

OPERATIONAL code: every function here receives only N and public parameters
(f, r, c, exponent offsets, seeds).  It never uses a factorisation of N.
Every inversion modulo N is gcd-checked; nonunits are returned, not swallowed.

Objects
-------
* sketch(N, f)        = F_N mod (f, N),  F_N = (1+X)^N - 1 - X^N
* sketch_G(N, f, k)   = (1+X)^(N-k) - X^(N-k) mod (f, N)   ("divided" quantities,
                        see notes/trackBD_analysis.md section D)
* lifted_derivative(N, r, c) : (1+X)^N mod (N^2, X^r - c), differentiated, /N

Readouts (each returns a list of proper divisors found, possibly empty)
* coef_gcds(coeffs, N)          gcd(N, a_j) for every coefficient
* content_gcd(coeffs, N)        gcd(N, a_0, ..., a_{r-1})
* euclid_readout(f, S, N)       Euclid / resultant of (f, S) over Z/N with
                                gcd-checked leading-coefficient inversions

Cost accounting: q_mul / q_sqr / q_mul_dimsum via src.common.metrics keys;
`gcd`, `inv` for the gcd-checked inversions; `euclid_steps`.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import FlintQuotient, NonUnit, inv_mod, power_of_linear

_NMOD_LIMIT = 1 << 62


# ---------------------------------------------------------------------------
# polynomial helpers over Z/N (lists low -> high, ints in [0, N))
# ---------------------------------------------------------------------------
def binomial_modulus(r: int, c: int) -> List[int]:
    """X^r - c as an integer coefficient list (low -> high)."""
    return [-c] + [0] * (r - 1) + [1]


def _mk(N: int, cs: Sequence[int]):
    if N < _NMOD_LIMIT:
        return flint.nmod_poly([int(x) % N for x in cs], N)
    ctx = flint.fmpz_mod_poly_ctx(N)
    return ctx([int(x) % N for x in cs])


def _coeffs(P, d: Optional[int] = None) -> List[int]:
    cs = [int(x) for x in P.coeffs()]
    if d is not None:
        cs = cs + [0] * (d - len(cs))
    return cs


def _strip(cs: List[int], N: int) -> List[int]:
    cs = [x % N for x in cs]
    while cs and cs[-1] == 0:
        cs.pop()
    return cs


# ---------------------------------------------------------------------------
# sketches
# ---------------------------------------------------------------------------
def power_sketch(N: int, f: Sequence[int], e: int, a: int = 1, b: int = 1,
                 m: Optional[Metrics] = None) -> List[int]:
    """Coefficients of (a + bX)^e mod (f, N), f monic integer list."""
    Q = FlintQuotient(N, list(f), m)
    if e == 0:
        return Q.coeffs(Q.one)
    lin = Q.reduce([a % N, b % N])
    return Q.coeffs(Q.pow(lin, e))


def sketch(N: int, f: Sequence[int], m: Optional[Metrics] = None) -> List[int]:
    """F_N mod (f, N) with F_N = (1+X)^N - 1 - X^N (repeated squaring)."""
    Q = FlintQuotient(N, list(f), m)
    x = Q.gen()
    S = Q.sub(Q.sub(power_of_linear(Q, N), Q.one), Q.pow(x, N))
    return Q.coeffs(S)


def sketch_G(N: int, f: Sequence[int], k: int = 1, m: Optional[Metrics] = None) -> List[int]:
    """G_k := (1+X)^(N-k) - X^(N-k) mod (f, N).

    D^{[k]} F_N = binom(N,k) G_k (Hasse derivative) and F_N' = N G_1, so the
    'divided' derivative quantities of the lift are unit multiples of G_k mod N
    (proved in notes/trackBD_analysis.md, Section D1)."""
    Q = FlintQuotient(N, list(f), m)
    S = Q.sub(power_of_linear(Q, N - k), Q.pow(Q.gen(), N - k))
    return Q.coeffs(S)


def scalar_F(N: int, c: int) -> int:
    """F_N(c) mod N = sketch modulo X - c (degree-1 quotient)."""
    return (pow(1 + c, N, N) - 1 - pow(c, N, N)) % N


def scalar_G(N: int, c: int, k: int = 1) -> int:
    return (pow(1 + c, N - k, N) - pow(c, N - k, N)) % N


# ---------------------------------------------------------------------------
# lifts (Track D)
# ---------------------------------------------------------------------------
def lifted_power(N: int, r: int, c: int, t: int = 2, m: Optional[Metrics] = None) -> List[int]:
    """(1+X)^N mod (N^t, X^r - c)."""
    M = N ** t
    Q = FlintQuotient(M, binomial_modulus(r, c), m)
    return Q.coeffs(power_of_linear(Q, N))


def lifted_derivative_divided(N: int, r: int, c: int, m: Optional[Metrics] = None) -> Dict[str, List[int]]:
    """Compute in the lift Z/N^2 [X] (no quotient ring) the derivative identity.

    Returns
      'G_from_lift' : (d/dX (1+X)^N  - d/dX X^N) / N  mod (N, X^r - c), where the
                      derivative is taken of the *unreduced* polynomial, realised as
                      N (1+X)^(N-1) - N X^(N-1) computed mod (N^2, X^r - c) and then
                      divided by N (exact division in Z/N^2 -> Z/N);
      'dS'          : derivative of the reduced sketch S = F_N mod (N^2, X^r-c),
                      as a polynomial of degree < r, then mod N;
      'quot'        : the sketch of the quotient Qt in F_N = Qt (X^r - c) + S,
                      recovered from the identity F_N' = Qt' (X^r-c) + r X^(r-1) Qt + S',
                      i.e. r X^(r-1) Qt == N G - S'  (mod X^r - c); requires r, c units.
    """
    M = N * N
    Q = FlintQuotient(M, binomial_modulus(r, c), m)
    x = Q.gen()
    Pn1 = power_of_linear(Q, N - 1)
    Xn1 = Q.pow(x, N - 1)
    deriv = Q.scale(N, Q.sub(Pn1, Xn1))            # N((1+X)^(N-1) - X^(N-1)) mod N^2
    dc = Q.coeffs(deriv)
    assert all(v % N == 0 for v in dc)
    G = [(v // N) % N for v in dc]
    # sketch modulo N^2 and its derivative inside the quotient
    S2 = Q.coeffs(Q.sub(Q.sub(Q.mul(Pn1, Q.reduce([1, 1])), Q.one), Q.mul(Xn1, x)))
    dS = [(j * S2[j]) % N for j in range(1, r)] + [0]
    return {"G_from_lift": G, "S_modN2": S2, "dS_modN": dS}


# ---------------------------------------------------------------------------
# readouts
# ---------------------------------------------------------------------------
def coef_gcds(coeffs: Sequence[int], N: int, m: Optional[Metrics] = None) -> List[Tuple[int, int]]:
    """[(j, d)] for every coefficient a_j with 1 < gcd(a_j, N) = d < N."""
    out = []
    for j, a in enumerate(coeffs):
        if m is not None:
            m.inc("gcd")
        d = gcd(a % N, N)
        if 1 < d < N:
            out.append((j, d))
    return out


def content_gcd(coeffs: Sequence[int], N: int) -> int:
    g = N
    for a in coeffs:
        g = gcd(g, a % N)
    return g


def euclid_readout(f: Sequence[int], S: Sequence[int], N: int,
                   m: Optional[Metrics] = None) -> Dict:
    """Euclid's algorithm on (f, S) in (Z/N)[X] with gcd-checked inversions.

    Returns dict with
      'factor'   : proper divisor found at a nonunit leading coefficient, or None
      'deg_gcd'  : degree of the computed gcd if Euclid completed (else None)
      'res'      : resultant Res(f, S) mod N if completed (else None)
      'steps'    : number of division steps
    A leading coefficient that is 0 mod N is a genuine degree drop (stripped);
    a nonzero nonunit leading coefficient is reported.  The remainder sequence
    over Z/N then coincides with the remainder sequences over every F_p | N
    up to the first nonunit (degrees agree), so a factor is reported iff the
    degree sequences mod different primes diverge or a leading coefficient
    vanishes modulo one prime by accident (classified in the harness).
    """
    m = m if m is not None else NullMetrics()
    a = _strip(list(f), N)
    b = _strip(list(S), N)
    res = 1
    steps = 0
    sign_flips = 0
    if not b:
        return {"factor": None, "deg_gcd": len(a) - 1, "res": 0, "steps": 0}
    # Res(a,b) for a monic: Res(a,b) = prod_{a(x)=0} b(x)
    # recurrence: Res(a, b) = (-1)^{da db} lc(b)^{da - dr} Res(b, r), r = a mod b
    while True:
        da, db = len(a) - 1, len(b) - 1
        if db == 0:
            res = res * pow(b[0], da, N) % N
            return {"factor": None, "deg_gcd": 0, "res": res, "steps": steps}
        lc = b[-1]
        m.inc("gcd")
        g = gcd(lc, N)
        if g != 1:
            return {"factor": g if 1 < g < N else None, "deg_gcd": None, "res": None,
                    "steps": steps, "where": "lc", "deg_b": db}
        m.inc("inv")
        inv = pow(lc, -1, N)
        bm = [(x * inv) % N for x in b]
        A, B = _mk(N, a), _mk(N, bm)
        R = A % B
        r = _strip(_coeffs(R), N)
        steps += 1
        m.inc("euclid_steps")
        if not r:
            return {"factor": None, "deg_gcd": db, "res": 0, "steps": steps}
        dr = len(r) - 1
        # Res(a, b) = (-1)^{da db} Res(b, a);  Res(b, a) = lc(b)^{da-dr} Res(b, r)
        if (da * db) % 2:
            res = -res
        res = res * pow(lc, da - dr, N) % N
        a, b = b, r


def try_inverse(f: Sequence[int], S: Sequence[int], N: int, m: Optional[Metrics] = None) -> Dict:
    """Attempt to invert S in (Z/N)[X]/(f): via Euclid; a nonunit is a factor."""
    out = euclid_readout(f, S, N, m)
    out["invertible"] = (out.get("deg_gcd") == 0)
    return out


def all_readouts(N: int, f: Sequence[int], S: Sequence[int], m: Optional[Metrics] = None) -> Dict:
    """Coefficient gcds, content gcd, and Euclid readout for one sketch."""
    cg = coef_gcds(S, N, m)
    ct = content_gcd(S, N)
    eu = euclid_readout(f, S, N, m)
    divs = sorted({d for _, d in cg} | ({ct} if 1 < ct < N else set())
                  | ({eu["factor"]} if eu.get("factor") else set()))
    return {"coef": cg, "content": ct, "euclid": eu, "divisors": divs}
