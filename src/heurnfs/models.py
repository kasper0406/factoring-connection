"""Heuristic L[1/3] cost programs for NFS variants (agent heurnfs; notes/heurnfs_results.md).

Units. x = log N. Every cost is an L_N[1/3, .] exponent. A norm of size exp(s x^{2/3} (log x)^{1/3}) has
"size coefficient" s. Degree d = delta (x / log x)^{1/3}. Sieve box E = L[1/3, e] per coordinate (search
entropy 2e), factor base B = L[1/3, beta].

HEURISTIC H (used by every program): a norm with size coefficient s is B-smooth with probability
L[1/3, -s/(3 beta)] (Canfield-Erdos-Pomerance applied to polynomial values), joint events multiply up to
L^{o(1)}, smoothness tests cost L^{o(1)} (ECM or batch), and sparse linear algebra on C columns costs
C^omega with omega = 2 (block Wiedemann).

Ledger (base-m, sieving dimension t = 2): rational size 1/delta, algebraic size delta*e + 1/delta.
A polynomial-height factor kappa multiplies the 1/delta terms: kappa = 2 is general N, kappa = 1 is SNFS.

All numeric optima use src/heurnfs/optim.py (floating point = analysis only). Closed forms are exact sympy.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import sqrt
from typing import Dict

import numpy as np
import sympy as sp

from src.heurnfs.optim import max_violation, minimise

# ----------------------------------------------------------------------------- closed forms (KNOWN)


def gnfs_closed(kappa=2) -> sp.Expr:
    """GNFS-type constant (32 kappa / 9)^{1/3}: kappa=2 -> (64/9)^{1/3}, kappa=1 -> SNFS (32/9)^{1/3}.

    Derivation: e = beta at balance, relation constraint 3 beta^2 >= delta beta + kappa/delta, minimised at
    delta = sqrt(kappa/beta): beta^3 = 4 kappa / 9, c = 2 beta.
    """
    k = sp.nsimplify(kappa)
    return sp.cbrt(sp.Rational(32, 9) * k)


def mnfs_closed() -> sp.Expr:
    """Coppersmith 1993: ((92 + 26 sqrt 13)/27)^{1/3}."""
    return sp.cbrt((92 + 26 * sp.sqrt(13)) / sp.Integer(27))


def tdim_closed(t: int) -> sp.Expr:
    """GNFS with sieving polynomials of degree t-1 (t coefficients, t = 2 is classical).

    Ledger: norm product size 2(t-1)/delta + 2 delta e / t (heights^{t-1}, A^{d+1}, A = E^{2/t}).
    Substituting delta = t delta'/2 gives kappa_t = 4(t-1)/t, so c(t) = (128 (t-1) / (9 t))^{1/3}.
    """
    return sp.cbrt(sp.Rational(128, 9) * sp.Rational(t - 1, t))


def init_split_closed() -> sp.Expr:
    """Plain initial split of a random target: R = u/v with u, v ~ N^{1/2}, both L[2/3, y]-smooth (ECM).

    Cost exponent 1/(3y) + sqrt(4y/3), minimised at y = 3^{-1/3}: value 3^{1/3} (Commeine-Semaev).
    """
    return sp.cbrt(3)


# ----------------------------------------------------------------------------- numeric programs


def _pack(x, names, val, cons) -> Dict:
    d = {n: float(v) for n, v in zip(names, x)}
    d["c"] = float(val)
    d["viol"] = max_violation(cons, x)
    return d


def gnfs(kappa: float = 2.0, omega: float = 2.0, seed: int = 0) -> Dict:
    """min max(2e, omega*beta) s.t. 2e - (delta e + kappa/delta)/(3 beta) >= beta."""
    def obj(x):
        e, b, d = x
        return max(2 * e, omega * b)

    def cons(x):
        e, b, d = x
        return [2 * e - (d * e + kappa / d) / (3 * b) - b]
    x, v = minimise(obj, cons, [.3, .3, .3], [1.5, 1.5, 3.], seed=seed)
    return _pack(x, ("e", "beta", "delta"), v, cons)


def tdim(t: int, seed: int = 0) -> Dict:
    """GNFS with t-dimensional sieving; search entropy 2e, norm size 2(t-1)/delta + 2 delta e / t."""
    def obj(x):
        e, b, d = x
        return max(2 * e, 2 * b)

    def cons(x):
        e, b, d = x
        return [2 * e - (2 * (t - 1) / d + 2 * d * e / t) / (3 * b) - b]
    x, v = minimise(obj, cons, [.3, .3, .3], [1.5, 1.5, 4.], seed=seed)
    return _pack(x, ("e", "beta", "delta"), v, cons)


def mnfs(omega: float = 2.0, enum: str = "box", la: str = "flat", kr: float = 1.0, ka: float = 1.0,
         omega_red: float = None, tests: bool = True, seed: int = 0) -> Dict:
    """Coppersmith MNFS: one rational side, V = L[v] algebraic sides f + (c1 x + c0) g.

    Variables (e, beta, beta_a, v, delta). r = kr/(3 beta delta) (rational smoothness penalty),
    a = (delta e + ka/delta)/(3 beta_a) (one algebraic side), s = max(beta, beta_a + v) (all columns).
    Relations: 2e + v - r - a >= s.
    Costs: sieve 2e ('box') or 2e - r ('output': HYPOTHETICAL output-sensitive enumeration of the
    rational-smooth pairs); algebraic tests 2e - r + v; LA omega*s ('flat') or, HYPOTHETICALLY,
    max(omega_red*beta, v + 2 beta_a) ('block': per-field elimination then a reduced rational system).
    kr, ka: height multipliers of the 1/delta terms (ka = 2 models a w ~ N^{1/d} multiplier in a_0).
    """
    om_r = omega if omega_red is None else omega_red

    def parts(x):
        e, b, ba, v, d = x
        r = kr / (3 * b * d)
        a = (d * e + ka / d) / (3 * ba)
        s = max(b, ba + v)
        return e, b, ba, v, d, r, a, s

    def obj(x):
        e, b, ba, v, d, r, a, s = parts(x)
        sieve = 2 * e if enum == "box" else 2 * e - r
        lac = omega * s if la == "flat" else max(om_r * b, v + 2 * ba)
        return max(sieve, (2 * e - r + v) if tests else -1e9, lac)

    def cons(x):
        e, b, ba, v, d, r, a, s = parts(x)
        return [2 * e + v - r - a - s]
    x, val = minimise(obj, cons, [.3, .3, .3, .01, .5], [1.3, 1.3, 1.3, .6, 2.5], restarts=60, seed=seed)
    out = _pack(x, ("e", "beta", "beta_a", "v", "delta"), val, cons)
    e, b, ba, v, d, r, a, s = parts(x)
    out.update(r=float(r), a=float(a), test_cost=float(2 * e - r + v))
    return out


def factory(seed: int = 0) -> Dict:
    """Coppersmith factorization factory: per-N cost max(2e - r, 2 max(beta, beta_a)) after a shared
    rational table (precomputation 2e). Second pass: least precomputation at the per-N optimum."""
    def parts(x):
        e, b, ba, d = x
        return e, b, ba, d, 1 / (3 * b * d), (d * e + 1 / d) / (3 * ba)

    def per_n(x):
        e, b, ba, d, r, a = parts(x)
        return max(2 * e - r, 2 * max(b, ba))

    def cons(x):
        e, b, ba, d, r, a = parts(x)
        return [2 * e - r - a - max(b, ba)]
    x, c = minimise(per_n, cons, [.5, .5, .5, .5], [1.5, 1.5, 1.5, 2.], seed=seed)

    def cons2(y):
        return cons(y) + [c + 1e-7 - per_n(y)]
    y, pre = minimise(lambda y: 2 * y[0], cons2, list(np.maximum(x - .05, .01)), list(x + .05), seed=seed)
    out = _pack(x, ("e", "beta", "beta_a", "delta"), c, cons)
    out["precomp"] = float(pre)
    out["precomp_viol"] = max_violation(cons2, y)
    return out


def mnfs_symmetric(seed: int = 0) -> Dict:
    """'Two-of-V' MNFS: V = L[v] polynomials of degree d sharing a root, heights N^{1/(d+1)} (lattice), a
    relation needs any two sides smooth. Each side size delta e + 1/delta; P1 = that/(3 beta).
    Relations 2e + 2v - 2 P1 >= v + beta (union approx valid while v <= P1); cost max(2e + v, 2(v + beta))."""
    def parts(x):
        e, b, v, d = x
        return e, b, v, d, (d * e + 1 / d) / (3 * b)

    def obj(x):
        e, b, v, d, p1 = parts(x)
        return max(2 * e + v, 2 * (v + b))

    def cons(x):
        e, b, v, d, p1 = parts(x)
        return [2 * e + 2 * v - 2 * p1 - v - b, p1 - v]
    x, val = minimise(obj, cons, [.3, .3, .01, .3], [1.5, 1.5, .6, 3.], restarts=60, seed=seed)
    return _pack(x, ("e", "beta", "v", "delta"), val, cons)


def descent_rows(c_init: float = None, seed: int = 0) -> Dict:
    """Relations produced as descent trees (rows = descended targets; Dixon descent).

    Bottom bound (model-PROVED): every row contains >= 1 bottom special-q relation, q ~ B; its q-lattice
    region must hold L[P] points, so 2e >= beta + P with P = (delta e + 2/delta)/(3 beta); B rows cost
    >= L[beta + P]. With c_init (cost exponent of the initial split of each target) the row cost is also
    >= beta + c_init. Objective max(2e, 2 beta, beta + c_init)."""
    ci = -1e9 if c_init is None else c_init

    def obj(x):
        e, b, d = x
        return max(2 * e, 2 * b, b + ci)

    def cons(x):
        e, b, d = x
        P = (d * e + 2 / d) / (3 * b)
        return [2 * e - b - P]
    x, v = minimise(obj, cons, [.3, .3, .3], [1.5, 1.5, 3.], seed=seed)
    return _pack(x, ("e", "beta", "delta"), v, cons)


def init_split_numeric() -> float:
    ys = np.linspace(0.05, 3, 200001)
    return float(np.min(1 / (3 * ys) + np.sqrt(4 * ys / 3)))


# ----------------------------------------------------------------------------- exact height LP


def height_lp(d1: int, d2: int) -> Fr:
    """min h1 + h2 (in units of x = log N) over polynomial pairs of degrees d1 <= d2 sharing a root mod N:
    counting (HEURISTIC G, Gaussian/lattice count):  (d1+1) h1 + (d2+1) h2 >= 1
    resultant (PROVED: N | Res(f, g) != 0, Res ~ H1^{d2} H2^{d1}):  d2 h1 + d1 h2 >= 1.
    Two-variable LP: the optimum is at a vertex; exact Fractions."""
    A = [(Fr(d1 + 1), Fr(d2 + 1)), (Fr(d2), Fr(d1))]
    cand = [(Fr(0), max(1 / A[0][1], 1 / A[1][1])), (max(1 / A[0][0], 1 / A[1][0]), Fr(0))]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    if det != 0:
        h1 = (A[1][1] - A[0][1]) / det
        h2 = (A[0][0] - A[1][0]) / det
        if h1 >= 0 and h2 >= 0:
            cand.append((h1, h2))
    feas = [c for c in cand if all(a * c[0] + b * c[1] >= 1 for a, b in A)]
    return min(c[0] + c[1] for c in feas)


def tower_ledger(eta: int, d: int) -> Dict[str, Fr]:
    """Tower (TNFS-style) variant for factoring: base ring Z[iota], deg h = eta, ideal n = (N, iota - r)
    of norm N; f over Z[iota] of degree d, g = x - m(iota). In integer-height units H' = eta * h:
    counting  (d+1) Hf' + Hm' >= 1,  resultant  Hf' + d Hm' >= 1  (norm of f(m) in n is >= N),
    norm product  Hf' + Hm' + (d+1) C / 2  for search entropy C (= 2 eta log A).
    eta cancels: identical to the 1-D ledger with a monic linear g. Returns the LP value (units of x).
    OPTIMISTIC: an ideal (N, iota - r) needs h with a root mod N of small height, which violates G
    (height >= N^{1/(eta+1)}); granting it free only helps the tower. eta does not enter this LP by
    construction. The genuine TNFS analogue is tnfs_genuine_coeff (strictly worse)."""
    A = [(Fr(d + 1), Fr(1)), (Fr(1), Fr(d))]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    h1 = (A[1][1] - A[0][1]) / det
    h2 = (A[0][0] - A[1][0]) / det
    cands = [(h1, h2), (Fr(0), Fr(1)), (Fr(1), Fr(0))]
    feas = [c for c in cands if c[0] >= 0 and c[1] >= 0 and all(a * c[0] + b * c[1] >= 1 for a, b in A)]
    return {"eta": Fr(eta), "heights": min(c[0] + c[1] for c in feas), "search_coeff": Fr(d + 1, 2)}


def tnfs_genuine_coeff(eta: int) -> Fr:
    """Genuine TNFS analogue for factoring: h of degree eta irreducible modulo the primes of N, f, g in Z[x]
    with a common root mod N, elements a(iota) - b(iota) x, search entropy C = 2 eta log A. Norm over Q is
    (||f|| A^d)^eta (||g|| A)^eta, so the ledger is eta * 2x/D + D C / 2: the coefficient of x/D is 2 eta
    (vs 2 for GNFS). Strictly worse for eta >= 2."""
    return Fr(2 * eta)


# ----------------------------------------------------------------------------- exponent barriers


def alpha_floor_generic() -> Fr:
    """Size-model barrier: log X >= 2 sqrt(x C) with C ~ x^alpha gives log X ~ x^{(1+alpha)/2};
    smoothing at B = L[alpha] costs exponent (1+alpha)/2 - alpha = (1-alpha)/2 <= alpha  =>  alpha >= 1/3."""
    # (1 - a)/2 <= a  <=>  a >= 1/3
    return Fr(1, 3)


def alpha_floor_one_level_descent() -> Fr:
    """One descent level with initial bound L[eta]: 1 - eta <= alpha and (1+eta)/2 <= 2 alpha => 2/5."""
    # eta >= 1 - alpha; (2 - alpha)/2 <= 2 alpha  <=>  alpha >= 2/5
    return Fr(2, 5)


def alpha_floor_aggregate(r_exp: Fr = Fr(1, 4), d_exp: Fr = Fr(1, 2)) -> Fr:
    """HYPOTHETICAL r20 #5 ledger log X << x/d + (d/r) C with d ~ x^{d_exp}, r ~ x^{r_exp}, C ~ x^alpha.
    log X exponent = max(1 - d_exp, d_exp - r_exp + alpha); smoothing at L[alpha] needs it <= 2 alpha.
    Returns the least alpha allowed for the given (d_exp, r_exp)."""
    # need 1 - d_exp <= 2 a and d_exp - r_exp + a <= 2 a
    return max((1 - d_exp) / 2, d_exp - r_exp)


def la_omega_curve(omegas=(1.0, 1.5, 2.0, 2.373)) -> Dict[float, Dict[str, float]]:
    """HYPOTHETICAL: constants if sparse linear algebra cost C^omega."""
    out = {}
    for w in omegas:
        g = gnfs(omega=w)["c"]
        # MNFS contains GNFS (v = 0, beta_a = beta), so MNFS <= GNFS; min() removes Nelder-Mead artifacts.
        out[w] = {"gnfs": g, "mnfs": min(mnfs(omega=w)["c"], g),
                  "mnfs_no_test_term": min(mnfs(omega=w, tests=False)["c"], g)}
    return out


def table() -> Dict[str, Dict]:
    """Every row of notes/heurnfs_results.md Table 1 (numeric)."""
    out = {
        "gnfs": gnfs(),
        "snfs_kappa1": gnfs(kappa=1.0),
        "gnfs_w_in_a0": gnfs(kappa=3.0),
        "mnfs": mnfs(),
        "mnfs_w_in_a0": mnfs(ka=2.0),
        "factory": factory(),
        "tdim3": tdim(3),
        "tdim4": tdim(4),
        "mnfs_symmetric": mnfs_symmetric(),
        "descent_bottom": descent_rows(),
        "descent_rows_plain_init": descent_rows(c_init=float(init_split_closed())),
        "mnfs_block_la": mnfs(la="block"),
        "mnfs_output_enum": mnfs(enum="output"),
    }
    return out


if __name__ == "__main__":  # pragma: no cover
    import json
    t = table()
    for k, v in t.items():
        print(f"{k:26s} c = {v['c']:.4f}   " + json.dumps({a: round(b, 4) for a, b in v.items() if a != 'c'}))
    print("closed: gnfs", float(gnfs_closed()), "mnfs", float(mnfs_closed()), "snfs", float(gnfs_closed(1)),
          "t3", float(tdim_closed(3)), "init", float(init_split_closed()))
    print("omega curve", la_omega_curve())
    print("height LP (1,d):", [str(height_lp(1, d)) for d in range(2, 8)])
    print("height LP (d,d):", [str(height_lp(d, d)) for d in range(2, 8)])
    print("height LP (2,10):", height_lp(2, 10), "tower:", tower_ledger(3, 5))
