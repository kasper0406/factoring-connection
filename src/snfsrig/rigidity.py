"""Diagnostic checkers for the rigidity statements of notes/snfsrig_results.md section 2.

  Lemma RES  (resultant floor):  f, g in Z[x] with a common root m mod N and Res(f, g) != 0  =>  N | Res(f, g),
                                 hence |Res| >= N and (Hadamard) ||f||_2^{deg g} ||g||_2^{deg f} >= N.
  Lemma HR   (height rigidity):  f0 irreducible, f0(m) = 0 mod N, f(m) = 0 mod N, f0 does not divide f
                                 =>  ||f||_2 >= (N / ||f0||_2^{deg f})^{1/deg f0}.
  Lemma RS   (rational-side rigidity): (g0, g1) with g0 = m g1 (mod N), max(|g0|, |g1|) < N/(2m + 2)
                                 => (g0, g1) = g1 (m, 1).
  Prop DI    (device invariance): the devices N*w, kN, root twists, f(x + t), u^d f(x/u), Moebius transforms map the
                                 value pair (F, G)(P) to (c1 F0(T P), c2 G0(T P)) with T in GL_2(Q).

These are checks of proved lemmas on small instances (exact integer arithmetic; sympy only for polynomial algebra
and resultants of explicit small polynomials -- never to factor N).  Not operational code.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction as Fr
from typing import Iterable, List, Sequence, Tuple

import sympy as sp

X, A, Bv = sp.symbols("x a b")


def poly_eval_mod(coeffs_low_to_high: Sequence[int], m: int, N: int) -> int:
    v = 0
    for c in reversed(coeffs_low_to_high):
        v = (v * m + c) % N
    return v


def l2(coeffs: Sequence[int]) -> float:  # only for reporting, never for decisions
    return math.sqrt(sum(c * c for c in coeffs))


def l2sq(coeffs: Sequence[int]) -> int:
    return sum(c * c for c in coeffs)


def resultant(f: Sequence[int], g: Sequence[int]) -> int:
    F = sp.Poly(list(reversed(list(f))), X)
    G = sp.Poly(list(reversed(list(g))), X)
    return int(sp.resultant(F, G))


def check_res_floor(f: Sequence[int], g: Sequence[int], m: int, N: int) -> dict:
    """Lemma RES on one instance: f(m) = g(m) = 0 mod N.  Returns the facts checked exactly."""
    assert poly_eval_mod(f, m, N) == 0 and poly_eval_mod(g, m, N) == 0
    R = resultant(f, g)
    df, dg = len(f) - 1, len(g) - 1
    had_sq = l2sq(f) ** dg * l2sq(g) ** df          # Hadamard: R^2 <= ||f||^{2 dg} ||g||^{2 df}
    return {"res": R, "divisible": R % N == 0, "nonzero": R != 0,
            "hadamard_ok": R * R <= had_sq, "abs_ge_N": (R == 0) or abs(R) >= N}


def snfs_pair(r: int, e: int, s: int, d: int) -> Tuple[List[int], int, int, int]:
    """N = r^e - s.  With q = ceil(e/d), m = r^q and j = d q - e:  r^j N = m^d - s r^j, so
    f0(x) = x^d - s r^j has f0(m) = r^j N = 0 (mod N) (when gcd(r, N) = 1 this is the root used).
    Returns (f0 low->high, m, N, j)."""
    N = r ** e - s
    q = -(-e // d)
    j = d * q - e
    m = r ** q
    f0 = [-s * r ** j] + [0] * (d - 1) + [1]
    assert poly_eval_mod(f0, m, N) == 0
    return f0, m, N, j


def is_irreducible(f: Sequence[int]) -> bool:
    return sp.Poly(list(reversed(list(f))), X).is_irreducible


def divides_over_Q(f0: Sequence[int], f: Sequence[int]) -> bool:
    F0 = sp.Poly(list(reversed(list(f0))), X, domain="QQ")
    F = sp.Poly(list(reversed(list(f))), X, domain="QQ")
    if F.is_zero:
        return True
    return F.rem(F0).is_zero


def hr_bound_sq_pow(f0: Sequence[int], deg_f: int, N: int) -> Tuple[int, int]:
    """Lemma HR in squared integer form: ||f||_2^{2 d0} >= N^2 / ||f0||_2^{2 deg_f}.  Returns (lhs_needed_num, den)
    so that the claim reads  l2sq(f)^{d0} * den >= num  with num = N^2, den = l2sq(f0)^{deg_f}."""
    return N * N, l2sq(f0) ** deg_f


def exhaustive_small_height(f0: Sequence[int], m: int, N: int, deg: int, H: int) -> dict:
    """Enumerate every f of degree <= deg with coefficients in [-H, H] and f(m) = 0 mod N; report those not
    divisible by f0 (over Q).  Lemma HR predicts none when l2sq(f)^{d0} * l2sq(f0)^{deg} < N^2 for all such f,
    i.e. in particular when ((deg+1) H^2)^{d0} * l2sq(f0)^{deg} < N^2."""
    d0 = len(f0) - 1
    found, nondiv = 0, []
    for cs in itertools.product(range(-H, H + 1), repeat=deg + 1):
        if all(c == 0 for c in cs):
            continue
        if poly_eval_mod(cs, m, N) != 0:
            continue
        found += 1
        if not divides_over_Q(f0, cs):
            nondiv.append(cs)
    guaranteed = ((deg + 1) * H * H) ** d0 * l2sq(f0) ** deg < N * N
    # every non-multiple must violate the guarantee region (consistency with Lemma HR)
    consistent = all(l2sq(c) ** d0 * l2sq(f0) ** deg >= N * N for c in nondiv)
    return {"found": found, "nondivisible": nondiv, "hr_guarantees_none": guaranteed, "consistent_with_HR": consistent}


def rational_side_vectors(m: int, N: int, Hg: int) -> List[Tuple[int, int]]:
    """All (g0, g1) with |g0|, |g1| <= Hg, (g0, g1) != 0, g0 = m g1 (mod N).  Enumerates g1 and the few g0."""
    out = []
    for g1 in range(-Hg, Hg + 1):
        base = (m * g1) % N
        for g0 in (base - N * t for t in range(-2, 3)):
            if abs(g0) <= Hg and (g0, g1) != (0, 0):
                out.append((g0, g1))
    return sorted(set(out))


def check_rs(m: int, N: int, Hg: int) -> dict:
    """Lemma RS: if Hg < N/(2m + 2) (and 0 <= m < N) every vector is g1 (m, 1)."""
    vecs = rational_side_vectors(m, N, Hg)
    bad = [(g0, g1) for (g0, g1) in vecs if g0 != m * g1]
    return {"n": len(vecs), "non_multiples": bad, "hypothesis": Hg * (2 * m + 2) < N}


# ---------------------------------------------------------------------------------------------------------
# Prop DI: device invariance (symbolic identities)
# ---------------------------------------------------------------------------------------------------------


def hom(f_low_high: Sequence[int], a, b):
    d = len(f_low_high) - 1
    return sum(sp.Integer(c) * a ** i * b ** (d - i) for i, c in enumerate(f_low_high))


def device_identities(f0: Sequence[int], m: int) -> dict:
    """For the devices of Prop DI, return sympy expressions that must simplify to 0.
    F0(a, b) = b^d f0(a/b), G0(a, b) = a - m b.  T acts on column vectors P = (a, b)."""
    d = len(f0) - 1
    F0 = hom(f0, A, Bv)
    G0 = A - m * Bv
    u, t, w, lam = sp.Integer(3), sp.Integer(5), sp.Integer(7), sp.Integer(2)
    al, be, ga, de = sp.Integer(2), sp.Integer(1), sp.Integer(1), sp.Integer(1)   # Moebius x -> (al x + be)/(ga x + de)
    fx = sp.Poly(list(reversed(list(f0))), X).as_expr()
    out = {}
    # (b) multiplier w: f = w f0, g = g0
    out["w_multiplier"] = sp.expand(hom([w * c for c in f0], A, Bv) - w * F0)
    # (a') root twist by lambda: root lambda m, f_l(x) = lambda^d f0(x/lambda), g = x - lambda m
    f_l = sp.Poly(sp.expand(lam ** d * fx.subs(X, X / lam)), X)
    F_l = hom(list(reversed(f_l.all_coeffs())), A, Bv)
    out["root_twist_F"] = sp.expand(F_l - F0.subs(Bv, lam * Bv))
    out["root_twist_G"] = sp.expand((A - lam * m * Bv) - G0.subs(Bv, lam * Bv))
    # (c) shift f(x + t): root m - t, F_t(a, b) = F0(a + t b, b)
    f_t = sp.Poly(sp.expand(fx.subs(X, X + t)), X)
    F_t = hom(list(reversed(f_t.all_coeffs())), A, Bv)
    out["shift_F"] = sp.expand(F_t - F0.subs(A, A + t * Bv))
    out["shift_G"] = sp.expand((A - (m - t) * Bv) - G0.subs(A, A + t * Bv))
    # (c') scaling u^d f(x/u): root u m, F(a, b) = F0(a, u b)
    f_u = sp.Poly(sp.expand(u ** d * fx.subs(X, X / u)), X)
    F_u = hom(list(reversed(f_u.all_coeffs())), A, Bv)
    out["scale_F"] = sp.expand(F_u - F0.subs(Bv, u * Bv))
    # (c'') Moebius: f_M(x) = (ga x + de)^d f0((al x + be)/(ga x + de)); F_M(a, b) = F0(al a + be b, ga a + de b)
    f_M = sp.Poly(sp.expand(sp.cancel((ga * X + de) ** d * fx.subs(X, (al * X + be) / (ga * X + de)))), X)
    F_M = hom(list(reversed(f_M.all_coeffs())), A, Bv)
    out["moebius_F"] = sp.expand(F_M - F0.subs({A: al * A + be * Bv, Bv: ga * A + de * Bv}, simultaneous=True))
    return out


# ---------------------------------------------------------------------------------------------------------
# Prop RB: size of the total value set of a fixed-family one-level design
# ---------------------------------------------------------------------------------------------------------


def value_set_bound(n_devices: int, n_targets: int, box_points: int) -> int:
    """|V| <= 2 * (#fields/devices) * (#special-q targets) * (#box points): every query is G(P)/l or F(P) for one
    of the listed (device, target, point) triples."""
    return 2 * n_devices * n_targets * box_points
