"""POSITIVE CONTROL: the small-characteristic systematic-equation relation family (BGJT/KW).

Setting (Kleinjung-Wesolowski 2019, Sect. 2-3.4; BGJT 2014):
  * base field k' = F_Q with Q = q^2, q = p^{e0};
  * target F_{Q^kk} = F_Q[X]/(I), I an irreducible factor of degree kk of h1 X^q - h0,
    deg h0, h1 <= 2, so that X^q == h0/h1 (mod I)   ("Frobenius is a small rational function");
  * V = span(X^{q+1}, X^q, X, 1).  For f = A X^{q+1} + B X^q + C X + D in V,
        h1 f == r_f := A X h0 + B h0 + C X h1 + D h1   (mod I),   deg r_f <= 3;
  * a RELATION is an f in P(V)(F_Q) that splits into linear factors over F_Q and whose
    r_f also splits.  The PGL_2(F_Q)-orbit of X^q - X (size q^3 + q) consists of split f.

Everything is exact.  `exhaustive_v_split` counts the split elements of P(V)(F_Q) exactly
(all ~Q^3 of them), using the fact that for f = t(x) + D the roots are the fibre t^{-1}(-D)
and that a repeated root can only sit at the unique zero of f' (= A x^q + C in char p).
"""
from __future__ import annotations

import random
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import numpy as np

from .ffield import GF, split_count_formula


def make_field(p: int, e0: int) -> Tuple[GF, int]:
    """F_{q^2} with q = p^e0. Returns (field, q)."""
    return GF(p, 2 * e0), p ** e0


def _mult_at(F: GF, f: List[int], x0: int) -> int:
    mu = 0
    while len(f) > 1 and F.peval(f, x0) == 0:
        f, rem = F.pdivmod(f, [int(F.neg[x0]), 1])
        assert not rem
        mu += 1
    return mu


def exhaustive_v_split(F: GF, q: int) -> Dict:
    """Exact enumeration of P(V)(F_Q). Returns counts and the list of split (A,B,C,D)."""
    Q = F.Q
    xs = F.elems
    xq = F.powvec(q)
    xq1 = F.mul[xq, xs]
    split: List[Tuple[int, int, int, int]] = []
    total = 0
    # case A = 1: f = x^{q+1} + B x^q + C x + D ; f' = x^q + C
    # x0 with x0^q = -C : Frobenius is a bijection, x0 = (-C)^{q^{?}}; find by table.
    inv_frob = np.zeros(Q, dtype=np.int64)
    inv_frob[xq] = xs  # xq is a permutation
    for B in range(Q):
        base = F.add[xq1, F.mul[B, xq]]
        for C in range(Q):
            t = F.add[base, F.mul[C, xs]]
            cnt = np.bincount(t, minlength=Q)
            x0 = int(inv_frob[int(F.neg[C])])
            D0 = int(F.neg[int(t[x0])])  # the only D for which x0 is a (repeated) root
            total += Q
            for val in np.nonzero(cnt >= 1)[0]:
                D = int(F.neg[int(val)])
                c = int(cnt[val])
                if D == D0:
                    f = [D, C] + [0] * (q - 2) + [B, 1]
                    mu = _mult_at(F, f, x0)
                    if c - 1 + mu == q + 1:
                        split.append((1, B, C, D))
                elif c == q + 1:
                    split.append((1, B, C, D))
    # case A = 0, B = 1: f = x^q + C x + D ; f' = C
    for C in range(Q):
        t = F.add[xq, F.mul[C, xs]]
        cnt = np.bincount(t, minlength=Q)
        total += Q
        for D in range(Q):
            if C == 0:
                split.append((0, 1, 0, D))  # (x + D^{1/q})^q
            elif cnt[int(F.neg[D])] == q:
                split.append((0, 1, C, D))
    # case A = B = 0, C = 1 (degree 1) and the constant form: always split
    for D in range(Q):
        split.append((0, 0, 1, D))
    split.append((0, 0, 0, 1))
    total += Q + 1
    assert total == (Q ** 4 - 1) // (Q - 1)
    return {"total": total, "split": split}


def orbit_size(q: int) -> int:
    """|PGL2(F_{q^2})| / |PGL2(F_q)| = q^3 + q."""
    return q ** 3 + q


def find_model(F: GF, q: int, seed: int, min_k: int = 3, tries: int = 4000):
    """Random h0, h1 of degree <= 2 with an irreducible factor I of h1 x^q - h0 of the largest
    degree kk >= min_k found (DDF piece of degree kk that is itself a single factor)."""
    rng = random.Random(f"systeq|bgjt-model|{F.Q}|{q}|{seed}")
    best = None
    for _ in range(tries):
        h0 = F.trim([rng.randrange(F.Q) for _ in range(3)])
        h1 = F.trim([rng.randrange(F.Q) for _ in range(2)] + [rng.randrange(1, F.Q)])
        if not h0:
            continue
        P = F.psub(F.pmul(h1, [0] * q + [1]), h0)
        dP = F.trim([F.m(i % F.p, c) for i, c in enumerate(P)][1:])
        if len(F.pgcd(P, dP)) > 1:
            continue
        for deg, g in F.irreducible_factors_ddf(P):
            if len(g) - 1 == deg and deg >= min_k and (best is None or deg > best[0]):
                best = (deg, h0, h1, g)
        if best is not None and best[0] >= q - 1:
            break
    if best is None:
        raise RuntimeError("no model found")
    kk, h0, h1, I = best
    return {"k": kk, "h0": h0, "h1": h1, "I": I}


def r_of(F: GF, f4, h0, h1):
    A, B, C, D = f4
    r = F.padd(F.pmul(F.trim([B, A]), h0), F.pmul(F.trim([D, C]), h1))
    return r


def f_poly(F: GF, q: int, f4):
    A, B, C, D = f4
    return F.trim([D, C] + [0] * (q - 2) + [B, A])


def relation_yield(F: GF, q: int, model: Dict, split_list) -> Dict:
    """Among split f in P(V), count those whose r_f is nonzero and splits.  Verifies
    h1 f == r_f (mod I) exactly for every relation."""
    h0, h1, I = model["h0"], model["h1"], model["I"]
    rel = 0
    checked = 0
    for f4 in split_list:
        r = r_of(F, f4, h0, h1)
        if not r:
            continue
        if F.splits(r):
            rel += 1
            f = f_poly(F, q, f4)
            lhs = F.pmod(F.psub(F.pmul(h1, f), r), I)
            assert not lhs, "congruence h1 f == r_f mod I failed"
            checked += 1
    return {"relations": rel, "verified_mod_I": checked}


def control_summary(p: int, e0: int, seed: int = 0) -> Dict:
    F, q = make_field(p, e0)
    Q = F.Q
    ex = exhaustive_v_split(F, q)
    S = len(ex["split"])
    model = find_model(F, q, seed)
    ry = relation_yield(F, q, model, ex["split"])
    sq, tq = split_count_formula(Q, q + 1)
    s3, t3 = split_count_formula(Q, 3)
    P_rand_q1 = Fraction(sq, tq)
    P_rand_3 = Fraction(s3, t3)
    frac_split = Fraction(S, ex["total"])
    Y = Fraction(ry["relations"], ex["total"])
    return {
        "p": p, "q": q, "Q": Q, "k": model["k"], "P_V": ex["total"], "split_in_V": S,
        "orbit_q3_plus_q": orbit_size(q),
        "frac_split_V": frac_split, "P_rand_split_deg_q+1": P_rand_q1,
        "G_V": frac_split / P_rand_q1,
        "relations": ry["relations"], "verified_mod_I": ry["verified_mod_I"],
        "yield_per_candidate": Y, "pred_size": P_rand_q1 * P_rand_3,
        "G_FF": Y / (P_rand_q1 * P_rand_3),
        "rel_over_factor_base": Fraction(ry["relations"], Q + 1),
    }
