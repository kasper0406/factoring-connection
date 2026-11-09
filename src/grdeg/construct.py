"""Construct frame-less type-D (subfield / Weil-restriction / CM) linked modules and measure box small-norm excess.

DIAGNOSTIC / harness code. Reuses src.rank3.forms and src.rankdegen.fh by import. No factoring; the geometry
(Delta_F, frames, box norm counts) is computed from the module alone and never needs p, q, or a root mod N.

A module is (f, P, g, Q): K = Q[x]/f, K' = Q[x]/g, generators u_j = P_j(alpha), u'_j = Q_j(beta). Linking Q_j = P_j
(mod N) is arranged by Q_j = P_j + N*(integer poly) so any common root m of f,g mod N links them; E1 needs only the
geometry. Norms W(a)=Res(f, sum a_j P_j), W'(a)=Res(g, sum a_j Q_j); F = W W'.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Sequence, Tuple

from src.rank3 import forms as RF
from src.rankdegen import fh as FH

Poly = List[int]


# --- number-field data -------------------------------------------------------

# K = Q(sqrt2, sqrt3), primitive element theta = sqrt2 + sqrt3, min poly x^4 - 10 x^2 + 1.
QUARTIC = [1, 0, -10, 0, 1]
# in terms of theta:  2*1 = 2 ; 2*sqrt2 = theta^3 - 9 theta ; 2*sqrt3 = 11 theta - theta^3 ; 2*sqrt6 = theta^2 - 5
TWO = [2]
TWO_S2 = [0, -9, 0, 1]      # theta^3 - 9 theta
TWO_S3 = [0, 11, 0, -1]     # 11 theta - theta^3
TWO_S6 = [-5, 0, 1]         # theta^2 - 5


def _shift(p: Poly, add: Poly, N: int) -> Poly:
    """return p + N*add (coefficientwise)."""
    n = max(len(p), len(add))
    out = [0] * n
    for k in range(len(p)):
        out[k] += p[k]
    for k in range(len(add)):
        out[k] += N * add[k]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def weil_restriction_q2q3(N: int) -> Dict[str, object]:
    """D1: rank-4 module over k=Q(sqrt2) inside K=K'=Q(sqrt2,sqrt3).
    u  = (2, 2 sqrt2, 2 sqrt3, 2 sqrt6)          [closed under * sqrt2: an O_{Q(sqrt2)}-module]
    u' = same, with the sqrt3 and sqrt6 generators N-shifted (type-D 'shifted generator', rank3/rankdegen).
    """
    f = QUARTIC
    P = [TWO, TWO_S2, TWO_S3, TWO_S6]
    g = QUARTIC
    # shift the z=sqrt3 generator by N (and sqrt6 = sqrt2*sqrt3 accordingly so the O_k-module structure is kept)
    Q = [TWO, TWO_S2, _shift(TWO_S3, TWO, N), _shift(TWO_S6, TWO_S2, N)]
    return {"name": "D1 Weil-restriction Q(sqrt2,sqrt3)/Q(sqrt2)", "f": f, "P": P, "g": g, "Q": Q, "r": 4}


def cyclotomic_z5(N: int) -> Dict[str, object]:
    """D2: rank-3 module on (1, zeta5, zeta5^2) over k=Q(zeta5), both sides equal up to an N-shift.
    K=K'=Q(zeta5), degree 4, min poly x^4+x^3+x^2+x+1. This is rank3's GP example; the N-shift of one
    generator is used to force a degeneracy. (Tests whether cyclotomic CM structure gives frame-less excess.)
    """
    f = [1, 1, 1, 1, 1]
    P = [[1], [0, 1], [0, 0, 1]]
    g = [1, 1, 1, 1, 1]
    Q = [[1], _shift([0, 1], [1], N), [0, 0, 1]]
    return {"name": "D2 cyclotomic Q(zeta5)", "f": f, "P": P, "g": g, "Q": Q, "r": 3}


def gp_control_cubics(N: int) -> Dict[str, object]:
    """General-position rank-3 control: two cubic power bases sharing a root mod N (rank3 E1 style)."""
    # two monic cubics; make them share root m mod N via g = f + (x - m) small; here just two fixed cubics.
    f = [-2, 0, 0, 1]   # x^3 - 2
    g = [-3, 1, 0, 1]   # x^3 + x - 3
    P = [[1], [0, 1], [0, 0, 1]]
    Q = [[1], [0, 1], [0, 0, 1]]
    return {"name": "GP control: cubics x^3-2, x^3+x-3 (power basis)", "f": f, "P": P, "g": g, "Q": Q, "r": 3}


# --- geometry ----------------------------------------------------------------

def geometry(mod: Dict[str, object]) -> Dict[str, object]:
    f, P, g, Q, r = mod["f"], mod["P"], mod["g"], mod["Q"], mod["r"]
    D = (len(f) - 1) + (len(g) - 1)
    dinfo = RF.delta_F(f, P, g, Q, r)
    delta_is_zero = (dinfo["Q2"] == 0)
    table = FH.minor_table(f, P, g, Q, r)
    fr = FH.frames(table)
    nmin = sum(1 for v in table["status"].values() if v == 0)
    ntot = len(table["status"])
    return {
        "D": D,
        "Q2_is_zero": delta_is_zero,
        "frame_free": (len(fr) == 0),
        "num_frames": len(fr),
        "num_zero_minors": nmin,
        "num_minors": ntot,
    }


def box_small_norm_counts(mod: Dict[str, object], E: int, thresholds: Sequence[int]) -> Dict[int, int]:
    """#{a in box, 0<|F|<=X} for each X in thresholds."""
    vals = RF.box_abs_values(mod["f"], mod["P"], mod["g"], mod["Q"], mod["r"], E)
    vals = [v for v in vals if v > 0]
    out = {}
    for X in thresholds:
        out[X] = sum(1 for v in vals if v <= X)
    return out


if __name__ == "__main__":
    N = 1001  # 7*11*13; geometry is N-independent except via the shift magnitude
    for ctor in (weil_restriction_q2q3, cyclotomic_z5, gp_control_cubics):
        mod = ctor(N)
        geo = geometry(mod)
        print(f"{mod['name']}: D={geo['D']} r={mod['r']} "
              f"Delta_F==0? {geo['Q2_is_zero']}  frame_free? {geo['frame_free']} "
              f"(frames={geo['num_frames']}, zero_minors={geo['num_zero_minors']})")
