"""Exact asymptotic bookkeeping for Hyp U' at the two-level 3/8 parameters (notes/uprime_results.md).

Every log-size is an `AExp`: a finite sum  sum_j c_j (log N)^{p_j} (log log N)^{q_j}  with Fraction coefficients
and exponents (class copied from src/lv5unif/params.py, AGENT_RULES rule 7).  X "beats" Y iff the leading term of
log X - log Y is positive and of positive order.  Leading-order comparisons only: this checks the notes' exponent
arithmetic; it is not an independent proof.

Design 2L (multidescent §3.2), D := log N/loglog N, lam := loglog N, Lambda := D^{3/8} lam:
  d ~ d0 D^{1/4} (odd), k ~ kappa D^{3/8}, nu = 1/(4k), mu = (1-nu)/d, h = 2 nu + eta0, eta0 = D^{-1/2};
  B  = exp(a Lambda)            (factor base, modulus smoothness y1)
  B2 = exp(b2 D^{1/2} lam)      (child bound, value smoothness y2 of the level-1 algebraic value)
  log Y = log N/(2k) = D^{5/8} lam/(2 kappa)      (special-q primes l = L[5/8])
  log X_r ~ mu log N ~ D^{3/4} lam/d0            (level-1 rational cofactor)
  log X_F ~ (mu + h + d nu) log N ~ (d0/(4 kappa)) D^{7/8} lam     (level-1 algebraic value)
  u_r  = log X_r/log B  ~ D^{3/8}/(a d0);   u_F2 = log X_F/log B2 ~ (d0/(4 kappa b2)) D^{3/8};  log u ~ (3/8) lam
  sigma1 = lim (u_r log u_r + u_F2 log u_F2)/Lambda = 3/(8 a d0) + 3 d0/(32 kappa b2)        (mixed scan, S1^2)
  level 2 (child q <= B2, V2 = sqrt(q) S2):  log X_r2 ~ D^{3/4} lam/d0,
    log(X_F2/q) <= (1/d0 + d0 b2/2) D^{3/4} lam,  sigma2 = (3/8)(2/d0 + d0 b2/2)/a          (child box W2 = S2^2)
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple

Key = Tuple[Fr, Fr]


class AExp:
    """sum c (log N)^p (loglog N)^q, exact (copied from src/lv5unif/params.py)."""

    def __init__(self, terms: Dict[Key, Fr] | None = None):
        self.t: Dict[Key, Fr] = {}
        for k, v in (terms or {}).items():
            kk = (Fr(k[0]), Fr(k[1]))
            self.t[kk] = self.t.get(kk, Fr(0)) + Fr(v)
        self.t = {k: v for k, v in self.t.items() if v != 0}

    @staticmethod
    def mono(c, p, q) -> "AExp":
        return AExp({(Fr(p), Fr(q)): Fr(c)})

    def __add__(self, o: "AExp") -> "AExp":
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = r.get(k, Fr(0)) + v
        return AExp(r)

    def __neg__(self) -> "AExp":
        return AExp({k: -v for k, v in self.t.items()})

    def __sub__(self, o: "AExp") -> "AExp":
        return self + (-o)

    def scale(self, c) -> "AExp":
        return AExp({k: Fr(c) * v for k, v in self.t.items()})

    def leading(self):
        if not self.t:
            return None
        k = max(self.t)
        return k, self.t[k]

    def tends_to_plus_infinity(self) -> bool:
        ld = self.leading()
        return ld is not None and ld[1] > 0 and ld[0] > (Fr(0), Fr(0))

    def tends_to_zero(self) -> bool:
        """For a *size* (not a log): every term has negative order (leading order < (0, 0))."""
        ld = self.leading()
        return ld is None or ld[0] < (Fr(0), Fr(0))

    def __repr__(self) -> str:
        return " + ".join(f"{v}*(logN)^{k[0]}*(loglogN)^{k[1]}" for k, v in sorted(self.t.items(), reverse=True))


def Dl(c, s, extra_lam=1) -> AExp:
    """c * D^s * lam^{extra_lam} = c (log N)^s (loglog N)^{extra_lam - s}."""
    return AExp.mono(c, Fr(s), Fr(extra_lam) - Fr(s))


LAM = AExp.mono(1, 0, 1)                      # loglog N


def sizes(a, d0, kappa, b2) -> Dict[str, AExp]:
    a, d0, kappa, b2 = map(Fr, (a, d0, kappa, b2))
    Lam = Dl(1, Fr(3, 8))
    return {
        "Lambda": Lam,
        "logB": Lam.scale(a),
        "logB2": Dl(b2, Fr(1, 2)),
        "logY": Dl(1 / (2 * kappa), Fr(5, 8)),
        "logH": Dl(1 / (2 * kappa), Fr(5, 8)),                  # (2 nu + eta0) log N, leading term
        "eta0_logN": Dl(1, Fr(1, 2)),
        "logXr": Dl(1 / d0, Fr(3, 4)),
        "logXF": Dl(d0 / (4 * kappa), Fr(7, 8)),
        "d": Dl(d0, Fr(1, 4), 0),                               # d itself (a size)
        "u_r": Dl(1 / (a * d0), Fr(3, 8), 0),
        "u_F2": Dl(d0 / (4 * kappa * b2), Fr(3, 8), 0),
        "logS1": Lam.scale(sigma1(a, d0, kappa, b2) / 2),
        "logXr2": Dl(1 / d0, Fr(3, 4)),
        "logXF2_over_q_max": Dl(1 / d0 + d0 * b2 / 2, Fr(3, 4)),
        "logS2": Lam.scale(sigma2(a, d0, kappa, b2) / 2),
        "omega_upper": Dl(1 / a, Fr(5, 8), 0),                  # omega(N) <= log N/log B = D^{5/8}/a
        "mbar": Dl(d0 / (4 * kappa * a), Fr(1, 2), 0),          # max #children = log X_F/log B
    }


def sigma1(a, d0, kappa, b2) -> Fr:
    """log(1/(rho_r(B) rho_F(B2))) = sigma1 Lambda (1+o(1)): u log u with log u ~ (3/8) lam on both sides."""
    a, d0, kappa, b2 = map(Fr, (a, d0, kappa, b2))
    return Fr(3, 8) * (1 / (a * d0) + d0 / (4 * kappa * b2))


def sigma1_single_y(a, d0, kappa) -> Fr:
    """What Lemma CS (both Rankin factors at B) would give: the algebraic u at B is ~ D^{1/2}, log u ~ (1/2) lam,
    so u_F log u_F ~ (d0/(8 kappa a)) D^{1/2} lam  =  L[1/2]-scale.  Returned: the coefficient of D^{1/2} lam."""
    a, d0, kappa = map(Fr, (a, d0, kappa))
    return Fr(1, 2) * d0 / (4 * kappa * a)


def sigma2(a, d0, kappa, b2) -> Fr:
    """Child box (worst window q ~ B2): u_r2 log u_r2 + u_F2' log u_F2' with B-smoothness on both sides."""
    a, d0, kappa, b2 = map(Fr, (a, d0, kappa, b2))
    return Fr(3, 8) * ((1 / d0) + (1 / d0 + d0 * b2 / 2)) / a


def sq_constraint(a, d0, kappa, b2) -> bool:
    """Item 3 (union bound): sum_{p > B} 1/p^2 ~ 1/(B log B) must be << rho_F(B2) (log N)^{-C}:
    a Lambda > u_F2 log u_F2 (1+o(1)) = (3 d0/(32 kappa b2)) Lambda."""
    a, d0, kappa, b2 = map(Fr, (a, d0, kappa, b2))
    return a > Fr(3) * d0 / (32 * kappa * b2)


def cm_level_constraint(d0, b2) -> bool:
    """Hyp CM proof plan, algebraic pair: modulus |G(P')| * Delta/q ~ M V1 S2^2 must sit below level 1/2 of the
    value F(P')/q (leading order D^{3/4} lam):  1/d0 < (1/2)(1/d0 + d0 b2/2)  <=>  d0^2 b2 > 2."""
    d0, b2 = Fr(d0), Fr(b2)
    return d0 * d0 * b2 > 2


def margins(a, d0, kappa, b2, C=Fr(100)) -> List[Tuple[str, bool, str]]:
    """Leading-order growth comparisons used in notes/uprime_results.md.  (name, holds, where used)."""
    s = sizes(a, d0, kappa, b2)
    C = Fr(C)
    out: List[Tuple[str, bool, str]] = []

    def chk(name, expr: AExp, where):
        out.append((name, expr.tends_to_plus_infinity(), where))

    def small(name, expr: AExp, where):
        out.append((name, expr.tends_to_zero(), where))

    ll = LAM
    # ---------------- hybrid restricted BV (Lemma 8.5H), y1 = B, y2 = B2, x in the grid around X_F
    chk("(H1) log B - C (loglog x)^2", s["logB"] - AExp.mono(C, 0, 2), "Steps A (both y), E6 (beta >= 3/4)")
    chk("(H1') log B - C loglog x", s["logB"] - ll.scale(C), "K (1 - beta(Q', y1)) <= delta_D/4")
    chk("(H2) u_F2 - (loglog x)^4", s["u_F2"] - AExp.mono(1, 0, 4), "Drappeau saving beats polylog prefactors")
    chk("(H2') log x - C log B2", s["logXF"] - s["logB2"].scale(C), "R1 = y2^K <= (x/d)^{eta_D}")
    # (H3): u log(4u log y) <= K0 log y2 on both sides; u log u ~ (3/8) u lam
    chk("(H3a) log B2 - C u_F2 lam", s["logB2"] - AExp({(k[0], k[1] + 1): v * C for k, v in s["u_F2"].t.items()}),
        "density floor rho(x, y2) >= y2^{-1-K0}")
    chk("(H3b) log B2 - C u_r lam", s["logB2"] - AExp({(k[0], k[1] + 1): v * C for k, v in s["u_r"].t.items()}),
        "density floor rho(Q', y1) >= y2^{-1-K0}")
    # (H4) the new transfer condition: K (log y2/log y1)(1 + C_alpha) lam <= c u_F2/log^2 u_F2
    # log y2/log y1 = (b2/a) D^{1/8};  times lam
    ratio_lam = Dl(Fr(b2) / Fr(a) * C, Fr(1, 8))                 # C (b2/a) D^{1/8} lam
    rhs = AExp({(k[0], k[1] - 2): v for k, v in s["u_F2"].t.items()})   # u_F2 / lam^2 (log u ~ lam)
    chk("(H4) u_F2/log^2 u - C (logB2/logB) loglog x", rhs - ratio_lam, "Step C(a): R1^{1-beta(Q',y1)} absorbed")
    chk("(H5a) log Y - C log B2", s["logY"] - s["logB2"].scale(C), "l > y2 (lift averaging), Part 2 large sieve")
    chk("(H5b) log(Q/Y) - C log B2", s["logXr"] - s["logB2"].scale(C), "level gap, lower side")
    chk("(H5c) log(sqrt x/Q) - C log B2", s["logXF"].scale(Fr(1, 2)) - s["logXr"] - s["logB2"].scale(C),
        "level gap, upper side")
    chk("saving beats polylog: u_F2/lam^2 - C lam", rhs - ll.scale(C), "FM bad mass: |X| eps (log N)^C -> 0")
    # ---------------- GM-U' / DS-U' / U1 at y2
    small("U1 at y2: (1 - alpha2) log varpi^2 -> 0",
          AExp({(Fr(3, 8) - Fr(1, 2), Fr(5, 8) + 1 - Fr(1, 2)): Fr(4) * Fr(a) + 2 * sigma1(a, d0, kappa, b2)}),
          "shifts by varpi^2 = B^4 S1^4 cost 1 + o(1) at y = B2")
    chk("B^C <= B2 (HT Thm 3 needs c <= y)", s["logB2"] - s["logB"].scale(C), "Lemma U1(ii) at y2 in one step")
    # ---------------- mixed scan
    chk("d log S1 << log X_F", s["logXF"] - AExp({(k[0] + Fr(1, 4), k[1] - Fr(1, 4)): v * C * Fr(d0)
                                                for k, v in s["logS1"].t.items()}), "S^d inside X_F (CS existence)")
    chk("2^{d+1} <= B (doubling S)", s["logB"] - s["d"].scale(C), "m <= 64 step of CS")
    # ---------------- squarefree / discriminant union bound
    out.append(("(SQ) a > 3 d0/(32 kappa b2)", sq_constraint(a, d0, kappa, b2), "Item 3: 1/(B log B) << rho_F2"))
    chk("B2^2 <= (2H+1)/d", s["logH"] - s["logB2"].scale(2) - ll.scale(C), "Pr_c[p^2 | F(P)] <= 1/p^2 + 1/(2H+1)")
    # ---------------- level-1 heights / Kern / Lemma W at eta0 = D^{-1/2}
    chk("eta0 logN/2 - C logS1 - C logB", s["eta0_logN"].scale(Fr(1, 2)) - s["logS1"].scale(C) - s["logB"].scale(C),
        "Kern-U H >= dV1^2 N^{eta0/2}; GM-U errors; Lemma W-U")
    chk("eta0 logN - C d", s["eta0_logN"] - s["d"].scale(C), "2^{2d+5} absorbed by N^{-eta0/4}")
    chk("log X_F - C eta0 logN", s["logXF"] - s["eta0_logN"].scale(C), "eta0 = o(theta)")
    # ---------------- recursion / IR / split
    ir = AExp({(Fr(5, 8) - Fr(1, 4), Fr(3, 8) + Fr(1, 4) - 1): Fr(1) / (2 * Fr(kappa)) / (13 * Fr(d0))})
    chk("Lemma IR: logH/(13 d loglogH) - C lam", ir - ll.scale(C), "Pr(rejection) = o(gamma_2)")
    # value factorisation by Theorem U on X <= X_F: log L_X[2/5] ~ (log X)^{2/5}(loglog)^{3/5} = O(D^{7/20} lam)
    vf = Dl(1, Fr(7, 20))
    chk("Lambda - C log L_{X_F}[2/5]", s["Lambda"] - vf.scale(C), "recursive value factorisations are L^{o(1)}")
    chk("log Y - omega(N) - C lam", s["logY"] - s["omega_upper"] - ll.scale(C), "Lemma S'' hypothesis")
    chk("log V1 - log B2 (b prime > B2)", s["logY"].scale(Fr(1, 2)) - s["logB2"].scale(C), "q never divides b")
    # ---------------- child boxes
    chk("log H - log B2 - 2 log S2 - C lam", s["logH"] - s["logB2"] - s["logS2"].scale(2) - ll.scale(C),
        "child pairs: Kern-U at level 2 (H >= d V2^2)")
    small("mbar/log N -> 0 (child cap polylog)", AExp({(k[0] - 1, k[1]): v for k, v in s["mbar"].t.items()}),
          "deterministic child cap O(D^{1/2})")
    return out


def all_margins_hold(a, d0, kappa, b2) -> bool:
    return all(ok for _, ok, _ in margins(a, d0, kappa, b2))


# --- cost constant ------------------------------------------------------------------------------------------
def two_pow_3_8_bracket(den: int = 10 ** 6) -> Tuple[Fr, Fr]:
    """Rationals lo < 2^{3/8} < hi, hi - lo = 1/den, decided exactly via x^8 vs 2^3."""
    lo_n = 1296839 * den // 10 ** 6
    while Fr(lo_n + 1, den) ** 8 < 8:
        lo_n += 1
    while Fr(lo_n, den) ** 8 > 8:
        lo_n -= 1
    return Fr(lo_n, den), Fr(lo_n + 1, den)


def cost_constant_2L(a, d0, kappa, b2, C_gamma) -> Tuple[Fr, Fr]:
    """Bracket for c = max(2a, a + s + max(a, sigma1 + sigma2), (a + s)/(1 - 2^{-3/8})), s = 2 kappa (1 + C_gamma).
    (Theorem 2L's Lemma F_alpha recursion with the special-q search W1 (1 + mbar W2) = L[3/8, sigma1 + sigma2].)"""
    a, d0, kappa, b2, C_gamma = map(Fr, (a, d0, kappa, b2, C_gamma))
    s = 2 * kappa * (1 + C_gamma)
    lo_t, hi_t = two_pow_3_8_bracket()
    flo, fhi = hi_t / (hi_t - 1), lo_t / (lo_t - 1)        # t/(t-1) decreasing in t
    p = max(2 * a, a + s + max(a, sigma1(a, d0, kappa, b2) + sigma2(a, d0, kappa, b2)))
    return max(p, (a + s) * flo), max(p, (a + s) * fhi)
