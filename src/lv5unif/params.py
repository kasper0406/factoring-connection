"""Exact asymptotic bookkeeping for Theorem U (rigorous L[2/5]; notes/lv5unif_results.md §§1, 10).

Every log-size is an `AExp`: a finite sum  sum_j c_j (log N)^{p_j} (log log N)^{q_j}  with Fraction
coefficients and exponents.  A quantity X "beats" Y (X/Y -> infinity) iff the leading term of
log X - log Y has positive coefficient, leading = lexicographically largest (p, q).  This is exactly how
the margins are compared in the notes (every margin there is a ratio exp(c (log N)^p (loglog N)^q)).

Parameters (D := log N / log log N, Lambda := Lambda_{2/5} = (log N)^{2/5} (log log N)^{3/5} = D^{2/5} log log N):
  d   ~ d0 * D^{1/5}   (odd),        k ~ kappa * D^{2/5},      nu = 1/(4k),
  B   = exp(a Lambda),               eta0 = D^{-1/2},
  theta * D^{1/5} -> theta0(d0, kappa) = 3/d0           if kappa >= d0^2/4   (level+window floors bind)
                                       = 2/d0 + d0/(4 kappa) otherwise     (Kern+window floors bind)
  log(1/(rho_r rho_F)) = sigma Lambda (1 + o(1)),  sigma = 2 theta0 / (5 a)        [CEP, u ~ D^{2/5}]
  S^2 = (log N)^{C3} / (rho_r rho_F)  (Lemma CS)  =>  log S^2 = sigma Lambda (1 + o(1)).
Split acceptance per attempt >= (gamma/(16 log N))^{2k+1}/(128 k^2), gamma >= (log N)^{-C_gamma}
  => log(1/p) = 2 kappa (1 + C_gamma) Lambda (1 + o(1)).
Cost (Prop select-U + Lemma F_alpha):  c = max( 2a,  a + s + max(a, sigma),  (a + s)/(1 - 2^{-2/5}) ),
  s := 2 kappa (1 + C_gamma).
No floating point is used for any decision; 2^{2/5} is bracketed by exact rationals.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple

Key = Tuple[Fr, Fr]


class AExp:
    """sum c (log N)^p (loglog N)^q, exact."""

    def __init__(self, terms: Dict[Key, Fr] | None = None):
        self.t: Dict[Key, Fr] = {}
        for k, v in (terms or {}).items():
            if v != 0:
                self.t[(Fr(k[0]), Fr(k[1]))] = self.t.get((Fr(k[0]), Fr(k[1])), Fr(0)) + Fr(v)
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

    def times_mono(self, p, q) -> "AExp":
        return AExp({(k[0] + Fr(p), k[1] + Fr(q)): v for k, v in self.t.items()})

    def leading(self) -> Tuple[Key, Fr] | None:
        if not self.t:
            return None
        k = max(self.t)
        return k, self.t[k]

    def tends_to_plus_infinity(self) -> bool:
        """True iff the leading term is positive and of positive order (p, q) > (0, 0)."""
        ld = self.leading()
        return ld is not None and ld[1] > 0 and ld[0] > (Fr(0), Fr(0))

    def __repr__(self) -> str:
        return " + ".join(f"{v}*(logN)^{k[0]}*(loglogN)^{k[1]}" for k, v in sorted(self.t.items(), reverse=True))


# --- basic sizes ------------------------------------------------------------------------------------------
LOGN = AExp.mono(1, 1, 0)
LOGLOGN = AExp.mono(1, 0, 1)
LAM = AExp.mono(1, Fr(2, 5), Fr(3, 5))          # Lambda_{2/5}


def D_pow(s) -> Key:
    """D^s = (log N)^s (log log N)^{-s}, as an exponent pair."""
    return (Fr(s), -Fr(s))


def theta0(d0: Fr, kappa: Fr) -> Fr:
    d0, kappa = Fr(d0), Fr(kappa)
    if kappa >= d0 * d0 / 4:
        return Fr(3) / d0
    return Fr(2) / d0 + d0 / (4 * kappa)


def sigma(a: Fr, d0: Fr, kappa: Fr) -> Fr:
    return 2 * theta0(d0, kappa) / (5 * Fr(a))


def sizes(a: Fr, d0: Fr, kappa: Fr) -> Dict[str, AExp]:
    a, d0, kappa = Fr(a), Fr(d0), Fr(kappa)
    th0 = theta0(d0, kappa)
    sg = sigma(a, d0, kappa)
    p, q = D_pow(Fr(-1, 5))
    out = {
        "logN": LOGN,
        "Lambda": LAM,
        "logB": LAM.scale(a),
        "log_ell": AExp.mono(Fr(1) / (2 * kappa), Fr(3, 5), Fr(2, 5)),      # 2 nu log N = log N/(2k)
        "nu_logN": AExp.mono(Fr(1) / (4 * kappa), Fr(3, 5), Fr(2, 5)),
        "theta_logN": AExp.mono(th0, 1 + p, q),                               # theta0 D^{-1/5} log N
        "eta0_logN": AExp.mono(1, Fr(1, 2), Fr(1, 2)),                        # D^{-1/2} log N
        "logS2": LAM.scale(sg),
        "d": AExp.mono(d0, *D_pow(Fr(1, 5))),                                 # d itself (a size, not a log)
        "k": AExp.mono(kappa, *D_pow(Fr(2, 5))),
        "d_logd": AExp.mono(d0 / 5, Fr(1, 5), Fr(4, 5)),                      # d log d ~ d0 D^{1/5} (1/5) loglogN
        "omega_log2_upper": AExp.mono(Fr(1) / a, Fr(3, 5), Fr(-3, 5)),        # omega(N) <= log N/log B (times log 2 <= 1)
        "u": AExp.mono(th0 / a, *D_pow(Fr(2, 5))),                            # total u_r + u_F
    }
    return out


def margins(a: Fr, d0: Fr, kappa: Fr, C: Fr = Fr(100)) -> List[Tuple[str, bool, str]]:
    """The growth comparisons listed in notes/lv5unif_results.md (Lemma P and the margins named in §§4-6, 11), with a
    generous constant C where the notes say 'B^{O(1)}' or '(log N)^{O(1)}'.  Returns (name, holds, where-used).
    Leading-order only: each check compares the leading terms of two log-sizes."""
    s = sizes(a, d0, kappa)
    C = Fr(C)
    th0 = theta0(d0, kappa)
    checks = []

    def chk(name, expr: AExp, where):
        checks.append((name, expr.tends_to_plus_infinity(), where))

    # Lemma 8.5'(R'): Y >= y^K for fixed K;  x^lambda >= y^{O(1)}  (lambda log x = eta0 log N / 4)
    chk("log Y - C log B", s["log_ell"] - s["logB"].scale(C), "Lemma 8.5'(R') Part 2, L2_fam class (C)")
    chk("eta0 logN/4 - C log B", s["eta0_logN"].scale(Fr(1, 4)) - s["logB"].scale(C), "8.5'(R'4), FM errors, Lemma W, window floor")
    # u >= (loglog x)^4 for Drappeau's saving to beat (log x)^{O(K)}
    chk("u - (loglogN)^4", s["u"].scale(Fr(1, 2)) - AExp.mono(1, 0, 4), "Lemma 8.5'(R'2) (u_F >= u/2 up to constants)")
    # CEP validity u <= (1-eps) log x / loglog x : log x/loglog x ~ theta0 D^{4/5}
    chk("log x/loglog x - 2u", AExp.mono(th0 / 2, Fr(4, 5), Fr(-4, 5)) - s["u"].scale(2), "CEP for rho asymptotics (cost only)")
    # characters: (2d)! bits cost exp(O(d log d)) = L^{o(1)}
    chk("Lambda - C d log d", LAM - s["d_logd"].scale(C), "LV6-U characters, Lemma SQ sign vectors 2^{d-1}")
    # box: S^d = N^{o(theta)}:  d log S << theta log N
    chk("theta logN - C d logS", s["theta_logN"] - AExp.mono(d0 * sigma(a, d0, kappa) / 2 * C, Fr(3, 5), Fr(2, 5)), "R4: norms N^{theta(1+o(1))}")
    # eta0 = o(theta): theta log N - C eta0 log N
    chk("theta logN - C eta0 logN", s["theta_logN"] - s["eta0_logN"].scale(C), "eta0 slack is o(theta) in the constant")
    # eta0 log N >> d (window floor with 2^{2d+6}; H >= d V^2 margins)
    chk("eta0 logN - C d", s["eta0_logN"] - s["d"].scale(C), "window floor constant 2^{2d+6}, Kern-U H >= dV^2")
    # Lemma S'': |D_A| >= 4k^3 2^omega:  log Y - omega log 2 - C loglogN
    chk("log Y - omega log2 - C loglogN", s["log_ell"] - s["omega_log2_upper"] - LOGLOGN.scale(C), "Lemma S'' hypothesis")
    # Lemma IR: log H/(13 d loglog H) >> C_gamma loglog N; log H >= 2 nu log N
    ir = AExp.mono(Fr(1) / (2 * kappa) / (13 * d0), Fr(3, 5) - Fr(1, 5), Fr(2, 5) + Fr(1, 5) - 1)
    chk("logH/(13 d loglogH) - C loglogN", ir - LOGLOGN.scale(C), "Lemma IR: Pr(rejection) = o(gamma)")
    # Theorem K: split cost exactly of order Lambda: 2k loglogN = 2 kappa Lambda
    two_k_ll = s["k"].scale(2).times_mono(0, 1)
    checks.append(("2k loglogN == 2 kappa Lambda", (two_k_ll - LAM.scale(2 * Fr(kappa))).leading() is None,
                   "split cost enters c (not o(1))"))
    # mu - nu >= eta0 : mu ~ 1/d ~ D^{-1/5}/d0 >> eta0 = D^{-1/2}
    chk("mu logN - C eta0 logN", AExp.mono(Fr(1) / d0, Fr(4, 5), Fr(1, 5)) - s["eta0_logN"].scale(C), "FM step 1, Lemma 8.5' q0 - 2nu >= eta0")
    # class (C) of L2_fam first term N^{-2nu} B^{O(1)}
    chk("2 nu logN - C logB", s["log_ell"] - s["logB"].scale(C), "L2_fam (C) first term")
    # Lemma 8.5'(R'): log y >> loglog x (so 1 - beta <= 1/4 etc.)
    chk("logB - C loglogN", s["logB"] - LOGLOGN.scale(C), "(alpha1), (alpha2), HT range")
    # --- added after codex lv5unif_2 #6 ---
    # level gap at the smallest grid point: Gamma >= (eta0/2) log N - O(Lambda) must beat (2K0+U0) log B + 4 loglog N
    chk("Gamma - C logB - C loglogN", s["eta0_logN"].scale(Fr(1, 2)) - LAM.scale(C) - s["logB"].scale(C) - LOGLOGN.scale(C),
        "Lemma 8.5'(R'4) level gap, every x in the FM grid")
    # (R'3) is a fixed-coefficient comparison: u log(4u log B)/log B -> (4/5) theta_F0/a^2 <= theta0/a^2 < K0
    K0 = 1 + (th0 / (Fr(a) * Fr(a))).__ceil__()
    checks.append(("(R'3): (4/5) theta0/a^2 < K0", Fr(4, 5) * th0 / (Fr(a) * Fr(a)) < K0, "Lemma 8.5'(R'3)"))
    # P10 heights: (h - 2 nu) log N >= eta0 log N must beat log d + log S^2 + (eta0/2) log N  and log(3 d varpi)
    chk("eta0 logN/2 - C logS2 - C logB", s["eta0_logN"].scale(Fr(1, 2)) - s["logS2"].scale(C) - s["logB"].scale(C),
        "P10: H >= d V^2 N^{eta0/2}, H >= 3 d varpi V (varpi = B^2 S^2)")
    # class (B) 2nd ratio N^{-h-(d+1)nu} B^{O(1)} d^{O(1)}: (h+(d+1)nu) log N >= 2 nu log N
    chk("(h+(d+1)nu)logN - C logB", s["log_ell"] - s["logB"].scale(C), "L2_fam class (B) second ratio")
    # class (C) 2nd ratio N^{(5-d)nu - h} B^{O(1)}: for d >= 5 the exponent is <= -h <= -2 nu
    chk("(h-(5-d)nu)logN - C logB  [d>=5]", s["log_ell"] - s["logB"].scale(C), "L2_fam class (C) second ratio")
    # primitivity (Moebius) error 2^{omega(l n)} <= exp(C mu log N/loglog N) against M/V = N^{mu - nu}
    chk("(mu-nu)logN - C mu logN/loglogN", AExp.mono(Fr(1) / d0, Fr(4, 5), Fr(1, 5)).scale(Fr(1, 2))
        - AExp.mono(C / d0, Fr(4, 5), Fr(1, 5) - 1), "FM fibres: Moebius error")
    # MF: with the extra floor log H >= 13 (sigma+1) d Lambda loglog N, rejection <= exp(-(sigma+1) Lambda (1+o(1)))
    sg = sigma(a, d0, kappa)
    chk("MF: logH/(13 d loglogH) - sigma Lambda - C loglogN", LAM.scale(sg + 1) - LAM.scale(sg) - LOGLOGN.scale(C),
        "Prop MF: Pr(rejected) <= p1/2 (with the extra height floor)")
    chk("MF extra floor is o(theta logN)", s["theta_logN"] - AExp.mono(13 * (sg + 1) * d0 * C, Fr(3, 5), Fr(7, 5)),
        "Prop MF: extra floor does not move the exponent")
    return checks


# --- 2^{2/5} bracket and the cost constant -------------------------------------------------------------------
def two_pow_2_5_bracket(den: int = 10 ** 6) -> Tuple[Fr, Fr]:
    """Rationals lo < 2^{2/5} < hi with hi - lo = 1/den, decided exactly via x^5 vs 4."""
    lo_n = 1319507 * den // 10 ** 6
    while Fr(lo_n + 1, den) ** 5 < 4:
        lo_n += 1
    while Fr(lo_n, den) ** 5 > 4:
        lo_n -= 1
    return Fr(lo_n, den), Fr(lo_n + 1, den)


def recursion_factor_bracket() -> Tuple[Fr, Fr]:
    """Bracket for 1/(1 - 2^{-2/5}) = t/(t - 1), t = 2^{2/5} (decreasing in t)."""
    lo_t, hi_t = two_pow_2_5_bracket()
    return hi_t / (hi_t - 1), lo_t / (lo_t - 1)


def cost_constant(a: Fr, d0: Fr, kappa: Fr, C_gamma: Fr) -> Tuple[Fr, Fr]:
    """Exact bracket [lo, hi] for c = max(2a, a + s + max(a, sigma), (a + s)/(1 - 2^{-2/5})), s = 2 kappa (1 + C_gamma)."""
    a, kappa, C_gamma = Fr(a), Fr(kappa), Fr(C_gamma)
    sg = sigma(a, d0, kappa)
    s = 2 * kappa * (1 + C_gamma)
    flo, fhi = recursion_factor_bracket()
    p = max(2 * a, a + s + max(a, sg))
    return max(p, (a + s) * flo), max(p, (a + s) * fhi)


def optimise(C_gamma: Fr, grid_a=None, grid_kappa=None) -> Dict:
    """Exact grid search (rational grid) for the smallest upper bracket of c; d0 chosen as the best rational
    in a small grid around sqrt(8 kappa) (Kern regime optimum) and 2 sqrt(kappa) (boundary)."""
    best = None
    grid_a = grid_a or [Fr(i, 20) for i in range(1, 61)]
    grid_kappa = grid_kappa or [Fr(1, 2 ** j) for j in range(0, 14)] + [Fr(i, 400) for i in range(1, 41)]
    for kappa in grid_kappa:
        cands = set()
        for num in range(1, 400):
            d0 = Fr(num, 20)
            # restrict to d0 near the optimum sqrt(8 kappa) (compare squares exactly)
            if Fr(1, 4) * 8 * kappa <= d0 * d0 <= 4 * 8 * kappa:
                cands.add(d0)
        for d0 in cands:
            for a in grid_a:
                lo, hi = cost_constant(a, d0, kappa, C_gamma)
                if best is None or hi < best["hi"]:
                    best = {"a": a, "kappa": kappa, "d0": d0, "lo": lo, "hi": hi,
                            "theta0": theta0(d0, kappa), "sigma": sigma(a, d0, kappa)}
    return best


# --- the tau table (data; the proofs are in the notes) ------------------------------------------------------
TAU_TABLE = [
    # ingredient, where d,k enter (as stated in the sources), exponential?, removable by, tau
    ("FM (first moment)", "#P_l ~ V^2/(d log V); one-point smoothing TV 3d^2V/H, |M0| <= dV^{d-1}; "
     "X_F/(t y') <= 2^{O(d)} varpi^2 (rho-shift costs exp(O(d loglogN/Lambda)) = 1+o(1)); C_d=2^{2d+5} only "
     "against N^{-eta0/2}", "no (after b-window + common V)", "b in (V(1-1/d),V], common V, X_F=(2d+1)MHV^d", 0),
    ("Lemma 8.5' (restricted BV)", "none: regime (R') hypotheses (u >= (loglog x)^4, rho >= y^{-K0}, Y >= y^K, "
     "x^lambda >= y^{2K0+6}) hold at y = L[2/5], x = N^{theta}, Y = L[3/5]", "no", "restate regime", 0),
    ("Lemma Kern", "C_d of the cube-section volume ((d+1)H)^{d-2} and (V^2/bb')^{d-2} <= 4^{d-2}",
     "yes as written (d^{O(d)}, 4^d)", "triangular fibre count: <= e/((2H+1)^2 (bb')^{d-2}) if bb' <= (2H+1)/d; b-window", 0),
    ("L2_fam", "(log N)^{O(C0)}, C0 = O(log(1/rho)/log B) fixed; tau(q); margins N^{-h}, N^{-2nu}, N^{(5-d)nu-h}",
     "no", "-", 0),
    ("LD* (second moment)", "C_d (Kern), (V_max/V_-)^{2d} <= 2^d, X_F constant 2^{2d+5} squared, <= 12 partners, "
     "C_W = 24/eta0", "yes as written (2^{O(d)})", "common V; X_F=(2d+1)MHV^d; Kern-U; eta0 = D^{-1/2} (C_W polylog)", 0),
    ("Lemma W", "C_W = 24/eta0", "no (polylog for eta0 = D^{-1/2})", "-", 0),
    ("Theorem K / Lemma S''", "(gamma/(16 log N))^{2k+1}/(128k^2): the (log N)^{2k} is the design's split cost, "
     "of order L[2/5] by choice of k; no C(2k,k), no 2^omega", "k only via the split cost (exact order Lambda)", "S'' replaces S'", 0),
    ("LV6", "Stark nu_1 ~ 1/((2d)! log Delta): character primes of (2d)!(log N)^{O(1)} bits; SQ 2^{d-1} sign vectors; "
     "t_chi = O(d^2 log N)", "cost only: exp(O(d log d)) = L[2/5]^{o(1)}", "large character primes", 0),
    ("Lemma IR / rejections", "exp(-log H/(13 d loglog H)) = exp(-(log N)^{2/5-o(1)}) = o(gamma)", "no", "-", 0),
]
