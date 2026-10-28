"""Exact asymptotic bookkeeping for Theorem U-MF (notes/mfU_theorem.md §§2, 10).

Design MF with integer multipliers: mu = 1/d, h = nu + eta0 + h_MF (window floor at mu = 1/d plus the MF height floor
log H >= 30 (sigmabar + 1) d Lambda loglog N), b in (V/2, V], one box radius V, S^2 ~ 1/rho_r, J_f ~ 1/rho_F.

Leading orders (D = log N/loglog N, Lambda = D^{2/5} loglog N):
  log X_r ~ (1/d0) D^{4/5} loglog N,   log X_F ~ (1/d0 + d0/(4 kappa)) D^{4/5} loglog N,
  sigma_r = 2/(5 a d0),  sigma_F = (2/(5a)) (1/d0 + d0/(4 kappa)),  sigmabar = sigma_r + sigma_F.
Cost (Prop select + Prop PP + Lemma F_alpha):
  c_MF = max(2(a + sigma_F), a + sigma_F + s + max(a, sigma_r, sigma_F), (a + sigma_F + s)/(1 - 2^{-2/5})),
  s = 2 kappa (1 + C_gamma).
Only exact rationals are used for decisions; 2^{2/5} is bracketed exactly (imported from lv5unif).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple

from src.lv5unif.params import AExp, LAM, LOGLOGN, LOGN, D_pow, recursion_factor_bracket
from src.lv5unif.params import cost_constant as cost_constant_U


def sigma_r(a, d0) -> Fr:
    return Fr(2) / (5 * Fr(a) * Fr(d0))


def sigma_F(a, d0, kappa) -> Fr:
    return Fr(2) / (5 * Fr(a)) * (Fr(1) / Fr(d0) + Fr(d0) / (4 * Fr(kappa)))


def sigma_bar(a, d0, kappa) -> Fr:
    return sigma_r(a, d0) + sigma_F(a, d0, kappa)


def sizes(a, d0, kappa) -> Dict[str, AExp]:
    a, d0, kappa = Fr(a), Fr(d0), Fr(kappa)
    sb = sigma_bar(a, d0, kappa)
    return {
        "logB": LAM.scale(a),
        "logY": AExp.mono(Fr(1) / (2 * kappa), Fr(3, 5), Fr(2, 5)),                 # log N/(2k)
        "nu_logN": AExp.mono(Fr(1) / (4 * kappa), Fr(3, 5), Fr(2, 5)),
        "mu_logN": AExp.mono(Fr(1) / d0, Fr(4, 5), Fr(1, 5)),                        # log N/d
        "eta0_logN": AExp.mono(1, Fr(1, 2), Fr(1, 2)),
        "logS2": LAM.scale(sigma_r(a, d0)),
        "logJf": LAM.scale(sigma_F(a, d0, kappa)),
        "log_rho_r_inv": LAM.scale(sigma_r(a, d0)),
        "log_rho_F_inv": LAM.scale(sigma_F(a, d0, kappa)),
        "logXr": AExp.mono(Fr(1) / d0, Fr(4, 5), Fr(1, 5)),
        "logXF": AExp.mono(Fr(1) / d0 + d0 / (4 * kappa), Fr(4, 5), Fr(1, 5)),
        "d": AExp.mono(d0, *D_pow(Fr(1, 5))),
        "d_logd": AExp.mono(d0 / 5, Fr(1, 5), Fr(4, 5)),
        "d_logS": AExp.mono(d0 * sigma_r(a, d0) / 2, Fr(3, 5), Fr(2, 5)),            # d * (sigma_r/2) Lambda
        "hMF_logN": AExp.mono(30 * (sb + 1) * d0, Fr(3, 5), Fr(7, 5)),               # 30(sb+1) d Lambda loglogN
        "omega_upper": AExp.mono(Fr(1) / a, Fr(3, 5), Fr(-3, 5)),                    # log N/log B
    }


def margins(a, d0, kappa, C=Fr(100)) -> List[Tuple[str, bool, str]]:
    """Leading-order checks of Lemma P (notes/mfU_theorem.md §2) and of the margins used in Lemmas OP, REJ, CS, FIB.
    Each entry: (name, holds, where used).  'holds' means: the log-ratio tends to +infinity."""
    s = sizes(a, d0, kappa)
    C = Fr(C)
    sb = sigma_bar(a, d0, kappa)
    out = []

    def chk(name, e: AExp, where):
        out.append((name, e.tends_to_plus_infinity(), where))

    chk("P1: log Y - C log B", s["logY"] - s["logB"].scale(C), "Lemma 4.4, 4.5, (P1)")
    chk("P1: log Y - C log S^2", s["logY"] - s["logS2"].scale(C), "2S^2 < Y; b > 4S^2 >= q in Prop R2")
    chk("P2: eta0 logN/4 - C log B", s["eta0_logN"].scale(Fr(1, 4)) - s["logB"].scale(C), "Lemma 4.3, Lemma OP")
    chk("P2: eta0 logN/4 - C log S^2", s["eta0_logN"].scale(Fr(1, 4)) - s["logS2"].scale(C), "W >= 2MV N^{eta0/3}")
    chk("P2: eta0 logN/4 - C d", s["eta0_logN"].scale(Fr(1, 4)) - s["d"].scale(C), "2^{d+5} <= N^{eta0/4}")
    chk("OP: eta0 logN/2 - log(1/rho_F) - C loglogN",
        s["eta0_logN"].scale(Fr(1, 2)) - s["log_rho_F_inv"] - LOGLOGN.scale(C), "Lemma OP: 3d^2V/H << rho_F")
    chk("P3: (mu - nu) logN - eta0 logN (M/V)", s["mu_logN"] - s["nu_logN"] - s["eta0_logN"].scale(C), "Lemma FIB")
    chk("FIB: Moebius 2^omega vs M/V", s["mu_logN"].scale(Fr(1, 2)) - AExp.mono(C / Fr(d0), Fr(4, 5), Fr(1, 5) - 1),
        "Lemma FIB: log 2^{omega(t)} = O(mu logN/loglogN)")
    chk("P4: Lambda - C d log d", LAM - s["d_logd"].scale(C), "characters, SQ sign vectors")
    chk("P4: theta logN - C d log S", s["logXF"] - s["d_logS"].scale(C), "S^d inside X_F is lower order")
    chk("P4: theta logN - C h_MF logN", s["logXF"] - s["hMF_logN"].scale(1), "MF height floor is o(theta log N)")
    chk("P5: log Y - omega log2 - C loglogN", s["logY"] - s["omega_upper"] - LOGLOGN.scale(C), "Lemma S'' hypothesis")
    # P6: log H/(24 d loglog H) >= 1.25 (sb+1) Lambda: log H >= 30(sb+1) d Lambda loglogN and loglog H <= loglog N
    ir_exp = AExp.mono(30 * (sb + 1) * Fr(d0) / (24 * Fr(d0)), Fr(2, 5), Fr(3, 5))   # = 1.25 (sb+1) Lambda
    chk("P6: IR exponent - log(1/rho_F) - C loglogN", ir_exp - s["log_rho_F_inv"] - LOGLOGN.scale(C),
        "Lemma REJ: rejection << rho_F (the MF height floor)")
    # without the floor: log H ~ nu log N, IR exponent ~ nu logN/(24 d loglogN) ~ D^{2/5}/(96 kappa d0): loses to rho_F
    no_floor = AExp.mono(Fr(1) / (96 * Fr(kappa) * Fr(d0)), Fr(2, 5), Fr(-2, 5))  # D^{2/5} = (logN)^{2/5}(loglogN)^{-2/5}
    out.append(("negative control: without the MF floor the IR exponent loses to log(1/rho_F)",
                not (no_floor - s["log_rho_F_inv"]).tends_to_plus_infinity(), "codex lv5unif_2 #3"))
    chk("P6: log N/(2d) - (log H/2) log 2", AExp.mono(Fr(1) / (2 * Fr(d0)), Fr(4, 5), Fr(1, 5))
        - s["hMF_logN"].scale(1) - s["nu_logN"], "Lemma REJ: 2^{|P_IR|} <= W N^{-1/(2d)}")
    chk("P7: CEP range for X_F", AExp.mono((Fr(1) / Fr(d0) + Fr(d0) / (4 * Fr(kappa))) / 2, Fr(4, 5), Fr(-4, 5))
        - AExp.mono((Fr(1) / Fr(d0) + Fr(d0) / (4 * Fr(kappa))) / Fr(a) * 2, Fr(2, 5), Fr(-2, 5)), "Lemma U1(iv)")
    chk("CS: J_f slack (logN)^{loglogloglogN} is o(Lambda)", LAM - AExp.mono(C, 0, 2), "Lemma CS (lambda_3 <= loglogN)")
    return out


def cost_constant(a, d0, kappa, C_gamma) -> Tuple[Fr, Fr]:
    """Exact bracket [lo, hi] for c_MF (notes/mfU_theorem.md Theorem U-MF)."""
    a, d0, kappa, C_gamma = Fr(a), Fr(d0), Fr(kappa), Fr(C_gamma)
    sr, sF = sigma_r(a, d0), sigma_F(a, d0, kappa)
    s = 2 * kappa * (1 + C_gamma)
    tq = max(a, sr, sF)
    p = max(2 * (a + sF), a + sF + s + tq)
    flo, fhi = recursion_factor_bracket()
    return max(p, (a + sF + s) * flo), max(p, (a + sF + s) * fhi)


def optimise(C_gamma, grid_a=None, grid_kappa=None, grid_d0=None) -> Dict:
    best = None
    grid_a = grid_a or [Fr(i, 20) for i in range(4, 61)]
    grid_kappa = grid_kappa or [Fr(i, 400) for i in range(1, 241)]
    grid_d0 = grid_d0 or [Fr(i, 8) for i in range(1, 49)]
    for kappa in grid_kappa:
        for d0 in grid_d0:
            for a in grid_a:
                lo, hi = cost_constant(a, d0, kappa, C_gamma)
                if best is None or hi < best["hi"]:
                    best = {"a": a, "kappa": kappa, "d0": d0, "lo": lo, "hi": hi,
                            "sigma_r": sigma_r(a, d0), "sigma_F": sigma_F(a, d0, kappa)}
    return best


if __name__ == "__main__":  # pragma: no cover
    import json
    import sys
    res = {}
    for name, ok, where in margins(Fr(3, 5), Fr(1), Fr(1, 8)):
        print(("OK  " if ok else "FAIL") + "  " + name + "   [" + where + "]")
    for Cg in (0, 1, 5, 20):
        b = optimise(Fr(Cg))
        res[str(Cg)] = {k: str(v) for k, v in b.items()}
        print(Cg, {k: (float(v) if isinstance(v, Fr) else v) for k, v in b.items()})
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as f:
            json.dump(res, f, indent=1)
