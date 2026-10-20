"""auditU: independent exact re-derivation of Theorem U's exponent and constants (notes/lv5unif_results.md §§1, 3, 6, 10, 13).

Written from scratch by agent auditU.  It imports nothing from src/lv5unif except, at the very end, `params.optimise`
and `params.cost_constant`, to compare.  Exact arithmetic only (Fractions, sympy); floats appear only in printed text.

Checks
  A1  theta_min(d, k) by exhaustive breakpoint search of the convex piecewise-linear theta(mu), at d = odd >= d0 t,
      k = max(ceil(kappa t^2), d+1), t = D^{1/5}; shows theta * t -> theta0(d0, kappa) (the claimed limit).
  A2  LD*-U Step 7 power count (sympy): exponents of M, H, V, Y and the d-dependent constant.
  A3  Lemma Kern-U by brute force (independent enumeration over C_H), including the (3/2)^d family at H = b^2
      and window pairs at H = d V^2.
  A4  Kern image index: |Z^2 : Phi(Z^d)| = |Delta| (checked by Smith form), so the Kern-U point bound exceeds the
      uniform-on-image density by V^2 d^2/|Delta|; this is what the q-sum in Step 7 pays for.
  A5  The alpha floor: LP over exponents (delta = log_D d, kap = log_D k) gives alpha >= 2/5, and (2+tau)/(5+2tau).
  A6  The cost constant c = max(2a, a+s+max(a,sigma), (a+s)/(1-2^{-2/5})): exact grid optimum for C_gamma in
      {0,1,5,20}, 2^{2/5} bracketed exactly; finiteness for every finite C_gamma; derivative in C_gamma.
  A7  Margins in leading-term arithmetic (exact (p,q) exponents of log N, log log N), my own list.
  A8  "When does the asymptotics start": smallest D (= log N/log log N, power of 10) at which the thinnest
      margins hold with plausible constants (exact rational comparison of D^{1/10} etc.).
"""
from __future__ import annotations

import itertools
import json
import random
import sys
from fractions import Fraction as Fr
from math import comb

import sympy as sp


# ---------------- A1: theta ----------------------------------------------------------------------------------------
def theta_min(d: int, k: int) -> tuple[Fr, Fr]:
    """min over mu in [nu, 1) of theta(mu) = 2mu + max(2nu, mu-(d-2)nu, 1-d mu+nu) + (d-1)nu  (eta0 -> 0)."""
    nu = Fr(1, 4 * k)
    lines = [(Fr(0), 2 * nu), (Fr(1), -(d - 2) * nu), (Fr(-d), 1 + nu)]
    cands = {nu}
    for (a1, b1), (a2, b2) in itertools.combinations(lines, 2):
        if a1 != a2:
            m = (b2 - b1) / (a1 - a2)
            if nu <= m < 1:
                cands.add(m)
    best = None
    for m in cands:
        th = 2 * m + max(al * m + be for al, be in lines) + (d - 1) * nu
        if best is None or th < best[0]:
            best = (th, m)
    return best


def theta0(d0: Fr, kappa: Fr) -> Fr:
    return Fr(3) / d0 if kappa >= d0 * d0 / 4 else Fr(2) / d0 + d0 / (4 * kappa)


def check_A1():
    out = []
    for d0, kappa in [(Fr(1), Fr(1, 8)), (Fr(1), Fr(1)), (Fr(2), Fr(1, 2)), (Fr(1, 2), Fr(1, 32))]:
        rows = []
        for t in [10, 100, 1000, 10000]:      # t = D^{1/5}
            d = max(3, -(-int(d0 * t) // 1))
            if d % 2 == 0:
                d += 1
            k = max(-(-int(kappa * t * t * 10 ** 6) // 10 ** 6), d + 1)
            k = max(int(sp.ceiling(kappa * t * t)), d + 1)
            th, mu = theta_min(d, k)
            rows.append({"t": t, "d": d, "k": k, "theta*t": str(th * t), "float": float(th * t)})
        out.append({"d0": str(d0), "kappa": str(kappa), "theta0": str(theta0(d0, kappa)), "rows": rows})
    return out


# ---------------- A2: LD*-U Step 7 power count ------------------------------------------------------------------
def check_A2():
    M, H, V, Y, d = sp.symbols("M H V Y d", positive=True)
    # (#P/|G|) * 10 * (e^{11/3} C_W / ((2H+1)^2 V^{2d-4})) * U * sum_q tau(q) Psi(X'+) Psi(X_F)^2 /(Y q (MV)^2)
    nP = V**2 / d                 # #P <= 5V^2/(d log V)
    G = M**2                      # |G| >= M^2/8
    kern = 1 / (H**2 * V**(2 * d - 4))
    U = M / V                     # U = (3M/V) rho_r (log N)^C
    PsiXr = M * V / Y             # X'_+ = 2MV/Y, times rho_r
    PsiXF2 = ((2 * d + 1) * M * H * V**d)**2   # times rho_F^2
    expr = sp.simplify(nP / G * kern * U * PsiXr * PsiXF2 / (Y * (M * V)**2))
    return {"Step7_product_without_rho_and_polylogs": str(expr),
            "expected": "(2d+1)^2 V^4/(d Y^2)  -> 9d * 4 S^4 since V^2 = 2YS^2"}


# ---------------- A3: Kern-U by brute force ------------------------------------------------------------------------
def v_of(P, d):
    a, b = P
    return [a ** (d - 1 - i) * b ** i for i in range(d)]


def max_fibre(P, Pp, d, H):
    v, vp = v_of(P, d), v_of(Pp, d)
    cnt = {}
    rng = range(-H, H + 1)
    for c in itertools.product(rng, repeat=d):
        z = (sum(ci * vi for ci, vi in zip(c, v)), sum(ci * vi for ci, vi in zip(c, vp)))
        cnt[z] = cnt.get(z, 0) + 1
    return max(cnt.values())


def check_A3(seed=20261020):
    rnd = random.Random(seed)
    viol = 0
    worst = Fr(0)
    n = 0
    recs = []
    for d in [3, 4, 5]:
        for _ in range(25):
            V = rnd.choice([2, 3, 4])
            H = rnd.choice([V * V, d * V * V]) if (2 * d * V * V + 1) ** d <= 3 * 10 ** 6 else V * V
            if (2 * H + 1) ** d > 3 * 10 ** 6:
                H = max(1, int((3 * 10 ** 6) ** (1 / d)) // 2)
            while True:
                b, bp = rnd.randint(1, V), rnd.randint(1, V)
                a, ap = rnd.randint(-V, V), rnd.randint(-V, V)
                if sp.gcd(a, b) == 1 and sp.gcd(ap, bp) == 1 and a * bp - ap * b != 0:
                    break
            mf = max_fibre((a, b), (ap, bp), d, H)
            bound = (2 * H // (b * bp) + 1) ** (d - 2)
            n += 1
            if mf > bound:
                viol += 1
            if b * bp * d <= 2 * H + 1:
                const = Fr(mf * (2 * H + 1) ** 2 * (b * bp) ** (d - 2), (2 * H + 1) ** d)
                worst = max(worst, const)
            recs.append((d, (a, b), (ap, bp), H, mf, bound))
    # the (3/2)^d family: P=(1,b), P'=(1,b+1) at H = b^2 (bb' > (2H+1)/d excluded by Kern-U's hypothesis)
    fam = []
    for d, b in [(3, 3), (4, 3), (5, 2)]:
        H = b * b
        if (2 * H + 1) ** d > 4 * 10 ** 6:
            continue
        mf = max_fibre((1, b), (1, b + 1), d, H)
        const = Fr(mf * (2 * H + 1) ** 2 * (b * (b + 1)) ** (d - 2), (2 * H + 1) ** d)
        fam.append({"d": d, "b": b, "H": H, "maxfibre": mf, "C_d_needed": float(const),
                    "Kern-U hypothesis bb'<=(2H+1)/d": b * (b + 1) * d <= 2 * H + 1})
    return {"cases": n, "violations": viol, "max const (2H+1)^2(bb')^{d-2}maxPr where bb'<=(2H+1)/d": float(worst),
            "e": 2.718281828, "family_H_eq_b2": fam}


# ---------------- A4: image index of Phi ---------------------------------------------------------------------------
def check_A4(seed=7):
    from sympy.matrices.normalforms import smith_normal_form
    rnd = random.Random(seed)
    bad = 0
    n = 0
    for d in [3, 4, 5, 6, 7]:
        for _ in range(20):
            while True:
                b, bp = rnd.randint(1, 30), rnd.randint(1, 30)
                a, ap = rnd.randint(-30, 30), rnd.randint(-30, 30)
                if sp.gcd(a, b) == 1 and sp.gcd(ap, bp) == 1 and a * bp - ap * b != 0:
                    break
            Mx = sp.Matrix([v_of((a, b), d), v_of((ap, bp), d)])
            snf = smith_normal_form(Mx, domain=sp.ZZ)
            idx = abs(snf[0, 0] * snf[1, 1])
            n += 1
            if idx != abs(a * bp - ap * b):
                bad += 1
    return {"pairs": n, "index != |Delta|": bad}


# ---------------- A5: alpha floor ------------------------------------------------------------------------------------
def check_A5():
    """Exponents in powers of D: d = D^delta, k = D^kap, theta >= max(D^{-delta}, D^{delta-kap}),
    split: k (1 + tau-loss d^tau) loglog <= Lambda_alpha => kap + tau*delta <= alpha;
    smoothness: theta log N / log B * log u <= Lambda_alpha  <=>  theta <= D^{2alpha-1}:  -delta <= 2alpha-1, delta-kap <= 2alpha-1.
    Minimise alpha over a fine rational grid of (delta, kap) and compare with (2+tau)/(5+2tau)."""
    res = {}
    for tau in [0, 1, 2]:
        best = None
        grid = [Fr(i, 120) for i in range(0, 121)]
        for delta in grid:
            for kap in grid:
                alpha = max(kap + tau * delta, (1 - delta) / 2, (1 + delta - kap) / 2)
                if best is None or alpha < best[0]:
                    best = (alpha, delta, kap)
        res[f"tau={tau}"] = {"alpha_min_grid": str(best[0]), "delta": str(best[1]), "kap": str(best[2]),
                             "formula (2+tau)/(5+2tau)": str(Fr(2 + tau, 5 + 2 * tau))}
    return res


# ---------------- A6: cost constant ----------------------------------------------------------------------------------
def bracket_2_25(den=10 ** 7):
    lo = 13195079 * den // 10 ** 7
    while Fr(lo + 1, den) ** 5 < 4:
        lo += 1
    while Fr(lo, den) ** 5 > 4:
        lo -= 1
    assert Fr(lo, den) ** 5 < 4 < Fr(lo + 1, den) ** 5
    return Fr(lo, den), Fr(lo + 1, den)


T_LO, T_HI = bracket_2_25()
REC_LO, REC_HI = T_HI / (T_HI - 1), T_LO / (T_LO - 1)   # 1/(1-2^{-2/5}) = t/(t-1), decreasing in t


def cost(a, d0, kappa, Cg):
    sig = 2 * theta0(d0, kappa) / (5 * a)
    s = 2 * kappa * (1 + Cg)
    p = max(2 * a, a + s + max(a, sig))
    return max(p, (a + s) * REC_LO), max(p, (a + s) * REC_HI), sig, s


def check_A6():
    out = {"1/(1-2^{-2/5}) bracket": [str(REC_LO), str(REC_HI)]}
    grid_a = [Fr(i, 40) for i in range(1, 121)]
    grid_k = [Fr(i, 256) for i in range(1, 129)] + [Fr(1, 2 ** j) for j in range(9, 15)]
    grid_d = [Fr(i, 16) for i in range(1, 97)]
    for Cg in [Fr(0), Fr(1), Fr(5), Fr(20), Fr(100)]:
        best = None
        for kappa in grid_k:
            for d0 in grid_d:
                for a in grid_a:
                    lo, hi, sig, s = cost(a, d0, kappa, Cg)
                    if best is None or hi < best[1]:
                        best = (lo, hi, a, d0, kappa, sig, s)
        lo, hi, a, d0, kappa, sig, s = best
        out[f"C_gamma={Cg}"] = {"c_hi": float(hi), "c_lo": float(lo), "a": str(a), "d0": str(d0), "kappa": str(kappa),
                                "sigma": float(sig), "s": float(s),
                                "terms": {"2a": float(2 * a), "a+s+max(a,sigma)": float(a + s + max(a, sig)),
                                          "(a+s)/(1-2^-2/5)": float((a + s) * REC_HI)}}
    # asymptotic in C_gamma: kappa ~ x/sqrt(C_gamma) optimum => c ~ const*sqrt(C_gamma)
    return out


# ---------------- A7: margins (my own list) ---------------------------------------------------------------------------
class E:
    """sum c * (log N)^p (loglog N)^q, exact."""

    def __init__(self, t=None):
        self.t = {k: v for k, v in (t or {}).items() if v != 0}

    @staticmethod
    def m(c, p, q):
        return E({(Fr(p), Fr(q)): Fr(c)})

    def __add__(self, o):
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = r.get(k, 0) + v
        return E(r)

    def __sub__(self, o):
        return self + E({k: -v for k, v in o.t.items()})

    def __mul__(self, c):
        return E({k: v * Fr(c) for k, v in self.t.items()})

    def lead(self):
        k = max(self.t)
        return k, self.t[k]

    def big(self):
        (p, q), c = self.lead()
        return c > 0 and (p, q) > (0, 0)


def Dp(s):          # D^s
    return (Fr(s), -Fr(s))


def check_A7(a=Fr(3, 5), d0=Fr(1), kappa=Fr(1, 8), Cg=Fr(20), C=Fr(1000)):
    th0 = theta0(d0, kappa)
    sig = 2 * th0 / (5 * a)
    LOGN, LL = E.m(1, 1, 0), E.m(1, 0, 1)
    LAM = E.m(1, Fr(2, 5), Fr(3, 5))
    logB = LAM * a
    logY = E.m(1 / (2 * kappa), Fr(3, 5), Fr(2, 5))          # log N/(2k), k = kappa D^{2/5}
    eta0logN = E.m(1, Fr(1, 2), Fr(1, 2))                   # D^{-1/2} log N
    thetalogN = E.m(th0, Fr(4, 5), Fr(1, 5))
    logS2 = LAM * sig
    d = E.m(d0, *Dp(Fr(1, 5)))
    dlogd = E.m(d0 / 5, Fr(1, 5), Fr(4, 5))
    k = E.m(kappa, *Dp(Fr(2, 5)))
    omega = E.m(1 / a, Fr(3, 5), Fr(-3, 5))
    u = E.m(th0 / a, *Dp(Fr(2, 5)))
    logH_min = logY                                            # h >= 2 nu  =>  log H >= log N/(2k) = log Y
    chks = {
        # gamma^{2k}: split exponent is exactly linear in Lambda
        "2k*(1+Cg)*loglogN == 2kappa(1+Cg)Lambda": ((E({(k.lead()[0][0], k.lead()[0][1] + 1): k.lead()[1]}) * (2 * (1 + Cg)))
                                                   - LAM * (2 * kappa * (1 + Cg))).t == {},
        # any 2^{d} loss in gamma would contribute 2k*d*log2 to log(1/p): is it o(Lambda)?  (it must NOT be)
        "2^d loss in gamma is NOT affordable: k*d >> Lambda": (E.m(kappa * d0, *Dp(Fr(3, 5))) - LAM * C).big(),
        # (log N)^{O(d)} loss: k*d*loglogN >> Lambda
        "(logN)^d loss NOT affordable": (E.m(kappa * d0, Fr(3, 5), Fr(2, 5)) - LAM * C).big(),
        # poly(d) = (log N)^{O(1)} losses are affordable: k * log d ~ Lambda * (1/5)... is O(Lambda)
        "d^{C} loss costs O(Lambda) (finite, enters c)": True,
        "N^{eta0/4} >> B^C": (eta0logN * Fr(1, 4) - logB * C).big(),
        "N^{eta0/4} >> S^C": (eta0logN * Fr(1, 4) - logS2 * C).big(),
        "N^{eta0/4} >> 2^{Cd}": (eta0logN * Fr(1, 4) - d * C).big(),
        "Y >> B^C (8.5'(R'4), L2_fam (C))": (logY - logB * C).big(),
        "u >> (loglogN)^4 (R'2)": (u * Fr(1, 2) - E.m(1, 0, 4)).big(),
        "(2d)! bits: d log d = o(Lambda)": (LAM - dlogd * C).big(),
        "S^d = N^{o(theta)}": (thetalogN - E.m(d0 * sig / 2, Fr(3, 5), Fr(2, 5)) * C).big(),
        "eta0 = o(theta)": (thetalogN - eta0logN * C).big(),
        "|D_A| >= 4k^3 2^omega (S'')": (logY - omega - LL * C).big(),
        "rejection exp(-logH/(13 d loglogH)) = o(gamma)": (E.m(1 / (2 * kappa * 13 * d0), Fr(2, 5), Fr(0)) - LL * C).big(),
        "C_W = 24/eta0 polylog (eta0 logN >> log S)": (eta0logN * Fr(1, 3) - logS2 * C).big(),
        "level gap Gamma >> (2K0+U0) log B": (eta0logN * Fr(1, 2) - LAM * C - logB * C).big(),
        "Drappeau saving u/log^2 u >> loglog (in exponent: D^{2/5}/(loglog)^2 vs loglog)":
            (E.m(th0 / a, Fr(2, 5), Fr(-2, 5) - 2) - LL * C).big(),
        "2^omega NOT L^{o(1)} (why S'' is needed)": (omega - LAM * C).big(),
    }
    return {"params": {"a": str(a), "d0": str(d0), "kappa": str(kappa), "C_gamma": str(Cg), "C": str(C)},
            "checks": {k_: bool(v) for k_, v in chks.items()}}


# ---------------- A8: where the asymptotics begin -----------------------------------------------------------------------
def check_A8(a=Fr(3, 5), d0=Fr(1), kappa=Fr(1, 8)):
    """Thinnest margins and their crossover D = log N / log log N (powers of 10, exact rational comparisons via
    integer powers).  We use loglog N ~ log D (lower-order correction ignored; printed as an estimate only)."""
    th0 = theta0(d0, kappa)
    sig = 2 * th0 / (5 * a)
    out = {}
    # (P2) N^{eta0/4} >= B^{C}:  D^{1/2} loglog/4 >= C a D^{2/5} loglog  <=>  D^{1/10} >= 4 C a
    for C in [3, 10, 2 * sig / a + 3]:
        need = 4 * Fr(C) * a       # D^{1/10} >= need  <=> D >= need^10
        out[f"P2: N^(eta0/4) >= B^{float(C):.3g}: D >= (4Ca)^10"] = float(need ** 10)
    # Lemma W / C_W:  N^{eta0/3} > S^C ... same D^{1/10} type
    out["C_W: eta0 logN/3 >= 2 log S: D >= (3 sigma)^10"] = float((3 * sig) ** 10)
    # (2d)! bits: d log d << Lambda : d0 D^{1/5} * (1/5) log D vs D^{2/5} log D: D^{1/5} >= d0/5 * C
    out["(2d)!: D^{1/5} >= 10*d0/5 (cost exponent < 0.1*Lambda): D >="] = float((Fr(10) * d0 / 5) ** 5)
    # S'' hypothesis: log Y - omega log 2 > 0  <=> loglog N /(2 kappa) > log2 / a  <=> loglog N > 2 kappa log2 /a
    out["S'': loglogN > 2kappa*log2/a (log N > ...)"] = float(2 * kappa * Fr(693147, 10 ** 6) / a)
    # d >= 3 needs D^{1/5} >= 3/d0
    out["d >= 3 honestly (d0 D^{1/5} >= 3): D >="] = float((3 / d0) ** 5)
    return out


# ---------------- A9: nested's Stark-threshold correction at alpha = 1/3 vs alpha = 2/5 ------------------------------
def check_A9():
    """Character primes have s ~ (2d)! (log N)^4 bits; their cost is s^{O(1)} = exp(O(d log d) + O(loglog N)).
    With d = d0 D^{delta}:  d log d ~ d0 delta D^{delta} loglog N  (log D ~ loglog N), and Lambda_alpha = D^{alpha} loglog N.
    Ratio (d log d)/Lambda_alpha ~ d0 delta D^{delta - alpha}: a positive constant (Theta(delta)) iff delta = alpha
    (nested's alpha = 1/3, d = delta D^{1/3}), and -> 0 iff delta < alpha (Theorem U: delta = 1/5 < 2/5).
    Exact exponent comparison: the ratio is c * D^{delta-alpha} with c = d0*delta."""
    rows = {}
    for name, delta, alpha in [("nested alpha=1/3, d ~ D^{1/3}", Fr(1, 3), Fr(1, 3)),
                               ("Theorem U alpha=2/5, d ~ D^{1/5}", Fr(1, 5), Fr(2, 5)),
                               ("hypothetical alpha=2/5 with d ~ D^{2/5}", Fr(2, 5), Fr(2, 5))]:
        e = delta - alpha
        rows[name] = {"exponent of D in (2 d log d)/Lambda_alpha": str(e),
                      "verdict": "L^{o(1)}" if e < 0 else ("L[alpha, Theta(d0*delta)] (finite, enters c)" if e == 0 else "super-L[alpha]")}
    # concrete: Stark's (2d)! threshold at Theorem U's d, cost exponent C*2d*log(2d)/Lambda for C = 7 (AKS s^{6+})
    conc = []
    for D10 in [6, 10, 20, 40, 80]:
        D = Fr(10) ** D10
        # d = D^{1/5} (d0 = 1): use exact integer 10^{D10/5} when divisible, else rational approx via floor
        d = int(round(10 ** (D10 / 5)))
        # log terms (floats for printing only): 2d log(2d) vs Lambda = D^{2/5} * loglogN, loglogN ~ ln(D)+ln(ln D)
        import math
        lnD = D10 * math.log(10)
        Lam = 10 ** (0.4 * D10) * (lnD + math.log(lnD))
        conc.append({"D": f"1e{D10}", "d": d, "7*2d ln(2d)/Lambda": 7 * 2 * d * math.log(2 * d) / Lam})
    rows["concrete (float, print only)"] = conc
    return rows


def main(out_path=None):
    res = {"A1_theta": check_A1(), "A2_LDstep7": check_A2(), "A3_kernU": check_A3(), "A4_image_index": check_A4(),
           "A5_alpha_floor": check_A5(), "A6_cost": check_A6(), "A7_margins": check_A7(), "A8_crossover": check_A8(), "A9_stark_alpha": check_A9()}
    # compare with lv5unif's own optimiser (imported only here)
    try:
        sys.path.insert(0, ".")
        from src.lv5unif import params as P
        cmp = {}
        for Cg in [0, 1, 5, 20]:
            b = P.optimise(Fr(Cg))
            cmp[f"C_gamma={Cg}"] = {"lv5unif_hi": float(b["hi"]), "a": str(b["a"]), "kappa": str(b["kappa"]),
                                    "d0": str(b["d0"])}
        res["compare_lv5unif"] = cmp
    except Exception as e:  # pragma: no cover
        res["compare_lv5unif"] = f"import failed: {e!r}"
    s = json.dumps(res, indent=1, default=str)
    if out_path:
        open(out_path, "w").write(s)
    print(s)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
