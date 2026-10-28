"""critic2: independent leading-order audit of the 2L_DW(alpha) design sizes, the window
(3/8, 2/5), the endpoint constraint a > 2 c2 and the cost constant c(alpha).

Written from scratch (exact rationals / sympy); imports nothing from the project's src/.

Conventions: lambda = log log N, D = log N / lambda, so log N = D*lambda and log D = lambda(1+o(1)).
A log-size is a list of terms (coef, p, q) meaning coef * D^p * lambda^q.  Leading term = the
lexicographically largest (p, q).  u log u for u = c D^e (e > 0) is c*e*D^e*lambda (1+o(1)).

Run:  .venv/bin/python -m src.audit2L.critic2_params
"""
from fractions import Fraction as Fr
import itertools
import sympy as sp

a, kap, d0, b2 = sp.symbols("a kappa d0 b2", positive=True)


# ----------------------------------------------------------------------------------------------
# Part A: sizes from the design definitions (c12alg section 1.1)
# ----------------------------------------------------------------------------------------------
def T(c, p, q=1):
    return [(sp.nsimplify(c), Fr(p), Fr(q))]


def add(*xs):
    out = []
    for x in xs:
        out += x
    return out


def scale(c, x):
    return [(sp.simplify(c * t[0]), t[1], t[2]) for t in x]


def mul_D(dp, c, x):
    """multiply by c * D^dp (no lambda)"""
    return [(sp.simplify(c * t[0]), t[1] + dp, t[2]) for t in x]


def lead(x):
    """leading term, summing coefficients of equal (p,q)"""
    acc = {}
    for c, p, q in x:
        acc[(p, q)] = sp.simplify(acc.get((p, q), 0) + c)
    keys = sorted([k for k, v in acc.items() if v != 0], reverse=True)
    k = keys[0]
    return acc[k], k[0], k[1]


def second(x):
    acc = {}
    for c, p, q in x:
        acc[(p, q)] = sp.simplify(acc.get((p, q), 0) + c)
    keys = sorted([k for k, v in acc.items() if v != 0], reverse=True)
    return keys[1] if len(keys) > 1 else None


def sizes(al, sig1=sp.Symbol("s1"), sig2=sp.Symbol("s2")):
    al = Fr(al)
    logN = T(1, 1)
    d_c, d_p = d0, 1 - 2 * al                       # d ~ d0 D^{1-2a}
    k_c, k_p = kap, al                              # k ~ kappa D^a
    nu_logN = T(1 / (4 * kap), 1 - al)              # nu log N = log N/(4k)
    logM = T(1 / d0, 2 * al)                        # mu log N = (1-nu) log N/d
    eta_logN = T(1, 3 * al - Fr(1, 2))              # eta0 = 1/ceil(D^{3/2-3a})
    logY = T(1 / (2 * kap), 1 - al)                 # Y = N^{1/(2k)}
    logS1 = T(sig1 / 2, al)
    logV1 = add(scale(sp.Rational(1, 2), logY), logS1)          # V1 = sqrt(2Y) S1
    logWm = add(scale(2, logM), scale(2, logV1), eta_logN)       # Wm = (2 M V1)^2 N^eta0
    logH = add(logWm, eta_logN, T(d0 * sp.log(2), 1 - 2 * al, 0))  # H = 2^{d+8} d^4 Wm N^eta0
    dlogV1 = mul_D(d_p, d_c, logV1)
    logXF = add(logM, logH, dlogV1)                 # XF = (2d+1) M H V1^d
    logXr = add(logM, logV1, scale(-1, logY))       # Xr = 2 M V1 / Y
    logB = T(a, al)
    logB2 = T(b2, 4 * al - 1)
    # child (conservative box, worst skew): |P'| <= V2 ~ 2 S2 q, q <= B2
    logS2 = T(sig2 / 2, al)
    logV2 = add(logS2, logB2)
    logXr2 = add(logM, logV2)
    logXF2 = add(logM, logH, mul_D(d_p, d_c, logV2), scale(-1, logB2))
    return dict(logN=logN, nu_logN=nu_logN, logM=logM, eta_logN=eta_logN, logY=logY,
                logV1=logV1, logWm=logWm, logH=logH, logXF=logXF, logXr=logXr, logB=logB,
                logB2=logB2, logXr2=logXr2, logXF2=logXF2, d=(d_c, d_p), k=(k_c, k_p))


def ulogu(logX, logy):
    cX, pX, qX = lead(logX)
    cy, py, qy = lead(logy)
    assert qX == qy == 1
    e = pX - py
    assert e > 0
    return sp.simplify(cX / cy * e), e     # coefficient of D^e * lambda


def sigmas(al):
    s = sizes(al)
    al = Fr(al)
    c_r, e_r = ulogu(s["logXr"], s["logB"])
    c_F2, e_F2 = ulogu(s["logXF"], s["logB2"])
    c_r2, e_r2 = ulogu(s["logXr2"], s["logB"])
    c_f2, e_f2 = ulogu(s["logXF2"], s["logB"])
    assert e_r == e_r2 == e_f2 == al
    sigma1 = c_r + (c_F2 if e_F2 == al else 0)
    sigma2 = sp.simplify(c_r2 + c_f2)
    return dict(sigma1=sp.simplify(sigma1), sigma2=sigma2, u2logu2=(c_F2, e_F2),
                leadXF=lead(s["logXF"]), secondXF=second(s["logXF"]), leadXr=lead(s["logXr"]),
                leadXF2=lead(s["logXF2"]), leadXr2=lead(s["logXr2"]), leadH=lead(s["logH"]))


# ----------------------------------------------------------------------------------------------
# Part B: every leading-order exponent comparison found in the excerpt, as f(alpha) > 0
# ----------------------------------------------------------------------------------------------
x = sp.Symbol("alpha")
theta = sp.Max(3 * x - 1, 1 - 2 * x, 2 * x - sp.Rational(1, 2))
theta_p = sp.Max(1 - 2 * x, 2 - 5 * x)
zeta = (theta + x) / 2
zeta_p = (theta_p + 4 * x - 1) / 2

CHECKS = [
    # (name, f, where used)
    ("theta < alpha (child loss vs child BDH saving)", x - theta, "AB2 margins, RAT"),
    ("theta_p < 3-7a (parent loss vs B2 saving)", (3 - 7 * x) - theta_p, "AB2(ii), O3"),
    ("theta < 1-... child gap 1/2-alpha", sp.Rational(1, 2) - x, "= alpha - (2alpha-1/2)"),
    ("theta_p gap 2-5a", 2 - 5 * x, "= (3-7a)-(1-2a)"),
    ("h_* >= B^{-5/12}: theta < alpha", x - theta, "BL, OPV-A"),
    ("h_*^p >= B2^{-5/12}: theta_p < 4a-1", (4 * x - 1) - theta_p, "par. blocks"),
    ("zeta < alpha", x - zeta, "EXC-A"), ("theta < zeta", zeta - theta, "EXC-A"),
    ("zeta_p < 4a-1", (4 * x - 1) - zeta_p, "EXC-A"), ("theta_p < zeta_p", zeta_p - theta_p, "EXC-A"),
    ("zeta > a - 7/20", zeta - (x - sp.Rational(7, 20)), "EXC-A"),
    ("zeta_p > 1-a-7/20", zeta_p - (1 - x - sp.Rational(7, 20)), "EXC-A"),
    ("omega(parent Delta') = o(P0=D^{7/10}): 1-a < 7/10", sp.Rational(7, 10) - (1 - x), "FO"),
    ("omega(child modulus 4Q2 S2^2) = o(P0): 4a-1 < 7/10", sp.Rational(7, 10) - (4 * x - 1), "FO"),
    ("Rankin parent tail: (1-a)-0.343 < 4a-1", (4 * x - 1) - (1 - x - sp.Rational(343, 1000)), "par FO"),
    ("child modulus x^o(1): 4a-1 < 2a", 2 * x - (4 * x - 1), "Fact D range"),
    ("parent alg modulus x^o(1): 1-a < 2-3a", (2 - 3 * x) - (1 - x), "Fact D range"),
    ("parent rat modulus x^o(1): 1-a < 2a", 2 * x - (1 - x), "Fact D range"),
    ("N^eta0 >> L[a]: 3a-1/2 > a", (3 * x - sp.Rational(1, 2)) - x, "TV vs densities"),
    ("N^eta0 >> B2: 3a-1/2 > 4a-1", (3 * x - sp.Rational(1, 2)) - (4 * x - 1), "TV vs q"),
    ("eta0 = o(mu): 3a-1/2 < 2a", 2 * x - (3 * x - sp.Rational(1, 2)), "design"),
    ("d log2 << eta0 log N: 1-2a < 3a-1/2", (3 * x - sp.Rational(1, 2)) - (1 - 2 * x), "h ~ 2mu+2nu+2eta"),
    ("filter spread (N^eta0)^gamma: 2a-1/2 < a", x - (2 * x - sp.Rational(1, 2)), "theta"),
    ("child |P'|^gamma: 3a-1 < a", x - (3 * x - 1), "theta"),
    ("PB with AB (first order): A/eta, -a+(3-7a) < 0", 8 * x - 3, "Prop PB window"),
    ("PB with AB: sqrt(AB)/eta: (3-8a)/2 < 0", (8 * x - 3) / 2, "Prop PB window"),
    ("PB with AB2: A/eta: 3-9a < 0", 9 * x - 3, "Prop PB 3/8"),
    ("PB with AB2: sqrt(AB)/eta: max(3-9a,-a) < 0", -sp.Max(3 - 9 * x, -x), "Prop PB 3/8"),
    ("AB2 -> AB: D^{-2a} <= D^{7a-3}", (7 * x - 3) - (-2 * x), "Cor of AB2"),
    ("AB2 -> AB: D^{-2(3-7a)} <= D^{7a-3}", (7 * x - 3) + 2 * (3 - 7 * x), "Cor of AB2"),
    ("(SQ)/subtraction/sigma1: u2 log u2 (3-7a) vs log B (a)", x - (3 - 7 * x), "FM-U', eta, sec11"),
    ("B2 << V1: 4a-1 < 1-a", (1 - x) - (4 * x - 1), "DOM(i), parents"),
    ("one-point parent term B2/#J: 4a-1 < 1-a", (1 - x) - (4 * x - 1), "VD/z"),
    ("u2 log u2 << log B2: 3-7a < 4a-1 (a>4/11)", (4 * x - 1) - (3 - 7 * x), "Lemma eta, Hildebrand"),
    ("recursion L_{XF}[2/5]: (2/5)(2-3a) < a (a>4/11)", x - sp.Rational(2, 5) * (2 - 3 * x), "uprime 6"),
    ("LV6 e^{O(d log d)}: 1-2a < a", x - (1 - 2 * x), "LV6-U"),
    ("2^d <= B: 1-2a < a", x - (1 - 2 * x), "CS overshoot"),
    ("d nu -> 0: 1-2a < a", x - (1 - 2 * x), "POSTPROCESS |D|"),
    ("d eta0 -> 0: 1-2a < 3/2-3a", (sp.Rational(3, 2) - 3 * x) - (1 - 2 * x), "POSTPROCESS |D|"),
    ("d log V1 dominates log MH: 2-3a > 2a", (2 - 3 * x) - 2 * x, "size of XF"),
    ("mu_opt branch k < d^2/4: a < 2-4a", (2 - 4 * x) - x, "mu=(1-nu)/d"),
    ("mu-(d-2)nu < 0 (h_* = 2nu): 2a-1 < 1-3a", (1 - 3 * x) - (2 * x - 1), "height floor"),
    ("N/M^d <= M: 1-a < 2a", 2 * x - (1 - x), "GB-DW"),
    ("Xr leading 2a > 1-a", 2 * x - (1 - x), "size of Xr"),
    ("children above B: a < 4a-1", (4 * x - 1) - x, "eta v0->0, DS-U' B^C<=B2"),
    ("DOM(iii) u_B >> u2: 2-4a > 3-7a", (2 - 4 * x) - (3 - 7 * x), "W >= mu1"),
    ("V2 << M: 4a-1 < 2a", 2 * x - (4 * x - 1), "CLB fibres"),
    ("B2^3 << H: 4a-1 < 2a", 2 * x - (4 * x - 1), "SQF CRT"),
    ("parent cell length >> Delta': 1-a < 2a", 2 * x - (1 - x), "BL(i) parents"),
    ("rejections e^{-D^{4a-1}}, e^{-D^{1-2a}} << polylog", sp.Min(4 * x - 1, 1 - 2 * x), "reject-DW"),
    ("u2 -> infinity: 3-7a > 0", 3 - 7 * x, "U1_2, DS-U'"),
    ("U1_2 shifts (1-a2)log w^2 -> 0: 1-3a < 0", 3 * x - 1, "U1_2"),
    ("CS' d log S = o(log XF): 1-a < 2-3a", (2 - 3 * x) - (1 - x), "CS'"),
    ("CS' XF-shift in u2: 2-5a < 3-7a", (3 - 7 * x) - (2 - 5 * x), "CS'"),
    ("eta small part: (1-2a)-3+(2-4a) < 0", 6 * x, "Lemma eta"),
    ("skew spread (S2 T0)^{d gamma}: 1/2-a < theta", (2 * x - sp.Rational(1, 2)) - (sp.Rational(1, 2) - x), "SPR"),
    ("beta-drift child: zeta < 2a", 2 * x - zeta, "OPV-A"),
    ("beta-drift parent: zeta_p < 2-3a", (2 - 3 * x) - zeta_p, "OPV-A par"),
    ("D^{2-8a} <= D^{-2a}", 6 * x - 2, "OPV-A par"),
    ("D^{14a-6} << D^{7a-3}", 3 - 7 * x, "RAT2"),
    ("range restatement needed iff 3-7a < 1/4", sp.Rational(1, 4) - (3 - 7 * x), "c12alg 2.7 (info)"),
]


def threshold(f):
    """return set of alpha in [1/3, 1/2] where f changes sign (piecewise linear)"""
    sols = set()
    for piece in sp.piecewise_fold(f).atoms(sp.Max, sp.Min) or []:
        pass
    # brute: solve on each linear branch by sampling breakpoints
    cands = set()
    for num in range(1000, 1500 + 1):
        pass
    try:
        for s in sp.solve(sp.Eq(f, 0), x):
            if s.is_real:
                cands.add(sp.Rational(s))
    except Exception:
        pass
    return sorted(cands)


def evalf(f, al):
    return sp.nsimplify(f.subs(x, sp.Rational(al)))


def run_checks():
    pts = [Fr(1, 3), Fr(4, 11), Fr(3, 8), Fr(7, 18), Fr(11, 28), Fr(79, 200), Fr(2, 5)]
    grid = [Fr(3, 8) + Fr(i, 4000) for i in range(0, 101)]   # [3/8, 2/5]
    print("\n== Part B: exponent comparisons f(alpha) > 0 (exact) ==")
    print(f"{'check':58s} {'3/8':>7s} {'7/18':>7s} {'2/5':>7s}  open-window  thresholds")
    fails = []
    for name, f, src in CHECKS:
        v38, v718, v25 = (evalf(f, p) for p in (Fr(3, 8), Fr(7, 18), Fr(2, 5)))
        inside = all(evalf(f, g) > 0 for g in grid[1:-1])
        th = []
        for g1, g2 in zip(grid, grid[1:]):
            s1, s2 = evalf(f, g1), evalf(f, g2)
            if s1 == 0:
                th.append(str(g1))
        # solve exactly for zero crossings in [1/3,1/2]
        try:
            zs = [s for s in sp.solve(f, x) if s.is_real and sp.Rational(1, 3) <= s <= sp.Rational(1, 2)]
        except Exception:
            zs = []
        tag = lambda v: ("ok" if v > 0 else ("TIE" if v == 0 else "FAIL"))
        print(f"{name[:58]:58s} {tag(v38):>7s} {tag(v718):>7s} {tag(v25):>7s}  {str(inside):>5s}  {zs}")
        if not inside and "info" not in src:
            fails.append(name)
    print("open-window failures (should be none):", fails)


# ----------------------------------------------------------------------------------------------
# Part C: cost constant
# ----------------------------------------------------------------------------------------------
def two_pow_neg_bracket(al, iters=60):
    """exact rational bracket [lo, hi] for 2^{-al}, al = p/q, via x^q vs 2^{-p}"""
    al = Fr(al)
    p, q = al.numerator, al.denominator
    lo, hi = Fr(0), Fr(1)
    target = Fr(1, 2 ** p)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if mid ** q <= target:
            lo = mid
        else:
            hi = mid
    return lo, hi


def cost(al, A, K, D0, B2, Cg=0, endpoint=None):
    al, A, K, D0, B2 = map(Fr, (al, A, K, D0, B2))
    c2 = Fr(3) * D0 / (32 * K * B2)
    sig1 = al / (A * D0) + (c2 if al == Fr(3, 8) else 0)
    sig2 = (al / A) * (Fr(4) / D0 + D0 * B2)
    s = 2 * K * (1 + Fr(Cg))
    lo, hi = two_pow_neg_bracket(al)
    rec_lo, rec_hi = (A + s) / (1 - lo), (A + s) / (1 - hi)
    t1, t2 = 2 * A, A + s + max(A, sig1 + sig2)
    return dict(c_lo=max(t1, t2, rec_lo), c_hi=max(t1, t2, rec_hi), terms=(t1, t2, (rec_lo, rec_hi)),
                sig1=sig1, sig2=sig2, c2=c2)


def grid38(constraint_factor):
    best = None
    As = [Fr(i, 20) for i in range(4, 41)]
    Ks = [Fr(1, 32), Fr(1, 16), Fr(1, 12), Fr(1, 8), Fr(3, 16), Fr(1, 4), Fr(3, 8), Fr(1, 2)]
    D0s = [Fr(1, 2), Fr(2, 3), Fr(3, 4), Fr(1), Fr(5, 4), Fr(3, 2), Fr(2), Fr(5, 2), Fr(3)]
    B2s = [Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1), Fr(3, 2), Fr(2), Fr(3), Fr(4), Fr(6), Fr(8)]
    for A, K, D0, B2 in itertools.product(As, Ks, D0s, B2s):
        c2 = Fr(3) * D0 / (32 * K * B2)
        if not A > constraint_factor * c2:
            continue
        r = cost(Fr(3, 8), A, K, D0, B2)
        if best is None or r["c_hi"] < best[0]:
            best = (r["c_hi"], (A, K, D0, B2))
    return best


def main():
    print("== Part A: sizes re-derived from the design definitions ==")
    for al in (Fr(3, 8), Fr(7, 18), Fr(79, 200)):
        sg = sigmas(al)
        print(f"alpha={al}: lead logXF={sg['leadXF']} (next exponent {sg['secondXF']}), "
              f"lead logXr={sg['leadXr']}, lead logH={sg['leadH']}")
        print(f"   lead logXF2={sg['leadXF2']}, lead logXr2={sg['leadXr2']}")
        print(f"   u2 log u2 = {sg['u2logu2'][0]} * D^{sg['u2logu2'][1]} lambda")
        print(f"   sigma1 = {sg['sigma1']},  sigma2 = {sg['sigma2']}")
        claimed1 = Fr(al) / (a * d0) + (sp.Rational(3, 1) * d0 / (32 * kap * b2) if al == Fr(3, 8) else 0)
        claimed2 = (sp.Rational(al) / a) * (4 / d0 + d0 * b2)
        print("   match claimed sigma1, sigma2:", sp.simplify(sg["sigma1"] - claimed1) == 0,
              sp.simplify(sg["sigma2"] - claimed2) == 0)

    run_checks()

    print("\n== Part B': endpoint 3/8, parent acceptability subtraction ==")
    print("D^C/(B rho2^2) = exp((-a + 2 c2) D^{3/8} lambda (1+o(1))) -> needs a > 2 c2;"
          " weight bridge / FM-U' / Lemma eta denominator need only a > c2 (SQ).")

    print("\n== Part C: cost constant ==")
    for al, P in ((Fr(7, 18), (Fr(3, 5), Fr(1, 8), 1, 2)), (Fr(3, 8), (Fr(4, 5), Fr(1, 8), 1, 2))):
        r = cost(al, *P)
        print(f"alpha={al}, (a,kappa,d0,b2)={P}, C_gamma=0: sigma1={r['sig1']}={float(r['sig1']):.5f}, "
              f"sigma2={r['sig2']}={float(r['sig2']):.5f}, c2={r['c2']}")
        t1, t2, (rl, rh) = r["terms"]
        print(f"   2a={float(t1):.5f}, a+s+max(a,s1+s2)={t2}={float(t2):.5f}, "
              f"rec in [{float(rl):.5f},{float(rh):.5f}]  => c in [{float(r['c_lo']):.5f},{float(r['c_hi']):.5f}]")
        if al == Fr(3, 8):
            print("   a > 2c2 ?", P[0] > 2 * r["c2"])
    lo, hi = two_pow_neg_bracket(Fr(3, 8))
    print(f"1/(1-2^(-3/8)) in [{float(1/(1-lo)):.6f}, {float(1/(1-hi)):.6f}]")
    lo, hi = two_pow_neg_bracket(Fr(7, 18))
    print(f"1/(1-2^(-7/18)) in [{float(1/(1-lo)):.6f}, {float(1/(1-hi)):.6f}]")
    for Cg in (0, 1, 5):
        r = cost(Fr(7, 18), Fr(3, 5), Fr(1, 8), 1, 2, Cg)
        print(f"alpha=7/18 standard point, C_gamma={Cg}: c_hi={float(r['c_hi']):.4f}")
    b2c = grid38(2)
    b1c = grid38(1)
    print(f"3/8 grid (C_gamma=0) under a > 2c2: c <= {float(b2c[0]):.4f} at {b2c[1]}")
    print(f"3/8 grid (C_gamma=0) under a >  c2 (SQ only): c <= {float(b1c[0]):.4f} at {b1c[1]}")

    print("\n== Part D: degeneracy in the open window ==")
    al = Fr(7, 18)
    for t in (1, 2, 4, 8, 16):
        A, K, D0, B2 = Fr(1, t), Fr(1, 8 * t), Fr(4 * t * t), Fr(1, 64 * t ** 5)
        r = cost(al, A, K, D0, B2)
        print(f"t={t}: (a,kappa,d0,b2)=({A},{K},{D0},{B2}) -> c_hi(C_gamma=0)={float(r['c_hi']):.4f}")
    print("Also: for alpha' < alpha both in the window, L[alpha', c'] = o(L[alpha, eps]) for every eps>0,")
    print("so the open-window theorem already implies time L_N[alpha, o(1)]; c(alpha) is cosmetic there.")


if __name__ == "__main__":
    main()
