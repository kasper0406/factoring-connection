"""audit2L: from-scratch exact checks of the exponent comparisons behind Theorem 2L-final.

Written independently of src/thmAB, src/c12rat, src/c12alg (no imports from them).  Every quantity is a
leading-order exponent e such that the corresponding factor is exp(Theta(D^e * polylog(lambda))), with
D = log N / log log N.  A comparison "loss << saving" holds asymptotically iff exponent(loss) < exponent(saving)
(polylog(lambda) factors never matter when the D-exponents differ).  All arithmetic is exact (Fraction).

Run:  .venv/bin/python -m src.audit2L.exponents [out.json]
"""
from fractions import Fraction as F
import json
import sys

HALF = F(1, 2)


def sizes(al):
    """Leading-order D-exponents of the design 2L_DW(alpha) (log of each quantity ~ D^e * lambda)."""
    return dict(
        logB=al,                    # B = L[alpha, a]
        logB2=4 * al - 1,           # B2 = L[4alpha-1, b2]
        logY=1 - al,                # special-q primes, k = kappa D^alpha
        logM=2 * al,                # M = N^mu, mu ~ 1/d, d ~ D^{1-2alpha}
        logXr=2 * al,               # level-1 rational cofactor
        logXF=2 - 3 * al,           # level-1 algebraic value (d log V1 term)
        logXr2=2 * al,              # child rational value
        logXF2=2 * al,              # child algebraic value / q
        eta0logN=3 * al - HALF,     # eta0 log N
        logV1sq=1 - al,             # V1^2 ~ Y S1^2
        logV2=4 * al - 1,           # V2 = 4 Q2 S2 (Q2 <= B2 dominates S2 = L[alpha])
        logd=0,                     # log d ~ lambda  (D-exponent 0)
    )


def exps(al):
    s = sizes(al)
    gam = -al                       # gamma = 1 - beta ~ lambda / log B  (children, rational values, y = B)
    gam2 = 1 - 4 * al               # gamma_2 ~ lambda / log B2 (parents' algebraic values, y = B2)
    u_child = s['logXF2'] - s['logB']           # = alpha
    u2 = s['logXF'] - s['logB2']                # = 3 - 7 alpha
    theta = max(3 * al - 1, 1 - 2 * al, 2 * al - HALF)
    theta_p = max(1 - 2 * al, 2 - 5 * al)
    zeta = (theta + al) / 2
    zeta_p = (theta_p + 4 * al - 1) / 2
    return dict(gam=gam, gam2=gam2, u_child=u_child, u2=u2, theta=theta, theta_p=theta_p,
                zeta=zeta, zeta_p=zeta_p, **s)


def checks(al):
    """List of (name, lhs_exponent, rhs_exponent, relation) with relation '<' or '<='; True if it holds."""
    e = exps(al)
    th, thp, z, zp = e['theta'], e['theta_p'], e['zeta'], e['zeta_p']
    C = []
    add = lambda name, lhs, rhs, rel='<': C.append((name, lhs, rhs, rel))
    # --- sizes / design
    add('u_child = alpha', e['u_child'], al, '<=')
    add('u_child = alpha (rev)', al, e['u_child'], '<=')
    add('u2 = 3-7a', e['u2'], 3 - 7 * al, '<=')
    add('level-1 alg value: d log V1 dominates (2-3a > 2a)', 2 * al, 2 - 3 * al)
    add('eta0 = o(mu): 3a-1/2 < 2a', 3 * al - HALF, 2 * al)
    add('N^eta0 beats L[alpha]^O(1)', al, 3 * al - HALF)
    add('N^eta0 beats B2^O(1)', 4 * al - 1, 3 * al - HALF)
    add('B2 << V1 (child box inside level-1 scale)', 4 * al - 1, 1 - al)
    add('B2 > B (children exist): 4a-1 > a (fails below 1/3)', al, 4 * al - 1)
    # --- algebraic level-1 cost is L^{o(1)} in the open window (u2 log u2 vs D^alpha); tie at 3/8
    add('alg level-1 smoothness cost u2 log u2 vs L[alpha]', e['u2'], al, '<=')
    # --- children: losses vs BDH saving exp(-c D^alpha / lambda^2)
    add('theta < alpha (child loss vs child saving)', th, al)
    add('theta == 2a-1/2 on window', th, 2 * al - HALF, '<=')
    add('thin child loss (V2^2)^gamma: 4a-1-a <= theta', e['logV2'] + e['gam'], th, '<=')
    add('filter spread (N^eta0/4)^gamma <= theta', e['eta0logN'] + e['gam'], th, '<=')
    add('SPR alg spread d*gamma_a*log S2: (1-2a)+(-a)+a <= theta', (1 - 2 * al) + e['gam'] + al, th, '<=')
    add('SPR alg spread d*gamma_a*log T0: (1-2a)-a+theta < theta', (1 - 2 * al) - al + th, th)
    add('h_* >= B^{-5/12}: theta < alpha', th, al)
    add('zeta in (theta, alpha): theta<zeta', th, z)
    add('zeta < alpha', z, al)
    add('EXC-A Rankin mid: zeta > alpha - 7/20', al - F(7, 20), z)
    add('child modulus small conductor: log Delta\' (alpha) < log x (2alpha)', al, 2 * al)
    add('omega(Delta\') = o(P0): alpha < 7/10', al, F(7, 10))
    add('beta-drift^2 children: 2zeta-4a <= -2a', 2 * z - 4 * al, -2 * al, '<=')
    add('first-order child error^2: -2a < budget -a', -2 * al, -al)
    # --- parents
    add('theta_p < 3-7a (parent loss vs B2 saving)', thp, 3 - 7 * al)
    add('theta_p < 4a-1 (h_*^p >= B2^{-5/12})', thp, 4 * al - 1)
    add('zeta_p > theta_p', thp, zp)
    add('zeta_p < 4a-1', zp, 4 * al - 1)
    add('EXC-A parents: zeta_p > 1-a-7/20', 1 - al - F(7, 20), zp)
    add('parent modulus small conductor: 1-a < 2-3a', 1 - al, 2 - 3 * al)
    add('omega(parent Delta\') = o(P0): 1-a < 7/10', 1 - al, F(7, 10))
    add('rational parent pair loss (V1^2)^gamma: 1-2a', (1 - al) + e['gam'], thp, '<=')
    add('alg parent modulus loss (Delta\')^gamma2: 2-5a', (1 - al) + e['gam2'], thp, '<=')
    add('RAT-2 rational parents: 1-2a < alpha', 1 - 2 * al, al)
    add('parent beta-drift^2 8a-4 <= -2a', 8 * al - 4, -2 * al, '<=')
    add('parent E_FO^2 (1-4a part): 2-8a <= -2a', 2 - 8 * al, -2 * al, '<=')
    add('u2 log u2 << log B2 (Hildebrand err << 1/u2): 3-7a < 4a-1', 3 - 7 * al, 4 * al - 1)
    add('recursion on parent values L_{XF}[2/5]: (2/5)(2-3a) < a', F(2, 5) * (2 - 3 * al), al)
    add('Hildebrand range at B2: 0 < 4a-1', 0, 4 * al - 1)
    # --- precision budget (Prop PB): eta ~ D^{7a-3}
    eta = 7 * al - 3
    add('PB first order: A/eta, -a < 7a-3 (needs a>3/8)', -al, eta)
    add('PB first order: sqrt(AB)/eta, (-a+7a-3)/2 < 7a-3', (-al + eta) / 2, eta)
    A2 = -2 * al
    B2e = max(-2 * al, -2 * (3 - 7 * al))
    add('PB second order: A2 < eta', A2, eta)
    add('PB second order: sqrt(A2 B2) < eta', (A2 + B2e) / 2, eta)
    add('acceptability subtraction: 3-7a < a (open window; tie at 3/8)', 3 - 7 * al, al)
    add('field rejections: IR exp(-D^{4a-1}), content 2^{-d}: exps > 0', 0, min(4 * al - 1, 1 - 2 * al))
    return [(n, l, r, rel, (l < r) if rel == '<' else (l <= r)) for (n, l, r, rel) in C]


def cost(al, a, kappa, d0, b2, Cg=0):
    """c(alpha) of Theorem 2L-final, with 1/(1-2^{-alpha}) bracketed exactly (returned as float bracket)."""
    al, a, kappa, d0, b2 = map(F, (al, a, kappa, d0, b2))
    s = 2 * kappa * (1 + Cg)
    sig1 = al / (a * d0)
    c2 = 3 * d0 / (32 * kappa * b2)
    if al == F(3, 8):
        sig1 += c2
    sig2 = (al / a) * (4 / d0 + d0 * b2)
    # 2^{-alpha} = t: bisection on t^{den} vs 2^{-num} exactly
    p, q = al.numerator, al.denominator
    lo, hi = F(0), F(1)
    for _ in range(60):
        mid = (lo + hi) / 2
        if mid ** q * 2 ** p < 1:
            lo = mid
        else:
            hi = mid
    rec_hi = (a + s) / (1 - hi)
    rec_lo = (a + s) / (1 - lo)
    terms = dict(lin=2 * a, rel=a + s + max(a, sig1 + sig2), rec=(float(rec_lo), float(rec_hi)))
    c = max(float(terms['lin']), float(terms['rel']), float(rec_hi))
    return dict(c=c, sigma1=float(sig1), sigma2=float(sig2), s=float(s), c2=float(c2),
                a_gt_2c2=bool(a > 2 * c2), terms={k: (float(v) if not isinstance(v, tuple) else v)
                                                   for k, v in terms.items()})


def main(out=None):
    res = {}
    grid = [F(3, 8)] + [F(3, 8) + F(k, 400) for k in range(1, 10)] + [F(7, 18)] + [F(2, 5)]
    for al in grid:
        cs = checks(al)
        res[str(al)] = dict(failed=[(n, str(l), str(r), rel) for (n, l, r, rel, ok) in cs if not ok],
                            n=len(cs), theta=str(exps(al)['theta']), theta_p=str(exps(al)['theta_p']),
                            child_gap=str(al - exps(al)['theta']),
                            parent_gap=str(3 - 7 * al - exps(al)['theta_p']))
    res['cost_7_18'] = cost(F(7, 18), F(3, 5), F(1, 8), 1, 2)
    res['cost_3_8'] = cost(F(3, 8), F(4, 5), F(1, 8), 1, 2)
    for k in ('3/8', '7/18', '2/5'):
        print(k, 'failed:', res[k]['failed'], 'gaps child/parent:', res[k]['child_gap'], res[k]['parent_gap'])
    print('interior grid all ok:', all(not res[str(F(3, 8) + F(k, 400))]['failed'] for k in range(1, 10)))
    print('cost 7/18:', res['cost_7_18'])
    print('cost 3/8 :', res['cost_3_8'])
    if out:
        json.dump(res, open(out, 'w'), indent=1)
    return res


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
