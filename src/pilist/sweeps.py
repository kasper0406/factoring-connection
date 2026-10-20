"""ANALYSIS (not an experiment): Prop JR numeric sweeps (unweighted and weighted) at asymptotic MNFS parameters.
Run: PYTHONPATH=. .venv/bin/python -m src.pilist.sweeps"""
from math import exp, log, sqrt
from src.pilist.decoder import boneh_threshold, lam_of, log_P, primorial_exps, rankin_log_bound, weighted_log_bound
def unweighted():
    pass
    pass
    for bits in (512, 1024, 4096, 2**16):
        x = bits*log(2); Ls = x**(1/3)*log(x)**(2/3)
        logB = 0.951*Ls; logE = 0.951*Ls; logU = 2*logE
        # rational value size: N^{1/d} E, d = delta (x/log x)^{1/3}, delta=1.402
        d = 1.402*(x/log(x))**(1/3); logy = x/d + logE; u = logy/logB
        lam = lam_of(u, logB)
        for logV, tag in ((logE, 'row'), (logU, 'box')):
            best = (-1e9, None)
            for z in [10,15,20,30,40,60,80,120,160,240,320,480,640,960,1280,2560,5120,10240]:
                for pw in (1, 2, 3, 6):
                    Pd = primorial_exps(z, pw)
                    lp = log_P(Pd)
                    thr = sqrt(logV*lp)
                    if thr > logy or thr > lp: continue   # nothing decodable (S > P or S > value)
                    lb, s = rankin_log_bound(logV, Pd, thr, lam)
                    if lb > best[0]: best = (lb, (z, pw, round(lp,1), round(thr,1), s))
            print(bits, tag, 'u=%.2f lam=%.4f logV=%.1f logy=%.1f' % (u, lam, logV, logy), 'max log JR bound = %.3f' % best[0], best[1], 'L-scale', round(Ls,1))

def weighted():
    pass
    pass
    pass
    res = []
    for bits in (1024, 65536):
        x = bits*log(2); Ls = x**(1/3)*log(x)**(2/3)
        logB = 0.951*Ls; logE = 0.951*Ls
        d = 1.402*(x/log(x))**(1/3); logy = x/d + logE; u = logy/logB; lam = lam_of(u, logB)
        for logV, tag in ((logE, 'row'), (2*logE, 'box')):
            best = (-1e9, None)
            for z2 in (40, 160, 640, 2560):
                for z1 in (5, 20, 80):
                    if z1 >= z2: continue
                    for w2 in (0.05, 0.2, 0.5, 1.0, 2.0):
                        for pw in (1, 3):
                            Pd = primorial_exps(z2, pw)
                            w = {q: (1.0 if q <= z1 else w2) for q in Pd}
                            T = sqrt(logV*sum(e*w[q]**2*log(q) for q, e in Pd.items()))
                            Wmax = min(sum(e*w[q]*log(q) for q, e in Pd.items()), max(w.values())*logy)
                            if T > Wmax: continue
                            lb, s = weighted_log_bound(logV, Pd, w, lam)
                            if lb > best[0]: best = (lb, ('z1', z1, 'z2', z2, 'w2', w2, 'pw', pw, 's', s))
                    # power-law family
            for z in (40, 160, 640):
                for a in (-1.0, -0.5, 0.5, 1.0):
                    Pd = primorial_exps(z, 1)
                    w = {q: log(q)**a for q in Pd}
                    T = sqrt(logV*sum(e*w[q]**2*log(q) for q, e in Pd.items()))
                    Wmax = min(sum(e*w[q]*log(q) for q, e in Pd.items()), max(w.values())*logy)
                    if T > Wmax: continue
                    lb, s = weighted_log_bound(logV, Pd, w, lam)
                    if lb > best[0]: best = (lb, ('powlaw z', z, 'a', a, 's', s))
            print(bits, tag, 'max weighted log bound %.3f' % best[0], best[1])

if __name__ == "__main__":  # pragma: no cover
    unweighted(); weighted()
