"""Pre-registered experiments (notes/coverage38_prereg.md).  HARNESS code (uses p, q).

  e2  dichotomy: exhaustive worst case of h_J over S_32 (and S_36) restricted to
      Lehman-hard pairs, vs all pairs, vs the Poisson prediction.
  e2s sampled tail at Y = 2^(n/6), n = 60, 72: uncovered fraction, Lehman-easy vs hard.
  e3  linear-statistics barrier realised: pointwise Kloosterman statistics of real
      configurations x_l = p l^{-1} mod q, covered vs uncovered, vs random-complement
      surrogates (n = 48, band (2^14, 2^15]).

Run:  PYTHONPATH=. .venv/bin/python -m src.coverage38.experiments e2|e3
"""
from __future__ import annotations

import json
import os
import sys
import time
from math import isqrt
from typing import Dict, List

import numpy as np
import sympy

from src.dchworst import coverage as cv
from src.coverage38 import lehman as lm
from src.coverage38 import oracle as orc

RAW = "experiments/raw"


def _write(name: str, rec: Dict) -> None:
    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, f"coverage38_{name}.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")


def h_pred(npairs: int, surv: List[float]) -> int:
    for H, s in enumerate(surv):
        if npairs * s < 1:
            return H
    return len(surv)


def sweep_snapshots(P, Q, primes, Hcap, signed, snaps):
    """Copy of dchworst.coverage.sweep keeping per-pair h_J at the snapshot J's."""
    run = np.full(len(P), Hcap + 1, dtype=np.int64)
    logsurv = np.zeros(Hcap + 1)
    out = {}
    for j, l in enumerate(primes, start=1):
        tab = cv.hmin_table(l, Hcap, signed)
        tab[1] = 1                      # the diagonal weight (1,1) (allowed by Theorem R)
        nu = cv.nu_counts(tab, Hcap)
        with np.errstate(divide="ignore"):
            logsurv += np.log(np.maximum(1e-300, 1.0 - nu / (l - 1)))
        t = cv.t_residues(P, Q, l)
        np.minimum(run, tab[t], out=run)
        if j in snaps:
            out[j] = (run.copy(), np.exp(logsurv).tolist())
    return out


def surv_curve(h: np.ndarray, Hcap: int) -> List[float]:
    hist = np.bincount(h, minlength=Hcap + 2)
    tail = np.cumsum(hist[::-1])[::-1]          # tail[H] = #{h >= H}
    return [int(x) for x in tail]


def e2(ns=(32, 36)):
    Y = 2 ** 13
    band = cv.band(Y)                            # primes = 3 mod 4 in (Y, 2Y], 440
    Hcap = 110
    for n in ns:
        t0 = time.time()
        P, Q = cv.all_pairs(n)
        Pf = 2.0 ** (n / 2)
        Xd = Pf ** 0.4
        score = lm.lehman_score_vec(P, Q, Bmax=int(2 * Xd ** 0.5) + 2)
        snapsJ = [16, 64, 128, len(band)] if n == 32 else [16, 128, len(band)]
        snaps = sweep_snapshots(P, Q, band, Hcap, True, set(snapsJ))
        for G in (4, 16):
            X0 = int(Xd / G)
            easy = lm.lehman_classify(P, Q, X0, Bmax=isqrt(X0) + 1, score=score)
            hard = ~easy
            for J in snapsJ:
                h, surv = snaps[J]
                hs_all = int(h.max())
                hs_hard = int(h[hard].max())
                pred_all = h_pred(len(P), surv)
                pred_hard = h_pred(int(hard.sum()), surv)
                # tail enrichment: top tail = h >= hq where #{h >= hq} >= 0.001 |S|
                tail = surv_curve(h, Hcap)
                hq = max(H for H in range(len(tail)) if tail[H] >= 0.001 * len(P))
                top = h >= hq
                enrich = float(easy[top].mean() / easy.mean())
                hard_top = int((hard & (h == hs_all)).sum())
                rec = {"exp": "e2", "n": n, "Y": Y, "J": J, "G": G, "X0": X0,
                       "npairs": int(len(P)), "n_easy": int(easy.sum()),
                       "frac_easy": float(easy.mean()),
                       "hstar_all": hs_all, "hstar_hard": hs_hard,
                       "pred_all": pred_all, "pred_hard": pred_hard,
                       "tail_H": int(hq), "tail_count": int(top.sum()),
                       "tail_enrich_easy": enrich,
                       "argmax_hard_at_hstar_all": hard_top,
                       "worst_pair_all": [int(P[np.argmax(h)]), int(Q[np.argmax(h)])],
                       "L_worst_all": float(score[np.argmax(h)]),
                       "max_L": float(score.max()), "P25": float(Xd),
                       "secs": round(time.time() - t0, 1)}
                _write("e2", rec)
                print(json.dumps(rec), flush=True)
        del P, Q, snaps, score


def e2s(plan=((60, 200000), (72, 200000)), Hcap=40):
    """Sampled version at the D'-type scaling Y = 2^(n/6) (P ~ Y^3): uncovered fraction
    among Lehman-easy vs Lehman-hard pairs (block 0 = primes 3 mod 4 in (Y, 2Y],
    signed weights with (1,1))."""
    from src.dchworst.experiments import sample_pairs
    from src.dcltheorem.params import pow2_frac
    for n, cnt in plan:
        t0 = time.time()
        P, Q = sample_pairs(n, cnt, seed0=3800 + n, procs=8)
        Y = pow2_frac(n, 1, 6)
        band = cv.band(Y)
        Pf = 2.0 ** (n / 2)
        Xd = Pf ** 0.4
        score = lm.lehman_score_vec(P, Q, Bmax=int(2 * Xd ** 0.5) + 2)
        snaps = sweep_snapshots(P, Q, band, Hcap, True, {len(band)})
        h, surv = snaps[len(band)]
        tail = surv_curve(h, Hcap)
        Href = min(H for H in range(Hcap + 1) if tail[H + 1] <= 0.05 * len(P))
        for G in (4, 16):
            X0 = int(Xd / G)
            easy = lm.lehman_classify(P, Q, X0, Bmax=isqrt(X0) + 1, score=score)
            unc = h > Href
            fe = float(unc[easy].mean())
            fh = float(unc[~easy].mean())
            ne, nh = int(easy.sum()), int((~easy).sum())
            se = (fe * (1 - fe) / max(ne, 1)) ** 0.5
            sh = (fh * (1 - fh) / max(nh, 1)) ** 0.5
            R = fe / fh if fh > 0 else float("nan")
            rec = {"exp": "e2s", "n": n, "Y": Y, "J": len(band), "G": G, "X0": X0,
                   "npairs": int(len(P)), "n_easy": ne, "frac_easy": ne / len(P),
                   "Href": int(Href), "unc_frac_all": float(unc.mean()),
                   "unc_frac_easy": fe, "unc_frac_hard": fh, "R": R,
                   "z": (fe - fh) / max((se * se + sh * sh) ** 0.5, 1e-12),
                   "hstar_all": int(h.max()), "hstar_easy": int(h[easy].max()),
                   "hstar_hard": int(h[~easy].max()),
                   "poisson_unc_frac": float(surv[Href]),
                   "secs": round(time.time() - t0, 1)}
            _write("e2s", rec)
            print(json.dumps(rec), flush=True)


def e3(n=48, nq=3, np_per_q=1000, nsur=100, seed=2026):
    rng = np.random.default_rng(seed)
    Y = 2 ** 14
    band = list(sympy.primerange(Y + 1, 2 * Y + 1))
    J = len(band)
    lo, hi = 2 ** (n - 1), 2 ** n
    qs = []
    q = sympy.nextprime(int(0.83 * 2 ** (n / 2)))
    while len(qs) < nq:
        qs.append(int(q))
        q = sympy.nextprime(q + 1000003)
    for q in qs:
        x1 = orc.points(1, q, band)
        ind = np.zeros(q)
        ind[x1] = 1
        absK1 = np.abs(np.fft.fft(ind))
        plo, phi = max(q // 2 + 1, -(-lo // q)), min(q - 1, (hi - 1) // q)
        ps = set()
        while len(ps) < np_per_q:
            ps.add(int(sympy.prevprime(int(rng.integers(plo + 2, phi)))))
        ps = sorted(p for p in ps if plo <= p <= phi)
        for H in (4, 16):
            # choose k0 so that lam_q = J |S| / q ~ 3
            k0 = 1
            for _ in range(40):
                S = orc.farey_arc_set(q, H, k0)
                lam = J * S.sum() / q
                if lam >= 3.0:
                    break
                k0 = max(k0 + 1, int(k0 * 3.0 / max(lam, 1e-9) * 0.98))
            S = orc.farey_arc_set(q, H, k0)
            lam = J * S.sum() / q
            absS = orc.fourier_abs(S)
            Q0 = int(np.ceil(q / k0))
            Epred = float(absS[1:].sum() / q * (np.sqrt(np.pi) / 2) * np.sqrt(J * (1 - J / q)))
            rows = []
            for p in ps:
                xs = (x1 * p) % q
                hits = int(S[xs].sum())
                st = orc.stats_for_p(absK1, absS, p, q, Q0)
                st.update({"p": p, "hits": hits})
                rows.append(st)
            sur = []
            for _ in range(nsur):
                ys = orc.random_complement(q, S, J, rng)
                ind = np.zeros(q)
                ind[ys] = 1
                aK = np.abs(np.fft.fft(ind))
                sur.append({"E_abs": float((absS[1:] * aK[1:]).sum() / q),
                            "Mlow": float(aK[1:Q0 + 1].max()), "Mall": float(aK[1:].max())})
            rec = {"exp": "e3", "n": n, "q": q, "H": H, "k0": k0, "J": J, "lam_q": float(lam),
                   "Q0": Q0, "E_pred": Epred, "sqrtJ": float(np.sqrt(J)),
                   "hoeff": float(2 * lam + orc.hoeffding_level(J, q - 1)),
                   "real": rows, "surrogate": sur}
            _write("e3", rec)
            unc = [r for r in rows if r["hits"] == 0]
            cov = [r for r in rows if r["hits"] > 0]
            med = lambda L, k: float(np.median([r[k] for r in L])) if L else float("nan")
            print(json.dumps({"q": q, "H": H, "k0": k0, "lam": lam, "n_unc": len(unc),
                              "n_cov": len(cov), "Epred": Epred,
                              "E_unc": med(unc, "E_abs"), "E_cov": med(cov, "E_abs"),
                              "E_sur": med(sur, "E_abs"),
                              "Mlow_unc": med(unc, "Mlow"), "Mlow_cov": med(cov, "Mlow"),
                              "Mlow_sur": med(sur, "Mlow"), "Mall_real": med(rows, "Mall"),
                              "Mall_sur": med(sur, "Mall")}), flush=True)


if __name__ == "__main__":
    {"e2": e2, "e2s": e2s, "e3": e3}[sys.argv[1]]()
