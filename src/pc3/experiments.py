"""pc3 experiments E1, E2 (notes/pc3_prereg.md).  DIAGNOSTIC.

E1: two-point structure of the child box (rational block of DC1): empirical pair covariance C_emp vs the
    principal (local x archimedean) model C_mod = (1 + C_loc(det))(1 + C_arch) - 1.
E2: the conditioning on the parent (the 3-point object of PC3) vs the unconditional 2-point object (Prop DC):
    Phi (incidence-weighted conditional relative variance of the child count) vs phi (unconditional), the
    parent-count relative variance sigma^2, for coefficient scales K.

Usage: python -m src.pc3.experiments {E1|E2} out.jsonl [n_labels]
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from . import toy

# frozen parameters (prereg)
PARAMS = dict(y=150, Y=2000, M=2000, S=20, qlo=500, Ks=(20, 80, 320), n_par=1500,
              n_g_E1=300_000, n_g_E2=60_000, X=90_000_000, seed=20261024)


def labels(n: int, rng):
    ps = toy.primes_upto(2 * PARAMS["qlo"])
    ps = ps[ps > PARAMS["qlo"]]
    out = []
    for _ in range(n):
        q = int(rng.choice(ps))
        r = int(rng.integers(0, q))
        out.append((q, r))
    return out


def density_profile(smooth: np.ndarray, lo: int, hi: int, nb: int = 200):
    """Exact smooth density in nb log-bins on [lo, hi)."""
    edges = np.unique(np.round(np.exp(np.linspace(np.log(lo), np.log(hi), nb + 1))).astype(np.int64))
    cs = np.concatenate([[0], np.cumsum(smooth[: edges[-1]], dtype=np.int64)])
    dens = (cs[edges[1:]] - cs[edges[:-1]]) / np.maximum(edges[1:] - edges[:-1], 1)
    return edges, dens


def phi_of(vals: np.ndarray, edges, dens):
    k = np.clip(np.searchsorted(edges, vals, side="right") - 1, 0, len(dens) - 1)
    return dens[k]


def run_E1(out, n_labels, smooth, rng):
    """Principal model evaluated on the SAME values: 1_S(n) -> f(n) = local(n) * phi(|n|) (local = small-prime
    valuation classes, exact empirical factors; phi = exact smooth density by size).  C_emp - C_mod measures the
    non-principal part (+ model error).  The primitive-g law and the pair lattice enter through the values."""
    P = PARAMS
    hi = 2 * P["M"] * P["S"] * 40
    tau_cls, taup = toy.principal_tables(smooth, hi // 16, hi // 2, P["y"])
    edges, dens = density_profile(smooth, 1, hi)
    for (q, r) in labels(n_labels, rng):
        t0 = time.time()
        e1, e2 = toy.label_basis(q, r)
        box = toy.child_box(e1, e2, P["S"])
        pts = np.array([[b[2], b[3]] for b in box], dtype=np.int64)
        g = toy.sample_g(P["M"], P["n_g_E1"], rng)
        npts = len(pts)
        acc = {k: np.zeros((npts, npts)) for k in ("emp", "mod", "loc", "arch")}
        sums = {k: np.zeros(npts) for k in acc}
        for c0 in range(0, len(g), 20_000):
            v = toy.values(g[c0:c0 + 20_000], pts)
            ph = phi_of(v, edges, dens)
            lo_ = toy.principal_proxy(v, tau_cls, taup)
            mats = dict(emp=toy.indicator(smooth, v).astype(np.float64), arch=ph, loc=lo_, mod=lo_ * ph)
            for k, A in mats.items():
                acc[k] += A.T @ A
                sums[k] += A.sum(0)
        ng = len(g)
        iu = np.triu_indices(npts, 1)
        C = {}
        for k in acc:
            mu = sums[k] / ng
            C[k] = ((acc[k] / ng) / np.outer(mu, mu) - 1)[iu]
        p1 = sums["emp"] / ng
        ce, cm = C["emp"], C["mod"]
        rec = dict(q=q, r=r, npts=npts, m=float(p1.sum()), mean_p=float(p1.mean()), n_pairs=int(len(ce)),
                   mean_Cemp=float(ce.mean()), mean_Cmod=float(cm.mean()), mean_Cloc=float(C["loc"].mean()),
                   mean_Carch=float(C["arch"].mean()), mean_resid=float((ce - cm).mean()),
                   mean_abs_resid=float(np.abs(ce - cm).mean()), cov_emp_mod=float(np.cov(ce, cm)[0, 1]),
                   var_mod=float(cm.var(ddof=1)), var_emp=float(ce.var(ddof=1)),
                   noise_var=float(np.mean(1.0 / (ng * np.outer(p1, p1)[iu]))), secs=time.time() - t0)
        qs = np.quantile(cm, np.linspace(0, 1, 11))
        bb = np.clip(np.searchsorted(qs, cm, side="right") - 1, 0, 9)
        rec["bin_mod"] = [float(cm[bb == k].mean()) for k in range(10)]
        rec["bin_emp"] = [float(ce[bb == k].mean()) for k in range(10)]
        out.write(json.dumps(rec) + "\n")
        out.flush()
        print("E1", q, r, npts, "m=%.1f Cemp=%.4f Cmod=%.4f resid=%.4f slope=%.3f (%.0fs)" % (
            rec["m"], rec["mean_Cemp"], rec["mean_Cmod"], rec["mean_resid"],
            rec["cov_emp_mod"] / rec["var_mod"], rec["secs"]), flush=True)


def run_E2(out, n_labels, smooth, parent, rng):
    P = PARAMS
    for (q, r) in labels(n_labels, rng):
        e1, e2 = toy.label_basis(q, r)
        box = toy.child_box(e1, e2, P["S"])
        pts = np.array([[b[2], b[3]] for b in box], dtype=np.int64)
        g = toy.sample_g(P["M"], P["n_g_E2"], rng)
        N = np.zeros(len(g))
        for c0 in range(0, len(g), 20_000):
            N[c0:c0 + 20_000] = toy.indicator(smooth, toy.values(g[c0:c0 + 20_000], pts)).sum(1)
        m = N.mean()
        phi = N.var() / m ** 2
        for K in P["Ks"]:
            t0 = time.time()
            par = toy.sample_parents(e1, e2, K, P["n_par"], rng)
            pp = np.array([[b[2], b[3]] for b in par], dtype=np.int64)
            pp = pp[P["M"] * (np.abs(pp[:, 0]) + np.abs(pp[:, 1])) <= P["X"]]  # |G| <= M(|a|+|b|) <= X
            Z = np.zeros(len(g))
            SY = np.zeros(len(pp))
            SYN = np.zeros(len(pp))
            SYN2 = np.zeros(len(pp))
            for c0 in range(0, len(g), 5_000):
                Yi = toy.indicator(parent, toy.values(g[c0:c0 + 5_000], pp)).astype(np.float64)
                n = N[c0:c0 + 5_000]
                Z[c0:c0 + 5_000] = Yi.sum(1)
                SY += Yi.sum(0)
                SYN += n @ Yi
                SYN2 += (n * n) @ Yi
            ok = SY >= 30
            pE = SY / len(g)
            mE = SYN[ok] / SY[ok]
            vE = SYN2[ok] / SY[ok] - mE ** 2
            w = pE[ok]
            Phi = float((w * vE / mE ** 2).sum() / w.sum())
            Phi_poi = float((w / mE).sum() / w.sum())
            shift = float((w * (mE / m - 1)).sum() / w.sum())
            z = Z.mean()
            sig2 = Z.var() / z ** 2
            # direct weighted object of Prop DC: E[Z (N-m)^2] / (E[Z] m^2)
            dc_direct = float((Z * (N - m) ** 2).mean() / (z * m * m))
            rec = dict(q=q, r=r, K=K, npts=len(pts), n_par=len(pp), n_par_used=int(ok.sum()), m=float(m),
                       phi=float(phi), phi_ex=float(phi - 1 / m), Phi=Phi, Phi_ex=Phi - Phi_poi, shift=shift,
                       z=float(z), sig2=float(sig2), sig2_ex=float(sig2 - 1 / z), dc_direct=dc_direct,
                       mean_pE=float(pE.mean()), secs=time.time() - t0)
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print("E2", q, r, K, "m=%.1f phi_ex=%.4f Phi_ex=%.4f shift=%.4f z=%.1f sig2=%.4f" % (
                m, rec["phi_ex"], rec["Phi_ex"], shift, z, sig2), flush=True)


def main():
    which, path = sys.argv[1], sys.argv[2]
    n_labels = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    rng = np.random.default_rng(PARAMS["seed"] + (0 if which == "E1" else 1))
    t0 = time.time()
    smooth, parent = toy.build_tables(PARAMS["X"], PARAMS["y"], PARAMS["Y"])
    print("tables built in %.0fs" % (time.time() - t0), flush=True)
    with open(path, "w") as out:
        if which == "E1":
            del parent
            run_E1(out, n_labels, smooth, rng)
        else:
            run_E2(out, n_labels, smooth, parent, rng)


if __name__ == "__main__":
    main()
