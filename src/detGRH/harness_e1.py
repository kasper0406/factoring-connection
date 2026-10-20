"""detGRH E1 harness (HARNESS: imports diagnostics that know the factors).

Usage:
  PYTHONPATH=. .venv/bin/python -m src.detGRH.harness_e1 experiments/configs/detGRH_e1.json \
      experiments/raw/detGRH_e1.jsonl
"""
from __future__ import annotations

import json
import sys
import time
from math import gcd
from typing import Dict, List

from src.common.metrics import append_jsonl

from . import core, diag


def e1a(cfg: Dict) -> Dict:
    Nmax, Mmax = cfg["e1a_Nmax"], cfg["e1a_mult_Nmax"]
    viol1, viol2, tested1, tested2 = [], [], 0, 0
    for N in range(3, Nmax, 2):
        tested1 += 1
        pred = 1 if diag.all_3mod4(N) else 0
        got = diag.count_sections(N, brute=True)
        if got != pred:
            viol1.append((N, got, pred))
        if N < Mmax and pred == 1:
            tested2 += 1
            if not diag.principal_root_multiplicative(N):
                viol2.append(N)
    return {"exp": "E1a", "tested_MR1": tested1, "viol_MR1": viol1, "tested_MR2": tested2, "viol_MR2": viol2}


def run_instance(cls: str, bits: int, seed: int, cfg: Dict) -> Dict:
    t0 = time.time()
    rec: Dict = {"exp": "E1bcd", "cls": cls, "bits": bits, "seed": seed}
    # instance: all factors must exceed B so that trial division does not end the run
    Bn = core.fb_bound(1 << (bits - 1))
    N, fs = diag.make_instance(cls, bits, seed, min_factor=core.fb_bound(1 << bits))
    rec["N"] = N
    rec["omega"] = len(fs)
    try:
        fb = core.FactorBase.make(N)
    except core.FactorExit as fe:
        rec["factor_exit"] = fe.where
        return rec
    n = fb.n
    rec["B"], rec["n"] = fb.B, n
    abound = core.grh_bound(N)
    rec["a_max"] = abound
    try:
        rel = core.collect(N, fb, K=cfg["K"], X=fb.B, T=cfg["T"], surplus=cfg["surplus"])
    except core.FactorExit as fe:
        rec["factor_exit"] = fe.where
        return rec
    rec["n_kernel"] = len(rel.kernel)
    rec["n_squares"] = len(rel.squares)
    rec["tested"] = rel.tested
    rk, D = core.lattice_index(rel.kernel, n)
    rec["qrank"] = rk
    rec["full_rank"] = rk == n
    rec["D_bits"] = D.bit_length()
    # ---------------- E1b: Z route
    if rk == n:
        g, a = core.miller(N, D, abound)
        rec["miller_factor_proper"] = g is not None and 1 < g < N and N % g == 0
        rec["miller_base"] = a
        g2, a2 = core.miller(N, D << n, abound)      # control: lattice 2*Lambda_c has index 2^n D
        rec["ctrl_miller_proper"] = g2 is not None and 1 < g2 < N and N % g2 == 0
    rec["gen_small_primes"] = diag.small_primes_generate(N, fs, abound + 1)
    # F2 route on Lambda_c and on the control 2*Lambda_c
    ndep, nfr, fac = core.f2_fruitful(fb, rel.kernel)
    rec["f2_deps"], rec["f2_fruitful_basis"] = ndep, nfr
    rec["f2_fruitful"] = nfr > 0
    dbl = [[2 * a for a in v] for v in rel.kernel]
    ndep2, nfr2, _ = core.f2_fruitful(fb, dbl)
    rec["ctrl_f2_deps"], rec["ctrl_f2_fruitful_basis"] = ndep2, nfr2
    # ---------------- E1c: full lattice
    Lam = diag.full_lattice(fb, fs)
    detL = diag.lattice_det(Lam)
    phi = 1
    for p in fs:
        phi *= p - 1
    rec["GB_full"] = detL == phi
    rec["dimWprime"] = diag.dim_W(Lam)
    rec["dimW"] = rec["dimWprime"] - 1
    if rk == n:
        rec["c"] = diag.two_rank_quotient(Lam, rel.kernel)
    # ---------------- E1d: canonical completion y, defects
    basis = core.F2Basis(rel.kernel, n)
    ys: List = []
    nchi, nondiag, viol = 0, 0, 0
    for (r, R, eR) in rel.squares:
        y = core.canonical_root(fb, basis, eR)
        ys.append(y)
        if y is None:
            continue
        if (y * y - R) % N:
            viol += 1
        nchi += 1
        chi = y * pow(r, -1, N) % N
        if chi not in (1, N - 1):
            nondiag += 1
    rec["chi_n"], rec["chi_nondiag"] = nchi, nondiag
    pairs = smooth_pairs = defect_nontriv = gcd_proper = 0
    m = len(rel.squares)
    done = False
    for i in range(m):
        if done:
            break
        if ys[i] is None:
            continue
        for j in range(i + 1, m):
            if ys[j] is None:
                continue
            pairs += 1
            if pairs > cfg["pair_cap"]:
                done = True
                break
            Rij = rel.squares[i][1] * rel.squares[j][1] % N
            e = fb.vector(Rij)
            if e is None:
                continue
            yij = core.canonical_root(fb, basis, e)
            if yij is None:
                continue
            smooth_pairs += 1
            if (yij * yij - Rij) % N:
                viol += 1
            prod_ = ys[i] * ys[j] % N
            delta = yij * pow(prod_, -1, N) % N
            if delta * delta % N != 1:
                viol += 1
            if delta not in (1, N - 1):
                defect_nontriv += 1
                g = gcd(yij - prod_, N)
                if 1 < g < N:
                    gcd_proper += 1
                else:
                    viol += 1
    rec.update({"pairs": min(pairs, cfg["pair_cap"]), "smooth_pairs": smooth_pairs,
                "defect_nontriv": defect_nontriv, "defect_gcd_proper": gcd_proper, "cy_violations": viol})
    rec["wall_s"] = round(time.time() - t0, 2)
    return rec


def main(cfg_path: str, out_path: str) -> None:
    cfg = json.load(open(cfg_path))
    open(out_path, "w").close()
    t = time.time()
    append_jsonl(out_path, e1a(cfg))
    print("E1a done", round(time.time() - t, 1), flush=True)
    for cls in cfg["classes"]:
        for bits in cfg["bits"]:
            b = bits
            if cls == "tri" and str(bits) in cfg.get("tri_bits_replacement", {}):
                b = cfg["tri_bits_replacement"][str(bits)]
            for seed in range(cfg["seeds"]):
                rec = run_instance(cls, b, seed, cfg)
                append_jsonl(out_path, rec)
            print(cls, b, "done", round(time.time() - t, 1), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
