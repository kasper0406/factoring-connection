"""Harness for experiment E1 of notes/pcest_prereg.md (toy PC for designs 1D and 2D).

Usage: .venv/bin/python -m src.pcest.e1 experiments/configs/pcest_e1.json \
           experiments/raw/pcest_e1.jsonl [workers]
       .venv/bin/python -m src.pcest.e1 --analyze experiments/raw/pcest_e1.jsonl
HARNESS: instance generation uses src/common/instances (factor labels are never read here).
"""
from __future__ import annotations

import json
import random
import sys
from multiprocessing import Pool
from typing import Dict, List

import gmpy2

from ..common.instances import balanced_semiprime, multi_prime
from ..nfsfruit.smooth import primes_upto
from . import geometry as gm
from . import toy


def instance(cls: str, bits: int, seed: int) -> int:
    if cls == "semi":
        return balanced_semiprime(bits, seed).N
    return multi_prime(bits, 3, seed).N


def job(args) -> Dict:
    cls, bits, seed, design, cfg = args
    N = instance(cls, bits, seed)
    d = cfg["d"]
    M = int(gmpy2.iroot(gmpy2.mpz(N), d + 1)[0])
    rng_l = random.Random(f"pcest-e1-ells|{cls}|{bits}|{seed}")
    lo = 1 << cfg["lbits"]
    J = [p for p in primes_upto(2 * lo) if p > lo]
    ells = sorted(rng_l.sample(J, cfg["n_ell"]))
    rng = random.Random(f"pcest-e1|{cls}|{bits}|{seed}|{design}")
    sm_r, sm_a = toy.Smooth(cfg["B_rat"]), toy.Smooth(cfg["B_alg"])
    fields, viol, lemmaQ_checks, lemmaQ_viol = [], 0, 0, 0
    for _ in range(cfg["n_fields"]):
        g, f = toy.sample_field(N, d, M, cfg["H"], design, rng)
        nrat, nboth = [], []
        for l in ells:
            c = counts_two(N, d, g, f, l, cfg["S"], sm_r, sm_a)
            nrat.append(c["nrat"]); nboth.append(c["nboth"]); viol += c["viol"]
            # Lemma Q exact check on the first two rational-smooth points (2D design only)
            if design == "2d" and len(c["pts"]) >= 2 and lemmaQ_checks < cfg["lemmaQ_checks"]:
                P, Pp = c["pts"][0], c["pts"][1]
                try:
                    gens, q, lam = gm.rational_pair_lattice(l, P, Pp)
                    gn = (gm.G_value(*g, *P) // l, gm.G_value(*g, *Pp) // l)
                    if not gm.in_gamma_q(gn[0], gn[1], q, lam):
                        lemmaQ_viol += 1
                except AssertionError:
                    lemmaQ_viol += 1
                lemmaQ_checks += 1
        fields.append(dict(g=list(g), nrat=nrat, nboth=nboth))
    return dict(exp="e1", cls=cls, bits=bits, seed=seed, design=design, N=N, M=M,
                ells=ells, cfg=cfg, fields=fields, viol=viol,
                lemmaQ_checks=lemmaQ_checks, lemmaQ_viol=lemmaQ_viol)


def counts_two(N, d, g, f, l, S, sm_r, sm_a) -> Dict:
    """As toy.counts but with separate rational / algebraic smoothness bounds."""
    g0, g1 = g
    out = dict(nrat=0, nboth=0, viol=0, pts=[])
    if g1 % l == 0:
        return out
    m_l = g0 * pow(g1, -1, l) % l
    from ..hypD.specialq import box_points
    for a, b in box_points(l, m_l, S):
        G = gm.G_value(g0, g1, a, b)
        if G == 0 or G % l:
            out["viol"] += 1
            continue
        if not sm_r(G // l):
            continue
        out["nrat"] += 1
        out["pts"].append((a, b))
        F = gm.hom_eval(f, a, b)
        if (F * g1 ** d - N * b ** d) % G:
            out["viol"] += 1
        if F != 0 and sm_a(F):
            out["nboth"] += 1
    return out


def moments(xs: List[int]):
    n = len(xs)
    m1 = sum(xs) / n
    m2 = sum(x * (x - 1) for x in xs) / n
    return m1, m2


def analyze(path: str) -> Dict:
    recs = [json.loads(l) for l in open(path)]
    by = {}
    for r in recs:
        cells_r = [x for fl in r["fields"] for x in fl["nrat"]]
        cells_b = [x for fl in r["fields"] for x in fl["nboth"]]
        mr1, mr2 = moments(cells_r)
        mb1, mb2 = moments(cells_b)
        pc_rat = mr2 / mr1 ** 2 if mr1 else None
        pc = mb2 / mb1 ** 2 if mb1 else None
        pc_alg = (mb2 / mr2) / (mb1 / mr1) ** 2 if (mr2 and mb1) else None
        gamma = sum(1 for x in cells_b if x > 0) / len(cells_b)
        key = (r["cls"], r["bits"], r["seed"])
        by.setdefault(key, {})[r["design"]] = dict(
            EN_rat=mr1, EN=mb1, PC_rat=pc_rat, PC=pc, PC_alg=pc_alg, gamma=gamma,
            viol=r["viol"], lemmaQ=(r["lemmaQ_checks"], r["lemmaQ_viol"]))
    rows = []
    for key in sorted(by):
        one, two = by[key].get("1d"), by[key].get("2d")
        rows.append(dict(inst=list(key), d1=one, d2=two))
    return dict(rows=rows)


def main():
    if sys.argv[1] == "--analyze":
        print(json.dumps(analyze(sys.argv[2]), indent=1))
        return
    cfg_all = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    jobs = [(c, b, s, des, cfg_all["params"]) for c in cfg_all["classes"] for b in cfg_all["bits"]
            for s in range(cfg_all["seeds"]) for des in ("1d", "2d")]
    with Pool(workers) as pool, open(out, "w") as fh:
        for rec in pool.imap_unordered(job, jobs):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
