"""Run midk E1 (ledger) and E2 (independence/construction) and write JSONL to experiments/raw/.

Usage: .venv/bin/python -m src.midk.experiments [E1|E2|all]
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp

from src.midk import ledger as L
from src.midk import independence as I

RAW = Path("experiments/raw")


def _fr(x):
    return str(x)


def run_E1(path: Path):
    recs = []
    # Family A exponents
    info = L.amgm_floor_exponents(1)
    recs.append({"exp": "E1", "family": "A", "item": "amgm_floor",
                 "logX_min": str(info["logX_min"]), "d_star": str(info["d_star"]),
                 "alpha_fallback": _fr(L.family_A_alpha("fallback")),
                 "alpha_Gr": _fr(L.family_A_alpha("Gr")),
                 "pred_P1a": L.family_A_alpha("fallback") == Fr(1, 4),
                 "pred_P1b": L.family_A_alpha("Gr") == Fr(1, 3)})
    sc = L.aggregate_optimum_scaling()
    recs.append({"exp": "E1", "family": "A", "item": "attaining_scaling",
                 "scaling": {k: (str(v) if not isinstance(v, tuple) else [str(v[0]), str(v[1])])
                             for k, v in sc.items()},
                 "u_x_exponent": str(sc["u"][0]), "alpha": _fr(sc["alpha"])})
    # monotonicity of G_r and KL floors in r (so r=2 optimal)
    for d in (6, 10, 20):
        gr = L.gr_floor_monotone(d)
        kl = L.kl_floor_monotone(d)
        gr_mono = all(gr[i][1] <= gr[i + 1][1] for i in range(len(gr) - 1))
        kl_mono = all(kl[i][1] <= kl[i + 1][1] for i in range(len(kl) - 1))
        recs.append({"exp": "E1", "family": "A", "item": "floor_monotone", "d": d,
                     "gr": [[r, _fr(v)] for r, v in gr], "kl": [[r, _fr(v)] for r, v in kl],
                     "gr_nondecreasing": gr_mono, "kl_nondecreasing": kl_mono})
    # Family B ladder
    tab = L.ladder_table(10)
    recs.append({"exp": "E1", "family": "B", "item": "ladder",
                 "table": [[j, _fr(a)] for j, a in tab],
                 "decreasing": all(tab[i][1] > tab[i + 1][1] for i in range(len(tab) - 1)),
                 "all_above_third": all(a > Fr(1, 3) for _, a in tab),
                 "inf": "1/3", "endpoint_j1": _fr(tab[0][1])})
    # Family C floor-sum deterministic
    recs.append({"exp": "E1", "family": "C", "item": "floorsum",
                 "table": [[j, _fr(L.floorsum_exponent(j))] for j in range(2, 7)],
                 "metric": "deterministic N^alpha (not L[1/3]); N^{1/(j+2)} >> L[1/3]"})
    with path.open("w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return recs


def run_E2(path: Path):
    recs = []
    # P2a independence baseline
    for p in (101, 1009, 10007):
        for r in (2, 3, 4):
            res = I.linking_independence_count(p, r, trials=200000, seed=p * 10 + r)
            res.update({"exp": "E2", "item": "P2a_independence"})
            recs.append(res)
    # P2b(II) / Lemma SPL-INT: dependence that lives only mod N splits N; scalar (global) does not.
    primes = [(1000003, 1000033), (99991, 99989), (1299709, 1299721)]
    for (p, q) in primes:
        N = p * q
        for m in (2, 3, 4):
            # structured mod-p-only symmetry -> should split N
            T, lam = I.modp_only_symmetry(p, q, m, lam_p=7 % p, seed=p + m)
            fac = I.split_from_intertwiner(N, T, lam)
            split_ok = fac is not None and 1 < fac < N and N % fac == 0
            recs.append({"exp": "E2", "item": "P2b_modN_dependence_splits", "N_bits": N.bit_length(),
                         "m": m, "split": split_ok, "factor_is_pq": split_ok and fac in (p, q)})
            # global scalar symmetry -> no split (consistent global eigenvector)
            Ts = I.scalar_intertwiner(N, m, lam=7)
            facs = I.split_from_intertwiner(N, Ts, lam=7)
            recs.append({"exp": "E2", "item": "P2b_global_scalar_no_split", "N_bits": N.bit_length(),
                         "m": m, "split": facs is not None})
    # P2b(I) structured constructions pay height / are useless
    for d in (6, 10):
        for r in (3, 4, 5):
            recs.append({"exp": "E2", "item": "P2b_ring_power_height",
                         **I.ring_power_height(d, base_height=1000, r=r)})
    recs.append({"exp": "E2", "item": "P2b_nnorm_fullrank",
                 **I.nnorm_fullrank_minima([1, 0, 0, 0, 0, 2], N=1000003 * 1000033)})
    with path.open("w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return recs


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("E1", "all"):
        run_E1(RAW / "midk_e1.jsonl")
        print("wrote", RAW / "midk_e1.jsonl")
    if which in ("E2", "all"):
        run_E2(RAW / "midk_e2.jsonl")
        print("wrote", RAW / "midk_e2.jsonl")


if __name__ == "__main__":
    main()
