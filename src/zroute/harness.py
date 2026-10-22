"""zroute experiments (HARNESS: may read instance labels for verification only, after the operational run).

E1: seedless deterministic NFS Z route (base-m, d = 3) + XK-Z + PIN + deterministic Miller.
E2: d = 1 rank growth: deterministic orbit targets (2^j) versus random-exponent seeds (Lemma HR control).
E3: pointwise supply uniformity across N: relation counts for a fixed budget of orbit targets versus uniformly
    random targets (control), over many N.

Usage: PYTHONPATH=. .venv/bin/python -m src.zroute.harness {e1,e2,e3} config.json out.jsonl
"""
from __future__ import annotations

import json
import random
import resource
import sys
import time
from typing import Dict, List

import flint

from src.common import instances
from src.detGRH.core import FactorBase, FactorExit, fb_bound, grh_bound
from src.zroute import core


def _inst(fam: str, bits: int, seed: int):
    if fam == "semi":
        return instances.balanced_semiprime(bits, seed)
    if fam == "tri":
        return instances.multi_prime(bits, 3, seed)
    raise ValueError(fam)


def _rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # macOS: bytes


def _check_factor(inst, g) -> bool:
    return g is not None and 1 < g < inst.N and inst.N % g == 0


# ----------------------------------------------------------------------------- E1

def run_e1(cfg: Dict, out) -> None:
    for fam, bits, seeds in cfg["cells"]:
        for seed in range(seeds):
            inst = _inst(fam, bits, seed)
            N = inst.N
            rec = {"exp": "e1", "fam": fam, "bits": bits, "seed": seed, "N": N}
            t0 = time.time()
            try:
                fb = FactorBase.make(N, fb_bound(N))
                fld = core.base_m_field(N, cfg["d"])
                ctr: Dict = {}
                allrows = core.nfs_rows(fld, fb, cfg["A"], cfg["Bmax"], 10 ** 9, ctr)
                rec.update(n=fb.n, rows_available=len(allrows), pairs_tested=ctr.get("pairs_tested"),
                           bad_prime_skips=ctr.get("bad_prime_skips", 0), disc=str(fld.disc))
                ids = set()
                prefix_needed = []
                for i, r in enumerate(allrows):
                    ids |= set(r.alg)
                    prefix_needed.append(len(ids))
                stages = []
                result = None
                for mult in cfg["surplus_mults"]:
                    # smallest prefix M with M >= n + #ideals(prefix) + mult*n
                    M = next((i + 1 for i in range(len(allrows))
                              if i + 1 >= fb.n + prefix_needed[i] + int(mult * fb.n)), None)
                    if M is None:
                        stages.append({"mult": mult, "M": None})
                        break
                    rows = allrows[:M]
                    c2: Dict = {}
                    t1 = time.time()
                    zs = core.xkz(fld, rows, c2)
                    rels = core.rational_relations(fb, rows, zs, c2)
                    res = core.pin_and_exponent(fb, rels, grh_bound(N))
                    raw = flint.fmpz_mat([r.rat for r in rows]).rank()
                    st = {"mult": mult, "M": M, "C_alg": c2.get("C_alg"), "dim_L1": c2.get("dim_Lambda1"),
                          "unit_rank_seen": c2.get("unit_rank_seen"), "torsion": c2.get("torsion_vectors"),
                          "undecided": c2.get("xk_undecided"), "nonunit": c2.get("NONUNIT", 0),
                          "verify_fail": c2.get("rel_verify_fail"), "rels": len(rels), "rank_out": res["rank"],
                          "rank_raw_rational": raw, "PIN": res["PIN"], "pinned": res["pinned"],
                          "n_small": res["n_small"], "unpinned": res["unpinned"][:20], "full": res["full"],
                          "secs": round(time.time() - t1, 1)}
                    g, base = core.miller_split(N, res["D"])
                    st["miller_base"] = base
                    st["miller_ok"] = _check_factor(inst, g)          # harness verification only
                    stages.append(st)
                    if res["PIN"]:
                        result = st
                        break
                rec["stages"] = stages
                rec["PIN"] = result is not None
                rec["factored"] = bool(result and result["miller_ok"])
                rec["factored_any_stage"] = any(s.get("miller_ok") for s in stages)
            except FactorExit as fe:
                rec["factor_exit"] = fe.where
                rec["factored"] = _check_factor(inst, fe.divisor)
            rec["secs"] = round(time.time() - t0, 1)
            rec["rss_mb"] = round(_rss_mb(), 1)
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(fam, bits, seed, rec.get("PIN"), rec.get("factored"), rec["secs"], rec["rss_mb"], flush=True)


# ----------------------------------------------------------------------------- E2

def core_modrank(n):
    from src.detGRH.core import ModRank
    return ModRank(n)


def run_e2(cfg: Dict, out) -> None:
    for bits, seeds in cfg["cells"]:
        for seed in range(seeds):
            inst = _inst("semi", bits, seed)
            N = inst.N
            fb = FactorBase.make(N, fb_bound(N))
            n = fb.n
            for design in ("orbit", "random"):
                t0 = time.time()
                ctr: Dict = {}
                cap = int(cfg[design + "_cap_mult"] * n)
                if design == "orbit":
                    j0 = N.bit_length() + 1
                    rows = core.orbit_rows(N, fb, 2, j0, cfg["J_cap"], cap, ctr)
                    tested = ctr["targets_tested"]
                else:
                    rng = random.Random(f"zroute-e2|{N}|{seed}")      # control: randomized seeds
                    X = N * N
                    rows, tested = [], 0
                    while len(rows) < cap and tested < cfg["J_cap"]:
                        tested += 1
                        x = [rng.randrange(X) for _ in range(n)]
                        lam = core.ev_fast(fb, x)
                        r = core.split_row(N, fb, lam, x, "random")
                        if r is not None:
                            rows.append(r)
                rats = [r.rat for r in rows]
                # rows needed to reach full rank (mod 2^61-1; = Q-rank lower bound), and rank at checkpoints
                mr = core_modrank(n)
                full_at = None
                for i, r in enumerate(rats):
                    mr.add(r)
                    if full_at is None and mr.rank == n:
                        full_at = i + 1
                        break
                res = core.pin_and_exponent(fb, rats, grh_bound(N)) if rats else {"PIN": False, "D": 0, "rank": 0}
                g, base = core.miller_split(N, res["D"]) if res.get("D") else (None, None)
                rec = {"exp": "e2", "bits": bits, "seed": seed, "N": N, "design": design, "n": n,
                       "targets_tested": tested, "rows": len(rows), "rows_to_full": full_at,
                       "rank": res["rank"], "PIN": res["PIN"], "miller_base": base,
                       "miller_ok": _check_factor(inst, g), "secs": round(time.time() - t0, 1),
                       "rss_mb": round(_rss_mb(), 1)}
                out.write(json.dumps(rec) + "\n")
                out.flush()
                print(bits, seed, design, len(rows), full_at, res["PIN"], rec["miller_ok"], rec["secs"], flush=True)


# ----------------------------------------------------------------------------- E3

def run_e3(cfg: Dict, out) -> None:
    bits, count, J = cfg["bits"], cfg["count"], cfg["J"]
    for seed in range(count):
        inst = _inst("semi", bits, seed)
        N = inst.N
        fb = FactorBase.make(N, fb_bound(N))
        j0 = N.bit_length() + 1
        lam = pow(2, j0, N)
        orbit_hits, seen, cycle = 0, set(), None
        for j in range(j0, j0 + J):
            if lam in seen and cycle is None:
                cycle = j
            seen.add(lam)
            (u, v), _ = core.gauss_reduce(N, lam)
            if v != 0 and fb.vector(u) is not None and fb.vector(v) is not None:
                orbit_hits += 1
            lam = lam * 2 % N
        rng = random.Random(f"zroute-e3|{N}")                          # control
        rand_hits = 0
        for _ in range(J):
            lam = rng.randrange(1, N)
            (u, v), _ = core.gauss_reduce(N, lam)
            if v != 0 and fb.vector(u) is not None and fb.vector(v) is not None:
                rand_hits += 1
        rec = {"exp": "e3", "bits": bits, "seed": seed, "N": N, "J": J, "orbit_hits": orbit_hits,
               "rand_hits": rand_hits, "orbit_distinct": len(seen), "cycle_at": cycle}
        out.write(json.dumps(rec) + "\n")
        if seed % 100 == 0:
            out.flush()
            print(seed, orbit_hits, rand_hits, flush=True)


def main() -> None:
    which, cfg_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    cfg = json.load(open(cfg_path))[which]
    with open(out_path, "w") as out:
        {"e1": run_e1, "e2": run_e2, "e3": run_e3}[which](cfg, out)


if __name__ == "__main__":
    main()
