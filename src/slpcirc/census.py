"""E2: census of tiny SLPs (pre-registered in notes/slpcirc_prereg.md).

Usage: python -m src.slpcirc.census SMAX OUT.jsonl
"""
from __future__ import annotations
import json, sys, time
from collections import defaultdict
import flint
from sympy import primerange
from .enum_slp import enumerate_tau
from .classify import classify_factor, TYPES

W = list(primerange(67, 252))
PRIO = {t: i for i, t in enumerate(TYPES)}


def main(smax: int, out: str):
    t0 = time.time()
    tau = enumerate_tau(smax)
    polys = [p for p in tau if len(p) >= 2]  # nonconstant
    fac_cache = {}   # factor coeff tuple -> type
    root_cache = {}  # (factor tuple, p) -> frozenset of roots or None if g == 0 mod p
    dense_pairs = []
    n_pairs = 0
    n_zero_mod_p = 0
    best_min_density = (-1.0, None)
    other_factors = {}
    for F in polys:
        f = flint.fmpz_poly(list(F))
        _, facs = f.factor()
        gs = []
        for g, _e in facs:
            key = tuple(int(c) for c in g.coeffs())
            if key not in fac_cache:
                fac_cache[key] = classify_factor(g)
                if fac_cache[key] == "OTHER" and g.degree() >= 8:
                    other_factors[key] = g
            gs.append((key, fac_cache[key]))
        gs.sort(key=lambda kt: PRIO[kt[1]])
        min_dens = 2.0
        for p in W:
            fp = flint.nmod_poly(list(F), p)
            if fp.is_zero():
                n_zero_mod_p += 1
                min_dens = min(min_dens, 1.0)
                continue
            n_pairs += 1
            seen = {}
            for key, ty in gs:
                rk = (key, p)
                if rk not in root_cache:
                    gp = flint.nmod_poly(list(key), p)
                    root_cache[rk] = None if gp.is_zero() else frozenset(int(r) for r, _ in gp.roots())
                rs = root_cache[rk]
                if rs is None:
                    continue
                for r in rs:
                    if r not in seen:
                        seen[r] = ty
            r = len(seen)
            min_dens = min(min_dens, r / p)
            if 8 * r >= p:
                cnt = defaultdict(int)
                for ty in seen.values():
                    cnt[ty] += 1
                dense_pairs.append({"F": list(F), "tau": tau[F], "p": p, "r": r, "types": dict(cnt)})
        if min_dens > best_min_density[0]:
            best_min_density = (min_dens, F)
    # H2
    m_stats = []
    for key, g in other_factors.items():
        disc = int(g.discriminant()) * int(g.leading_coefficient())
        tot = 0; cnt = 0
        for p in W:
            if disc % p == 0:
                continue
            gp = flint.nmod_poly(list(key), p)
            tot += len(gp.roots()); cnt += 1
        if cnt:
            m_stats.append({"g": list(key), "deg": g.degree(), "m": tot / cnt, "nprimes": cnt})
    good = sum(1 for d in dense_pairs if (d["r"] - d["types"].get("OTHER", 0)) * 10 >= 9 * d["r"])
    Fbest = best_min_density[1]
    res = {
        "exp": "E2", "smax": smax, "npolys": len(polys), "n_pairs": n_pairs, "n_zero_mod_p": n_zero_mod_p,
        "n_dense": len(dense_pairs), "n_dense_group_or_explicit_90": good,
        "H1_frac": good / len(dense_pairs) if dense_pairs else None,
        "dense_failures": [d for d in dense_pairs if (d["r"] - d["types"].get("OTHER", 0)) * 10 < 9 * d["r"]],
        "n_other_deg8": len(m_stats),
        "H2_mean_m": sum(s["m"] for s in m_stats) / len(m_stats) if m_stats else None,
        "H2_max_m": max((s["m"] for s in m_stats), default=None),
        "H2_over2": sorted([s for s in m_stats if s["m"] > 2], key=lambda s: -s["m"]),
        "H3_best_min_density": best_min_density[0], "H3_best_F": list(Fbest) if Fbest else None,
        "type_counts": {t: sum(1 for v in fac_cache.values() if v == t) for t in TYPES},
        "seconds": round(time.time() - t0, 1),
    }
    # H3: types of roots of best F at each p
    if Fbest:
        f = flint.fmpz_poly(list(Fbest))
        res["H3_best_factors"] = [[[int(c) for c in g.coeffs()], fac_cache[tuple(int(c) for c in g.coeffs())]] for g, _ in f.factor()[1]]
    top = sorted(dense_pairs, key=lambda d: -d["r"] / d["p"])[:25]
    res["top_dense"] = top
    with open(out, "a") as fh:
        fh.write(json.dumps(res) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("dense_failures", "H2_over2", "top_dense")}))


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2])
