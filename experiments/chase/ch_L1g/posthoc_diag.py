"""Disclosed post-hoc diagnostic for P2 (honest sublattice setup accounting + last-prime windowing)
and P4 at penultimate depth m-1. Does NOT modify pre-registered results.jsonl or summary.json.
"""
import json
import math
import random
from pathlib import Path
from time import perf_counter

import experiment as exp

OUT_DIR = Path(__file__).resolve().parent
W_BUDGET = exp.W_BUDGET


def run_rnd_window(form, K, good, primorial, rng):
    """Random branch with last prime chosen from (2K^2/mod, 4K^2/mod] to avoid trivial >B overshoot."""
    target = 2 * K * K
    good_sorted = sorted(good, key=lambda x: x[0])
    ps = [p for p, _ in good_sorted]
    work = 0
    seen = set()
    hits = set()
    t0 = perf_counter()
    while work < W_BUDGET:
        mod = 1
        r = 0
        used = set()
        # Multiply random primes while 4*K^2 // mod > max_p, then pick last prime in (target//mod, 2*target//mod]
        while mod <= target and work < W_BUDGET:
            need_lo = target // mod + 1
            need_hi = (2 * target) // mod
            if need_lo <= ps[-1]:
                cands = [item for item in good_sorted if need_lo <= item[0] <= need_hi and item[0] not in used]
                if cands:
                    p, rts = rng.choice(cands)
                    z = rts[rng.randrange(2)]
                    r = z if mod == 1 else r + mod * (((z - r) * pow(mod, -1, p)) % p)
                    mod *= p
                    work += 1
                    break
            p, rts = rng.choice(good_sorted)
            if p in used:
                continue
            used.add(p)
            z = rts[rng.randrange(2)]
            r = z if mod == 1 else r + mod * (((z - r) * pow(mod, -1, p)) % p)
            mod *= p
            work += 1
        if mod <= target:
            break
        pt = exp.reconstruct_uv(r, mod, K)
        work += 1
        if pt is not None and pt not in seen and work < W_BUDGET:
            seen.add(pt)
            work += 1
            if exp.is_smooth(exp.eval_q(form, pt[0], pt[1]), primorial):
                hits.add(pt)
    ms = (perf_counter() - t0) * 1000.0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {"work": work, "tested": len(seen), "hits": len(hits), "main_hits": len(main_hits), "ms": round(ms, 2)}


def run_tree_sieve_honest(form, B, K, good, primes, rng):
    """Sublattice sieve charging 1 op per prime root projection onto Lambda_{d',r'} + 1 op per box point + 1 per incidence."""
    lo_d = (2 * K * K) // B
    hi_d = 2 * K * K
    work = 0
    seen = set()
    hits = set()
    t0 = perf_counter()
    while work < W_BUDGET:
        mod = 1
        r = 0
        used = set()
        while mod <= lo_d and work < W_BUDGET:
            p, rts = rng.choice(good)
            if p in used:
                continue
            if mod * p > hi_d:
                if mod > lo_d:
                    break
                if len(used) > 8:
                    break
                continue
            used.add(p)
            z = rts[rng.randrange(2)]
            mod, r = (p, z) if mod == 1 else (mod * p, r + mod * (((z - r) * pow(mod, -1, p)) % p))
            work += 1
        if not (lo_d < mod <= hi_d):
            continue
        pts = exp.enumerate_sublattice_box(r, mod, K)
        work += 1 + len(pts)
        if not pts or work >= W_BUDGET:
            continue
        rems = [exp.eval_q(form, u, v) // mod for u, v in pts]
        for p in primes:
            if work >= W_BUDGET:
                break
            work += 1  # charge 1 modular progression setup per prime on this sublattice
            for idx in range(len(pts)):
                if rems[idx] % p == 0:
                    work += 1
                    while rems[idx] % p == 0:
                        rems[idx] //= p
        if work <= W_BUDGET:
            for pt, rem in zip(pts, rems):
                if pt not in seen:
                    seen.add(pt)
                    if rem == 1:
                        hits.add(pt)
    ms = (perf_counter() - t0) * 1000.0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {"work": work, "tested": len(seen), "hits": len(hits), "main_hits": len(main_hits), "ms": round(ms, 2)}


def measure_penultimate_corr(form, N, K, good, primorial, rng, max_samples=160):
    target = 2 * K * K
    pairs = []
    attempts = 0
    while len(pairs) < max_samples and attempts < 4000:
        attempts += 1
        mod = 1
        r = 0
        used = set()
        prev_mod, prev_r = 1, 0
        while mod <= target:
            p, rts = rng.choice(good)
            if p in used:
                continue
            used.add(p)
            prev_mod, prev_r = mod, r
            z = rts[rng.randrange(2)]
            mod, r = (p, z) if mod == 1 else (mod * p, r + mod * (((z - r) * pow(mod, -1, p)) % p))
        pt = exp.reconstruct_uv(r, mod, K)
        if pt is None or max(pt) <= K // 4 or prev_mod <= 1:
            continue
        qmin_pen, b1 = exp.min_qnorm_on_lattice(form, N, prev_r, prev_mod)
        # Also record whether penultimate shortest Q-vector equals terminal (u,v)
        same_vec = (abs(b1[0]), abs(b1[1])) == pt
        cof = exp.eval_q(form, pt[0], pt[1]) // mod
        sp, _ = exp.smooth_part_and_flag(cof, primorial) if cof > 1 else (1, True)
        frac = 1.0 if cof <= 1 else math.log2(sp) / math.log2(cof)
        pairs.append((-math.log2(qmin_pen), frac, same_vec, prev_mod < (3 * (K // 4) ** 2) // 4))
    return pairs


def main():
    schedule = [(72, 512, 512), (84, 1024, 768), (96, 2048, 1536)]
    pooled = {k: {"hits": 0, "main_hits": 0, "tested": 0, "work": 0, "ms": 0.0}
              for k in ("RND_WINDOW", "TREE_SIEVE_HONEST", "BEAM_SV", "BEAM_QNORM")}
    pen_all = []
    pen_nonvis = []
    for bits, B, K in schedule:
        primes = exp.primes_upto(B)
        primorial = math.prod(primes)
        for rep in range(6):
            seed = 920000 + 1000 * bits + rep
            rng = random.Random(seed)
            N, form = exp.make_fixture(bits, primes, rng)
            good = exp.good_split_primes(form, N, primes)
            rw = run_rnd_window(form, K, good, primorial, random.Random(seed + 21))
            ts = run_tree_sieve_honest(form, B, K, good, primes, random.Random(seed + 22))
            bsv = exp.run_beam(form, N, K, good, primorial, random.Random(seed + 4), False)
            bqn = exp.run_beam(form, N, K, good, primorial, random.Random(seed + 6), True)
            for name, d in (("RND_WINDOW", rw), ("TREE_SIEVE_HONEST", ts), ("BEAM_SV", bsv), ("BEAM_QNORM", bqn)):
                for fld in ("hits", "main_hits", "tested", "work", "ms"):
                    pooled[name][fld] += d[fld]
            for sc, fr, same, below in measure_penultimate_corr(form, N, K, good, primorial, random.Random(seed + 23)):
                pen_all.append((sc, fr, same))
                if below:
                    pen_nonvis.append((sc, fr, same))

    out = {
        "pooled_72_84_96_posthoc": {k: {kk: round(vv, 2) if isinstance(vv, float) else vv for kk, vv in v.items()}
                                    for k, v in pooled.items()},
        "penultimate_corr_all_r": round(exp.pearson_r([(x, y) for x, y, _ in pen_all]), 5),
        "penultimate_corr_all_n": len(pen_all),
        "penultimate_same_vec_all": sum(1 for _, _, s in pen_all if s),
        "penultimate_below_nonvis_r": round(exp.pearson_r([(x, y) for x, y, _ in pen_nonvis]), 5),
        "penultimate_below_nonvis_n": len(pen_nonvis),
        "penultimate_below_nonvis_same_vec": sum(1 for _, _, s in pen_nonvis if s),
    }
    (OUT_DIR / "posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
