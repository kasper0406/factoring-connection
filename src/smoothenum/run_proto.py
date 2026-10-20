"""Pre-registered prototype runs (notes/smoothenum_prereg.md). Streams JSONL; one row / one lattice in memory.

  .venv/bin/python -m src.smoothenum.run_proto calib     # P_s curves (random integers), side tests
  .venv/bin/python -m src.smoothenum.run_proto e1        # box line-sieve baseline (rational + doubly)
  .venv/bin/python -m src.smoothenum.run_proto e2        # smooth-modulus lattice sweep (rational side)
  .venv/bin/python -m src.smoothenum.run_proto e3        # bottom-up generation (toy instance)
  .venv/bin/python -m src.smoothenum.run_proto e4ours    # special-(q1,q2) to K distinct relations
  (las is run by the shell commands in the prereg; e4check parses its relations)

Floating point is used only for sieve thresholds and for the *predictions* (interpolated P_s curves);
every counted relation is confirmed by exact integer tests.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from math import gcd, log2
from typing import Dict, List

import numpy as np

from src.smoothenum.mech import (Counters, bottom_up_hit, box_line_sieve, crt_root, lattice_count,
                                 lattice_points, lattice_sample)
from src.smoothenum.nfsdata import (PolyPair, SideTest, algebraic_fb, base_m_pair, load_cado_poly,
                                    primes_upto, rational_fb)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw", "smoothenum")
POLY = os.path.join(RAW, "c61.poly")
# CADO params.c60
LIM0, LIM1, LPB0, LPB1, MFB0, MFB1 = 78682, 111342, 18, 19, 17, 38
A_BOX, BM_BOX = 2 ** 21, 2 ** 9          # box |a| < 2^21, 1 <= b <= 2^9; U = 2^31 points
U_BOX = 2 * A_BOX * BM_BOX
SEED = 20261020


def setup():
    P = load_cado_poly(POLY)
    return P, rational_fb(P, LIM0), algebraic_fb(P, LIM1), SideTest(LIM0, LPB0, MFB0), SideTest(LIM1, LPB1, MFB1)


def ps_curve(test: SideTest, bits: List[float], n: int, seed: int) -> Dict[float, Dict[str, float]]:
    """Fraction of uniform random integers in [2^(k-1/2), 2^(k+1/2)) passing the side test (prediction input)."""
    rng = random.Random(seed)
    out = {}
    for k in bits:
        lo, hi = int(2 ** (k - 0.5)), int(2 ** (k + 0.5))
        hits = sum(1 for _ in range(n) if test.ok(rng.randrange(lo, hi)))
        out[k] = {"n": n, "hits": hits, "p": hits / n}
    return out


def load_curve(name: str):
    d = json.load(open(os.path.join(RAW, name)))
    items = sorted((float(k), v["p"]) for k, v in d.items())
    ks = np.array([k for k, _ in items])
    ps = np.array([max(p, 1e-12) for _, p in items])
    return ks, np.log(ps)


def ps_interp(curve, logbits: np.ndarray) -> np.ndarray:
    ks, lps = curve
    v = np.exp(np.interp(np.clip(logbits, ks[0], ks[-1]), ks, lps))
    # below the first calibration size (2^14 < lim0) every integer is factor-base smooth
    return np.where(logbits < ks[0], 1.0, v)


def cmd_calib():
    P, rfb, afb, rt, at = setup()
    t0 = time.time()
    c0 = ps_curve(rt, [float(k) for k in range(14, 62, 2)], 20000, SEED)
    json.dump({str(k): v for k, v in c0.items()}, open(os.path.join(RAW, "ps_rat.json"), "w"), indent=1)
    c1 = ps_curve(at, [float(k) for k in range(20, 116, 4)], 60000, SEED + 1)
    json.dump({str(k): v for k, v in c1.items()}, open(os.path.join(RAW, "ps_alg.json"), "w"), indent=1)
    print("calib done", time.time() - t0)


def e1_rows():
    rng = random.Random(SEED + 10)
    rows = []
    for k in range(32):
        b = 1 + 16 * k
        a0 = rng.randrange(-A_BOX, A_BOX - 2 ** 16)
        rows.append((b, a0, a0 + 2 ** 16))
    return rows


def cmd_e1():
    P, rfb, afb, rt, at = setup()
    c0 = load_curve("ps_rat.json")
    out = open(os.path.join(RAW, "e1_box.jsonl"), "w")
    tot = {"points": 0, "coprime": 0, "rat_found": 0, "rat_pred": 0.0, "both_found": 0, "t_rat": 0.0, "t_both": 0.0}
    for (b, lo, hi) in e1_rows():
        t = time.time()
        rels, cnt = box_line_sieve(P, A_BOX, b, b + 1, rfb, rt, a_window=(lo, hi))
        t_r = time.time() - t
        t = time.time()
        rels2, cnt2 = box_line_sieve(P, A_BOX, b, b + 1, rfb, rt, afb, at, a_window=(lo, hi))
        t_b = time.time() - t
        a = np.arange(lo, hi, dtype=np.int64)
        cop = np.gcd(a, b) == 1
        gl = np.log2(np.abs(P.Y1 * a.astype(np.float64) + P.Y0 * float(b)))
        pred = float(np.sum(ps_interp(c0, gl[cop])))
        rec = {"b": b, "a_lo": lo, "a_hi": hi, "points": hi - lo, "coprime": int(cop.sum()),
               "rat_found": len(rels), "rat_pred": pred, "both_found": len(rels2),
               "mean_log2G": float(gl.mean()), "t_rat": t_r, "t_both": t_b,
               "rat_counters": cnt.as_dict(), "both_counters": cnt2.as_dict(), "both_rels": rels2}
        out.write(json.dumps(rec) + "\n")
        out.flush()
        for k in ("points", "coprime", "rat_found", "rat_pred", "both_found", "t_rat", "t_both"):
            tot[k] += rec[k]
    # recall check on the first two rows: exact test of every coprime point
    rec_ex, rec_sv = 0, 0
    for (b, lo, hi) in e1_rows()[:2]:
        rels, _ = box_line_sieve(P, A_BOX, b, b + 1, rfb, rt, a_window=(lo, hi))
        rec_sv += len(rels)
        rec_ex += sum(1 for a in range(lo, hi) if gcd(a, b) == 1 and rt.ok(P.G(a, b)))
    tot.update(recall_rat=rec_sv / rec_ex, recall_exact=rec_ex, recall_sieve=rec_sv,
               ratio_obs_pred=tot["rat_found"] / tot["rat_pred"],
               Ybox_rat=tot["rat_found"] / tot["points"], Ybox_both=tot["both_found"] / tot["points"])
    json.dump(tot, open(os.path.join(RAW, "e1_summary.json"), "w"), indent=1)
    print(json.dumps(tot, indent=1))


def cmd_e1exact():
    """Post-review check (codex smoothenum_2): exact two-sided test of every coprime point of the E1 windows,
    to measure the sieve's two-sided recall (the E4 gain is quoted against this exact baseline)."""
    P, rfb, afb, rt, at = setup()
    tot = {"points": 0, "rat_exact": 0, "both_exact": 0}
    for (b, lo, hi) in e1_rows():
        for a in range(lo, hi):
            if gcd(a, b) != 1:
                continue
            if rt.ok(P.G(a, b)):
                tot["rat_exact"] += 1
                if at.ok(P.F(a, b)):
                    tot["both_exact"] += 1
        tot["points"] += hi - lo
    tot["Ybox_both_exact"] = tot["both_exact"] / tot["points"]
    tot["Ybox_rat_exact"] = tot["rat_exact"] / tot["points"]
    json.dump(tot, open(os.path.join(RAW, "e1_exact.json"), "w"), indent=1)
    print(json.dumps(tot, indent=1))


def random_modulus(rng: random.Random, primes: List[int], roots: Dict[int, int], t: int):
    """Squarefree M in [2^t, 2^(t+1)) from distinct rational FB primes >= 2^8 (rejection sampling)."""
    while True:
        M, used = 1, []
        while M < 2 ** t:
            p = rng.choice(primes)
            if p in used:
                continue
            used.append(p)
            M *= p
        if M < 2 ** (t + 1):
            return M, used


E2_BUCKETS = {8: 60, 12: 60, 16: 60, 20: 60, 24: 60, 28: 60, 31: 2000, 34: 4000, 37: 20000, 40: 40000}


def cmd_e2():
    P, rfb, afb, rt, at = setup()
    c0 = load_curve("ps_rat.json")
    rng = random.Random(SEED + 20)
    primes = [int(p) for p in rfb.p.tolist() if p >= 256]
    roots = {int(p): int(r) for p, r in zip(rfb.p.tolist(), rfb.r.tolist())}
    out = open(os.path.join(RAW, "e2_lattice.jsonl"), "w")
    summ = {}
    for t, K in E2_BUCKETS.items():
        agg = {"lattices": 0, "points": 0, "cost": 0, "rel_est": 0.0, "pred_pos": 0.0, "found": 0,
               "sampled": 0, "geo_points_pred": 0.0, "t_sec": 0.0}
        t0 = time.time()
        for _ in range(K):
            M, used = random_modulus(rng, primes, roots, t)
            R, M2 = crt_root([(roots[p], p) for p in used])
            assert M2 == M
            n = lattice_count(R, M, A_BOX, BM_BOX)
            if n <= 2000:
                cnt = Counters()
                pts = list(lattice_points(R, M, A_BOX, BM_BOX, cnt))
            else:
                pts = lattice_sample(R, M, A_BOX, BM_BOX, 2000, rng)
            found, pred = 0, 0.0
            for a, b in pts:
                if gcd(a, b) != 1:
                    continue
                G = P.G(a, b)
                assert G % M == 0
                pred += float(ps_interp(c0, np.array([log2(abs(G) // M)]))[0]) if abs(G) >= M else 1.0
                if rt.ok(G):
                    found += 1
            w = n / len(pts) if pts else 0.0
            agg["lattices"] += 1
            agg["points"] += n
            agg["cost"] += n + 1
            agg["rel_est"] += found * w
            agg["pred_pos"] += pred * w
            agg["found"] += found
            agg["sampled"] += len(pts)
            agg["geo_points_pred"] += U_BOX / M
            if t <= 28 or found:
                out.write(json.dumps({"t": t, "M": M, "primes": used, "n": n, "sampled": len(pts),
                                      "found": found, "pred": pred}) + "\n")
        agg["t_sec"] = time.time() - t0
        agg["Y"] = agg["rel_est"] / agg["cost"]
        agg["Y_pred_pos"] = agg["pred_pos"] / agg["cost"]
        agg["count_ratio"] = agg["points"] / agg["geo_points_pred"]
        summ[t] = agg
        print(t, json.dumps({k: round(v, 5) if isinstance(v, float) else v for k, v in agg.items()}), flush=True)
    json.dump(summ, open(os.path.join(RAW, "e2_summary.json"), "w"), indent=1)


def geo_prediction(t_list=tuple(E2_BUCKETS)) -> Dict[int, float]:
    """P4 prediction from geometry + the E1 box |G| sample only (no E2 data): for M uniform-in-log over
    [2^t, 2^(t+1)), Yhat(t) = E_M[ (U/M) * c6 * E_box[P_s(|G|/M)] / (U/M + 1) ], c6 = coprime fraction."""
    c0 = load_curve("ps_rat.json")
    P = load_cado_poly(POLY)
    gls = []
    for (b, lo, hi) in e1_rows():
        a = np.arange(lo, hi, 64, dtype=np.float64)
        gls.append(np.log2(np.abs(P.Y1 * a + P.Y0 * float(b))))
    gl = np.concatenate(gls)
    c6 = 6 / np.pi ** 2
    out = {}
    for t in t_list:
        vals = []
        for f in np.linspace(0.02, 0.98, 25):
            lm = t + f
            pts = U_BOX / 2 ** lm
            ps = float(np.mean(ps_interp(c0, gl - lm)))
            vals.append(pts * c6 * ps / (pts + 1))
        out[t] = float(np.mean(vals))
    box = float(np.mean(ps_interp(c0, gl))) * c6
    out["box"] = box
    return out


def cmd_e3():
    from src.common.instances import balanced_semiprime
    N = balanced_semiprime(64, 1503).N
    P = base_m_pair(N, 3)
    m = -P.Y0
    A, Bm = 2 ** 9, 2 ** 8
    assert 2 * A < m
    rng = random.Random(SEED + 30)
    ps = primes_upto(2 ** 10).tolist()
    n_gen, hits, t0 = 200000, 0, time.time()
    for _ in range(n_gen):
        while True:
            s = 1
            while s < 2 ** 20:
                s *= rng.choice(ps)
            if s < 2 ** 23:
                break
        if bottom_up_hit(P, s, A, Bm) or bottom_up_hit(P, -s, A, Bm):
            hits += 1
    pred = 2 * A / m
    res = {"N": N, "m": m, "A": A, "Bm": Bm, "generated": n_gen, "hits": hits, "rate": hits / n_gen,
           "pred_rate": pred, "ratio": hits / n_gen / pred, "t_sec": time.time() - t0,
           "c61_pred_rate": 2 * A_BOX / abs(load_cado_poly(POLY).Y0)}
    json.dump(res, open(os.path.join(RAW, "e3_bottomup.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


def cmd_e4ours(K: int, tmax: float):
    """Special-(q1, q2): q1 | G (rational FB prime >= 2^13), q2 | F (algebraic FB ideal >= 2^13), lattice of index
    q1 q2 in the E1 box; enumerate all box points, test both sides exactly; stop at K distinct relations."""
    P, rfb, afb, rt, at = setup()
    rng = random.Random(SEED + 40)
    rq = [(int(r), int(p)) for p, r in zip(rfb.p.tolist(), rfb.r.tolist()) if p >= 2 ** 13]
    aq = [(int(r), int(p)) for p, r in zip(afb.p.tolist(), afb.r.tolist()) if p >= 2 ** 13]
    seen, cnt, t0 = set(), Counters(), time.time()
    tests = 0
    tag = "" if K == 4436 else f"_K{K}"
    log = open(os.path.join(RAW, f"e4_ours{tag}.jsonl"), "w")
    while len(seen) < K and time.time() - t0 < tmax:
        r1, q1 = rng.choice(rq)
        r2, q2 = rng.choice(aq)
        if q1 == q2:
            continue
        R, M = crt_root([(r1, q1), (r2, q2)])
        for a, b in lattice_points(R, M, A_BOX, BM_BOX, cnt):
            if gcd(a, b) != 1:
                continue
            tests += 1
            if rt.ok(P.G(a, b)) and at.ok(P.F(a, b)):
                if (a, b) not in seen:
                    seen.add((a, b))
                    log.write(json.dumps([a, b, q1, q2]) + "\n")
    wall = time.time() - t0
    res = {"K_target": K, "distinct": len(seen), "wall": wall, **cnt.as_dict(), "tests": tests,
           "rels_per_cost": len(seen) / cnt.cost(), "sec_per_rel": wall / max(1, len(seen))}
    json.dump(res, open(os.path.join(RAW, f"e4_ours{tag}_summary.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


def cmd_e4check(relfile: str):
    """Parse las relations 'a,b:...' and re-test them exactly with our side tests (P6a)."""
    P, rfb, afb, rt, at = setup()
    ab = set()
    for line in open(relfile):
        if line.startswith("#") or ":" not in line:
            continue
        a, b = line.split(":")[0].split(",")
        ab.add((int(a), int(b)))
    ok = sum(1 for a, b in ab if rt.ok(P.G(a, b)) and at.ok(P.F(a, b)))
    inbox = sum(1 for a, b in ab if abs(a) < A_BOX and 1 <= b <= BM_BOX)
    res = {"las_unique": len(ab), "pass_ours": ok, "frac": ok / max(1, len(ab)), "in_our_box": inbox}
    json.dump(res, open(os.path.join(RAW, os.path.basename(relfile) + ".check.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":  # pragma: no cover
    c = sys.argv[1]
    if c == "calib":
        cmd_calib()
    elif c == "e1":
        cmd_e1()
    elif c == "e1exact":
        cmd_e1exact()
    elif c == "e2":
        cmd_e2()
    elif c == "geo":
        print(json.dumps(geo_prediction(), indent=1))
    elif c == "e3":
        cmd_e3()
    elif c == "e4ours":
        cmd_e4ours(int(sys.argv[2]), float(sys.argv[3]))
    elif c == "e4check":
        cmd_e4check(sys.argv[2])
