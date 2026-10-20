"""Pre-registered experiments of agent pilist (notes/pilist_prereg.md). Usage:
    PYTHONPATH=. .venv/bin/python -m src.pilist.run_experiment E1|E2|E3
Raw results: experiments/raw/pilist_E{1,2,3}.jsonl; config experiments/configs/pilist.json.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from itertools import combinations
from math import exp, floor, gcd, log, prod
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CFG = json.loads((ROOT / "experiments/configs/pilist.json").read_text())
RAW = ROOT / "experiments/raw"


def _dump(name, rows):
    with open(RAW / f"pilist_{name}.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")


# ------------------------------------------------------------------------------------------- E1
def run_E1():
    from src.pilist.pencil import enumerate_split, find_lines, line_params, line_stats, small_c_triples, tl_class
    rows = []
    for M, B in CFG["E1"]["cells"]:
        t0 = time.time()
        ms, arr = enumerate_split(M, B)
        for kproj in (0, 1):
            if kproj == 0:
                lines = find_lines(arr)
                arr_l = arr
            else:
                # DEVIATION (disclosed): collapse equal projections (same e_1..e_{M-1}) into classes; lines are
                # found among distinct projections and a member triple must use three distinct classes.
                proj = arr[:, : M - 1]
                uniq, inv = np.unique(proj, axis=0, return_inverse=True)
                inv = inv.ravel()
                cls_members = {}
                for i, u in enumerate(inv):
                    cls_members.setdefault(int(u), []).append(i)
                lines = []
                for memu, D in find_lines(uniq):
                    mem = tuple(sorted(i for u in memu for i in cls_members[u]))
                    lines.append((mem, D))
                arr_l = arr
            by = Counter()
            viol_disjoint = viol_nonsqf = 0
            max_m_e3 = 0
            e_full_smallc = 0           # lines with e = M (no fixed part), >= 3 members, a |c| <= 2 triple
            tl = Counter()              # TL classes of small-c triples
            inregime_C = []             # class-C small-c triples inside the Lemma R-type regime
            examples = []
            for mem, D in lines:
                st = line_stats(ms, mem)
                by[(st["e"], min(st["m"], 9))] += 1
                if kproj == 0:
                    if not st["disjoint"]:
                        viol_disjoint += 1
                    if st["e"] >= 1 and st["nonsqf"] > max(2 * st["e"] - 2, 0):
                        viol_nonsqf += 1
                if st["e"] >= 3:
                    max_m_e3 = max(max_m_e3, st["m"])
                ts = line_params(arr_l, mem, D, kproj=kproj)
                nsc = small_c_triples(ts, CFG["E1"]["cmax"])
                if nsc and st["e"] == M:
                    e_full_smallc += 1
                if not nsc:
                    continue
                # classify the small-c triples
                for (i0, t0_), (i1, t1), (i2, t2) in combinations(list(zip(mem, ts)), 3):
                    if len({t0_, t1, t2}) < 3:
                        continue
                    v = (t1 - t2, t2 - t0_, t0_ - t1)
                    g = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
                    c = tuple(x // g for x in v)
                    if max(abs(x) for x in c) > CFG["E1"]["cmax"]:
                        continue
                    # G = sum c_i P_i: zero (kproj = 0) or a constant (kproj = 1)
                    Gvec = c[0] * arr[i0] + c[1] * arr[i1] + c[2] * arr[i2]
                    Gc = int(Gvec[-1])
                    assert all(int(x) == 0 for x in Gvec[: M - kproj])
                    h = abs(Gc)
                    L = int(floor(log(h, 2))) if h >= 1 else -1
                    cls = tl_class(ms, i0, i1, i2, max(L, 0))
                    tl[cls] += 1
                    # in regime (Lemma R shape with K = 1 for a constant G, j = 0): M >= 1 + 2L + 3; G = 0 needs e >= 6
                    inreg = (h >= 1 and M >= 1 + 2 * L + 3) or (h == 0 and st["e"] >= 6)
                    if cls == "C" and inreg:
                        inregime_C.append({"roots": [ms[i0], ms[i1], ms[i2]], "c": c, "G": Gc})
                    if len(examples) < 3 and cls == "C":
                        examples.append({"roots": [ms[i0], ms[i1], ms[i2]], "c": c, "G": Gc, "e": st["e"]})
            rows.append({"M": M, "B": B, "kproj": kproj, "n_split": len(ms), "n_lines": len(lines),
                         "by_e_m": {f"{e},{m}": v for (e, m), v in sorted(by.items())},
                         "viol_disjoint": viol_disjoint, "viol_nonsqf": viol_nonsqf, "max_members_e3": max_m_e3,
                         "e_full_smallc_lines": e_full_smallc, "tl_classes": dict(tl),
                         "inregime_C": inregime_C[:20], "n_inregime_C": len(inregime_C), "examples_C": examples,
                         "secs": round(time.time() - t0, 2)})
            print(json.dumps({k: rows[-1][k] for k in ("M", "B", "kproj", "n_lines", "viol_disjoint", "viol_nonsqf",
                                                        "max_members_e3", "e_full_smallc_lines", "tl_classes",
                                                        "n_inregime_C", "secs")}), flush=True)
    _dump("E1", rows)


# ------------------------------------------------------------------------------------------- E2
def run_E2():
    from src.pilist.decoder import (boneh_threshold, count_model_fraction, hg_row_decoder, lam_of, log_P,
                                    primes_upto, rankin_log_bound, sieve_row_smooth)
    c = CFG["E2"]
    meta = json.loads((ROOT / c["testbed"] / "meta.json").read_text())
    Y0, Y1 = int(meta["f0_rational_low_to_high"][0]), int(meta["f0_rational_low_to_high"][1])
    Bfb = c["B_smooth"]
    fb = primes_upto(Bfb)
    rows_b = [b for b in range(1, 100000) if gcd(b, Y1) == 1]
    out = []
    for X, nrows in c["X_rows"]:
        bs = rows_b[:nrows]
        t0 = time.time()
        smooth = {b: set(sieve_row_smooth(Y1, Y0, b, X, fb)) for b in bs}
        npts = sum(sum(1 for a in range(-X, X + 1) if gcd(a, b) == 1) for b in bs)
        nsm = sum(len(s) for s in smooth.values())
        ybox = nsm / npts
        print(f"X={X}: rows={len(bs)} box smooth rate {ybox:.5f} ({nsm}/{npts}) {time.time() - t0:.1f}s", flush=True)
        for z in c["z"]:
            Pd = {q: 1 for q in primes_upto(z) if Y1 % q}
            P = prod(Pd)
            lV = log(2 * X + 1)
            thr = boneh_threshold(lV, Pd)
            if thr > log_P(Pd):
                out.append({"X": X, "z": z, "skipped": "threshold above log P"})
                continue
            lo, hi = count_model_fraction(Pd, thr)
            pred = hi * (2 * X + 1)
            truth_tot = smooth_out = 0
            dec_rows = dec_miss = 0
            cert = []
            for ib, b in enumerate(bs):
                cc = (Y0 * b * pow(Y1, -1, P)) % P
                arr_a = np.arange(-X, X + 1)
                # exact d(a) = gcd(G(a,b), P) only where the log-sieve says it is large
                lg = np.zeros(2 * X + 1)
                for q in Pd:
                    r = (-cc) % q
                    lg[(r + X) % q::q] += log(q)
                cand = np.nonzero(lg >= thr - 1e-9)[0]
                tr = [int(arr_a[i]) for i in cand if log(gcd(Y1 * int(arr_a[i]) + Y0 * b, P)) >= thr - 1e-12]
                truth_tot += len(tr)
                smooth_out += sum(1 for a in tr if a in smooth[b])
                if ib < c["decoder_rows"]:
                    # smallest m in the ladder whose certified threshold is within 0.5 nats of thr
                    for m in c["m_ladder"]:
                        roots, norms = hg_row_decoder(cc, P, X, m, 2 * m)
                        sc = log(norms[0]) / m
                        if sc <= thr + 0.5:
                            break
                    cert.append(sc - thr)
                    dec_rows += 1
                    need = [int(arr_a[i]) for i in np.nonzero(lg >= sc - 1e-9)[0]
                            if log(gcd(Y1 * int(arr_a[i]) + Y0 * b, P)) > sc + 1e-12]
                    dec_miss += sum(1 for a in need if a not in roots)
            u = log(abs(Y0) * bs[-1] + abs(Y1) * X) / log(Bfb)
            lam = lam_of(u, log(Bfb))
            jr = exp(rankin_log_bound(lV, Pd, thr, lam)[0])
            per_call = smooth_out / len(bs)
            out.append({"X": X, "z": z, "rows": len(bs), "logP": log_P(Pd), "thr": thr,
                        "pred_out_per_call": pred, "pred_lo": lo * (2 * X + 1), "obs_out_per_call": truth_tot / len(bs),
                        "truth_total": truth_tot, "pred_total": pred * len(bs),
                        "smooth_outputs": smooth_out, "ybox": ybox, "gamma_hat": per_call / ybox,
                        "jr_bound": jr, "decoder_rows": dec_rows, "decoder_misses": dec_miss,
                        "cert_minus_thr_max": max(cert) if cert else None})
            print(json.dumps(out[-1]), flush=True)
    _dump("E2", out)


# ------------------------------------------------------------------------------------------- E3
def run_E3():
    from src.pilist.tails import r_threshold, simulate
    c = CFG["E3"]
    out = []
    for n in c["n"]:
        for k, j in c["kj"]:
            Rlb = r_threshold(k, j, 1 / n, target=c["target"])
            strategies = ["nonadaptive", "greedy", "oracle"] if j > 0 else ["nonadaptive"]
            for strat in strategies:
                rs = []
                for seed in range(c["seeds"]):
                    jj = j if strat != "nonadaptive" else 0
                    r = simulate(n, k, jj, strat, seed)
                    rs.append(r["R_star"])
                rr = sorted(x for x in rs if x is not None)
                med = rr[len(rr) // 2] if rr else None
                below = sum(1 for x in rs if x is not None and x <= Rlb)
                out.append({"n": n, "k": k, "j": j, "strategy": strat, "R_lb": Rlb, "R_star": rs,
                            "median": med, "median_over_n": med / n if med else None, "n_below_Rlb": below})
                print(json.dumps({kk: out[-1][kk] for kk in ("n", "k", "j", "strategy", "R_lb", "median",
                                                              "median_over_n", "n_below_Rlb")}), flush=True)
    _dump("E3", out)


if __name__ == "__main__":
    {"E1": run_E1, "E2": run_E2, "E3": run_E3}[sys.argv[1]]()
