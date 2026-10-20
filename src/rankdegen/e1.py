"""Experiment E1 (pre-registered, notes/rankdegen_prereg.md): Theorem FH and the linking exponent, exact, r = 3.

Harness code: N is generated here (primality tests used by the rank3 generator); the checked quantities use only the
module, N and its prime list (for labelling). No factoring method is implemented.
Run: PYTHONPATH=. .venv/bin/python -m src.rankdegen.e1   (writes experiments/raw/rankdegen_E1.jsonl)
"""
from __future__ import annotations

import bisect
import json
import math
import random
import time
from fractions import Fraction as Fr

from flint import arb, ctx

from src.rank3.e1 import inst_A
from src.rankdegen import fh

OUT = "experiments/raw/rankdegen_E1.jsonl"
SEED = 20261020
E = 30
R = 3


def crt(residues, moduli):
    x, M = 0, 1
    for a, m in zip(residues, moduli):
        t = ((a - x) * pow(M, -1, m)) % m
        x += M * t
        M *= m
    return x % M


def instances():
    out = []
    roots = {31: 2, 41: 13, 47: 16}
    for ps in ([31], [31, 41], [31, 41, 47]):
        N = math.prod(ps)
        f = [-2, 0, 0, 0, 0, 0, 1]
        out.append({"kind": "C2", "N": N, "primes": ps, "f": f, "g": f,
                    "P": [[1], [0, 1], [0, 0, 0, 1]], "Q": [[1], [N, 1], [0, 0, 0, 1]],
                    "mu": crt([roots[p] for p in ps], ps)})
    N = 2759
    f = [1, 0, 1, 0, 0, 0, 1]
    out.append({"kind": "B", "N": N, "primes": [31, 89], "f": f, "g": f,
                "P": [[1], [0, 0, 1], [0, 0, 0, 0, 1]], "Q": [[1], [N, 0, 1], [0, 0, 0, 0, 1]], "mu": 14})
    phi5 = [1, 1, 1, 1, 1]
    out.append({"kind": "A", "N": 31 * 41, "primes": [31, 41], "f": phi5, "g": phi5,
                "P": [[1], [0, 1], [0, 0, 1]], "Q": [[1], [0, 1], [0, 0, 1]], "mu": None})
    for s in range(3):
        rnd = random.Random(f"{SEED}-rankdegen-E1-{s}")
        I = inst_A(rnd, 3, 3)
        N, m = I["N"], I["m"]
        mm = m - N if m > N // 2 else m
        base = {"N": N, "primes": I["primes"], "f": I["f"], "g": I["g"], "mu": m, "seed": s}
        out.append(dict(base, kind="GP", P=[[1], [0, 1], [0, 0, 1]], Q=[[1], [0, 1], [0, 0, 1]]))
        out.append(dict(base, kind="C1", P=[[1], [0, 1], [0]], Q=[[1], [0, 1], [0, -mm, 1]]))
        out.append(dict(base, kind="FP", P=[[1], [0, 1], [0]], Q=[[1], [0, 1], [-mm, 1]]))
    return out


def run_one(I):
    f, g, P, Q, N = I["f"], I["g"], I["P"], I["Q"], I["N"]
    T = fh.minor_table(f, P, g, Q, R)
    d, D = T["d"], T["D"]
    fr = fh.frames(T)
    dists = []
    for s in range(R + 2):
        w = fh.pure_s(T, fr, s)
        if sum(w) > 0:
            dists.append((f"pure{s}", w, True))
    lpval, lpw = fh.maxmin_coverage(T, fr, [Fr(R + 1, D)] * D) if fr else (None, None)
    if lpval is not None:
        dists.append(("LPuniform", lpw, False))
    pairs = fh.box_norm_pairs(f, P, g, Q, R, E)
    err = R * D * (2 * E + 1) ** (R - 1)
    # 2D grid of thresholds
    k1max = max(a.bit_length() for a, _ in pairs)
    k2max = max(b.bit_length() for _, b in pairs)
    by1 = sorted(pairs)
    recs = []
    ctx.prec = 256
    logN = arb(N).log()
    viol = checks = informative = 0
    p2viol = p2checks = 0
    for name, w, invariant in dists:
        st = fh.distribution_stats(T, fr, w)
        yK = Fr(R, R + 1) * st["marg"][0]
        yK2 = Fr(R, R + 1) * st["marg"][d]
        t = st["mincov"] if name != "LPuniform" else lpval
        logJ_N = arb(2).log() * R + fh.bean_thunder_C(R).log() - arb(R - 1) / (R + 1) * (arb(t.numerator) / t.denominator) * logN
        bounds = [("Nonly", logJ_N)]
        p2 = None
        if invariant:
            logJ_m = fh.fh_log_bound(T, fr, w)
            bounds.append(("minors", logJ_m))
            if t > 0:
                lhs = arb(0)
                for S, ws in zip(fr, w):
                    if ws:
                        lhs += arb(ws.numerator) / ws.denominator * fh.log_phi(T, S)
                rhs = arb(R - 1) * (arb(t.numerator) / t.denominator) * logN
                p2 = bool(lhs.lower() >= rhs.upper())
                p2checks += 1
                p2viol += not p2
        cnt_ok = True
        ninf = 0
        nchk = 0
        for k1 in range(0, k1max + 1):
            Y1 = 1 << k1
            sub = [b for a, b in by1[:bisect.bisect_right(by1, (Y1, 1 << 4000))]]
            sub.sort()
            for k2 in range(0, k2max + 1):
                Y2 = 1 << k2
                c = bisect.bisect_right(sub, Y2)
                for bname, lj in bounds:
                    main = (lj + arb(yK.numerator) / yK.denominator * arb(Y1).log()
                            + arb(yK2.numerator) / yK2.denominator * arb(Y2).log()).exp()
                    ok = bool(c <= main.lower() + err)  # certified: compare with the lower end of the ball
                    nchk += 1
                    cnt_ok &= ok
                    viol += not ok
                    if c > err:
                        ninf += 1
        checks += nchk
        informative += ninf
        recs.append({"dist": name, "yK": str(yK), "yK2": str(yK2), "mincov": str(t),
                     "logJ_Nonly": float(logJ_N.mid()), "logJ_minors": float(bounds[1][1].mid()) if len(bounds) > 1 else None,
                     "P2": p2, "count_ok": cnt_ok, "checks": nchk, "informative": ninf})
    return {"kind": I["kind"], "seed": I.get("seed"), "N": N, "primes": I["primes"], "f": f, "g": g, "P": P, "Q": Q,
            "D": D, "zero_minors": sum(1 for v in T["status"].values() if v == 0), "n_minors": len(T["status"]),
            "frames": len(fr), "LP_value": str(lpval) if lpval is not None else None, "box_points": len(pairs),
            "err": err, "dists": recs, "viol": viol, "checks": checks, "informative": informative,
            "P2_viol": p2viol, "P2_checks": p2checks}


def main():
    tot = {"viol": 0, "checks": 0, "inf": 0, "p2v": 0, "p2c": 0}
    with open(OUT, "w") as fhd:
        for I in instances():
            t0 = time.time()
            res = run_one(I)
            res["secs"] = round(time.time() - t0, 1)
            fhd.write(json.dumps(res) + "\n")
            fhd.flush()
            tot["viol"] += res["viol"]; tot["checks"] += res["checks"]; tot["inf"] += res["informative"]
            tot["p2v"] += res["P2_viol"]; tot["p2c"] += res["P2_checks"]
            print(res["kind"], res["seed"], "N", res["N"], "zeros", res["zero_minors"], "/", res["n_minors"],
                  "frames", res["frames"], "LP", res["LP_value"], "viol", res["viol"], "/", res["checks"],
                  "inf", res["informative"], "P2", res["P2_viol"], "/", res["P2_checks"],
                  [(d_["dist"], d_["mincov"]) for d_ in res["dists"]], res["secs"], "s", flush=True)
    print("TOTAL", tot)


if __name__ == "__main__":
    main()
