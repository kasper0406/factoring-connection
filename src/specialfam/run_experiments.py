"""HARNESS driver for the pre-registered experiments E1-E5 (notes/specialfam_prereg.md).

Usage:  .venv/bin/python -m src.specialfam.run_experiments e1|e2c|e3|e4|e5|summarize
Raw JSONL -> experiments/raw/specialfam_<exp>.jsonl; summary -> experiments/analysis/specialfam_summary.json
Operational algorithms receive only N and the public family parameters; the harness
then checks the output against the hidden factors and classifies the hit.
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

from src.common.metrics import Metrics, append_jsonl, read_jsonl
from src.specialfam import harness as H
from src.specialfam.normbox import normbox_factor, disc_box_factor, near_square_factor
from src.specialfam.smoothcopp import smoothcopp_factor, common_factor_factor, balanced_band
from src.specialfam.copper import best_params

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "experiments" / "raw"
CONFIG = json.loads((ROOT / "experiments" / "configs" / "specialfam.json").read_text())


def out(name):
    p = RAW / f"specialfam_{name}.jsonl"
    if p.exists():
        p.unlink()
    return p


def e1():
    """E1/E2: Algorithm A on the HL-type family; correctness, classification, cost."""
    c = CONFIG["e1"]
    eta = Fraction(c["eta"])
    path = out("e1")
    for beta_s in c["betas"]:
        beta = Fraction(beta_s)
        for n in c["sizes"]:
            for seed in range(c["per_cell"]):
                inst = H.normbox_instance(n, beta, eta, seed)
                N = inst["N"]
                kmin = max(4, math.floor(beta * N.bit_length()) - c["kmin_back"])
                mt = Metrics()
                t0 = time.perf_counter()
                res = normbox_factor(N, eta, kmin=kmin, metrics=mt)
                wall = time.perf_counter() - t0
                f = res.get("factor")
                ok = bool(f) and N % f == 0 and 1 < f < N
                rec = {"exp": "e1", "beta": beta_s, "n": n, "nbits": N.bit_length(), "seed": seed,
                       "kp": inst["kp"], "b": inst["b"], "ok": ok, "band": res.get("band"),
                       "how": res.get("how"), "class": H.classify_normbox(inst, res, eta),
                       "units_band": res.get("units"), "units_total": res.get("units_total"),
                       "m": res.get("m"), "K": res.get("K"), "wall_s": round(wall, 4),
                       "counts": mt.counts}
                append_jsonl(path, rec)
            print("e1", beta_s, n, flush=True)


def e2c():
    """E2-control: unstructured balanced semiprimes, same algorithm, bands up to n/2."""
    c = CONFIG["e2c"]
    eta = Fraction(c["eta"])
    path = out("e2c")
    for n in c["sizes"]:
        for seed in range(c["per_cell"]):
            inst = H.control_instance(n, seed)
            N = inst["N"]
            kp = inst["p"].bit_length()
            res = normbox_factor(N, eta, kmin=max(4, N.bit_length() // 2 - 3), kmax=N.bit_length() // 2 + 1)
            f = res.get("factor")
            ok = bool(f) and N % f == 0 and 1 < f < N
            wrong = bool(f) and not ok
            rep = H.sumsq_rep(inst["p"])
            promised = bool(rep) and rep[1] ** 1000 <= inst["p"] ** 119
            append_jsonl(path, {"exp": "e2c", "n": n, "seed": seed, "found": ok, "wrong": wrong,
                                "how": res.get("how"), "band": res.get("band"),
                                "p_has_small_b": promised, "rep": rep})
        print("e2c", n, flush=True)


def e3():
    """E3: other Theorem-G families (HKL-type box with public centre; near-square)."""
    c = CONFIG["e3"]
    path = out("e3")
    import random
    for n in c["sizes"]:
        for seed in range(c["per_cell"]):
            rng = random.Random(f"specialfam|e3|{n}|{seed}")
            kp = n // 2
            # HKL-type: centre a0 = b0 = isqrt(2^{kp-1}/2)... public; side L = 2^{ceil(0.265 kp)}
            L = 1 << math.ceil(0.265 * kp)
            a0 = b0 = math.isqrt(1 << (kp - 2))
            while True:
                a, b = a0 + rng.randrange(L), b0 + rng.randrange(L)
                p = a * a + b * b
                if H.is_prime(p):
                    break
            q = H.rand_prime_in(rng, p + 1, 2 * p - 1)
            res = disc_box_factor(p * q, a0, b0, L)
            append_jsonl(path, {"exp": "e3", "fam": "hkl_box", "n": n, "seed": seed,
                                "ok": res.get("factor") in (p, q), "units": res.get("units")})
            # near square: p = a^2 + k, k <= 2^{0.238 kp}
            Kmax = 1 << math.ceil(0.238 * kp)
            Amax = math.isqrt(1 << kp)
            while True:
                a = rng.randrange(math.isqrt(1 << (kp - 1)), Amax)
                k = rng.randint(1, Kmax)
                p = a * a + k
                if H.is_prime(p):
                    break
            q = H.rand_prime_in(rng, p + 1, 2 * p - 1)
            res = near_square_factor(p * q, Amax, Kmax)
            append_jsonl(path, {"exp": "e3", "fam": "near_square", "n": n, "seed": seed,
                                "ok": res.get("factor") in (p, q), "units": res.get("units")})
        print("e3", n, flush=True)


def _smooth_gfac(gbits: int):
    """g = 2^a * 3 * 5 with about gbits bits (P+(g) = 5; 8 square roots of 1 mod g)."""
    a = max(1, gbits - 4)  # 15 ~ 2^3.9
    return {2: a, 3: 1, 5: 1}


def e4():
    """E4: Algorithm B, smooth g | gcd(p^2-1, q^2-1), mixed CRT signs."""
    c = CONFIG["e4"]
    path = out("e4")
    import random
    for gam in c["gammas"]:
        for n in c["sizes"]:
            gbits = round(gam * n)
            for seed in range(c["per_cell"]):
                rng = random.Random(f"specialfam|e4|{gam}|{n}|{seed}")
                fac = _smooth_gfac(gbits)
                sp = "".join(rng.choice("+-") for _ in fac)
                sq = "".join(rng.choice("+-") for _ in fac)
                inst = H.smooth_gcd_instance(n, fac, sp, sq, seed)
                N = inst["N"]
                Dmin = 1 << (inst["g"].bit_length() - 1)  # public promise g >= Dmin
                mt = Metrics()
                t0 = time.perf_counter()
                res = smoothcopp_factor(N, c["B"], Dmin, balanced_band(N), Dmax_lat=c["Dmax_lat"], metrics=mt)
                wall = time.perf_counter() - t0
                D, m, Y = best_params(N, balanced_band(N)[0], c["Dmax_lat"])
                lo, hi = balanced_band(N)
                pred = -(-((hi - lo) // Dmin + 1) // (2 * Y + 1))
                append_jsonl(path, {"exp": "e4", "gamma": gam, "n": n, "nbits": N.bit_length(), "seed": seed,
                                    "g": inst["g"], "signs": [sp, sq], "ok": res.get("factor") in (inst["p"], inst["q"]),
                                    "chunks": res.get("chunks"), "pairs": res.get("pairs"),
                                    "n_divisors": res.get("n_divisors"), "pred_chunks_per_pair": pred,
                                    "lattice": [D, m], "Y": Y, "wall_s": round(wall, 4),
                                    "ps_units": mt.counts.get("ps_units", 0)})
            print("e4", gam, n, flush=True)


def e5():
    """E5: Algorithm D, arbitrary (non-smooth) g = gcd(p-1, q-1) ~ N^gamma; WM comparison."""
    c = CONFIG["e5"]
    path = out("e5")
    for gam in c["gammas"]:
        for n in c["sizes"]:
            gbits = max(4, round(Fraction(gam) * n))
            for seed in range(c["per_cell"]):
                inst = H.common_gcd_instance(n, gbits, seed)
                N = inst["N"]
                mt = Metrics()
                t0 = time.perf_counter()
                res = common_factor_factor(N, balanced_band(N), Dmax_lat=c["Dmax_lat"], metrics=mt)
                wall = time.perf_counter() - t0
                append_jsonl(path, {"exp": "e5", "gamma": gam, "n": n, "nbits": N.bit_length(), "seed": seed,
                                    "g": inst["g"], "gcd_pm1": inst["gcd_pm1"],
                                    "ok": res.get("factor") in (inst["p"], inst["q"]),
                                    "d": res.get("d"), "t": res.get("t"), "chunks": res.get("chunks"),
                                    "ps_units": res.get("ps_units"), "lll_calls": mt.counts.get("lll_calls", 0),
                                    "wall_s": round(wall, 4)})
            print("e5", gam, n, flush=True)


def _fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def summarize():
    S = {}
    def med(v):
        v = sorted(v)
        return v[len(v) // 2]
    if (RAW / "specialfam_e1.jsonl").exists():
        R = read_jsonl(RAW / "specialfam_e1.jsonl")
        for beta in sorted({r["beta"] for r in R}):
            rs = [r for r in R if r["beta"] == beta]
            ns = sorted({r["n"] for r in rs})
            ub = [math.log2(med([r["units_band"] for r in rs if r["n"] == n])) for n in ns]
            ut = [math.log2(med([r["units_total"] for r in rs if r["n"] == n])) for n in ns]
            wl = [math.log2(med([r["wall_s"] for r in rs if r["n"] == n])) for n in ns]
            half = [i for i, n in enumerate(ns) if n >= ns[len(ns) // 2]]
            S[f"e1_beta{beta}"] = {
                "instances": len(rs), "ok": sum(r["ok"] for r in rs),
                "classes": {c: sum(r["class"] == c for r in rs) for c in {r["class"] for r in rs}},
                "slope_units_band": round(_fit(ns, ub), 4), "slope_units_total": round(_fit(ns, ut), 4),
                "slope_wall": round(_fit(ns, wl), 4),
                "slope_wall_upper_half": round(_fit([ns[i] for i in half], [wl[i] for i in half]), 4),
                "pred": float(Fraction(beta) * (Fraction(1, 4) + Fraction(CONFIG["e1"]["eta"]) / 2)),
                "full_gcd_blocks": sum(r["counts"].get("full_gcd_blocks", 0) for r in rs)}
    if (RAW / "specialfam_e2c.jsonl").exists():
        R = read_jsonl(RAW / "specialfam_e2c.jsonl")
        S["e2c"] = {"instances": len(R), "found": sum(r["found"] for r in R), "wrong": sum(r["wrong"] for r in R),
                    "p_has_small_b": sum(r["p_has_small_b"] for r in R)}
    if (RAW / "specialfam_e3.jsonl").exists():
        R = read_jsonl(RAW / "specialfam_e3.jsonl")
        S["e3"] = {f: {"instances": sum(r["fam"] == f for r in R), "ok": sum(r["ok"] for r in R if r["fam"] == f)}
                   for f in {r["fam"] for r in R}}
    if (RAW / "specialfam_e4.jsonl").exists():
        R = read_jsonl(RAW / "specialfam_e4.jsonl")
        for gam in sorted({r["gamma"] for r in R}):
            rs = [r for r in R if r["gamma"] == gam]
            ns = sorted({r["n"] for r in rs})
            ch = [math.log2(max(1, med([r["chunks"] for r in rs if r["n"] == n]))) for n in ns]
            pc = [math.log2(max(1, med([r["pred_chunks_per_pair"] for r in rs if r["n"] == n]))) for n in ns]
            S[f"e4_gamma{gam}"] = {"instances": len(rs), "ok": sum(r["ok"] for r in rs),
                                   "slope_chunks": round(_fit(ns, ch), 4), "slope_pred_chunks_per_pair": round(_fit(ns, pc), 4),
                                   "max_pairs": max(r["pairs"] for r in rs),
                                   "asymptotic": 0.25 - gam}
    if (RAW / "specialfam_e5.jsonl").exists():
        R = read_jsonl(RAW / "specialfam_e5.jsonl")
        for gam in sorted({r["gamma"] for r in R}):
            rs = [r for r in R if r["gamma"] == gam]
            ns = sorted({r["n"] for r in rs})
            ps = [math.log2(max(1, med([r["ps_units"] for r in rs if r["n"] == n]))) for n in ns]
            ch = [math.log2(max(1, med([r["chunks"] for r in rs if r["n"] == n]))) for n in ns]
            # P7 statistic (pre-registered; added to the summarizer after the freeze, analysis only)
            mx = [math.log2(max(1, med([max(r["ps_units"], r["chunks"]) for r in rs if r["n"] == n]))) for n in ns]
            g = float(Fraction(gam))
            wm = (0.25 - g if g <= 1 / 12 else 1 / 8 + g / 2 if g <= 1 / 8 else 0.25 - g / 2 if g <= 1 / 6
                  else 0.5 - 2 * g if g <= 0.2 else g / 2 if g <= 0.25 else 0.25 - g / 2)
            S[f"e5_gamma{gam}"] = {"instances": len(rs), "ok": sum(r["ok"] for r in rs),
                                   "slope_ps_units": round(_fit(ns, ps), 4), "slope_chunks": round(_fit(ns, ch), 4),
                                   "slope_max_ps_chunks": round(_fit(ns, mx), 4),
                                   "pred_ps": g / 2, "pred_chunks_asym": max(0.0, 0.25 - g),
                                   "ours_asym": max(g / 2, 0.25 - g), "WM_E": wm}
    p = ROOT / "experiments" / "analysis" / "specialfam_summary.json"
    p.write_text(json.dumps(S, indent=1, sort_keys=True))
    print(json.dumps(S, indent=1, sort_keys=True))


if __name__ == "__main__":
    {"e1": e1, "e2c": e2c, "e3": e3, "e4": e4, "e5": e5, "summarize": summarize}[sys.argv[1]]()
