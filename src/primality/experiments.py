"""HARNESS experiments for the primality-leakage study.

Subcommands (run from the project root):
  .venv/bin/python -m src.primality.experiments leakrates   # E1 leak rates 32-128 bits + diagnostic families
  .venv/bin/python -m src.primality.experiments kangaroo    # E2 interval-DLP factoring from residues, cost vs n
  .venv/bin/python -m src.primality.experiments iterate     # E3 iterated / combined residues
  .venv/bin/python -m src.primality.experiments epscert     # E4a frequency of one-sided eps certificates vs n
  .venv/bin/python -m src.primality.experiments epsexact    # E4b exact conditional transcript laws (small n)
  .venv/bin/python -m src.primality.experiments battery     # E4c pre-registered transcript/eps battery
  .venv/bin/python -m src.primality.experiments aks         # E5 AKS leak census
  .venv/bin/python -m src.primality.experiments menu        # E6 exponent x group menu
Configs: experiments/configs/primality_*.json.  Raw: experiments/raw/primality_*.jsonl.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
import time
from fractions import Fraction as Fr
from math import gcd, isqrt
from multiprocessing import Pool

import gmpy2
import sympy

from src.common.instances import balanced_semiprime, carmichael_chernick, _rng, random_prime, is_prime
from src.common.metrics import Metrics, append_jsonl
from src.primality import dlog, ops
from src.primality import formulas as F

RAW = "experiments/raw/"


def log2f(x):
    x = Fr(x)
    if x == 0:
        return None
    return math.log2(x.numerator) - math.log2(x.denominator)


# ------------------------------------------------------------------ E1
def diag_shared_gcd(bits, seed, gap=10):
    """DIAGNOSTIC family: p = 2Gk+1, q = 2Gl+1 with a shared prime G of about bits/2 - gap bits,
    so g = gcd(p-1, q-1) >= 2G and the leak rates are ~2^-gap (measurable)."""
    rng = _rng("diag_shared_gcd", {"bits": bits, "gap": gap}, seed)
    hb = bits // 2
    while True:
        G = random_prime(rng, hb - gap)
        ps = []
        for _ in range(200):
            k = rng.getrandbits(gap - 1) | (1 << (gap - 2))
            c = 2 * G * k + 1
            if c.bit_length() == hb and is_prime(c):
                ps.append(c)
                if len(ps) == 2 and ps[0] != ps[1]:
                    break
        if len(ps) >= 2 and ps[0] != ps[1]:
            p, q = sorted(ps[:2])
            return p, q, G


def e1_instance(args):
    fam, bits, seed, nb = args
    if fam == "balanced":
        inst = balanced_semiprime(bits, seed)
        p, q = inst.labels["factors"]
        N = inst.N
    elif fam == "shared":
        p, q, G = diag_shared_gcd(bits, seed)
        N = p * q
    rng = random.Random(f"e1|{fam}|{bits}|{seed}")
    cm = F.closed_mr_semiprime(p, q)
    ss = F.closed_ss_semiprime(p, q)
    sel = ops.selfridge(N)
    rec = {"fam": fam, "bits": bits, "seed": seed, "N": N,
           "g": gcd(p - 1, q - 1), "u": F.v2(p - 1), "v": F.v2(q - 1),
           "pred": {"fermat": float(cm["single"]), "mr_prefix": float(cm["prefix"]), "mr_sqrt": float(cm["sqrt"]),
                    "ss": float(ss)}}
    if "P" in sel:
        D = sel["D"]
        dp, dq = int(sympy.legendre_symbol(D % p, p)), int(sympy.legendre_symbol(D % q, q))
        lu = F.closed_ladder(p - dp, q - dq, N + 1)
        rec["pred"]["lucas_prefix_uniform"] = float(lu["prefix"])
        rec["lucas_G"] = [gcd(N + 1, p - dp), gcd(N + 1, q - dq)]
    obs = {"fermat": 0, "mr_prefix": 0, "mr_pm": 0, "mr_sqrt": 0, "ss": 0, "lucas_pm": 0, "nonunit": 0}
    for _ in range(nb):
        a = rng.randrange(2, N - 1)
        if gcd(a, N) != 1:
            obs["nonunit"] += 1
            continue
        if ops.fermat(N, a)["leaks"]:
            obs["fermat"] += 1
        mr = ops.miller_rabin(N, a)
        st = [l for l in mr["leaks"] if l["class"] == "structural"]
        obs["mr_prefix"] += any(l["kind"] == "minus1" for l in st)
        obs["mr_pm"] += bool(st)
        obs["mr_sqrt"] += mr["sqrt_split"] is not None and 1 < mr["sqrt_split"] < N
        obs["ss"] += any(l["class"] == "structural" for l in ops.solovay_strassen(N, a)["leaks"])
        if "P" in sel:
            P = rng.randrange(N)
            Q = (P * P - sel["D"]) * pow(4, -1, N) % N
            if gcd(Q, N) == 1:
                obs["lucas_pm"] += any(l["class"] == "structural" for l in ops.strong_lucas(N, P, Q)["leaks"])
    rec["obs"] = obs
    rec["bases"] = nb
    b = ops.bpsw(N)
    rec["bpsw_leak"] = [l for l in b["leaks"]]
    return rec


def run_leakrates(a):
    cfg = {"balanced": {32: (20000, 5), 48: (5000, 5), 64: (5000, 5), 96: (2000, 3), 128: (2000, 3)},
           "shared": {64: (200, 200), 96: (100, 200), 128: (100, 200)}, "chernick_kbits": [10, 14, 18],
           "chernick_seeds": 20, "chernick_bases": 200}
    json.dump(cfg, open("experiments/configs/primality_leakrates.json", "w"), indent=1, default=str)
    out = RAW + "primality_leakrates.jsonl"
    open(out, "w").close()
    jobs = []
    for fam in ("balanced", "shared"):
        for bits, (cnt, nb) in cfg[fam].items():
            jobs += [(fam, bits, s, nb) for s in range(cnt)]
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(e1_instance, jobs, chunksize=16):
            append_jsonl(out, rec)
    # Carmichael (Chernick): empirical MR factor rate vs exact formula
    for kb in cfg["chernick_kbits"]:
        for sd in range(cfg["chernick_seeds"]):
            inst = carmichael_chernick(kb, sd)
            N = inst.N
            ps = inst.labels["factors"]
            pred = F.mr_exact(ps)["prefix"]
            rng = random.Random(f"carm|{kb}|{sd}")
            hit = tot = 0
            for _ in range(cfg["chernick_bases"]):
                x = rng.randrange(2, N - 1)
                if gcd(x, N) != 1:
                    continue
                tot += 1
                hit += any(l["kind"] == "minus1" and l["class"] == "structural"
                           for l in ops.miller_rabin(N, x)["leaks"])
            append_jsonl(out, {"fam": "chernick", "kbits": kb, "seed": sd, "N": N, "bits": N.bit_length(),
                               "pred_mr_prefix": float(pred), "obs": hit, "bases": tot})


# ------------------------------------------------------------------ E2
def e2_job(args):
    kind, bits, seed, method = args
    inst = balanced_semiprime(bits, seed)
    m = Metrics()
    t0 = time.perf_counter()
    if kind == "fermat":
        r = dlog.factor_via_fermat_residue(inst.N, seed, method, m)
    else:
        r = dlog.factor_via_lucas_residue(inst.N, seed, method, m)
    ok = r.get("factors") is not None and tuple(sorted(r["factors"])) == tuple(inst.labels["factors"])
    p, q = inst.labels["factors"]
    return {"kind": kind, "method": method, "bits": bits, "seed": seed, "ok": ok, "class": r.get("class"),
            "attempts": r.get("attempts"), "grp_mul": m.counts.get("grp_mul", 0), "stored": m.maxima.get("stored", 0),
            "wall_s": time.perf_counter() - t0, "N4": inst.N ** 0.25, "s_minus_2sqrtN": (p + q) - 2 * math.isqrt(inst.N)}


def run_kangaroo(a):
    cfg = {"fermat_kangaroo": list(range(16, 69, 4)), "fermat_bsgs": list(range(16, 61, 4)),
           "lucas_kangaroo": list(range(16, 57, 8)), "per_size": 30, "seed0": 0}
    json.dump(cfg, open("experiments/configs/primality_kangaroo.json", "w"), indent=1)
    jobs = []
    for key, sizes in cfg.items():
        if not isinstance(sizes, list):
            continue
        kind, method = key.split("_")
        jobs += [(kind, b, s, method) for b in sizes for s in range(cfg["per_size"])]
    out = RAW + "primality_kangaroo.jsonl"
    open(out, "w").close()
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(e2_job, jobs):
            append_jsonl(out, rec)


# ------------------------------------------------------------------ E3
def e3_job(args):
    """Iterated residue y2 = y^(N-1) = a^((s-2)^2).  (i) kangaroo over the exponent interval of (s-2)^2;
    (ii) scan over s with finite differences (a^((x+1)^2) = a^(x^2) a^(2x+1)).  Compare with plain y."""
    bits, seed = args
    inst = balanced_semiprime(bits, seed)
    N = inst.N
    rng = random.Random(seed)
    lo, hi = dlog.fermat_interval(N)
    out = {"bits": bits, "seed": seed}
    a = rng.randrange(2, N - 1)
    y = pow(a, N - 1, N)
    y2 = pow(y, N - 1, N)
    assert y2 == pow(a, (N - 1) ** 2, N)
    # plain
    m = Metrics(); G = dlog.Mult(N, None, m)
    x = dlog.kangaroo(G, a, y, lo, hi, rng)
    out["plain_ok"] = x is not None and dlog.factor_from_s(N, x + 2) is not None
    out["plain_mul"] = m.counts.get("grp_mul", 0)
    # (i) kangaroo on the squared exponent range
    m = Metrics(); G = dlog.Mult(N, None, m)
    x2 = dlog.kangaroo(G, a, y2, lo * lo, hi * hi, rng, max_factor=60)
    ok = False
    if x2 is not None:
        r = isqrt(x2)
        ok = r * r == x2 and dlog.factor_from_s(N, r + 2) is not None
    out["sq_kangaroo_found"] = x2 is not None
    out["sq_kangaroo_sqrt_ok"] = ok
    # a collision x != (N-1)^2 still gives a multiple of ord(a): split with the single-base ladder
    ladder = None
    if x2 is not None and not ok:
        ladder = dlog.split_with_order_multiple(N, a, (N - 1) ** 2 - x2)
    out["sq_kangaroo_ladder_split"] = ladder is not None
    out["sq_kangaroo_ok"] = ok or ladder is not None
    out["sq_kangaroo_mul"] = m.counts.get("grp_mul", 0)
    # (ii) finite-difference scan over x in [lo, hi]
    m = Metrics(); G = dlog.Mult(N, None, m)
    cur = G.pow(a, lo * lo)
    inc = G.pow(a, 2 * lo + 1)
    a2 = G.mul(a, a)
    xx = lo
    found = None
    while xx <= hi:
        if cur == y2 and dlog.factor_from_s(N, xx + 2):
            found = xx
            break
        cur = G.mul(cur, inc)
        inc = G.mul(inc, a2)
        xx += 1
    out["sq_scan_ok"] = found is not None
    out["sq_scan_mul"] = m.counts.get("grp_mul", 0)
    # (iii) joint: Fermat + Lucas residues, BSGS on each alone (the joint unknown is 1-dimensional: s fixes q-p)
    m1 = Metrics(); r1 = dlog.factor_via_fermat_residue(N, seed, "bsgs", m1)
    m2 = Metrics(); r2 = dlog.factor_via_lucas_residue(N, seed, "bsgs", m2)
    out["bsgs_fermat_mul"] = m1.counts.get("grp_mul", 0)
    out["bsgs_lucas_mul"] = m2.counts.get("grp_mul", 0)
    out["bsgs_fermat_ok"] = r1.get("factors") is not None
    out["bsgs_lucas_ok"] = r2.get("factors") is not None
    out["W_fermat"] = hi - lo
    out["W_lucas"] = 2 * (isqrt(5 * N) + 1)
    return out


def run_iterate(a):
    cfg = {"bits": [16, 20, 24, 28, 32, 36], "per_size": 20}
    json.dump(cfg, open("experiments/configs/primality_iterate.json", "w"), indent=1)
    out = RAW + "primality_iterate.jsonl"
    open(out, "w").close()
    jobs = [(b, s) for b in cfg["bits"] for s in range(cfg["per_size"])]
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(e3_job, jobs):
            append_jsonl(out, rec)


# ------------------------------------------------------------------ E4a
def cf_controls(N, T):
    """Continued fraction of sqrt(N): Q_i via the standard recurrence.  PC1: period odd and <= T
    (Q_l = 1 at odd l).  PC2: some odd i <= T with Q_i a perfect square k^2, gcd(k, N) = 1.  Both certify eps = +1.
    Operational (uses N only).  Copied minimal recurrence (not imported from src/eps)."""
    a0 = isqrt(N)
    P, Q = 0, 1
    a_ = a0
    pc1 = pc2 = False
    first_sq = None
    for i in range(1, T + 1):
        P = a_ * Q - P
        Q = (N - P * P) // Q
        a_ = (a0 + P) // Q
        if i % 2 == 1:
            if Q == 1:
                pc1 = True
            if gmpy2.is_square(Q) and gcd(isqrt(Q), N) == 1:
                pc2 = True
                if first_sq is None:
                    first_sq = i
        if Q == 1:
            break
    return pc1, pc2, first_sq


def e4a_job(args):
    bits, seed, T = args
    inst = balanced_semiprime(bits, seed)
    p, q = inst.labels["factors"]
    N = inst.N
    rec = {"bits": bits, "seed": seed, "Nmod4": N % 4, "eps": int(sympy.legendre_symbol(p - 1, p)) if N % 4 == 1 else None}
    if N % 4 != 1:
        return rec
    k = min(F.v2(p - 1), F.v2(q - 1))
    u, v = F.v2(p - 1), F.v2(q - 1)
    g = gcd(p - 1, q - 1)
    Fq = Fr(g * g, (p - 1) * (q - 1))
    rec["mr_cert"] = log2f(Fq * Fr(4 ** k - 4, 3 * 4 ** k))
    rec["ss_cert"] = log2f(Fq / 4) if u != v else None
    pc1, pc2, fs = cf_controls(N, T)
    rec.update({"pc1": pc1, "pc2": pc2, "first_sq": fs, "T": T})
    return rec


def run_epscert(a):
    cfg = {"bits": [16, 24, 32, 48, 64, 96, 128], "per_size": 4000, "T": 1024}
    json.dump(cfg, open("experiments/configs/primality_epscert.json", "w"), indent=1)
    out = RAW + "primality_epscert.jsonl"
    open(out, "w").close()
    jobs = [(b, s, cfg["T"]) for b in cfg["bits"] for s in range(cfg["per_size"])]
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(e4a_job, jobs, chunksize=32):
            append_jsonl(out, rec)


# ------------------------------------------------------------------ E4b
def mr_pattern_law(p, q):
    """Exact law of the MR +-1 pattern (S0, L1, L2, ..., NONE) and of the leak flag, per base,
    from the tau laws (Thm 2/3, verified exhaustively)."""
    N = p * q
    ex = F.mr_exact([p, q])
    law = {}
    for j, pr in ex["liar_level"].items():
        law["S0" if j == 0 else f"L{j}"] = pr
    law["NONE"] = 1 - ex["liar"]
    return law, ex["pm"]


def ss_pattern_law(p, q):
    """Exact law of the Euler value class: b = +1, b = -1, other (mod N)."""
    N = p * q
    E = (N - 1) // 2
    dists = [F.indicator_dist(r - 1, [(E, 1), (E, -1)]) for r in (p, q)]
    pp = {"+1": Fr(0), "-1": Fr(0)}
    for v, pr in dists[0].items():
        if v in dists[1]:
            if v == (1, 0):
                pp["+1"] += pr * dists[1][v]
            if v == (0, 1):
                pp["-1"] += pr * dists[1][v]
    pp["other"] = 1 - pp["+1"] - pp["-1"]
    return pp


def run_epsexact(a):
    """Exact conditional transcript laws given (eps, N mod 8), aggregated over ALL pairs of primes p < q in
    [2^(b-1), 2^b) with pq = 1 mod 4 (each pair weight 1, each base uniform)."""
    res = []
    for b in (7, 8, 9, 10, 11):
        ps = list(sympy.primerange(2 ** (b - 1), 2 ** b))
        agg = {}
        for i, p in enumerate(ps):
            for q in ps[i + 1:]:
                N = p * q
                if N % 4 != 1:
                    continue
                eps = 1 if p % 4 == 1 else -1
                key = (eps, N % 8)
                law, leak = mr_pattern_law(p, q)
                ssl = ss_pattern_law(p, q)
                A = agg.setdefault(key, {"pairs": 0, "mr": {}, "ss": {}, "leak": Fr(0)})
                A["pairs"] += 1
                A["leak"] += leak
                for k_, v_ in law.items():
                    A["mr"][k_] = A["mr"].get(k_, 0) + v_
                for k_, v_ in ssl.items():
                    A["ss"][k_] = A["ss"].get(k_, 0) + v_
        for key, A in agg.items():
            n = A["pairs"]
            A["mr"] = {k_: v_ / n for k_, v_ in A["mr"].items()}
            A["ss"] = {k_: v_ / n for k_, v_ in A["ss"].items()}
            A["leak"] = A["leak"] / n
        for c in (1, 5):
            P_, M_ = agg.get((1, c)), agg.get((-1, c))
            if not P_ or not M_:
                continue
            keys = set(P_["mr"]) | set(M_["mr"])
            tv_mr = sum(abs(P_["mr"].get(k_, 0) - M_["mr"].get(k_, 0)) for k_ in keys) / 2
            keys = set(P_["ss"]) | set(M_["ss"])
            tv_ss = sum(abs(P_["ss"].get(k_, 0) - M_["ss"].get(k_, 0)) for k_ in keys) / 2
            row = {"b": b, "N_mod8": c, "pairs_eps+": P_["pairs"], "pairs_eps-": M_["pairs"],
                   "tv_mr": float(tv_mr), "tv_ss": float(tv_ss),
                   "mr_eps+": {k_: float(v_) for k_, v_ in sorted(P_["mr"].items())},
                   "mr_eps-": {k_: float(v_) for k_, v_ in sorted(M_["mr"].items())},
                   "ss_eps+": {k_: float(v_) for k_, v_ in P_["ss"].items()},
                   "ss_eps-": {k_: float(v_) for k_, v_ in M_["ss"].items()},
                   "leak_eps+": float(P_["leak"]), "leak_eps-": float(M_["leak"])}
            res.append(row)
            print(b, c, f"tv_mr={float(tv_mr):.3e} tv_ss={float(tv_ss):.3e}", "L>=2 eps-:",
                  sum(v for k_, v in M_["mr"].items() if k_.startswith("L") and int(k_[1:]) >= 2), flush=True)
    with open(RAW + "primality_epsexact.json", "w") as f:
        json.dump(res, f, indent=1)


# ------------------------------------------------------------------ E4c
BATTERY_FEATURES = [
    "mr_pattern_none", "mr_cert_sqrtm1", "ss_b_minus1", "ss_euler_class_pm1",
    "fermat_y_bit0", "fermat_y_bit1", "mr_b0_bit0",
    "jac_y_minus1", "jac_y_plus1", "jac_b0_plus1", "jac_euler_plus_jac",
    "lucas_pattern_none", "lucas_U_bit0", "jac_lucasV",
]


def e4c_job(args):
    bits, seed, nb = args[:3]
    inst = balanced_semiprime(bits, seed)
    N = inst.N
    if N % 4 != 1:
        return None
    p = inst.labels["factors"][0]
    eps = 1 if p % 4 == 1 else -1
    rng = random.Random(f"bat|{bits}|{seed}")
    feats = {f: 0 for f in BATTERY_FEATURES}
    sel = ops.selfridge(N)
    for _ in range(nb):
        a = rng.randrange(2, N - 1)
        mr = ops.miller_rabin(N, a)
        seq = mr["seq"]
        y = seq[-1]
        feats["mr_pattern_none"] += mr["pattern"] == "NONE"
        feats["mr_cert_sqrtm1"] += any(seq[i] == N - 1 for i in range(1, len(seq) - 1))
        ss = ops.solovay_strassen(N, a)
        feats["ss_b_minus1"] += ss["euler"] == N - 1
        feats["ss_euler_class_pm1"] += ss["euler_class"] != "other"
        feats["fermat_y_bit0"] += y & 1
        feats["fermat_y_bit1"] += (y >> 1) & 1
        feats["mr_b0_bit0"] += seq[0] & 1
        feats["jac_y_minus1"] += int(gmpy2.jacobi(y - 1, N)) == 1
        feats["jac_y_plus1"] += int(gmpy2.jacobi(y + 1, N)) == 1
        feats["jac_b0_plus1"] += int(gmpy2.jacobi(seq[0] + 1, N)) == 1
        feats["jac_euler_plus_jac"] += int(gmpy2.jacobi((ss["euler"] + ss["jacobi"]) % N, N)) == 1
        if "P" in sel:
            P = rng.randrange(N)
            Q = (P * P - sel["D"]) * pow(4, -1, N) % N
            lu = ops.strong_lucas(N, P, Q)
            feats["lucas_pattern_none"] += lu["pattern"] == "NONE"
            feats["lucas_U_bit0"] += lu["U"][-1] & 1
            feats["jac_lucasV"] += int(gmpy2.jacobi(lu["V"][-1], N)) == 1
    # controls (frozen in experiments/configs/primality_battery.json)
    pc1, pc2, _ = cf_controls(N, 256)
    feats["CONTROL_cf_pc2_T256"] = int(pc2)
    prng = random.Random(f"planted|{bits}|{seed}")
    feats["CONTROL_planted_rho0.06"] = int(eps == 1) if prng.random() < 0.53 else int(eps != 1)
    return {"bits": bits, "seed": seed, "eps": eps, "N_mod8": N % 8, "nb": nb, "feats": feats}


def run_battery(a):
    cfg = json.load(open("experiments/configs/primality_battery.json"))
    out = RAW + "primality_battery.jsonl"
    open(out, "w").close()
    jobs = [(b, s, cfg["bases_per_N"]) for b in cfg["bits"] for s in range(cfg["seed0"], cfg["seed0"] + cfg["instances_per_bits"])]
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(e4c_job, jobs, chunksize=64):
            if rec:
                append_jsonl(out, rec)


# ------------------------------------------------------------------ E5
def e5_job(args):
    bits, seed, n_a = args
    inst = balanced_semiprime(bits, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    m = Metrics()
    r, ell = ops.aks_parameters(N)
    t0 = time.perf_counter()
    res = ops.aks(N, r=r, a_values=range(1, n_a + 1), metrics=m)
    struct = [l for l in res["leaks"] if l["class"] == "structural"]
    trial = [l for l in res["leaks"] if l["kind"] == "aks_trial"]
    return {"bits": bits, "seed": seed, "r": r, "ell": ell, "a_checked": n_a, "verdict": res["verdict"],
            "structural": len(struct), "trial": len(trial), "coeff_hits": res["n_coeff_hits"],
            "null_coeff": n_a * r * (1 / p + 1 / q), "wall_s": time.perf_counter() - t0,
            "p_gt_r": p > r, "q_minus_p": q - p, "sparse_regime": 2 * (q - p + 1) < r,
            "coeff_divisors": sorted({l["divisor"] == p for l in res["leaks"] if l["class"] == "coefficient"})}


def run_aks(a):
    cfg = {"bits": [20, 24, 28, 32, 40], "per_size": 40, "a_checked": 64}
    json.dump(cfg, open("experiments/configs/primality_aks.json", "w"), indent=1)
    out = RAW + "primality_aks.jsonl"
    open(out, "w").close()
    jobs = []
    for b in cfg["bits"]:
        for s in range(cfg["per_size"]):
            jobs.append((b, s, 64))
    with Pool(a.workers) as pool:
        for rec in pool.imap_unordered(e5_job, jobs):
            append_jsonl(out, rec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("--workers", type=int, default=30)
    a = ap.parse_args()
    fn = {"leakrates": run_leakrates, "kangaroo": run_kangaroo, "iterate": run_iterate, "epscert": run_epscert,
          "epsexact": run_epsexact, "battery": run_battery, "aks": run_aks}
    if a.cmd == "menu":
        from src.primality import menu
        menu.run(a)
        return
    fn[a.cmd](a)


if __name__ == "__main__":
    main()
