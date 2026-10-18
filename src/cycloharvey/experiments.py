"""Step-4 experiments (HARNESS: uses instance labels for classification only).

  .venv/bin/python -m src.cycloharvey.experiments harvey
  .venv/bin/python -m src.cycloharvey.experiments cyclo
  .venv/bin/python -m src.cycloharvey.experiments onesided

E1 harvey  : operational Harvey-style search (src/cycloharvey/operational.py)
             on balanced semiprimes 32..56 bits; full-list counts vs the model
             count (calibration.cost_harvey) and early-exit work; hit class.
E2 cyclo   : the same targets with cyclotomic Frobenius handles, all phi(m)
             guesses; R-multiplications x D vs Z/N multiplications.
E3 onesided: T_{r,F} == alpha^{F(p,q)} mod p (identity) and whether it also
             holds mod q, for asymmetric / symmetric F, p == q or not mod m.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from math import gcd

from src.common.instances import balanced_semiprime
from src.common.metrics import append_jsonl, read_jsonl
from .calibration import best_harvey, cost_harvey
from .operational import Cyclo, cyclo_harvey_search, handle, harvey_search

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
CFG = os.path.join(ROOT, "experiments", "configs", "cycloharvey_experiments.json")

CONFIG = {
    "E1": {"bits": [32, 36, 40, 44, 48, 52, 56], "seeds": list(range(16))},
    "E2": {"bits": [32, 36, 40, 44, 48], "seeds": list(range(8)), "m": [5, 8, 12]},
    "E3": {"bits": [32, 40, 48], "seeds": list(range(40)), "m": [5, 8, 12],
           "F": {"asym_3q+2p": {"(0, 1)": 3, "(1, 0)": 2},
                 "sym_p+q": {"(0, 1)": 1, "(1, 0)": 1},
                 "sym_p2+q2": {"(0, 2)": 1, "(2, 0)": 1},
                 "asym_cubic_q3-9Kq": {"(0, 3)": 1, "(1, 0)": 2}}},
}


def _classify_harvey(N, g, hit, m):
    """STRUCTURAL iff the colliding giant's target equation holds for the
    prime found: i + j m + shift == a (N/g) + b g."""
    if hit is None:
        return "STRUCTURAL"                      # exact match: u verified
    i = hit["baby"]
    a, b, j, shift = hit["target"]
    return "STRUCTURAL" if i + j * m + shift == a * (N // g) + b * g else "ACCIDENTAL"


def e1():
    out = os.path.join(RAW, "cycloharvey_e1_harvey.jsonl")
    if os.path.exists(out):
        os.remove(out)
    for bits in CONFIG["E1"]["bits"]:
        for seed in CONFIG["E1"]["seeds"]:
            inst = balanced_semiprime(bits, seed)
            N = inst.N
            b = best_harvey(N)                   # parameters from N only
            r, m = b["r"], b["m"]
            t = time.time()
            full = harvey_search(N, r, m, seed, early_exit=False)
            t_full = time.time() - t
            early = harvey_search(N, r, m, seed, early_exit=True)
            p = inst.labels["p"]
            ok = full["factor"] in (p, N // p)
            cls = _classify_harvey(N, full["factor"], full.get("hit"), m) if ok else None
            rec = {"bits": bits, "seed": seed, "r": r, "m": m,
                   "model_giant": b["giant"], "model_cost": b["cost"],
                   "impl_giant": full["met"].get("giant"),
                   "status": full["status"], "success": ok, "class": cls,
                   "counts_full": full["met"]["counts"],
                   "counts_early": early["met"]["counts"], "status_early": early["status"],
                   "wall_full_s": round(t_full, 3),
                   "sqrtN_bits": (N.bit_length()) / 2}
            append_jsonl(out, rec)
            print(bits, seed, rec["status"], ok, cls, rec["impl_giant"], rec["model_giant"], flush=True)


def e2():
    out = os.path.join(RAW, "cycloharvey_e2_cyclo.jsonl")
    if os.path.exists(out):
        os.remove(out)
    for bits in CONFIG["E2"]["bits"]:
        for seed in CONFIG["E2"]["seeds"]:
            inst = balanced_semiprime(bits, seed)
            N = inst.N
            p = inst.labels["p"]
            b = best_harvey(N)
            r, mb = max(2, b["r"]), b["m"]
            base = harvey_search(N, r, mb, seed, early_exit=False)
            for mm in CONFIG["E2"]["m"]:
                t = time.time()
                res = cyclo_harvey_search(N, mm, r, mb, seed)
                g = res["factor"]
                ok = g in (p, N // p)
                cls = None
                if ok:
                    if res["status"] == "exact_match" or res.get("hit") is None:
                        cls = "STRUCTURAL" if res["status"] == "exact_match" else "OTHER_GCD"
                    else:
                        s, a, bb, j, shift = res["hit"]["target"]
                        i = res["hit"]["baby"]
                        right_guess = (s % mm) == (g % mm)
                        eq = i + j * mb + shift == a * (N // g) + bb * g
                        cls = "STRUCTURAL" if (right_guess and eq) else "ACCIDENTAL"
                c = res["met"]["counts"]
                rec = {"bits": bits, "seed": seed, "m": mm, "D": res["met"].get("D"),
                       "r": r, "mb": mb, "status": res["status"], "success": ok,
                       "class": cls, "R_mul": c.get("R_mul", 0),
                       "R_base_ops_lower": c.get("R_mul_D", 0),
                       "handles": c.get("handles", 0), "giant": res["met"].get("giant"),
                       "harvey_rn_mul": base["met"]["counts"].get("rn_mul", 0)
                       + base["met"]["counts"].get("rn_mul_pow", 0),
                       "harvey_giant": base["met"].get("giant"),
                       "wall_s": round(time.time() - t, 3),
                       "p_eq_q_mod_m": (p % mm) == ((N // p) % mm)}
                append_jsonl(out, rec)
                print(bits, seed, mm, rec["status"], ok, cls, rec["R_base_ops_lower"],
                      rec["harvey_rn_mul"], flush=True)


def _alpha_power_mod(R, alpha_coeffs, e, prime, m):
    """HARNESS: alpha^e in (Z/prime)[X]/Phi_m, e reduced mod prime^f - 1."""
    Rp = Cyclo(prime, m)
    f = 1
    while pow(prime, f, m) != 1:
        f += 1
    return Rp, Rp.coeffs(Rp.pow(Rp.elt(alpha_coeffs), e % (prime ** f - 1)))


def e3():
    out = os.path.join(RAW, "cycloharvey_e3_onesided.jsonl")
    if os.path.exists(out):
        os.remove(out)
    for bits in CONFIG["E3"]["bits"]:
        for seed in CONFIG["E3"]["seeds"]:
            inst = balanced_semiprime(bits, seed)
            N = inst.N
            p = inst.labels["p"]
            q = N // p
            for mm in CONFIG["E3"]["m"]:
                if gcd(mm, N) != 1:
                    continue
                R = Cyclo(N, mm)
                rng = random.Random(f"e3|{bits}|{seed}|{mm}")
                ac = [rng.randrange(N) for _ in range(R.D)]
                al = R.elt(ac)
                for name, Fs in CONFIG["E3"]["F"].items():
                    F = {tuple(int(x) for x in k.strip("()").split(",")): v for k, v in Fs.items()}
                    T = handle(R, al, p % mm, F)
                    e = sum(c * p ** i * q ** j for (i, j), c in F.items())
                    Rp, lhs_p = _alpha_power_mod(R, ac, e, p, mm)
                    Rq, lhs_q = _alpha_power_mod(R, ac, e, q, mm)
                    Tc = R.coeffs(T)
                    holds_p = lhs_p == [x % p for x in Tc]
                    holds_q = lhs_q == [x % q for x in Tc]
                    rec = {"bits": bits, "seed": seed, "m": mm, "F": name,
                           "p_eq_q_mod_m": p % mm == q % mm,
                           "identity_mod_p": holds_p, "also_mod_q": holds_q}
                    append_jsonl(out, rec)
    rows = read_jsonl(out)
    summ = {}
    for r in rows:
        k = (r["F"], r["p_eq_q_mod_m"])
        s = summ.setdefault(k, [0, 0, 0])
        s[0] += 1
        s[1] += r["identity_mod_p"]
        s[2] += r["also_mod_q"]
    for k, v in sorted(summ.items()):
        print(k, "n=%d identity_mod_p=%d two_sided=%d" % tuple(v))


if __name__ == "__main__":
    os.makedirs(os.path.dirname(CFG), exist_ok=True)
    with open(CFG, "w") as f:
        json.dump(CONFIG, f, indent=1)
    {"harvey": e1, "cyclo": e2, "onesided": e3}[sys.argv[1]]()
