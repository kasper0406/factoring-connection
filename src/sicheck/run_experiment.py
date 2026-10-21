"""sicheck experiments E1 (H_S and residue distribution), E3 (fruitfulness / rank), E4 (end-to-end toy).

HARNESS. Operational pieces come from src.sicheck.core (they see only N, kN and public parameters). This file also
uses p, q (instance labels) for stratification (Blum / non-Blum) and to verify factors; such uses are marked _diag.

Usage:
  .venv/bin/python -m src.sicheck.run_experiment E1 [cell ...]
  .venv/bin/python -m src.sicheck.run_experiment E3
  .venv/bin/python -m src.sicheck.run_experiment E4
Raw output: experiments/raw/sicheck_E{1,3,4}.jsonl (appended).
"""
from __future__ import annotations

import json
import os
import random
import resource
import sys
import time
from math import gcd

from src.common.instances import balanced_semiprime
from src.sicheck import core
from src.sicheck.diag import NullDensity, semiprimes_by_type_diag

RAW = "experiments/raw"

# E1 cells: name -> (nbits, y, h, K intervals). Fixed in notes/sicheck_prereg.md.
E1_CELLS = {
    "A": (40, 100, 6_300_000, 6000),
    "B": (48, 200, 8_000_000, 6000),
    "B10": (48, 200, 27_000_000, 3000),
    "C": (56, 600, 3_400_000, 6000),
    "D": (64, 1500, 2_800_000, 4000),
    "E": (80, 6000, 6_600_000, 2000),
}
E1_SEEDS = (0, 1)
CONS = ("next", "last", "rand", "cofactor")

# E3 cells: nbits -> (y, h); 6 Blum + 6 non-Blum N each.
E3_CELLS = {40: (150, 1_070_000), 48: (300, 1_120_000), 56: (600, 3_400_000)}
E3_QS_YFACTOR = 2.2   # qsfb uses the primes l <= 2.2 y with J_N(l) = +1 (about pi(y) of them)
E3_PER_TYPE = 6
E3_EXCESS = 20
E3_KCAP = 15_000
E3_RANDOM_DEPS = 64

# E4: nbits -> (y, h); seeds 0..2; constructors rand, next.
E4_CELLS = {40: (150, 1_070_000), 48: (300, 1_120_000), 56: (600, 3_400_000), 64: (1500, 2_800_000),
            72: (2500, 5_000_000), 80: (6000, 6_600_000)}
E4_NEXT_MAXBITS = 64  # constructor "next" also run up to 64 bits; "rand" at every size
E4_SEEDS = (0, 1, 2)


def rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


def dump(name: str, rec: dict) -> None:
    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, f"sicheck_{name}.jsonl"), "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")


# ----------------------------------------------------------------------------------------------------------------------
def run_E1_cell(cell: str) -> None:
    nbits, y, h, K = E1_CELLS[cell]
    primes = core.primes_upto(y)
    null = NullDensity(h, primes)
    for seed in E1_SEEDS:
        t0 = time.time()
        N = balanced_semiprime(nbits, seed).N
        counts = []
        per = {c: {"k": [], "v": [], "vs": [], "pn": [], "pe": [], "po": [], "ue": [], "g": [], "dup": [], "u0": []}
               for c in CONS}
        seen = {c: set() for c in CONS}
        for k in range(1, K + 1):
            X = k * N
            sm = core.smooth_in_interval(X, h + 1, primes)
            counts.append(len(sm))
            for c in CONS:
                if c == "cofactor":
                    ch = core.cofactor_choice(X, h, sm, primes)
                    u, u0 = (ch[1], ch[0]) if ch else (None, None)
                else:
                    u = core.CONSTRUCTORS[c](X, h, sm, primes, core.interval_rng(seed, X))
                    u0 = None
                if u is None:
                    continue
                v = u - X
                vs = core.is_smooth(v, primes)
                g = gcd(u, v)
                key = (u // g, v // g)
                dup = bool(vs and key in seen[c])
                if vs:
                    seen[c].add(key)
                d = per[c]
                d["k"].append(k); d["v"].append(v); d["vs"].append(int(vs)); d["pn"].append(null.p(v))
                d["pe"].append(null.p_parity(v, 0)); d["po"].append(null.p_parity(v, 1)); d["ue"].append(1 - u % 2)
                d["g"].append(g); d["dup"].append(int(dup)); d["u0"].append(u0)
        for c in CONS:
            dump("E1", {"cell": cell, "seed": seed, "nbits": nbits, "N": N, "y": y, "h": h, "K": K,
                        "constructor": c, "counts": counts if c == "next" else None, **per[c],
                        "secs": round(time.time() - t0, 1), "rss_mb": round(rss_mb(), 1)})
        print(f"E1 {cell} seed {seed}: {time.time() - t0:.0f}s rss {rss_mb():.0f}MB", flush=True)


# ----------------------------------------------------------------------------------------------------------------------
class RelState:
    def __init__(self, primes, need):
        self.primes = primes
        self.need = need
        self.rels = []
        self.seen = set()
        self.deps = core.Dependencies()
        self.phi_seq = []          # Phi of each new basis dependency, in order of appearance (prefix statistic)
        self.dups = 0
        self.vsmooth = 0
        self.outputs = 0

    def done(self):
        return len(self.rels) >= self.need

    def offer(self, k, u, X, N):
        self.outputs += 1
        v = u - X
        if v <= 0:
            return
        ev = core.exponent_vector(v, self.primes)
        if ev is None:
            return
        self.vsmooth += 1
        r = core.Relation(k, u, v, core.exponent_vector(u, self.primes), ev)
        if r.key in self.seen:
            self.dups += 1
            return
        self.seen.add(r.key)
        self.rels.append(r)
        mask = self.deps.add(r.parity())
        if mask is not None:
            Z, gg = core.square_root_of_unity(self.rels, mask, self.primes, N)
            self.phi_seq.append(None if gg != 1 else int(core.classify_root(Z, N)[0]))


def beta0_choice(X, h, sm, primes, N):
    """Adversarial but p-oblivious constructor (uses N via Jacobi symbols only): avoid outputs that would be a
    relation with beta = 1; otherwise the smallest smooth u."""
    ok = []
    for u in sm:
        v = u - X
        if v > 0 and core.is_smooth(v, primes) and core.beta_blum(u, v, primes, N) == 1:
            continue
        ok.append(u)
    return ok[0] if ok else None


def finish_state(st: RelState, N: int, blum: bool, seed: int) -> dict:
    primes = st.primes
    rng = random.Random(f"sicheck-E3|{seed}|{N}")
    basis = st.deps.basis
    phi_basis = []
    beta_mism = 0
    beta_checked = 0
    beta = [core.beta_blum(r.u, r.v, primes, N) for r in st.rels]

    def phi_and_check(mask):
        nonlocal beta_mism, beta_checked
        Z, gg = core.square_root_of_unity(st.rels, mask, primes, N)
        if gg != 1:
            return None
        fr = int(core.classify_root(Z, N)[0])
        pred = 0
        i, m = 0, mask
        while m:
            if m & 1:
                pred ^= beta[i]
            m >>= 1
            i += 1
        beta_checked += 1
        beta_mism += int(pred != fr)
        return fr

    for mask in basis:
        phi_basis.append(phi_and_check(mask))
    fruit = 0
    tried = 0
    for _ in range(E3_RANDOM_DEPS if basis else 0):
        mask = 0
        while mask == 0:
            for b in basis:
                if rng.getrandbits(1):
                    mask ^= b
        fr = phi_and_check(mask)
        if fr is not None:
            tried += 1
            fruit += fr
    return {"fb": len(primes), "relations": len(st.rels), "rank": st.deps.rank, "dim_deps": len(basis),
            "phi_basis": phi_basis, "phi_nonzero": any(x == 1 for x in phi_basis), "phi_seq": st.phi_seq,
            "random_tried": tried, "random_fruitful": fruit, "beta_checked": beta_checked,
            "beta_mismatch": beta_mism, "beta_ones": sum(beta), "dups": st.dups, "vsmooth": st.vsmooth,
            "outputs": st.outputs, "done": st.done()}


def run_E3() -> None:
    for nbits, (y, h) in E3_CELLS.items():
        primes = core.primes_upto(y)
        for blum in (True, False):
            for (seed, N, p, q) in semiprimes_by_type_diag(nbits, E3_PER_TYPE, blum):
                t0 = time.time()
                qfb = core.qs_factor_base(core.primes_upto(int(E3_QS_YFACTOR * y)), N)
                # J+ base ("qsfb"): primes l <= 2.2y with J_N(l) = +1; run only when J_N(2) = +1 so that 2 is in it
                names = list(CONS) + (["qsfb"] if core.jacobi(2, N) == 1 else []) + (["beta0"] if blum else [])
                st = {c: RelState(qfb if c == "qsfb" else primes, (len(qfb) if c == "qsfb" else len(primes))
                                  + E3_EXCESS) for c in names}
                k = 0
                while k < E3_KCAP and not all(s.done() for s in st.values()):
                    k += 1
                    X = k * N
                    sm = core.smooth_in_interval(X, h + 1, primes)
                    for c in names:
                        s = st[c]
                        if s.done():
                            continue
                        if c in CONS:
                            u = core.CONSTRUCTORS[c](X, h, sm, primes, core.interval_rng(seed, X))
                        elif c == "qsfb":
                            smq = core.smooth_in_interval(X, h + 1, qfb)
                            u = smq[0] if smq else None
                        else:
                            u = beta0_choice(X, h, sm, primes, N)
                        if u is not None:
                            s.offer(k, u, X, N)
                for c in names:
                    rec = finish_state(st[c], N, blum, seed)
                    dump("E3", {"nbits": nbits, "y": y, "h": h, "seed": seed, "N": N, "blum_diag": blum,
                                "p_mod4_diag": p % 4, "q_mod4_diag": q % 4, "constructor": c, "k_last": k, **rec})
                print(f"E3 {nbits} blum={blum} seed {seed}: k={k} {time.time() - t0:.0f}s rss {rss_mb():.0f}MB",
                      flush=True)


# ----------------------------------------------------------------------------------------------------------------------
E3B_CELLS = {32: (300, 300_000), 36: (600, 1_000_000)}
E3B_PER_TYPE = 8


def run_E3b() -> None:
    """Addendum A1: J+ base degeneracy on Blum N (Cor. B(a)), with full-base controls and non-Blum contrasts."""
    for nbits, (Y, h) in E3B_CELLS.items():
        full = core.primes_upto(Y)
        for blum in (True, False):
            found = 0
            s0 = 0
            while found < E3B_PER_TYPE:
                (seed, N, p, q), = semiprimes_by_type_diag(nbits, 1, blum, start_seed=s0)
                s0 = seed + 1
                if core.jacobi(2, N) != 1:
                    continue
                found += 1
                qfb = core.qs_factor_base(full, N)
                for base_name, base in (("Jplus", qfb), ("full", full)) if blum else (("Jplus", qfb),):
                    t0 = time.time()
                    res = core.factor_si(N, Y, h, constructor="next", seed=seed, excess=10, factor_base=base,
                                         k_max=20_000, n_random_deps=32)
                    st = res["stats"]
                    dump("E3b", {"nbits": nbits, "Y": Y, "h": h, "seed": seed, "N": N, "blum_diag": blum,
                                 "p_mod4_diag": p % 4, "q_mod4_diag": q % 4, "base": base_name, "fb": len(base),
                                 "relations": st["relations"], "dim_deps": st["dim_deps"],
                                 "phi_nonzero": st["phi_nonzero"], "random_deps": st["random_deps"],
                                 "random_fruitful": st["random_fruitful"], "intervals": st["intervals"],
                                 "factor_ok": bool(res["factor"] and N % res["factor"] == 0 and 1 < res["factor"] < N),
                                 "secs": round(time.time() - t0, 1)})
                print(f"E3b {nbits} blum={blum} seed {seed} done", flush=True)


# ----------------------------------------------------------------------------------------------------------------------
def run_E4() -> None:
    for nbits, (y, h) in E4_CELLS.items():
        primes = core.primes_upto(y)
        null = NullDensity(h, primes)
        for seed in E4_SEEDS:
            inst = balanced_semiprime(nbits, seed)
            N = inst.N
            for c in (("rand", "next") if nbits <= E4_NEXT_MAXBITS else ("rand",)):
                t0 = time.time()
                trace = []
                res = core.factor_si(N, y, h, constructor=c, seed=seed, excess=20, trace=trace)
                f = res["factor"]
                ok = bool(f and 1 < f < N and N % f == 0)  # harness check
                st = res.get("stats", {})
                pred = sum(null.p(v) for (_, v) in trace)
                dump("E4", {"nbits": nbits, "y": y, "h": h, "seed": seed, "N": N, "constructor": c,
                            "factor_ok": ok, "kind": res["kind"],
                            **{kk: vv for kk, vv in st.items() if kk != "phi_basis"},
                            "phi_basis": st.get("phi_basis"), "pred_vsmooth_null": pred,
                            "secs": round(time.time() - t0, 1), "rss_mb": round(rss_mb(), 1)})
                print(f"E4 {nbits} seed {seed} {c}: ok={ok} intervals={st.get('intervals')} "
                      f"rels={st.get('relations')} {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "E1":
        for cell in (sys.argv[2:] or list(E1_CELLS)):
            run_E1_cell(cell)
    elif what == "E3":
        run_E3()
    elif what == "E3b":
        run_E3b()
    elif what == "E4":
        run_E4()
    else:
        raise SystemExit("usage: E1 [cells] | E3 | E4")
