"""HARNESS for the ptecons pre-registered experiments (notes/ptecons_prereg.md).

  .venv/bin/python -m src.ptecons.run_experiment E1      # AVG check on N-independent structured families
  .venv/bin/python -m src.ptecons.run_experiment E2      # aligned-cluster sieve variants on general N
  .venv/bin/python -m src.ptecons.run_experiment pilot   # small pilot (parameter sanity only)

Outputs JSONL under experiments/raw/ptecons_*.jsonl.  Smoothness of residuals v is tested with a
product-of-primes remainder test (v is a residual, never N).  Instances come from
src/common/instances.balanced_semiprime (harness); operational generators see only N.
"""
from __future__ import annotations

import itertools
import json
import math
import os
import random
import sys
import time
from collections import Counter
from typing import Dict, Iterable, List

import gmpy2

from src.common.instances import balanced_semiprime
from . import acs
from .polys import centered, dickson, evaluate, from_roots

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def primes_upto(B: int) -> List[int]:
    sieve = bytearray([1]) * (B + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(B ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(B + 1) if sieve[i]]


class SmoothTester:
    def __init__(self, B: int):
        self.B = B
        self.P = gmpy2.mpz(1)
        for p in primes_upto(B):
            self.P *= p

    def is_smooth(self, v: int) -> bool:
        v = abs(v)
        if v <= 1:
            return True
        x = self.P % v
        e = 1
        while (1 << e) < v.bit_length() + 1:
            x = (x * x) % v
            e += 1
        x = (x * x) % v
        return x == 0


def lam(v: int, N: int) -> float:
    return math.log2(max(abs(v), 1)) / math.log2(N)


def summarize(rels, N: int, st: SmoothTester, rng: random.Random) -> Dict:
    vs, us, pieces = [], set(), set()
    for pc, u, v in rels:
        if u == v:            # exact (non-genuine) relation: u < N/2; excluded
            continue
        us.add(u)
        vs.append(v)
        pieces.update(pc)
    if not vs:
        return {"members": 0}
    lams = sorted(lam(v, N) for v in vs)
    sm = sum(1 for v in vs if st.is_smooth(v))
    ctrl = 0
    for v in vs:
        b = max(abs(v), 2).bit_length()
        c = rng.randrange(1 << (b - 1), 1 << b)
        ctrl += st.is_smooth(c)
    q = lambda f: lams[min(len(lams) - 1, int(f * len(lams)))]
    srt = sorted(abs(v) for v in vs)
    lam_at = {str(r): (lam(srt[r - 1], N) if len(srt) >= r else None) for r in (100, 1000, 10000)}
    return {"members": len(vs), "distinct_u": len(us), "distinct_pieces": len(pieces),
            "lam_min": lams[0], "lam_q10": q(0.10), "lam_median": q(0.5), "lam_at_rank": lam_at,
            "smooth": sm, "smooth_ctrl": ctrl,
            "smooth_ratio": (sm / ctrl) if ctrl else None}


def run_E2(cfg: Dict, out: str) -> None:
    st = SmoothTester(cfg["B"])
    with open(out, "w") as fh:
        for bits, seeds, variants in cfg["sizes"]:
            for seed in seeds:
                N = balanced_semiprime(bits, seed).N
                for name, kw in variants.items():
                    t0 = time.time()
                    gen = getattr(acs, kw["fn"])(N, **{k: v for k, v in kw.items() if k != "fn"},
                                                cap=cfg["cap"])
                    rec = summarize(gen, N, st, random.Random(f"ptecons|E2|{bits}|{seed}|{name}"))
                    rec.update({"exp": "E2", "bits": bits, "seed": seed, "variant": name,
                                "pred_cos_lam": 0.5 + math.log2(kw.get("A", 1)) / math.log2(N)
                                if kw["fn"] == "cos" else None,
                                "secs": round(time.time() - t0, 2)})
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                    print(bits, seed, name, {k: rec.get(k) for k in
                                             ("members", "lam_median", "lam_at_rank", "smooth_ratio",
                                              "secs")}, flush=True)


# ---------------- E1: N-independent structured families (AVG) -----------------

def fam_generic(M: int, A: int, m: int) -> List[int]:
    return sorted({math.prod(m - a for a in c)
                   for c in itertools.combinations_with_replacement(range(-A, A + 1), M)})


def fam_symmetric(J: int, A: int, cs: Iterable[int], m: int) -> List[int]:
    out = set()
    for c in cs:
        for rs in itertools.combinations_with_replacement(range(0, A + 1), J):
            out.add(math.prod((m - c) ** 2 - r * r for r in rs))
    return sorted(out)


def fam_pte_classes(M: int, A: int, t: int, m: int) -> List[int]:
    """Members of nontrivial classes with equal e_1..e_t (multisets in [-A, A])."""
    cls: Dict[tuple, list] = {}
    for c in itertools.combinations_with_replacement(range(-A, A + 1), M):
        P = from_roots(c)
        key = tuple(P[M - 1 - i] for i in range(t))
        cls.setdefault(key, []).append(c)
    out = set()
    for mem in cls.values():
        if len(mem) >= 2:
            for c in mem:
                out.add(math.prod(m - a for a in c))
    return sorted(out)


def fam_dickson(ms: Iterable[int]) -> List[int]:
    out = set()
    for n, ts in ((3, (-2, -1, 1, 2)), (4, (-2, 0, 2)), (6, (-2, -1, 1, 2))):
        D = dickson(n)
        for t in ts:
            c = list(D)
            c[0] -= evaluate(D, t)
            for m in ms:
                out.add(evaluate(c, m))
    return sorted(x for x in out if x)


def fam_composition(cs: Iterable[int], m0: int, m1: int) -> List[int]:
    """Double-quadratic composition: S(X) = (X^2 - c)^2, S(X) - S(a) splits iff a^2 + b^2 = 2c;
    P = prod over chosen pairs (X-a)(X+a)(X-b)(X+b), here one or two pairs per c."""
    out = set()
    for c in cs:
        pairs = []
        a = 0
        while a * a <= c:
            b2 = 2 * c - a * a
            b = math.isqrt(b2)
            if b * b == b2 and a <= b:
                pairs.append((a, b))
            a += 1
        for m in range(m0, m1):
            for (a, b) in pairs:
                out.add((m * m - a * a) * (m * m - b * b))
            for (a, b), (a2, b2) in itertools.combinations(pairs, 2):
                out.add((m * m - a * a) * (m * m - b * b) * (m * m - a2 * a2) * (m * m - b2 * b2))
    return sorted(x for x in out if x)


def run_E1(cfg: Dict, out: str) -> None:
    fams = {
        "generic_M6": fam_generic(6, cfg["A6"], cfg["m6"]),
        "symmetric_J3": fam_symmetric(3, cfg["Asym"], range(-cfg["Csym"], cfg["Csym"] + 1), cfg["msym"]),
        "pte_M5_t3": fam_pte_classes(5, cfg["Apte"], 3, cfg["mpte"]),
        "dickson": fam_dickson(range(cfg["dm0"], cfg["dm1"])),
        "composition": fam_composition(range(1, cfg["ccmax"]), cfg["cm0"], cfg["cm1"]),
    }
    Ns = [balanced_semiprime(cfg["bits"], s).N for s in range(cfg["nN"])]
    V = cfg["V"]
    with open(out, "w") as fh:
        for name, us in fams.items():
            obs = exp = 0.0
            n_used = 0
            for N in Ns:
                for u in us:
                    if abs(u) <= N:      # need |u| > N for a genuine (non-exact) congruence
                        continue
                    n_used += 1
                    exp += (2 * V + 1) / N
                    if abs(centered(u, N)) <= V:
                        obs += 1
            G = obs / exp if exp else None
            rec = {"exp": "E1", "family": name, "family_size": len(us), "pairs": n_used,
                   "observed": obs, "expected": exp, "G": G, "V": V, "bits": cfg["bits"], "nN": cfg["nN"]}
            fh.write(json.dumps(rec) + "\n")
            print(rec, flush=True)


def main() -> None:
    what = sys.argv[1]
    cfg = json.load(open(os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "configs",
                                      "ptecons.json")))
    os.makedirs(RAW, exist_ok=True)
    if what == "E1":
        run_E1(cfg["E1"], os.path.join(RAW, "ptecons_E1.jsonl"))
    elif what == "E2":
        run_E2(cfg["E2"], os.path.join(RAW, "ptecons_E2.jsonl"))
    elif what == "pilot":
        run_E2(cfg["pilot"], os.path.join(RAW, "ptecons_pilot.jsonl"))
    else:
        raise SystemExit("unknown")


if __name__ == "__main__":
    main()
