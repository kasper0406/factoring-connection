"""factorbits experiment E1 (pre-registered in notes/factorbits_prereg.md).  DIAGNOSTIC / HARNESS.

E1: Theorem PT on concrete pure generic-ring trees, exact over prime ranges.
  (a) partition lemma against real instrumented execution mod N (all pairs at n = 28; seeded random
      pairs at n = 32, 36);
  (b) Theorem PT bound for iid factors, targets bit1(min), bit2(min), bit(b-2)(min), bit(b-3)(min);
  (c) eps cross-population bound.
Usage: .venv/bin/python -m src.factorbits.run_experiments e1 [bits ...]
Writes experiments/raw/factorbits_e1.jsonl (one line per (bits, tree)).
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

import numpy as np

from . import pure as P

ROOT = Path(__file__).resolve().parents[2]
import os
RAW = Path(os.environ.get("FACTORBITS_RAW", ROOT / "experiments" / "raw" / "factorbits_e1.jsonl"))


def trees_for(bits: int):
    b = bits
    M = lambda B: P.smooth_exponent(B, b)
    out = {
        "pm1_nested": P.nonadaptive([P.PowTest(2, M(16)), P.PowTest(2, M(64)), P.PowTest(2, M(256))]),
        "mixed_adaptive": P.node(P.PowTest(2, M(64)),
                                 P.node(P.LucasTest(3, M(64)),
                                        P.node(P.ECMTest(7, 2, M(64)), P.leaf("none"), P.leaf("ecm")),
                                        P.leaf("p+1")),
                                 P.node(P.PowTest(3, M(16)), P.leaf("p-1 big"), P.leaf("p-1 small"))),
        "ecm3": P.nonadaptive([P.ECMTest(7, 2, M(128)), P.ECMTest(11, 3, M(128)), P.ECMTest(13, 5, M(128))]),
    }
    if b <= 16:
        lo = 3 << (b - 2)
        out["lipton_interval"] = P.nonadaptive([P.ProductTest(lo, lo + (1 << (b - 3)), "x")])
        out["qr_sniff_full"] = P.nonadaptive([P.ProductTest(1, 1 << (b - 1), "x2")])
        out["qr_sniff_partial"] = P.nonadaptive([P.ProductTest(1, 1 << (b // 2 + 2), "x2")])
    return out


def targets_for(bits: int):
    return {
        "bit1_min": lambda X, Y: (np.minimum(X, Y) >> 1) & 1,
        "bit2_min": lambda X, Y: (np.minimum(X, Y) >> 2) & 1,
        f"bit{bits-3}_min": lambda X, Y: (np.minimum(X, Y) >> (bits - 3)) & 1,
        f"bit{bits-4}_min": lambda X, Y: (np.minimum(X, Y) >> (bits - 4)) & 1,
        "x_mod4": lambda X, Y: X & 3,
    }


def pt_vectorised(primes: np.ndarray, cls_ids: np.ndarray, F) -> dict:
    """Same quantities as pure.pt_bound_check (iid uniform ordered pairs), vectorised."""
    n = len(primes)
    T = np.broadcast_to(F(primes[:, None], primes[None, :]), (n, n)).astype(np.int64)
    vals, Tc = np.unique(T, return_inverse=True)
    Tc = Tc.reshape(T.shape)
    baseline = np.bincount(Tc.ravel()).max() / (n * n)
    good = 0
    sizes = []
    for c in np.unique(cls_ids):
        S = np.nonzero(cls_ids == c)[0]
        sizes.append(len(S))
        good += np.bincount(Tc[np.ix_(S, S)].ravel()).max()
    sizes = np.array(sizes, dtype=np.int64)
    split = n * n - int((sizes ** 2).sum())
    eta = 1 - sizes.max() / n
    ps = split / (n * n)
    return {"success": (good + split) / (n * n), "baseline": float(baseline), "p_split": ps,
            "eta": float(eta), "bound": float(baseline + ps + eta ** 2)}


def run_e1(bits_list):
    RAW.parent.mkdir(parents=True, exist_ok=True)
    with open(RAW, "a") as fh:
        for bits in bits_list:
            primes = P.primes_top2(bits)
            arr = np.array(primes, dtype=np.int64)
            P1 = [x for x in primes if x % 4 == 1]
            P3 = [x for x in primes if x % 4 == 3]
            for name, tree in trees_for(bits).items():
                t0 = time.time()
                cls = P.classes(tree, primes)
                labels = sorted(set(map(repr, cls.values())))
                ids = np.array([labels.index(repr(cls[x])) for x in primes])
                # (a) partition lemma vs real execution
                rng = random.Random(f"factorbits|e1|{bits}|{name}")
                if bits <= 14:
                    pairs = [(primes[i], primes[j]) for i in range(len(primes)) for j in range(i + 1, len(primes))]
                else:
                    k = 20000 if not name.startswith(("lipton", "qr_")) else 2000
                    pairs = []
                    while len(pairs) < k:
                        x, y = rng.sample(primes, 2)
                        pairs.append((x, y))
                mism = 0
                nsplit = 0
                for x, y in pairs:
                    out = P.run_on_N(tree, x * y)
                    if cls[x] == cls[y]:
                        mism += out != ("leaf", cls[x])
                    else:
                        nsplit += 1
                        mism += not (out[0] == "split" and out[1] in (x, y))
                # (b) PT bound
                pt = {tn: pt_vectorised(arr, ids, F) for tn, F in targets_for(bits).items()}
                viol_pt = sum(r["success"] > r["bound"] + 1e-12 for r in pt.values())
                # (c) eps cross bound
                ec = P.eps_cross_check(cls, P1, P3)
                viol_eps = int(ec["success"] > ec["bound"] + 1e-12)
                rec = {"bits_prime": bits, "n": 2 * bits, "tree": name, "n_primes": len(primes),
                       "n_classes": len(labels), "pairs_checked": len(pairs), "pairs_split": nsplit,
                       "partition_mismatches": mism, "pt": pt, "pt_violations": int(viol_pt),
                       "eps_cross": ec, "eps_violations": viol_eps, "seconds": round(time.time() - t0, 1)}
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                print(bits, name, "classes", len(labels), "mism", mism, "viol", viol_pt, viol_eps,
                      "eps", {k: round(v, 4) for k, v in ec.items()}, f"{rec['seconds']}s", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "e1":
        run_e1([int(b) for b in sys.argv[2:]] or [14, 16, 18])
