"""crtdecode experiments E1, E2 (pre-registered: notes/crtdecode_prereg.md, commit 7d07797).

Usage: .venv/bin/python -m src.crtdecode.run_exp E1|E2
HARNESS: the ground-truth sieve (core.window_profile) is computed here and never passed to a decoder.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time

from .core import prefix_incidences, primes_upto, smooth_in_range, window_profile
from .decoders import divisor_trick, heavy_hitter, information_decode_and_certify, random_tests

CONFIGS = {
    "C1": {"bits": 44, "y": 1 << 9, "h": 1 << 14, "B": 3000},
    "C2": {"bits": 52, "y": 1 << 10, "h": 1 << 16, "B": 6000},
}
N_WINDOWS = 300
N_E2 = 100
RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def windows(cfg_name: str):
    c = CONFIGS[cfg_name]
    rng = random.Random(f"crtdecode-E1-{cfg_name}")
    return [rng.randrange(1 << (c["bits"] - 1), 1 << c["bits"]) for _ in range(N_WINDOWS)]


def repertoires(cfg_name: str):
    c = CONFIGS[cfg_name]
    h, y = c["h"], c["y"]
    r1 = smooth_in_range(h // 2, h, y)
    r2 = smooth_in_range(h // 2, 8 * h, y)
    rng = random.Random(f"crtdecode-rep-{cfg_name}")
    rng.shuffle(r1)
    rng.shuffle(r2)
    return r1, r2


def run_e1(out_path: str) -> None:
    with open(out_path, "w") as f:
        for name, c in CONFIGS.items():
            h, y, B = c["h"], c["y"], c["B"]
            r1, r2 = repertoires(name)
            print(f"{name}: |R1| = {len(r1)}, |R2| = {len(r2)}", flush=True)
            for i, X in enumerate(windows(name)):
                t0 = time.time()
                prof = window_profile(X, h, y)
                n_smooth = int(prof["smooth"].sum())
                first = prof["first"]
                rec = {"cfg": name, "i": i, "X": X, "n_smooth": n_smooth, "first": first}
                rec["SV"] = prefix_incidences(X, first, y) if first >= 0 else None
                order = list(range(h + 1))
                random.Random(f"crtdecode-RT-{name}-{i}").shuffle(order)
                decs = {
                    "DT": divisor_trick(X, h, y, r1, B),
                    "DT8": divisor_trick(X, h, y, r2, B),
                    "PH-1.25": heavy_hitter(X, h, y, r2, B, 1.25),
                    "PH-1.5": heavy_hitter(X, h, y, r2, B, 1.5),
                    "RT": random_tests(X, h, y, order, B),
                }
                for k, r in decs.items():
                    if r.found:
                        assert prof["smooth"][r.z - X], "decoder output not smooth"
                    agg = {}
                    for (a, b) in r.tested_log:
                        t = agg.setdefault(str(a), [0, 0])
                        t[0] += 1
                        t[1] += int(b)
                    rec[k] = {"found": r.found, "cost": r.cost, "tests": r.tests, "pointers": r.pointers,
                              "by_hits": agg}
                rec["secs"] = round(time.time() - t0, 3)
                f.write(json.dumps(rec) + "\n")
                if i % 50 == 0:
                    print(name, i, rec["secs"], flush=True)


def run_e2(out_path: str) -> None:
    with open(out_path, "w") as f:
        for name, c in CONFIGS.items():
            h, y, bits = c["h"], c["y"], c["bits"]
            for i, X in enumerate(windows(name)[:N_E2]):
                prof = window_profile(X, h, y)
                truth = [int(j) for j in (prof["smooth"]).nonzero()[0]]
                r = information_decode_and_certify(X, h, y, bits)
                rec = {"cfg": name, "i": i, "X": X, "pi_y": len(primes_upto(y)), "depth": r["depth"],
                       "groups": r["groups"], "exact": r["smooth_offsets"] == truth, "n_smooth": len(truth),
                       "output": r["output"], "leaf_incidences": r["incidences"],
                       "sieve_incidences": prof["incidences"], "h": h}
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    which = sys.argv[1]
    os.makedirs(RAW, exist_ok=True)
    if which == "E1":
        run_e1(os.path.join(RAW, "crtdecode_E1.jsonl"))
    elif which == "E2":
        run_e2(os.path.join(RAW, "crtdecode_E2.jsonl"))
    else:
        raise SystemExit("usage: run_exp E1|E2")
