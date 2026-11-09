"""grdeg experiments E1 (frame-less box excess) and E2 (linking obstruction). Writes exact JSONL raw data.

Usage: .venv/bin/python -m src.grdeg.run_experiment [E1|E2|all]
"""
from __future__ import annotations

import json
import os
from typing import Dict, List

from sympy import nextprime

from src.grdeg import construct as C
from src.grdeg import obstruction as O
from src.rank3 import forms as RF

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")
SEED = 20261109


def _semiprime(bits: int, offset: int) -> int:
    import random
    rng = random.Random(SEED + offset)
    half = bits // 2
    p = nextprime(rng.randrange(2 ** (half - 1), 2 ** half))
    q = nextprime(rng.randrange(2 ** (half - 1), 2 ** half))
    return int(p * q)


def e1() -> List[Dict]:
    """For each family and N: geometry + sublattice-confinement of small norms (P1,P2,P3 operationalised)."""
    rows = []
    Ns = [_semiprime(b, i) for i, b in enumerate((24, 34, 44))]
    ctors = {"D1_weil_q2q3": C.weil_restriction_q2q3,
             "D2_cyclo_z5": C.cyclotomic_z5,
             "GP_control_cubics": C.gp_control_cubics}
    E = 5
    for name, ctor in ctors.items():
        for N in Ns:
            mod = ctor(N)
            geo = C.geometry(mod)
            r = mod["r"]
            # sublattice vs genuine: 'genuine' = any coordinate beyond the unshifted Q-rank active.
            # For D1 the shifted generators are indices 2,3 (sqrt3, sqrt6); for others, index >=1 beyond first two.
            shifted = [j for j in range(r) if any(c != 0 for c in
                       [mod["Q"][j][k] - (mod["P"][j][k] if k < len(mod["P"][j]) else 0)
                        for k in range(len(mod["Q"][j]))])]
            sub_min = None
            gen_min = None
            n_sq = 0
            n_tot = 0
            for a in RF.box_points(r, E):
                w, w2 = RF.F_value(mod["f"], mod["P"], mod["g"], mod["Q"], a)
                if w == 0 or w2 == 0:
                    continue
                n_tot += 1
                v = abs(w * w2)
                active = any(a[j] != 0 for j in shifted)
                if active:
                    gen_min = v if gen_min is None else min(gen_min, v)
                else:
                    sub_min = v if sub_min is None else min(sub_min, v)
                import math
                rt = math.isqrt(abs(w))
                if rt * rt == abs(w):
                    n_sq += 1
            rows.append({
                "exp": "E1", "family": name, "N_bits": int(N).bit_length(), "D": geo["D"], "r": r,
                "Q2_zero(DeltaF=0)": bool(geo["Q2_is_zero"]), "frame_free": bool(geo["frame_free"]),
                "num_frames": geo["num_frames"], "zero_minors": geo["num_zero_minors"],
                "num_minors": geo["num_minors"], "box_E": E, "box_pts": n_tot,
                "min|F|_sublattice(shift inactive)": sub_min,
                "min|F|_genuine(shift active)": gen_min,
                "min|F|_genuine_over_N": (float(gen_min) / N if gen_min else None),
                "W_perfect_square_count": n_sq,
            })
    return rows


def e2() -> List[Dict]:
    """Linking obstruction: requirement statement, honest brute search cost, and Rabin split classification."""
    rows = []
    # k = Q(sqrt2): mu = x^2 - 2
    mu = [-2, 0, 1]
    for bits in (24, 34, 44, 60):
        # pick primes with 2 a QR (==+-1 mod 8) so k=Q(sqrt2) splits at both (linking map exists diagonally)
        import random
        rng = random.Random(SEED * 7 + bits)
        half = bits // 2
        def pick(start):
            p = nextprime(start)
            while p % 8 not in (1, 7):
                p = nextprime(p)
            return int(p)
        p = pick(rng.randrange(2 ** (half - 1), 2 ** half))
        q = pick(rng.randrange(2 ** (half - 1), 2 ** half))
        N = p * q
        req = O.link_requirement(N, mu)
        brute = O.naive_link_attempt(N, mu, search_bound=min(10 ** 5, p + 1))
        rab = O.rabin_split_harness(2, p, q)
        structural = [x for x in rab.get("results", []) if x["class"].startswith("STRUCTURAL")]
        rows.append({
            "exp": "E2", "subfield": "Q(sqrt2)", "N_bits": N.bit_length(),
            "e": req["subfield_degree_e"], "charged_subroutine": req["charged_subroutine"],
            "hardness": req["hardness"],
            "brute_root_found_below_1e5": brute["found"], "brute_trials": brute["trials"],
            "rabin_structural_splits": len(structural),
            "rabin_results": rab.get("results", []),
            "null_random_gcd_rate(1/p)": O.null_gcd_rate(p),
        })
    return rows


def _write(name: str, rows: List[Dict]):
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, name)
    with open(path, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"wrote {len(rows)} rows -> experiments/raw/{name}")


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("E1", "all"):
        r1 = e1()
        _write("grdeg_E1.jsonl", r1)
        for r in r1:
            print(r["family"], r["N_bits"], "DeltaF=0" if r["Q2_zero(DeltaF=0)"] else "GP",
                  "frame_free" if r["frame_free"] else "framed",
                  "sub_min|F|=", r["min|F|_sublattice(shift inactive)"],
                  "gen_min|F|/N=", r["min|F|_genuine_over_N"])
    if which in ("E2", "all"):
        r2 = e2()
        _write("grdeg_E2.jsonl", r2)
        for r in r2:
            print(r["N_bits"], "bit:", r["charged_subroutine"], "| brute<1e5:", r["brute_root_found_below_1e5"],
                  "| structural splits:", r["rabin_structural_splits"])
