"""Export real NFS testbeds for agent smoothenum (format documented in notes/cadobench_results.md §7).

experiments/raw/cadobench_testbed_c{D}/
  meta.json               N, polynomials, CADO parameters, alpha, conventions, CADO commit
  qlattices.jsonl         the fixed special-q set (V4 range: nq consecutive special-q from qmin)
  las_relations.jsonl.gz  every relation las output on that set (production parameters, I = I_CADO)
  las_summary.json        las CPU/elapsed/funnel/-T categories for that run (from V4, 3 repeats)
  truth_windows.json      (c60-c80) exhaustive windows: definition and exact counts (V2/V3)
  truth_hits.jsonl.gz     (c60-c80) every S_B point of those windows with flags SCADO / las

Usage: python -m cadobench.export --digits 90
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys

from .relstats import QHDR, parse_rel
from .smoothcount import QLattice, murphy_alpha
from .validate import RAW, WORKDIR_RAW, setup


def cado_commit() -> str:
    try:
        return subprocess.run(["git", "-C", os.path.expanduser("~/src/cado-nfs"), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def parse_raw_out(path: str):
    lats = []
    rels = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("# Sieving side-"):
                g = list(map(int, QHDR.match(line).groups()))
                lats.append(QLattice(q=g[1], rho=g[2], a0=g[3], b0=g[4], a1=g[5], b1=g[6], I=g[7], J=g[8]))
            elif line[0] != "#" and line.strip():
                rels.append(line.strip())
    return lats, rels


def export(digits: int) -> str:
    tb, P, side = setup(digits)
    outdir = os.path.join(RAW, f"cadobench_testbed_c{digits}")
    os.makedirs(outdir, exist_ok=True)
    from .numbers import NUMBERS
    N = NUMBERS[digits][0]
    meta = {
        "digits": digits, "N": str(N), "generator": f"common.instances.balanced_semiprime({NUMBERS[digits][4]}, 0)",
        "f0_rational_low_to_high": [str(c) for c in P["f0"]],
        "f1_algebraic_low_to_high": [str(c) for c in P["f1"]],
        "skew": P.get("skew"), "poly_file": P["raw"],
        "cado_params": tb.params, "qmin": tb.qmin,
        "alpha": [murphy_alpha(P["f0"]), murphy_alpha(P["f1"])],
        "cado_commit": cado_commit(),
        "conventions": {
            "homogenisation": "F(a,b) = sum_k c_k a^k b^(d-k)",
            "special_q_side": 1,
            "lattice": "a = i*a0 + j*a1, b = i*b0 + j*b1, then (a,b) -> (-a,-b) if b < 0; "
                       "i in [-I/2, I/2), j in [0, J); points with gcd(a,b) > 1 are not relations",
            "norms": "n0 = |G(a,b)|, n1 = |F(a,b)|/q",
            "S_B": "n0 is 2^lpb0-smooth and n1 is 2^lpb1-smooth",
            "S_CADO": "S_B and, per side, the cofactor above lim_i (n_i divided by its lim_i-smooth part) is < 2^mfb_i",
            "relation_line": "CADO format a,b:hex primes of G(a,b):hex primes of F(a,b) (with multiplicity; q included)",
        },
    }
    with open(os.path.join(outdir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    v4path = os.path.join(RAW, f"cadobench_v4_c{digits}.json")
    I0 = int(tb.params["I"])
    raw = os.path.join(WORKDIR_RAW, f"c{digits}_v4_I{I0}_t10.las.out")
    if not os.path.exists(raw):
        raw = os.path.join(WORKDIR_RAW, f"c{digits}_v4_I{I0}.las.out")
    if os.path.exists(raw):
        lats, rels = parse_raw_out(raw)
        uniq = {(l.q, l.rho): l for l in lats}
        per = {k: 0 for k in uniq}
        with gzip.open(os.path.join(outdir, "las_relations.jsonl.gz"), "wt") as fh:
            seen = set()
            for line in rels:
                a, b, p0, p1 = parse_rel(line)
                for q in sorted(set(p1)):
                    for (qq, r), lat in ((k, v) for k, v in uniq.items() if k[0] == q):
                        if (a - r * b) % q == 0 and (q, r, a, b) not in seen:
                            seen.add((q, r, a, b))
                            ij = lat.ij(a, b)
                            per[(qq, r)] += 1
                            fh.write(json.dumps({"q": q, "rho": r, "a": a, "b": b, "i": ij[0] if ij else None,
                                                 "j": ij[1] if ij else None, "line": line}) + "\n")
        with open(os.path.join(outdir, "qlattices.jsonl"), "w") as fh:
            for (q, r), l in uniq.items():
                fh.write(json.dumps({"q": q, "rho": r, "a0": l.a0, "b0": l.b0, "a1": l.a1, "b1": l.b1, "I": l.I,
                                     "J": l.J, "las_relations": per[(q, r)]}) + "\n")
    if os.path.exists(v4path):
        v4 = json.load(open(v4path))
        keep = [r for r in v4["runs"] if r["I"] == I0]
        summ = {"I": I0, "nq": v4["nq"], "q0": v4["q0"], "repeats": len(keep),
                "cpu_s": [r["cpu_s"] for r in keep], "elapsed_s": [r["elapsed_s"] for r in keep],
                "relations": keep[0]["relations"], "unique_ab": keep[0]["unique_ab"], "funnel": keep[0]["funnel"],
                "categories": keep[0]["categories"], "threads": keep[0]["threads"], "cmd": keep[0]["cmd"]}
        with open(os.path.join(outdir, "las_summary.json"), "w") as fh:
            json.dump(summ, fh, indent=1)
    v23 = os.path.join(RAW, f"cadobench_v23_c{digits}.json")
    if os.path.exists(v23):
        d = json.load(open(v23))
        with open(os.path.join(outdir, "truth_windows.json"), "w") as fh:
            json.dump({"rows": d["rows"], "totals": d["totals"]}, fh, indent=1)
        hits = v23.replace(".json", "_hits.jsonl.gz")
        if os.path.exists(hits):
            shutil.copy(hits, os.path.join(outdir, "truth_hits.jsonl.gz"))
    return outdir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--digits", type=int, nargs="+", required=True)
    a = ap.parse_args(argv)
    for d in a.digits:
        print(export(d))
    return 0


if __name__ == "__main__":
    sys.exit(main())
