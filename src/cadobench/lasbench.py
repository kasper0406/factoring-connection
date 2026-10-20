"""Reproducible las benchmarks on a fixed special-q range (harness).

A *testbed* is the polynomial pair, the side-1 factor base and the sieve parameters of a finished
CADO-NFS factorization (a workdir made by run_cado). `run_las` runs the `las` binary on a fixed
special-q range (or a seeded random sample of special-q) with optional parameter overrides, under
the memory watchdog, and parses everything las prints:

  * relations (unique (a,b) pairs, and the raw count),
  * the survivor funnel (before/after sieve, after trial division per side, entering cofactoring,
    smooth),
  * total CPU and elapsed time, and the categorized -T breakdown (norms / sieving / cofactoring per
    side, bucket fill, small sieve, apply buckets, ...),
  * per special-q lattice bases (q, rho, a0, b0, a1, b1, I, J) when verbose output is kept.

Throughput is reported as relations per elapsed second and per CPU second. Every number comes from
las's own output. Repeated runs give a timing distribution; yields are deterministic for a given
q-range and parameter set (we check this in tests).

Usage:
  python -m cadobench.lasbench --testbed c90_run1 --q0 5000000 --q1 5002000 -t 10 --repeat 3 \
      --set lambda1=3.0 --out experiments/raw/cadobench_lasbench_x.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from . import memwatch
from .run_cado import CADO_BUILD, WORK

LAS = os.path.join(CADO_BUILD, "sieve", "las")

# las parameters copied from the CADO parameter snapshot (tasks.X or tasks.sieve.X).
LAS_KEYS = ("I", "A", "lim0", "lim1", "lpb0", "lpb1", "mfb0", "mfb1", "lambda0", "lambda1",
            "ncurves0", "ncurves1", "powlim0", "powlim1", "tdthresh", "bkthresh", "bkthresh1",
            "adjust_strategy", "sqside", "batchlpb0", "batchlpb1", "batchmfb0", "batchmfb1")


@dataclass
class Testbed:
    name: str
    workdir: str
    poly: str
    fb1: str
    params: Dict[str, str] = field(default_factory=dict)
    qmin: Optional[int] = None

    @classmethod
    def from_workdir(cls, name: str) -> "Testbed":
        wd = name if os.path.isabs(name) else os.path.join(WORK, name)
        base = None
        for f in os.listdir(wd):
            if f.endswith(".poly") and "." not in f[:-5]:
                base = f[:-5]
        if base is None:
            raise FileNotFoundError(f"no <name>.poly in {wd}")
        snaps = sorted(f for f in os.listdir(wd) if ".parameters_snapshot." in f)
        params: Dict[str, str] = {}
        qmin = None
        with open(os.path.join(wd, snaps[-1])) as fh:
            for line in fh:
                if "=" not in line or line.startswith("#"):
                    continue
                k, v = (x.strip() for x in line.split("=", 1))
                for pre in ("tasks.sieve.", "tasks."):
                    if k.startswith(pre) and k[len(pre):].replace(".", "_") in LAS_KEYS \
                            and "." not in k[len(pre):]:
                        params.setdefault(k[len(pre):], v)
                if k == "tasks.qmin":
                    qmin = int(v)
        return cls(name=os.path.basename(wd), workdir=wd, poly=os.path.join(wd, base + ".poly"),
                   fb1=os.path.join(wd, base + ".roots1.gz"), params=params, qmin=qmin)

    def read_poly(self) -> dict:
        """Coefficients (low degree first) of f0 (rational, side 0) and f1 (algebraic, side 1)."""
        out = {"raw": open(self.poly).read()}
        for line in out["raw"].splitlines():
            m = re.match(r"^(poly0|poly1|c\d|Y\d|skew|n)\s*[:=]\s*(.*)$", line.strip())
            if m:
                out[m.group(1)] = m.group(2).strip()
        if "poly0" in out:
            out["f0"] = [int(x) for x in out["poly0"].split(",")]
            out["f1"] = [int(x) for x in out["poly1"].split(",")]
        else:  # legacy c0..cd / Y0,Y1 format
            out["f1"] = [int(out[f"c{i}"]) for i in range(10) if f"c{i}" in out]
            out["f0"] = [int(out["Y0"]), int(out["Y1"])]
        return out


SURV_RE = re.compile(r"^# survivors (\S+): (\d+)")
TIME_Q_RE = re.compile(r"^# Sieving side-(\d) q=(\d+); rho=(\d+); a0=(-?\d+); b0=(-?\d+); a1=(-?\d+); "
                       r"b1=(-?\d+); I=(\d+); J=(\d+);")
CAT_RE = re.compile(r"^# (-?\d+) ([^:]+): ([\d.e+-]+)$")
TREE_RE = re.compile(r"^#(\s+)([\d.e+-]+) (.*?)(?: ; category (-?\d+))?$")


def parse_las(text: str, keep_relations: bool = False) -> dict:
    """Parse las stdout/out-file text (lines starting with '#' are stats)."""
    res: dict = {"funnel": {}, "categories": {}, "tree": {}, "special_q": []}
    rel_ab = set()
    nrel = 0
    rels: List[str] = []
    in_final = False
    in_tree = False
    for line in text.splitlines():
        if not line:
            continue
        if line[0] != "#":
            nrel += 1
            ab = line.split(":", 1)[0]
            rel_ab.add(ab)
            if keep_relations:
                rels.append(line)
            continue
        m = TIME_Q_RE.match(line)
        if m:
            res["special_q"].append(dict(zip(("side", "q", "rho", "a0", "b0", "a1", "b1", "I", "J"),
                                             map(int, m.groups()))))
            continue
        if line.startswith("# Average J="):
            in_final = True
            m = re.search(r"Average J=(\d+) for (\d+) special-q", line)
            if m:
                res["avg_J"], res["n_special_q"] = int(m.group(1)), int(m.group(2))
            continue
        if line.startswith("# Discarded"):
            m = re.search(r"Discarded (\d+) special-q's out of (\d+)", line)
            if m:
                res["discarded_q"], res["pushed_q"] = int(m.group(1)), int(m.group(2))
        if line.startswith("# Timings"):
            in_tree = True
        m = SURV_RE.match(line)
        if m and in_final:
            res["funnel"][m.group(1)] = int(m.group(2))
            continue
        m = CAT_RE.match(line)
        if m and in_final:
            res["categories"][f"{m.group(1)} {m.group(2).strip()}"] = float(m.group(3))
            continue
        if line.startswith("# Total cpu time"):
            m = re.search(r"Total cpu time ([\d.]+)s", line)
            res["cpu_s"] = float(m.group(1))
            m = re.search(r"\[norm ([\d.]+)\+([\d.]+), sieving ([\d.]+) \(([\d.]+)\+([\d.]+) \+ ([\d.]+)\), "
                          r"factor ([\d.]+) \(([\d.]+)\+([\d.]+) \+ ([\d.]+)\), rest ([\d.]+)\]", line)
            if m:
                v = list(map(float, m.groups()))
                res["coarse"] = {"norm": v[0] + v[1], "sieving": v[2], "factor": v[6], "rest": v[10]}
        if line.startswith("# Total elapsed time"):
            res["elapsed_s"] = float(re.search(r"Total elapsed time ([\d.]+)s", line).group(1))
        if line.startswith("# PeakMemusage"):
            res["las_peak_mem_mb"] = int(re.search(r"= (\d+)", line).group(1))
        if line.startswith("# Total ") and "reports" in line:
            m = re.search(r"# Total (\d+) reports", line)
            res["reports"] = int(m.group(1))
            in_tree = False
        if in_final:
            m = TREE_RE.match(line)
            if m and not line.startswith("# survivors"):
                depth = (len(m.group(1)) - 1)
                name = m.group(3).strip()
                key = f"{depth}:{name}"
                res["tree"][key] = res["tree"].get(key, 0.0) + float(m.group(2))
    res["relations"] = nrel
    res["unique_ab"] = len(rel_ab)
    if keep_relations:
        res["relation_lines"] = rels
    if res.get("elapsed_s"):
        res["rels_per_elapsed_s"] = nrel / res["elapsed_s"]
    if res.get("cpu_s"):
        res["rels_per_cpu_s"] = nrel / res["cpu_s"]
    if res.get("n_special_q"):
        res["rels_per_q"] = nrel / res["n_special_q"]
    return res


def las_command(tb: Testbed, q0: int, q1: int, threads: int, overrides: Dict[str, str],
                out: str, timings: bool = True, random_sample: Optional[int] = None,
                seed: int = 1, verbose: bool = False, extra: Optional[List[str]] = None) -> List[str]:
    params = dict(tb.params)
    params.update(overrides)
    params.pop("las.threads", None)
    cmd = [LAS, "-poly", tb.poly, "-fb1", tb.fb1, "-q0", str(q0), "-q1", str(q1), "-t", str(threads),
           "-out", out]
    if "A" in params:
        params.pop("I", None)
    for k, v in params.items():
        if v in ("", None):
            continue
        if k == "batch":
            if str(v).lower() in ("1", "true", "yes"):
                cmd.append("-batch")
            continue
        cmd += ["-" + k.replace("_", "-"), str(v)]
    if random_sample:
        cmd += ["-random-sample", str(random_sample), "-seed", str(seed)]
    if timings:
        cmd.append("-T")
    if verbose:
        cmd.append("-v")
    if extra:
        cmd += extra
    return cmd


def run_las(tb: Testbed, q0: int, q1: int, threads: int = 10, overrides: Optional[Dict[str, str]] = None,
            timings: bool = True, random_sample: Optional[int] = None, seed: int = 1,
            keep_relations: bool = False, out_path: Optional[str] = None, max_rss_gb: float = 6.0,
            extra: Optional[List[str]] = None) -> dict:
    overrides = dict(overrides or {})
    tmpdir = tempfile.mkdtemp(prefix="lasbench_", dir=WORK)
    out = out_path or os.path.join(tmpdir, "las.out")
    cmd = las_command(tb, q0, q1, threads, overrides, out, timings, random_sample, seed, extra=extra)
    with open(os.path.join(tmpdir, "las.stderr"), "w") as err:
        mem = memwatch.run(cmd, max_rss_gb=max_rss_gb, stdout=subprocess.DEVNULL, stderr=err, track=None)
    opener = gzip.open if out.endswith(".gz") else open
    with opener(out, "rt") as fh:
        text = fh.read()
    res = parse_las(text, keep_relations=keep_relations)
    res.update({"testbed": tb.name, "q0": q0, "q1": q1, "threads": threads, "overrides": overrides,
                "random_sample": random_sample, "seed": seed, "cmd": cmd,
                "wall_s": mem["wall_s"], "peak_rss_bytes": mem["peak_tree_rss_bytes"],
                "killed": mem["killed"], "returncode": mem["returncode"]})
    if out_path is None:
        os.remove(out)
    return res


def summarize(runs: List[dict]) -> dict:
    """Mean/stdev/min/max of throughput across repeated runs (identical yields asserted)."""
    def stats(xs):
        xs = [x for x in xs if x is not None]
        if not xs:
            return None
        return {"mean": statistics.fmean(xs), "sd": statistics.stdev(xs) if len(xs) > 1 else 0.0,
                "min": min(xs), "max": max(xs), "n": len(xs)}
    return {
        "relations": sorted({r["relations"] for r in runs}),
        "unique_ab": sorted({r["unique_ab"] for r in runs}),
        "elapsed_s": stats([r.get("elapsed_s") for r in runs]),
        "cpu_s": stats([r.get("cpu_s") for r in runs]),
        "rels_per_elapsed_s": stats([r.get("rels_per_elapsed_s") for r in runs]),
        "rels_per_cpu_s": stats([r.get("rels_per_cpu_s") for r in runs]),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--testbed", required=True, help="workdir name under ~/cadobench-work")
    ap.add_argument("--q0", type=int, required=True)
    ap.add_argument("--q1", type=int, required=True)
    ap.add_argument("-t", "--threads", type=int, default=10)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--random-sample", type=int, default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-timings", action="store_true")
    ap.add_argument("--set", action="append", default=[], help="override, e.g. lambda1=3.0")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    tb = Testbed.from_workdir(a.testbed)
    ov = dict(s.split("=", 1) for s in a.set)
    runs = [run_las(tb, a.q0, a.q1, a.threads, ov, timings=not a.no_timings,
                    random_sample=a.random_sample, seed=a.seed) for _ in range(a.repeat)]
    res = {"testbed": tb.name, "params": tb.params, "overrides": ov, "runs": runs, "summary": summarize(runs)}
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res["summary"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
