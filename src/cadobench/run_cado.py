"""Run a full CADO-NFS factorization under the memory watchdog and extract per-phase timings.

Harness code (it knows N's factors only to check the output; nothing here is operational).

Phase wall times are reconstructed from the timestamped `<name>.log` that cado-nfs.py writes in its
workdir. Each log line carries a phase name ("Lattice Sieving", "Filtering - Merging", ...). The time
between two consecutive lines is charged to
  * the previous phase, if the new line is "<phase>: Starting";
  * otherwise to the new line's phase.
Book-keeping "phases" (Client Launcher, API server, HTTP server, Command, root) are transparent.
Phases are grouped as polyselect / sieve / filter / linalg / sqrt / other, matching the task.

Usage:
  python -m cadobench.run_cado --digits 70 --threads 10 --out experiments/raw/cadobench_c70_run1.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional

from . import memwatch
from .numbers import NUMBERS, factors_ok

CADO_SRC = os.path.expanduser("~/src/cado-nfs")
CADO_BUILD = os.path.join(CADO_SRC, "build", "kmini")
CADO_VENV_BIN = os.path.join(CADO_SRC, "cado-nfs.venv", "bin")
WORK = os.path.expanduser("~/cadobench-work")

TRANSPARENT = {"Client Launcher", "API server", "HTTP server", "Command", "root", "HTTP Server",
               "Complete Factorization / Discrete logarithm / Class group"}

GROUPS = [
    ("polyselect", ("Polynomial Selection",)),
    ("sieve", ("Lattice Sieving",)),
    ("filter", ("Filtering",)),
    ("linalg", ("Linear Algebra", "Quadratic Characters")),
    ("sqrt", ("Square Root",)),
]

LINE_RE = re.compile(r"^PID\d+ (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d,\d{3}) (\w+):([^:]+): (.*)$")


def group_of(phase: str) -> str:
    for g, prefixes in GROUPS:
        if any(phase.startswith(p) for p in prefixes):
            return g
    return "other"


def parse_log(path: str) -> dict:
    by_phase: Dict[str, float] = {}
    reported: Dict[str, str] = {}
    cur: Optional[str] = None
    t_prev = t_first = t_last = None
    factors: List[int] = []
    with open(path, errors="replace") as fh:
        for line in fh:
            m = LINE_RE.match(line.rstrip("\n"))
            if not m:
                continue
            ts = dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S,%f").timestamp()
            level, phase, msg = m.group(2), m.group(3).strip(), m.group(4)
            if t_first is None:
                t_first = t_prev = ts
            t_last = ts
            if level == "Info" and ("Total" in msg or "time" in msg.lower()):
                reported[f"{phase}: {msg.split(':')[0]}"] = msg
            if phase == "Square Root" and msg.startswith("Factors:"):
                factors = [int(x) for x in msg.split(":", 1)[1].split()]
            if phase in TRANSPARENT:
                continue
            dtsec = ts - t_prev
            if msg == "Starting":
                if cur is not None:
                    by_phase[cur] = by_phase.get(cur, 0.0) + dtsec
                cur = phase
            else:
                by_phase[phase] = by_phase.get(phase, 0.0) + dtsec
                cur = phase
            t_prev = ts
    groups: Dict[str, float] = {}
    for ph, s in by_phase.items():
        groups[group_of(ph)] = groups.get(group_of(ph), 0.0) + s
    return {
        "wall_total_s": (t_last - t_first) if t_first is not None else None,
        "wall_by_phase_s": by_phase,
        "wall_by_group_s": groups,
        "cado_reported": reported,
        "factors": factors,
    }


def sieve_params(workdir: str) -> dict:
    """Key sieve parameters from the parameters snapshot."""
    snaps = sorted(f for f in os.listdir(workdir) if ".parameters_snapshot." in f)
    out = {}
    if not snaps:
        return out
    keys = re.compile(r"^tasks\.(I|A|lim0|lim1|lpb0|lpb1|qmin|sieve\.(mfb0|mfb1|lambda0|lambda1|ncurves0|ncurves1|"
                      r"qrange|rels_wanted|batch|las\.threads)|polyselect\.(degree|admax|P|nq)|linalg\..*|threads)\s*=")
    with open(os.path.join(workdir, snaps[-1])) as fh:
        for line in fh:
            if keys.match(line):
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
    return out


def cado_env() -> dict:
    env = dict(os.environ)
    env["PATH"] = CADO_VENV_BIN + ":" + env.get("PATH", "")
    return env


def run(digits: int, threads: int, tag: str, extra: List[str], max_rss_gb: float) -> dict:
    N, p, q, seed, bits = NUMBERS[digits]
    workdir = os.path.join(WORK, f"c{digits}_{tag}")
    os.makedirs(WORK, exist_ok=True)
    cmd = [os.path.join(CADO_BUILD, "cado-nfs.py"), str(N), "-t", str(threads), "--workdir", workdir] + extra
    with open(workdir + ".stdout", "w") as out:
        mem = memwatch.run(cmd, max_rss_gb=max_rss_gb, env=cado_env(), stdout=out, stderr=subprocess.STDOUT,
                           log=sys.stderr)
    logs = [f for f in os.listdir(workdir) if f.endswith(".log")]
    parsed = parse_log(os.path.join(workdir, logs[0])) if logs else {}
    res = {
        "digits": digits, "N": str(N), "seed": seed, "bits": bits, "threads": threads, "tag": tag,
        "extra_args": extra, "workdir": workdir, "memwatch": mem, "params": sieve_params(workdir),
        **parsed,
    }
    res["factored_ok"] = factors_ok(digits, parsed.get("factors", []))  # harness check only
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--digits", type=int, required=True)
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--tag", default="run1")
    ap.add_argument("--max-rss-gb", type=float, default=8.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("extra", nargs="*", help="extra cado-nfs.py args, e.g. tasks.sieve.lambda0=2.0")
    a = ap.parse_args(argv)
    res = run(a.digits, a.threads, a.tag, a.extra, a.max_rss_gb)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res.get(k) for k in ("digits", "wall_total_s", "wall_by_group_s", "factored_ok")},
                     indent=1))
    return 0 if res["factored_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
