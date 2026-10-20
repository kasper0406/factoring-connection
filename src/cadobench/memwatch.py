"""Memory watchdog for heavy CADO-NFS jobs (harness; not operational code).

Runs a command in its own process group and samples once per second:
  * the summed RSS of the whole process tree, plus the peak RSS per command name
    (CADO runs one binary per phase: polyselect, las, dup1/dup2/purge/merge, bwc's
    prep/secure/krylov/mksol/lingen, sqrt), so "peak RSS by command" is the per-phase peak;
  * the system memory-pressure level (kern.memorystatus_vm_pressure_level:
    1 normal, 2 warn, 4 critical) and the free percentage from `memory_pressure -Q`.

It kills the whole process group (SIGTERM, then SIGKILL) when
  * the tree RSS exceeds --max-rss-gb (default 8), or
  * the "available" memory (free + inactive + speculative + purgeable pages) drops below
    --min-avail-gb (default 2), or
  * the pressure level reaches warn (2) or critical (4) for --pressure-grace samples in a row; if the
    system is already at warn when the job starts (baseline load from other processes), only critical
    or a rise above the starting level counts.

Usage:
  python -m cadobench.memwatch --json out.json [--max-rss-gb 8] -- cmd args...
The JSON summary records wall time, exit status, peak tree RSS, per-command peaks and the
reason for any kill.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import signal
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

PAGE = 16384


def _sysctl_int(name: str) -> Optional[int]:
    try:
        out = subprocess.run(["sysctl", "-n", name], capture_output=True, text=True, timeout=5)
        return int(out.stdout.strip())
    except Exception:
        return None


def available_bytes() -> Optional[int]:
    """free + inactive + speculative + purgeable pages, in bytes (macOS vm_stat)."""
    try:
        out = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=5).stdout
    except Exception:
        return None
    m = re.search(r"page size of (\d+) bytes", out)
    page = int(m.group(1)) if m else PAGE
    vals = {}
    for key in ("Pages free", "Pages inactive", "Pages speculative", "Pages purgeable"):
        mm = re.search(re.escape(key) + r":\s+(\d+)\.", out)
        vals[key] = int(mm.group(1)) if mm else 0
    return sum(vals.values()) * page


# CADO's clients run with --daemon and detach from the tree, so processes whose command line
# matches this pattern are tracked (and killed) too.
DEFAULT_TRACK = r"cadobench-work/|cado-nfs/build/"


def _ps() -> List[Tuple[int, int, int, str]]:
    out = subprocess.run(["ps", "-axww", "-o", "pid=,ppid=,rss=,command="], capture_output=True,
                         text=True).stdout
    rows = []
    for line in out.splitlines():
        parts = line.split(None, 3)
        if len(parts) < 4:
            continue
        rows.append((int(parts[0]), int(parts[1]), int(parts[2]) * 1024, parts[3]))
    return rows


def _cmdname(command: str) -> str:
    """Basename of argv[0]; for python scripts, the script basename."""
    toks = command.split()
    name = os.path.basename(toks[0]) if toks else "?"
    if name.lower().startswith("python") and len(toks) > 1 and not toks[1].startswith("-"):
        name = os.path.basename(toks[1])
    return name


def tracked_pids(root_pid: int, track: Optional[str] = DEFAULT_TRACK) -> Dict[int, Tuple[int, str]]:
    """pid -> (rss bytes, command name) for the tree of root_pid plus processes matching `track`."""
    rows = _ps()
    children: Dict[int, List[int]] = {}
    info: Dict[int, Tuple[int, str]] = {}
    for pid, ppid, rss, command in rows:
        info[pid] = (rss, command)
        children.setdefault(ppid, []).append(pid)
    me = os.getpid()
    sel = set()
    stack = [root_pid]
    while stack:
        p = stack.pop()
        if p in sel:
            continue
        sel.add(p)
        stack.extend(children.get(p, []))
    if track:
        rx = re.compile(track)
        for pid, (_, command) in info.items():
            if pid != me and rx.search(command) and "memwatch" not in command:
                sel.add(pid)
    return {p: (info[p][0], _cmdname(info[p][1])) for p in sel if p in info}


def tree_rss(root_pid: int, track: Optional[str] = DEFAULT_TRACK) -> Tuple[int, Dict[str, int]]:
    """Summed RSS (bytes) of the tracked processes, and RSS by command name."""
    tp = tracked_pids(root_pid, track)
    total = 0
    by_cmd: Dict[str, int] = {}
    for rss, cmd in tp.values():
        total += rss
        by_cmd[cmd] = by_cmd.get(cmd, 0) + rss
    return total, by_cmd


def run(cmd: List[str], max_rss_gb: float = 8.0, min_avail_gb: float = 2.0, pressure_grace: int = 5,
        interval: float = 1.0, log=None, cwd: Optional[str] = None, env=None,
        stdout=None, stderr=None, track: Optional[str] = DEFAULT_TRACK) -> dict:
    t0 = time.time()
    proc = subprocess.Popen(cmd, start_new_session=True, cwd=cwd, env=env, stdout=stdout, stderr=stderr)
    peak_total = 0
    peak_by_cmd: Dict[str, int] = {}
    min_avail = None
    max_level = 1
    bad_pressure = 0
    killed = None
    samples = 0
    # Baseline-aware pressure rule: when the system is already at "warn" before the job starts (other
    # agents / Docker), only a rise above that level, or "critical", counts.
    start_level = _sysctl_int("kern.memorystatus_vm_pressure_level") or 1
    while proc.poll() is None:
        time.sleep(interval)
        samples += 1
        total, by_cmd = tree_rss(proc.pid, track)
        peak_total = max(peak_total, total)
        for k, v in by_cmd.items():
            peak_by_cmd[k] = max(peak_by_cmd.get(k, 0), v)
        avail = available_bytes()
        if avail is not None:
            min_avail = avail if min_avail is None else min(min_avail, avail)
        level = _sysctl_int("kern.memorystatus_vm_pressure_level") or 1
        max_level = max(max_level, level)
        bad_pressure = bad_pressure + 1 if (level >= 4 or (level >= 2 and level > start_level)) else 0
        reason = None
        if total > max_rss_gb * 2**30:
            reason = f"tree RSS {total / 2**30:.2f} GB > {max_rss_gb} GB"
        elif avail is not None and avail < min_avail_gb * 2**30:
            reason = f"available memory {avail / 2**30:.2f} GB < {min_avail_gb} GB"
        elif bad_pressure >= pressure_grace:
            reason = f"memory pressure level {level} for {bad_pressure} samples"
        if log is not None and samples % 30 == 0:
            print(f"[memwatch] t={time.time() - t0:.0f}s rss={total / 2**30:.2f}GB "
                  f"avail={(avail or 0) / 2**30:.2f}GB level={level}", file=log, flush=True)
        if reason:
            killed = reason
            print(f"[memwatch] KILL: {reason}", file=sys.stderr, flush=True)
            victims = list(tracked_pids(proc.pid, track))
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(proc.pid, sig)
                except OSError:  # group already gone (macOS reports EPERM for a reaped group)
                    pass
                for v in victims:
                    try:
                        os.kill(v, sig)
                    except OSError:
                        pass
                time.sleep(5)
            break
    rc = proc.wait()
    return {
        "cmd": cmd,
        "returncode": rc,
        "wall_s": time.time() - t0,
        "peak_tree_rss_bytes": peak_total,
        "peak_rss_by_cmd_bytes": dict(sorted(peak_by_cmd.items(), key=lambda kv: -kv[1])),
        "min_available_bytes": min_avail,
        "max_pressure_level": max_level,
        "start_pressure_level": start_level,
        "killed": killed,
        "samples": samples,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--max-rss-gb", type=float, default=8.0)
    ap.add_argument("--min-avail-gb", type=float, default=2.0)
    ap.add_argument("--pressure-grace", type=int, default=5)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args(argv)
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    res = run(cmd, a.max_rss_gb, a.min_avail_gb, a.pressure_grace, log=sys.stderr)
    with open(a.json, "w") as fh:
        json.dump(res, fh, indent=1)
    return 0 if res["returncode"] == 0 and not res["killed"] else 1


if __name__ == "__main__":
    sys.exit(main())
