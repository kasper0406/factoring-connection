"""Relation-set statistics from a finished CADO run (diagnostic harness).

* special-q lattices and per-q yields, parsed from the raw las outputs in <workdir>/<name>.upload/;
* unique relations (dedup by (a, b)), large-ideal incidence histograms, and an independent
  singleton-removal pass (V6: assumption T of heurnfs §3.1);
* filtering / matrix numbers parsed from the CADO log (V5).

Ideals: side 0 (rational) -> p; side 1 (algebraic) -> (p, r) with r = a/b mod p (r = p when p | b).
"Large" means p > lim_side. Our purge removes singleton relations by large ideals only, iterated to a
fixed point; CADO's purge also removes cliques to tune the excess, so our counts are an upper bound
on CADO's post-purge set (disclosed in the notes).
"""
from __future__ import annotations

import gzip
import os
import re
from collections import Counter
from typing import Dict, Iterator, List, Optional, Tuple

from .smoothcount import QLattice

QHDR = re.compile(r"^# Sieving side-(\d) q=(\d+); rho=(\d+); a0=(-?\d+); b0=(-?\d+); a1=(-?\d+); b1=(-?\d+); "
                  r"I=(\d+); J=(\d+);")
NREL = re.compile(r"^# (\d+) relation\(s\)")


def upload_files(workdir: str) -> List[str]:
    """The relation files CADO accepted ("Found k relations in '<path>'"), in acceptance order.

    Workunits that finish after the relation target is reached are left in the upload directory but
    not used by CADO; they are excluded so that our counts match CADO's."""
    log = [f for f in os.listdir(workdir) if f.endswith(".log")][0]
    out, seen = [], set()
    for m in re.finditer(r"Found \d+ relations in '([^']+)'", open(os.path.join(workdir, log), errors="replace").read()):
        if m.group(1) not in seen and os.path.exists(m.group(1)):
            seen.add(m.group(1))
            out.append(m.group(1))
    return out


def iter_raw(workdir: str) -> Iterator[Tuple[str, str]]:
    """Yields ('q', header line) / ('n', '# k relation(s)') / ('r', relation line), file by file."""
    for path in upload_files(workdir):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    if line.startswith("# Sieving side-"):
                        yield "q", line
                    elif NREL.match(line):
                        yield "n", line
                else:
                    yield "r", line.rstrip("\n")


def qlattices(workdir: str, max_q: Optional[int] = None) -> List[Tuple[QLattice, int]]:
    """(lattice, #relations reported for that q) for the special-q sieved in the run.

    las prints '# Sieving side-1 q=...' then later '# k relation(s)' for that q; with several threads
    the per-q blocks are still emitted whole (las buffers per special-q), which we rely on.
    """
    out: List[Tuple[QLattice, int]] = []
    cur = None
    for kind, line in iter_raw(workdir):
        if kind == "q":
            m = QHDR.match(line)
            g = list(map(int, m.groups()))
            cur = QLattice(q=g[1], rho=g[2], a0=g[3], b0=g[4], a1=g[5], b1=g[6], I=g[7], J=g[8])
        elif kind == "n" and cur is not None:
            out.append((cur, int(NREL.match(line).group(1))))
            cur = None
            if max_q and len(out) >= max_q:
                break
    return out


def parse_rel(line: str) -> Tuple[int, int, List[int], List[int]]:
    ab, s0, s1 = line.split(":")
    a, b = ab.split(",")
    p0 = [int(x, 16) for x in s0.split(",")] if s0 else []
    p1 = [int(x, 16) for x in s1.split(",")] if s1 else []
    return int(a), int(b), p0, p1


def unique_relations(workdir: str) -> Tuple[Dict[Tuple[int, int], Tuple[List[int], List[int]]], int]:
    rels: Dict[Tuple[int, int], Tuple[List[int], List[int]]] = {}
    raw = 0
    for kind, line in iter_raw(workdir):
        if kind != "r":
            continue
        raw += 1
        a, b, p0, p1 = parse_rel(line)
        if (a, b) not in rels:
            rels[(a, b)] = (p0, p1)
    return rels, raw


def large_ideals(a: int, b: int, p0: List[int], p1: List[int], lim0: int, lim1: int):
    out = set()
    for p in p0:
        if p > lim0:
            out.add((0, p, 0))
    for p in p1:
        if p > lim1:
            r = p if b % p == 0 else (a * pow(b, -1, p)) % p
            out.add((1, p, r))
    return out


def histogram(counter: Counter) -> Dict[int, int]:
    h = Counter(counter.values())
    return dict(sorted(h.items()))


def _enc(side: int, p: int, r: int) -> int:
    return (r << 41) | (p << 1) | side


def _side(code: int) -> int:
    return code & 1


def incidence_stats(workdir: str, lim0: int, lim1: int) -> dict:
    """Large-ideal incidence before/after singleton removal (numpy; ~16 bytes per incidence).

    Ideal code: side 0 -> p << 1; side 1 -> ((p << 28 | r) << 1) | 1 (p, r < 2^28, fits int64).
    """
    import numpy as np
    ab_a: List[int] = []
    ab_b: List[int] = []
    inc_rel: List[int] = []
    inc_code: List[int] = []
    raw = 0
    for kind, line in iter_raw(workdir):
        if kind != "r":
            continue
        raw += 1
        a, b, p0, p1 = parse_rel(line)
        k = len(ab_a)
        ab_a.append(a)
        ab_b.append(b)
        codes = set()
        for p in p0:
            if p > lim0:
                codes.add(p << 1)
        for p in p1:
            if p > lim1:
                r = p if b % p == 0 else (a * pow(b, -1, p)) % p
                codes.add((((p << 28) | r) << 1) | 1)
        for c in codes:
            inc_rel.append(k)
            inc_code.append(c)
    A = np.array(ab_a, dtype=np.int64)
    Bv = np.array(ab_b, dtype=np.int64)
    del ab_a, ab_b
    _, first = np.unique(np.stack([A, Bv], axis=1), axis=0, return_index=True)
    keep_rel = np.zeros(len(A), dtype=bool)
    keep_rel[first] = True
    nuniq = int(keep_rel.sum())
    R = np.array(inc_rel, dtype=np.int64)
    C = np.array(inc_code, dtype=np.int64)
    del inc_rel, inc_code
    m = keep_rel[R]
    R, C = R[m], C[m]
    codes, cidx = np.unique(C, return_inverse=True)
    side = (codes & 1).astype(np.int8)
    nlp = np.bincount(np.bincount(R, minlength=len(A))[keep_rel])

    def summary(alive_inc: "np.ndarray") -> dict:
        cnt = np.bincount(cidx[alive_inc], minlength=len(codes))
        present = cnt > 0
        h = np.bincount(cnt[present])
        tot = int(present.sum())
        occ = int(cnt.sum())
        hist = {int(k): int(v) for k, v in enumerate(h) if v and k > 0}
        by_side = {}
        for sd in (0, 1):
            cs = cnt[present & (side == sd)]
            n = len(cs)
            by_side[sd] = {"ideals": n, "frac_exactly_2": float((cs == 2).sum() / n) if n else None,
                           "mean_occ": float(cs.mean()) if n else None}
        return {"ideals": tot, "incidences": occ, "mean_occ": occ / tot if tot else None,
                "frac_exactly_1": hist.get(1, 0) / tot if tot else None,
                "frac_exactly_2": hist.get(2, 0) / tot if tot else None,
                "frac_ge_3": sum(v for k, v in hist.items() if k >= 3) / tot if tot else None,
                "hist_head": {k: hist[k] for k in sorted(hist)[:12]}, "max_occ": max(hist) if hist else 0,
                "by_side": by_side}

    alive_inc = np.ones(len(R), dtype=bool)
    before = summary(alive_inc)
    rel_alive = keep_rel.copy()
    it = 0
    while True:
        it += 1
        cnt = np.bincount(cidx[alive_inc], minlength=len(codes))
        single_inc = alive_inc & (cnt[cidx] == 1)
        dead = np.unique(R[single_inc])
        if len(dead) == 0:
            break
        rel_alive[dead] = False
        alive_inc &= rel_alive[R]
    after = summary(alive_inc)
    return {"raw_relations": raw, "unique_relations": nuniq, "dup_ratio": nuniq / raw if raw else None,
            "large_primes_per_relation": {int(k): int(v) for k, v in enumerate(nlp) if v},
            "all_unique": before, "after_singleton_removal": {"relations": int(rel_alive.sum()), "iterations": it, **after}}


# ---------------------------------------------------------------- filtering/matrix numbers from the log

def filtering_numbers(workdir: str) -> dict:
    log = [f for f in os.listdir(workdir) if f.endswith(".log")][0]
    text = open(os.path.join(workdir, log), errors="replace").read()
    out: dict = {}
    m = re.findall(r"(\d+) unique relations remain in total", text)
    if m:
        out["unique_relations"] = int(m[-1])
    m = re.findall(r"Total number of relations: (\d+)", text)
    if m:
        out["raw_relations"] = int(m[-1])
    m = re.findall(r"After purge, (\d+) relations with (\d+) primes remain with excess (-?\d+)", text)
    if m:
        out["purge_rows"], out["purge_cols"], out["purge_excess"] = map(int, m[-1])
    out["purge_rounds"] = len(m)
    m = re.findall(r"Found (\d+) free relations", text)
    if m:
        out["free_relations"] = int(m[-1])
    base = log[:-4]
    mp = os.path.join(workdir, base + ".merge.stdout.1")
    if os.path.exists(mp):
        mm = re.findall(r"Final matrix has N=(\d+) nc=(\d+) \((\d+)\) W=(\d+)", open(mp).read())
        if mm:
            out["matrix_N"], out["matrix_nc"], _, out["matrix_W"] = map(int, mm[-1])
    m = re.search(r"Krylov: CPU time ([\d.]+), WCT time ([\d.]+).*?\((\d+) iterations\)", text)
    if m:
        out["krylov_cpu"], out["krylov_wct"], out["krylov_iter"] = float(m.group(1)), float(m.group(2)), int(m.group(3))
    m = re.search(r"Mksol: CPU time ([\d.]+),\s+WCT time ([\d.]+)", text)
    if m:
        out["mksol_cpu"], out["mksol_wct"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"Lingen CPU time ([\d.]+), WCT time ([\d.]+)", text)
    if m:
        out["lingen_cpu"], out["lingen_wct"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"Lattice Sieving: Total time: ([\d.e+]+)s", text)
    if m:
        out["las_cpu_total"] = float(m.group(1))
    m = re.search(r"Average J: ([\d.]+) for (\d+) special-q", text)
    if m:
        out["avg_J"], out["n_special_q"] = float(m.group(1)), int(m.group(2))
    return out
