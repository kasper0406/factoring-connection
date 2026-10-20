"""Pre-registered runs for ptebarrier (notes/ptebarrier_prereg.md).

    .venv/bin/python -m src.ptebarrier.run_experiment all | E1 | E1b | E3
E1 also produces the E2 fiber statistics (same enumeration).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

from .enum import split_polys, group_by_top
from .rigidity import h_R, injection_bound, mul, from_roots, trim, evalp, int_roots

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CFG = os.path.join(ROOT, "experiments", "configs", "ptebarrier.json")
RAW = os.path.join(ROOT, "experiments", "raw")


def load_cfg():
    with open(CFG) as f:
        return json.load(f)


def _write(name, rows):
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, f"ptebarrier_{name}.jsonl")
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print("wrote", path, len(rows), "rows")


def run_E1(cfg):
    rows = []
    for B, Mmax in cfg["E1"]["grid"]:
        q = 2 * B + 1
        for M in range(2, Mmax + 1):
            t0 = time.time()
            entries = list(split_polys(B, M))
            for K in cfg["E1"]["K"]:
                if K >= M:
                    continue
                groups = group_by_top(entries, K)
                hmin = None
                nontriv = 0
                maxf = 0
                sizes_sum = 0
                disjoint_ok = True
                for key, mem in groups.items():
                    n = len(mem)
                    maxf = max(maxf, n)
                    sizes_sum += n
                    if n < 2:
                        continue
                    nontriv += 1
                    for i in range(n):
                        li = mem[i][1]
                        for k2 in range(i + 1, n):
                            lk = mem[k2][1]
                            hgt = max(abs(x - y) for x, y in zip(li, lk))
                            if hmin is None or hgt < hmin:
                                hmin = hgt
                            if K == 1 and set(mem[i][0]) & set(mem[k2][0]):
                                disjoint_ok = False
                hr = h_R(B, M, K, 0)
                s = M - K
                row = dict(B=B, M=M, K=K, s=s, h_min=hmin, h_R=hr,
                           P1_pass=(hmin is None or hmin > hr),
                           n_multisets=len(entries), n_fibers=len(groups),
                           max_fiber=maxf, mean_fiber=sizes_sum / len(groups),
                           n_nontrivial_fibers=nontriv,
                           injection_bound=injection_bound(q, K),
                           P3_pass=maxf <= injection_bound(q, K),
                           P4_pass=(disjoint_ok if K == 1 else None),
                           s_ge_sqrt2M_plus2=(s >= math.sqrt(2 * M) + 2))
                rows.append(row)
            print(f"E1 B={B} M={M} n={len(entries)} {time.time()-t0:.1f}s", flush=True)
            del entries
    _write("E1", rows)
    return rows


def _Qlist(c):
    out = []
    for c1 in range(-c, c + 1):
        for c0 in range(-c, c + 1):
            if (c0, c1) != (0, 0):
                out.append(trim([c0, c1]))
    return out


def run_E1b(cfg):
    rows = []
    Qs = _Qlist(cfg["E1b"]["Qcoef"])
    for B, Mmax in cfg["E1b"]["grid"]:
        for M in range(2, Mmax + 1):
            t0 = time.time()
            polys = [list(co) + [1] for _, co in split_polys(B, M)]
            for K in cfg["E1b"]["K"]:
                groups: dict = {}
                for ip, P in enumerate(polys):
                    for iq, Q in enumerate(Qs):
                        R = mul(Q, P)
                        key = (len(R),) + tuple(R[K:])
                        groups.setdefault(key, []).append(tuple(R[:K]))
                hmin = None
                for mem in groups.values():
                    if len(mem) < 2:
                        continue
                    uniq = list(set(mem))
                    for i in range(len(uniq)):
                        for k2 in range(i + 1, len(uniq)):
                            hgt = max(abs(x - y) for x, y in zip(uniq[i], uniq[k2]))
                            if hmin is None or hgt < hmin:
                                hmin = hgt
                hr = h_R(B, M, K, 1)
                rows.append(dict(B=B, M=M, K=K, j=1, h_min=hmin, h_R=hr,
                                 P1b_pass=(hmin is None or hmin > hr), n_Q=len(Qs)))
                del groups
            print(f"E1b B={B} M={M} {time.time()-t0:.1f}s", flush=True)
    _write("E1b", rows)
    return rows


def split_over_Z(p):
    """Integer roots with multiplicity if p (monic) splits over Z, else None."""
    p = trim(p)
    roots = []
    while len(p) > 1:
        if p[0] == 0:
            roots.append(0)
            p = p[1:]
            continue
        found = None
        for r in sorted(int_roots(p), key=abs):
            found = r
            break
        if found is None:
            return None
        # deflate by (X - found)
        n = len(p) - 1
        qd = [0] * n
        acc = p[-1]
        for i in range(n - 1, -1, -1):
            qd[i] = acc
            acc = p[i] + acc * found
        assert acc == 0
        p = trim(qd)
        roots.append(found)
    return sorted(roots)


def run_E3(cfg):
    F = cfg["E3"]["F"]  # low->high
    R = cfg["E3"]["range"]
    rows = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            p = list(F)
            p[1] += a
            p[0] += b
            rts = split_over_Z(p)
            if rts is not None:
                rows.append(dict(a=a, b=b, roots=rts))
    _write("E3", rows)
    print("E3 split members:", rows)
    return rows


def main(argv):
    cfg = load_cfg()
    what = argv[1] if len(argv) > 1 else "all"
    if what in ("all", "E1"):
        run_E1(cfg)
    if what in ("all", "E1b"):
        run_E1b(cfg)
    if what in ("all", "E3"):
        run_E3(cfg)


if __name__ == "__main__":
    main(sys.argv)


def run_E1c_posthoc(cfg):
    """POST-HOC (disclosed): E1 cells were mostly vacuous (no PTE pairs at all for B <= 5),
    so repeat the R0 check at larger B and small M, where PTE pairs exist."""
    rows = []
    for B, Mmax in cfg["E1c_posthoc"]["grid"]:
        for M in range(3, Mmax + 1):
            entries = list(split_polys(B, M))
            for K in cfg["E1c_posthoc"]["K"]:
                if K >= M:
                    continue
                groups = group_by_top(entries, K)
                hmin, npairs = None, 0
                for mem in groups.values():
                    n = len(mem)
                    for i in range(n):
                        for k2 in range(i + 1, n):
                            npairs += 1
                            hgt = max(abs(x - y) for x, y in zip(mem[i][1], mem[k2][1]))
                            if hmin is None or hgt < hmin:
                                hmin = hgt
                hr = h_R(B, M, K, 0)
                rows.append(dict(B=B, M=M, K=K, h_min=hmin, h_R=hr, n_pairs=npairs,
                                 P1_pass=(hmin is None or hmin > hr)))
                print(rows[-1], flush=True)
            del entries
    _write("E1c_posthoc", rows)
    return rows


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "E1c":
    run_E1c_posthoc(load_cfg())
