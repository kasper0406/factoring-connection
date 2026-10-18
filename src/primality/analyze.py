"""HARNESS: turn experiments/raw/primality_* into experiments/analysis/primality_tables.md.
Run: .venv/bin/python -m src.primality.analyze
"""
from __future__ import annotations

import collections
import json
import math

import numpy as np

from src.common.instances import balanced_semiprime
from src.primality import formulas as F
from src.primality import ops

RAW = "experiments/raw/"
OUT = "experiments/analysis/primality_tables.md"


def jl(name):
    import gzip
    import os
    path = RAW + name
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        return [json.loads(l) for l in gzip.open(path + ".gz", "rt")]
    return [json.loads(l) for l in open(path)]


def leakrates(lines):
    R = jl("primality_leakrates.jsonl")
    g = collections.defaultdict(list)
    for r in R:
        if r["fam"] != "chernick":
            g[(r["fam"], r["bits"])].append(r)
    lines += ["## E1 leak rates (per unit base): expected (exact formula, summed over instances x bases) vs observed", "",
              "| family | bits | instances x bases | Fermat exp/obs | MR prefix exp/obs | MR sqrt exp/obs | SS exp/obs | median log2 MR-prefix rate | mean g |",
              "|---|---:|---:|---|---|---|---|---:|---:|"]
    for k in sorted(g):
        v = g[k]
        nb = v[0]["bases"]
        cells = []
        for key, okey in [("fermat", "fermat"), ("mr_prefix", "mr_prefix"), ("mr_sqrt", "mr_sqrt"), ("ss", "ss")]:
            e = sum(x["pred"][key] for x in v) * nb
            o = sum(x["obs"][okey] for x in v)
            cells.append(f"{e:.3g} / {o}")
        med = np.median([math.log2(x["pred"]["mr_prefix"]) for x in v])
        lines.append(f"| {k[0]} | {k[1]} | {len(v)} x {nb} | " + " | ".join(cells) + f" | {med:.1f} | {np.mean([x['g'] for x in v]):.3g} |")
    # exact pm (MR incl. b_t + 1) and P-sampled strong Lucas at 32 bits
    v = g[("balanced", 32)]
    e_pm = e_lu = 0.0
    e_bpsw = 0.0
    for x in v:
        inst = balanced_semiprime(32, x["seed"])
        p, q = inst.labels["factors"]
        e_pm += float(F.mr_exact([p, q])["pm"])
        sel = ops.selfridge(inst.N)
        if "P" in sel:
            lu = float(F.lucas_exact([p, q], sel["D"], True)["pm"])
            e_lu += lu
    nb = v[0]["bases"]
    lines += ["", f"32-bit balanced, all gcd(b_i +- 1, N) checks (MR 'pm'): expected {e_pm * nb:.1f}, observed "
              f"{sum(x['obs']['mr_pm'] for x in v)}; strong Lucas with U and V checks, P-sampled: expected {e_lu * nb:.1f}, "
              f"observed {sum(x['obs']['lucas_pm'] for x in v)}. BPSW (deterministic, one run per N): leaks on "
              f"{sum(len(x['bpsw_leak']) > 0 for x in v)} of {len(v)} N; heuristic expectation (MR-pm + Lucas-pm rate, "
              f"treating base 2 and P=1 as random) {e_pm + e_lu:.1f}."]
    C = [r for r in R if r["fam"] == "chernick"]
    lines += ["", "Chernick Carmichael numbers (3 primes), MR prefix factor rate per unit base:", "",
              "| kbits | median bits | expected | observed | bases | min exact rate |", "|---:|---:|---:|---:|---:|---:|"]
    for kb in sorted({r["kbits"] for r in C}):
        w = [r for r in C if r["kbits"] == kb]
        lines.append(f"| {kb} | {np.median([x['bits'] for x in w]):.0f} | {sum(x['pred_mr_prefix'] * x['bases'] for x in w):.1f} | "
                     f"{sum(x['obs'] for x in w)} | {sum(x['bases'] for x in w)} | {min(x['pred_mr_prefix'] for x in w):.4f} |")
    lines.append("")


def kangaroo(lines):
    R = jl("primality_kangaroo.jsonl")
    g = collections.defaultdict(list)
    for r in R:
        g[(r["kind"], r["method"], r["bits"])].append(r)
    lines += ["## E2 factoring from the residue alone (interval DLP), group multiplications", "",
              "| residue | method | bits | success | median mults | mean mults | N^(1/4) |", "|---|---|---:|---:|---:|---:|---:|"]
    for k in sorted(g):
        v = g[k]
        lines.append(f"| {k[0]} | {k[1]} | {k[2]} | {sum(x['ok'] for x in v)}/{len(v)} | {np.median([x['grp_mul'] for x in v]):.0f} | "
                     f"{np.mean([x['grp_mul'] for x in v]):.0f} | {np.median([x['N4'] for x in v]):.0f} |")
    lines.append("")
    for kind, meth in [("fermat", "kangaroo"), ("fermat", "bsgs"), ("lucas", "kangaroo")]:
        xs, ys = [], []
        for k in g:
            if k[:2] == (kind, meth) and k[2] >= 24:
                xs.append(k[2])
                ys.append(math.log2(np.mean([x["grp_mul"] for x in g[k]])))
        sl = np.polyfit(xs, ys, 1)
        lines.append(f"* fit log2(mean mults) vs n, n >= 24, {kind}/{meth}: slope {sl[0]:.3f} (N^(1/4) predicts 0.25)")
    lines.append("")


def iterate(lines):
    R = jl("primality_iterate.jsonl")
    lines += ["## E3 iterated residue a^((s-2)^2) vs plain a^(s-2)", "",
              "| bits | plain kangaroo median | sq-kangaroo found / split / of | sq-kangaroo median mults | Fermat-scan median | BSGS Fermat | BSGS Lucas |",
              "|---:|---:|---|---:|---:|---:|---:|"]
    xs, ys = [], []
    for b in sorted({r["bits"] for r in R}):
        v = [r for r in R if r["bits"] == b]
        lines.append(f"| {b} | {np.median([r['plain_mul'] for r in v]):.0f} | {sum(r['sq_kangaroo_found'] for r in v)} / "
                     f"{sum(r['sq_kangaroo_ok'] for r in v)} / {len(v)} | {np.median([r['sq_kangaroo_mul'] for r in v]):.0f} | "
                     f"{np.median([r['sq_scan_mul'] for r in v]):.0f} | {np.median([r['bsgs_fermat_mul'] for r in v]):.0f} | "
                     f"{np.median([r['bsgs_lucas_mul'] for r in v]):.0f} |")
        xs.append(b)
        ys.append(math.log2(np.median([r["sq_kangaroo_mul"] for r in v])))
    lines += ["", f"* slope of log2(median sq-kangaroo mults) vs n: {np.polyfit(xs, ys, 1)[0]:.3f}; no run recovered s "
              f"(sqrt check) -- every success came from a collision a^x = a^((N-1)^2), x != (N-1)^2, i.e. order finding.", ""]


def aks(lines):
    R = jl("primality_aks.jsonl")
    lines += ["## E5 AKS (64 values of a per N)", "",
              "| bits | median r | N | structural (whole poly) | trial-division | Kummer regime N | hits | predicted 64(r - cov) | generic N | hits | null r*64*(1/p+1/q) |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for b in sorted({r["bits"] for r in R}):
        v = [r for r in R if r["bits"] == b]
        rows = []
        for r in v:
            inst = balanced_semiprime(b, r["seed"])
            p, q = inst.labels["factors"]
            rr = r["r"]
            d = q - p
            cov = len({i % rr for i in range(min(d + 1, rr))} | {(p + i) % rr for i in range(min(d + 1, rr))}) if q < 2 * p else rr
            rows.append((r["coeff_hits"], r["a_checked"] * (rr - cov), r["null_coeff"], cov < rr))
        k = [x for x in rows if x[3]]
        nk = [x for x in rows if not x[3]]
        lines.append(f"| {b} | {np.median([r['r'] for r in v]):.0f} | {len(v)} | {sum(r['structural'] for r in v)} | "
                     f"{sum(r['trial'] for r in v)} | {len(k)} | {sum(x[0] for x in k)} | {sum(x[1] for x in k)} | {len(nk)} | "
                     f"{sum(x[0] for x in nk)} | {sum(x[2] for x in nk):.1f} |")
    lines.append("")


def epscert(lines):
    R = jl("primality_epscert.jsonl")
    lines += ["## E4a one-sided eps = +1 certificates on D_n (N = 1 mod 4), per instance with eps = +1", "",
              "| bits | eps=+1 N | MR sqrt(-1) mean rate/base | median log2 | SS b=-1 mean rate/base | CF PC1 (T=1024) | CF PC2 (T=1024) | false certs on eps=-1 |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for b in sorted({r["bits"] for r in R}):
        v = [r for r in R if r["bits"] == b and r.get("Nmod4") == 1]
        vp = [r for r in v if r["eps"] == 1]
        vm = [r for r in v if r["eps"] == -1]
        mr = [2 ** r["mr_cert"] if r["mr_cert"] is not None else 0 for r in vp]
        ss = [2 ** r["ss_cert"] if r["ss_cert"] is not None else 0 for r in vp]
        med = np.median([r["mr_cert"] for r in vp if r["mr_cert"] is not None])
        lines.append(f"| {b} | {len(vp)} | {np.mean(mr):.2e} | {med:.1f} | {np.mean(ss):.2e} | {sum(r['pc1'] for r in vp)} | "
                     f"{sum(r['pc2'] for r in vp)} | {sum(r['pc1'] or r['pc2'] for r in vm)} |")
    lines.append("")


def epsexact(lines):
    R = json.load(open(RAW + "primality_epsexact.json"))
    lines += ["## E4b exact transcript laws given (eps, N mod 8), all prime pairs in [2^(b-1), 2^b)", "",
              "| b | N mod 8 | pairs eps+/eps- | TV(MR pattern) | TV(SS class) | P(MR L>=2 \\| eps-) | P(MR L>=2 \\| eps+) | leak eps+ | leak eps- |",
              "|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        lp = sum(v for k, v in r["mr_eps+"].items() if k.startswith("L") and int(k[1:]) >= 2)
        lm = sum(v for k, v in r["mr_eps-"].items() if k.startswith("L") and int(k[1:]) >= 2)
        lines.append(f"| {r['b']} | {r['N_mod8']} | {r['pairs_eps+']}/{r['pairs_eps-']} | {r['tv_mr']:.2e} | {r['tv_ss']:.2e} | "
                     f"{lm:.1e} | {lp:.2e} | {r['leak_eps+']:.3f} | {r['leak_eps-']:.3f} |")
    lines.append("")


def battery(lines):
    cfg = json.load(open("experiments/configs/primality_battery.json"))
    R = jl("primality_battery.jsonl")
    feats = cfg["features"] + list(cfg["controls"])
    lines += ["## E4c pre-registered transcript battery (config frozen in commit before the run)", "",
              "| feature | " + " | ".join(f"z (n={b})" for b in cfg["bits"]) + " |", "|---|" + "---:|" * len(cfg["bits"])]
    usable = {}
    for f in feats:
        row = []
        for b in cfg["bits"]:
            v = [r for r in R if r["bits"] == b]
            usable[b] = len(v)
            num = den = 0.0
            degenerate = True
            for c in (1, 5):
                xp = np.array([r["feats"][f] for r in v if r["N_mod8"] == c and r["eps"] == 1], float)
                xm = np.array([r["feats"][f] for r in v if r["N_mod8"] == c and r["eps"] == -1], float)
                if len(xp) < 2 or len(xm) < 2:
                    continue
                var = xp.var(ddof=1) / len(xp) + xm.var(ddof=1) / len(xm)
                if var == 0:
                    continue
                degenerate = False
                diff = xp.mean() - xm.mean()
                num += diff / var
                den += 1 / var
            if degenerate:
                row.append("DEGENERATE")
            else:
                z = num / math.sqrt(den)
                flag = " **HIT**" if abs(z) >= 3.32 and not f.startswith("CONTROL") else ""
                row.append(f"{z:+.2f}{flag}")
        lines.append(f"| {f} | " + " | ".join(row) + " |")
    lines += ["", "Usable instances (N = 1 mod 4): " + ", ".join(f"n={b}: {usable[b]}" for b in cfg["bits"]) +
              ". Threshold |z| >= 3.32 (Bonferroni, 42 tests).", ""]


def menu(lines):
    R = jl("primality_menu.jsonl")
    keys = sorted(R[0]["leak"])
    lines += ["## E6 exponent x group menu: leak probability A_p + A_q - 2 A_p A_q (uniform element)", "",
              "| group | exponent | n=32 median log2 | n=48 median log2 | n=64 median log2 | n=64 P(leak >= 1/2) |", "|---|---|---:|---:|---:|---:|"]

    def ml(v):
        xs = [math.log2(x) if x > 0 else -999 for x in v]
        m = np.median(xs)
        return "-inf" if m <= -999 else f"{m:.1f}"
    order = ["Zn*", "torus_D", "unit_Z[sqrtD]", "T2_ideal", "F2_ideal", "T3_ideal", "F3_ideal", "T4_ideal", "F4_ideal", "T6_ideal", "F6_ideal"]
    enames = ["N-1", "N+1", "N", "(N-1)/2", "N^2-1", "N^3-1", "N^4-1", "N^6-1", "L100", "L1000", "L10000", "L100000",
              "(N-1)L1000", "(N^2-1)L1000", "(N^6-1)L1000", "(N^2-1)(N^6-1)L10000"]
    for gname in order:
        for e in enames:
            key = f"{gname}|{e}"
            if key not in R[0]["leak"]:
                continue
            cols = []
            for b in (32, 48, 64):
                v = [r["leak"][key] for r in R if r["bits"] == b]
                cols.append(ml(v))
            v64 = [r["leak"][key] for r in R if r["bits"] == 64]
            lines.append(f"| {gname} | {e} | " + " | ".join(cols) + f" | {np.mean([x >= 0.5 for x in v64]):.3f} |")
    ec = [row for r in R if "ec" in r for row in r["ec"]]
    lines += ["", f"Random elliptic curves, n = 32 ({len(ec)} curves; cyclic approximation gcd(E,#E)/#E):", "",
              "| exponent | median log2 leak | P(leak >= 1/2) |", "|---|---:|---:|"]
    for e in ec[0]:
        v = [row[e] for row in ec]
        lines.append(f"| {e} | {ml(v)} | {np.mean([x >= 0.5 for x in v]):.3f} |")
    lines.append("")


def main():
    lines = ["# primality: experiment tables (generated by src/primality/analyze.py; do not edit by hand)", ""]
    for fn in (leakrates, kangaroo, iterate, aks, epscert, epsexact, battery, menu):
        try:
            fn(lines)
        except FileNotFoundError as e:
            lines += [f"(missing {e.filename})", ""]
    open(OUT, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
