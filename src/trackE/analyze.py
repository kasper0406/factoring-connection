"""HARNESS: aggregate Track E raw data into experiments/analysis/trackE_tables.md.

Run: .venv/bin/python -m src.trackE.analyze
"""
from __future__ import annotations

import json
import math
import os
import statistics as st
from collections import defaultdict

from src.common.instances import balanced_semiprime
from src.common.metrics import read_jsonl
from src.trackE import oracles
from src.trackE.coppersmith import AffineHint, balanced_promise, best_m, log2_maxY, rigorous_ok

RAW = "experiments/raw/trackE_{}.jsonl"
OUT = "experiments/analysis/trackE_tables.md"


def _load(name):
    p = RAW.format(name)
    return read_jsonl(p) if os.path.exists(p) else []


def _Y_for(inst, fmt, k):
    N = inst.N
    L, U = balanced_promise(N)
    if fmt == "low":
        hint = AffineHint(1 << k, inst.labels["p"] % (1 << k), L, U)
    elif fmt == "oddM":
        import random
        rng = random.Random(f"trackE-oddM|{N.bit_length()}|{k}")
        M = rng.getrandbits(k) | (1 << (k - 1)) | 1
        hint = AffineHint(M, inst.labels["p"] % M, L, U)
    else:
        h = oracles.oracle_high_bits(inst, k)
        hint = AffineHint(1, 0, h["A"], h["A"] + h["W"])
    lo, hi = hint.x_range(); xm = (lo + hi) // 2
    return max(1, hi - xm, xm - lo), hint.L


def k_theory(inst, fmt, d, mode):
    """Smallest k satisfying the condition for dimension d (m near the optimum).
    mode: 'rig' exact HG+LLL(2^{1/4}) test; 'heur' hermite 1.02; 'svp' hermite 1.
    The condition is monotone in Y, hence in k: binary search."""
    N = inst.N
    n = N.bit_length()

    def ok(k):
        Y, L = _Y_for(inst, fmt, k)
        her = {"rig": 2 ** 0.25, "heur": 1.02, "svp": 1.0}[mode]
        m0 = best_m(n, math.log2(L), d, hermite=her)[0]
        ms = [m for m in range(m0 - 2, m0 + 3) if 1 <= m <= d]
        if mode == "rig":
            return any(rigorous_ok(N, L, Y, d, m) for m in ms)
        ly = max(log2_maxY(n, math.log2(L), d, m, hermite=her) for m in ms)
        return math.log2(Y) < ly

    lo, hi = 1, n // 2
    if not ok(hi):
        return None
    while lo < hi:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def threshold_table(lines):
    R = _load("threshold")
    if not R:
        lines.append("_no threshold data yet_\n")
        return
    by = defaultdict(list)
    for r in R:
        by[(r["n"], r["seed"], r["fmt"], r["d"])].append(r)
    kmin = {}
    tl = {}
    for key, rs in by.items():
        rs.sort(key=lambda r: r["k"])
        ok = [r for r in rs if r["success"]]
        if not ok:
            continue
        # smallest k such that every attempted k' >= k succeeded
        km = None
        for r in reversed(rs):
            if r["success"]:
                km = r["k"]
            else:
                break
        kmin[key] = km
        rr = [r for r in rs if r["k"] == km][0]
        tl[key] = (rr.get("time_lll", 0.0), rr["maxima"].get("max_operand_bits"))
    insts = {}
    lines.append("### Threshold: minimal known bits k vs lattice dimension d\n")
    lines.append("k = known bits (low: p mod 2^k; oddM: p mod M, M a random odd k-bit integer; "
                 "high: top k bits of p). Asymptote (d→∞, Theorem 1 of [May10]): "
                 "low/oddM k* = n/4 − 2 (unknown range √N/2 / 2^k ≤ 2N^{1/4}), high h* = n/4 − 1. "
                 "Columns: empirical k_min over 5 seeds (min/median/max); theory k for the same d with best m: "
                 "`rig` = exact Howgrave-Graham + worst-case LLL 2^{(d-1)/4} test (guaranteed); "
                 "`heur` = root-Hermite 1.02; `svp` = Hermite factor 1. LLL time (s, median) at k_min, max entry bits.\n")
    for fmt in ("low", "oddM", "high"):
        lines.append(f"\n**format {fmt}**\n")
        lines.append("| n | d | k_min min/med/max | excess over asymptote (med) | k rig | k heur | k svp | LLL s | entry bits |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        ns = sorted({k[0] for k in kmin})
        for n in ns:
            ds = sorted({k[3] for k in kmin if k[0] == n and k[2] == fmt})
            for d in ds:
                ks = [kmin[(n, s, fmt, d)] for s in range(5) if (n, s, fmt, d) in kmin]
                if not ks:
                    continue
                if (n, 0) not in insts:
                    insts[(n, 0)] = balanced_semiprime(n, 0)
                inst = insts[(n, 0)]
                asym = n / 4 - 2 if fmt != "high" else n / 4 - 1
                th = {mode: k_theory(inst, fmt, d, mode) for mode in ("rig", "heur", "svp")}
                tls = [tl[(n, s, fmt, d)][0] for s in range(5) if (n, s, fmt, d) in tl]
                bits = max(tl[(n, s, fmt, d)][1] or 0 for s in range(5) if (n, s, fmt, d) in tl)
                lines.append(f"| {n} | {d} | {min(ks)}/{st.median(ks):g}/{max(ks)} | {st.median(ks) - asym:+.1f} | "
                             f"{th['rig']} | {th['heur']} | {th['svp']} | {st.median(tls):.3g} | {bits} |")
    lines.append("")


def simple_tables(lines):
    # cost
    R = _load("cost")
    if R:
        lines.append("### Decode cost (format a, low bits), price per list element\n")
        lines.append("| n | k = n/4+δ | δ | d (empirical min) | decode s (med) | wrong-candidate s (med) | wrong found | d rigorous | rigorous decode s | entry bits (rig) |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        by = defaultdict(list)
        for r in R:
            by[(r["n"], r["delta"], r["mode"])].append(r)
        for (n, dl) in sorted({(k[0], k[1]) for k in by}):
            e = by.get((n, dl, "empirical_min_d"), [])
            g = by.get((n, dl, "rigorous_d"), [])
            ed = st.median([r["d"] for r in e]) if e else None
            es = st.median([r["decode_s"] for r in e]) if e else None
            ws = st.median([r["wrong_decode_s"] for r in e]) if e else None
            wf = sum(r["wrong_found"] for r in e) if e else None
            gd = [r["d"] for r in g if r.get("d")]
            gs = [r["decode_s"] for r in g if r.get("d")]
            gb = [r.get("max_operand_bits") for r in g if r.get("d") and r.get("max_operand_bits")]
            lines.append(f"| {n} | {n // 4 + dl} | {dl} | {ed} | {es:.4g} | {ws:.4g} | {wf} | "
                         f"{st.median(gd) if gd else '>80'} | {(f'{st.median(gs):.4g}') if gs else '-'} | {max(gb) if gb else '-'} |")
        lines.append("")
    # flips
    R = _load("flips")
    if R:
        lines.append("### Bit-flip noise in a low-bits hint\n")
        rs = [r for r in R if not r.get("sweep")]
        agg = defaultdict(lambda: [0, 0, 0, 0])
        for r in rs:
            a = agg[(r["n"], r["d"], r["extra"], r["e"])]
            a[0] += 1; a[1] += r["plain_ok"]; a[2] += r["trunc_ok"]; a[3] += (r["trunc_ok"] == r["predicted_trunc_ok"])
        lines.append("| n | d | k_trunc | extra bits | flips e | trials | plain ok | truncation ok | truncation = prediction (lowest flip ≥ k_trunc) |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        ktr = {(r["n"], r["d"]): r["k_trunc"] for r in rs}
        for key in sorted(agg):
            a = agg[key]
            lines.append(f"| {key[0]} | {key[1]} | {ktr[(key[0], key[1])]} | {key[2]} | {key[3]} | {a[0]} | {a[1]} | {a[2]} | {a[3]}/{a[0]} |")
        sw = [r for r in R if r.get("sweep")]
        if sw:
            lines.append("\nSingle-flip sweep (n=256, d=15): flip position j vs decoders "
                         "(plain = all k bits; trunc = low k_trunc bits; at-j = low j bits).\n")
            lines.append("| flip j | plain ok | trunc ok | at-j ok |")
            lines.append("|---|---|---|---|")
            byj = defaultdict(lambda: [0, 0, 0, 0])
            for r in sw:
                a = byj[r["flip"]]; a[0] += 1; a[1] += r["plain_ok"]; a[2] += r["trunc_ok"]; a[3] += r["trunc_at_flip_ok"]
            for j in sorted(byj):
                a = byj[j]
                lines.append(f"| {j} | {a[1]}/{a[0]} | {a[2]}/{a[0]} | {a[3]}/{a[0]} |")
            lines.append(f"\n(k = {sw[0]['k']}, k_trunc = {sw[0]['k_trunc']})\n")
    R = _load("highnoise")
    if R:
        lines.append("### Additive noise on the interval hint (n=256, d=15)\n")
        lines.append("Oracle interval [A+noise, A+noise+W] (noise unknown to the decoder); decoder widens to "
                     "[A'-tol, A'+W+tol]. Rows grouped by the width actually decoded.\n")
        lines.append("| log2 width decoded | ok when tol ≥ |noise| | ok when tol < |noise| |")
        lines.append("|---|---|---|")
        agg = defaultdict(lambda: [0, 0, 0, 0])
        for r in R:
            nz, tl = (r["noise_over_W"] or 0), (r["tol_over_W"] or 0)
            a = agg[r["log2_width"]]
            if tl >= nz:
                a[0] += 1; a[1] += r["ok"]
            else:
                a[2] += 1; a[3] += r["ok"]
        for w in sorted(agg):
            a = agg[w]
            lines.append(f"| {w} | {a[1]}/{a[0]} | {a[3]}/{a[2]} |")
        lines.append("")
    R = _load("crtnoise")
    if R:
        lines.append("### Format (iii) with wrong residues: random-subset decoder (n=128, d=11)\n")
        lines.append("| K residues | wrong e | need bits | total bits | ok | decodes (mean) | decodes (max) |")
        lines.append("|---|---|---|---|---|---|---|")
        agg = defaultdict(list)
        for r in R:
            agg[(r["K"], r["e"])].append(r)
        for key in sorted(agg):
            rs = agg[key]
            lines.append(f"| {key[0]} | {key[1]} | {rs[0]['need_bits']} | {rs[0]['total_bits']} | "
                         f"{sum(r['ok'] for r in rs)}/{len(rs)} | {st.mean(r['decodes'] for r in rs):.1f} | {max(r['decodes'] for r in rs)} |")
        lines.append("")
    R = _load("sumM")
    if R:
        lines.append("### Format (ii): s = p+q mod M, M = product of ω primes of a given size (n=128, d=7)\n")
        lines.append("| prime bits | labels | ω (med) | ambiguous ω' (med) | list size 2^{ω'} (med) | decodes to success (mean) | ok |")
        lines.append("|---|---|---|---|---|---|---|")
        agg = defaultdict(list)
        for r in R:
            if r.get("status") == "too_few_primes":
                continue
            agg[(r["lbits"], r["labels"])].append(r)
        for key in sorted(agg, key=lambda k: (k[1], -k[0])):
            rs = agg[key]
            lines.append(f"| {key[0]} | {key[1]} | {st.median(r['omega'] for r in rs):g} | {st.median(r['omega_ambiguous'] for r in rs):g} | "
                         f"{st.median(r['list_size'] for r in rs):g} | {st.mean(r['decodes'] for r in rs):.1f} | {sum(r['ok'] for r in rs)}/{len(rs)} |")
        lines.append("")
    R = _load("sum2k")
    if R:
        lines.append("### Format (i): s = p+q mod 2^k (v = v_2(q−p), harness label)\n")
        lines.append("| n | d | k − kmin_guess | v | trials | coset modulus bits = k − v ? | ok | decodes (mean) |")
        lines.append("|---|---|---|---|---|---|---|---|")
        agg = defaultdict(list)
        for r in R:
            agg[(r["n"], r["d"], r["k"] - r["kmin_guess"], min(r["v"], 4))].append(r)
        for key in sorted(agg):
            rs = agg[key]
            cb = all(set(r["coset_bits"]) == {r["k"] - r["v"]} for r in rs if r["k"] > 2 * r["v"])
            lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {key[3] if key[3] < 4 else '≥4'} | {len(rs)} | {cb} | "
                         f"{sum(r['ok'] for r in rs)}/{len(rs)} | {st.mean(r['decodes'] for r in rs):.2f} |")
        lines.append("")
    R = _load("hsbits")
    if R:
        lines.append("### Random known bits of p and q (fraction δ of the low t bits of each) → branch-and-prune → (a)\n")
        lines.append("| n | d | depth t | δ | known bits (mean) | max width (med) | final width (med) | decodes (mean) | ok |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        agg = defaultdict(list)
        for r in R:
            agg[(r["n"], r["d"], r["delta"])].append(r)
        for key in sorted(agg):
            rs = agg[key]
            ok = sum(r["status"] == "ok" for r in rs)
            mw = [r["max_width"] or 0 for r in rs]
            fw = [r["final_width"] or 0 for r in rs if r["final_width"] is not None]
            dc = [r["decodes"] or 0 for r in rs]
            lines.append(f"| {key[0]} | {key[1]} | {rs[0]['depth']} | {key[2]} | {st.mean(r['known_bits'] for r in rs):.1f} | "
                         f"{st.median(mw):g} | {st.median(fw) if fw else '-'} | {st.mean(dc):.1f} | {ok}/{len(rs)} |")
        lines.append("")
    R = _load("pm1")
    if R:
        agg = defaultdict(lambda: [0, 0])
        for r in R:
            agg[r["n"]][0] += 1; agg[r["n"]][1] += r["ok"]
        lines.append("### Format (v): e = N mod (p−1) → Pollard exponent N − e\n")
        lines.append("| n | ok / trials |")
        lines.append("|---|---|")
        for n in sorted(agg):
            lines.append(f"| {n} | {agg[n][1]}/{agg[n][0]} |")
        lines.append("")
    R = _load("legendre")
    if R:
        lines.append("### Format (iv): Legendre bits (p/l) only\n")
        rs = [r for r in R if not r.get("coset_route")]
        lines.append("Toy n = 36: survivors among all odd integers of the promise interval after K Legendre constraints "
                     "(exhaustive scan — exponential; shows the information rate only).\n")
        lines.append("| K | log2 survivors (mean) | log2 interval | true p always survives |")
        lines.append("|---|---|---|---|")
        agg = defaultdict(list)
        for r in rs:
            agg[r["K"]].append(r)
        for K in sorted(agg):
            a = agg[K]
            lines.append(f"| {K} | {st.mean(r['log2_survivors'] for r in a):.2f} | {a[0]['log2_interval']} | {all(r['divisors_among'] >= 1 for r in a)} |")
        cr = [r for r in R if r.get("coset_route")]
        lines.append("\nCoppersmith route: smallest primes with ∏l ≥ 2^need; coset list = ∏(l−1)/2.\n")
        lines.append("| n | d | need bits | #primes (= info bits) | largest l | log2 list size |")
        lines.append("|---|---|---|---|---|---|")
        for r in cr:
            lines.append(f"| {r['n']} | {r['d']} | {r['need_bits']} | {r['n_primes']} | {r['B']} | {r['log2_list']} |")
        lines.append("")


def info_floor_table(lines):
    R = _load("info_floor")
    if not R:
        return
    rows = [r for r in R if r.get("kind") != "uniformity"]
    lines.append("### Information floor: H(quantity | N mod m) in bits (exact fiber enumeration)\n")
    lines.append("`det` = determined by N mod m for every residue (zero new information). `sym` = invariant under p↔q. "
                 "`f(s)` = a function of s = p+q mod m given N mod m. Entropy is min/mean/max over N mod m.\n")
    ms = [13, 41, 16, 64, 12, 9]
    for m in ms:
        rs = [r for r in rows if r["m"] == m]
        if not rs:
            continue
        lines.append(f"\n**m = {m}** (log2 φ(m) = {rs[0]['log2_units']})\n")
        lines.append("| quantity | H min/mean/max | det | sym | f(s) |")
        lines.append("|---|---|---|---|---|")
        for r in rs:
            lines.append(f"| {r['quantity']} | {r['H_min']}/{r['H_mean']}/{r['H_max']} | {'Y' if r['determined_by_N_mod_m'] else '-'} | "
                         f"{'Y' if r['symmetric'] else '-'} | {'Y' if r['function_of_s_given_c'] else '-'} |")
    U = [r for r in R if r.get("kind") == "uniformity"]
    if U:
        mx = max(U, key=lambda r: r["chi2"] / r["dof"])
        lines.append(f"\nUniformity of p mod l given N mod l over 4000 random balanced 64-bit semiprimes "
                     f"(l ∈ {{3,5,7,11,13}}, {len(U)} cells): largest χ²/dof = {mx['chi2'] / mx['dof']:.2f} "
                     f"(l={mx['l']}, c={mx['c']}, n={mx['n']}, dof={mx['dof']}); "
                     f"mean χ²/dof = {st.mean(r['chi2'] / r['dof'] for r in U):.2f}.\n")


def main():
    lines = ["# Track E tables (generated by src/trackE/analyze.py)\n"]
    threshold_table(lines)
    simple_tables(lines)
    info_floor_table(lines)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w").write("\n".join(lines) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
