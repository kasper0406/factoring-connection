"""Score the pre-registered predictions (notes/sicheck_prereg.md) from experiments/raw/sicheck_E*.jsonl.

HARNESS. Floating point is used only for statistics (ratios, KS distances, Dickman rho for display).
Usage: .venv/bin/python -m src.sicheck.score [E1|E3|E4]
"""
from __future__ import annotations

import json
import math
import os
import sys
from collections import defaultdict

RAW = "experiments/raw"


def load(name):
    path = os.path.join(RAW, f"sicheck_{name}.jsonl")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [json.loads(l) for l in f if l.strip()]


def dickman_rho(u: float) -> float:
    """Dickman rho by integrating u rho'(u) = -rho(u - 1) on a grid (display only)."""
    if u <= 1:
        return 1.0
    step = 1e-3
    n1 = int(round(1 / step))
    vals = [1.0] * (n1 + 1)  # rho on [0, 1]
    m = int(math.ceil(u / step))
    for i in range(n1 + 1, m + 1):
        t = i * step
        # trapezoid on rho'(t) = -rho(t-1)/t
        d0 = -vals[i - 1 - n1] / (t - step)
        d1 = -vals[i - n1] / t
        vals.append(vals[-1] + step * (d0 + d1) / 2)
    return vals[m]


def ks_uniform(xs):
    xs = sorted(xs)
    n = len(xs)
    return max(max((i + 1) / n - x, x - i / n) for i, x in enumerate(xs)) if n else float("nan")


def ks_trunc_exp(xs, rate, h):
    xs = sorted(xs)
    n = len(xs)
    Z = 1 - math.exp(-rate * h)
    D = 0.0
    for i, x in enumerate(xs):
        F = (1 - math.exp(-rate * x)) / Z
        D = max(D, (i + 1) / n - F, F - i / n)
    return D


def ratio_line(O, E, V):
    se = math.sqrt(V) if V > 0 else float("nan")
    return f"O={O} E={E:.1f} O/E={O / E:.3f} z={(O - E) / se:+.2f}" if E > 0 else f"O={O} E=0"


# ----------------------------------------------------------------------------------------------------------------------
def score_E1():
    rows = load("E1")
    if not rows:
        print("no E1 data")
        return
    print("== E1 ==")
    cells = sorted({r["cell"] for r in rows}, key=lambda c: ["A", "B", "B10", "C", "D", "E"].index(c))
    cons = ["next", "last", "rand", "cofactor"]
    # P1
    print("-- P1: O/E (size-matched null), pooled over seeds; [0.85,1.15] per cell (E>=50), [0.92,1.08] all cells")
    allpool = defaultdict(lambda: [0, 0.0, 0.0])
    p1_ok = True
    for c in cons:
        for cell in cells:
            O = E = V = 0
            n = 0
            for r in rows:
                if r["cell"] == cell and r["constructor"] == c:
                    O += sum(r["vs"])
                    E += sum(r["pn"])
                    V += sum(p * (1 - p) for p in r["pn"])
                    n += len(r["v"])
            a = allpool[c]
            a[0] += O; a[1] += E; a[2] += V
            ok = (E < 50) or (0.85 <= O / E <= 1.15)
            p1_ok &= ok
            print(f"  {c:8s} {cell:4s} outputs={n:5d} {ratio_line(O, E, V)} {'ok' if ok else 'FAIL'}")
    for c in cons:
        O, E, V = allpool[c]
        ok = 0.92 <= O / E <= 1.08
        p1_ok &= ok
        print(f"  {c:8s} ALL  {ratio_line(O, E, V)} {'ok' if ok else 'FAIL'}")
    print(f"  P1 {'HELD' if p1_ok else 'FAILED'}")
    # P2
    print("-- P2: k-parity strata (next, rand)")
    p2a = p2b = True
    for c in ("next", "rand"):
        acc = {0: [0, 0.0, 0.0, 0.0, 0.0], 1: [0, 0.0, 0.0, 0.0, 0.0]}  # O, E_naive, V_naive, E_par, V_par
        for r in rows:
            if r["constructor"] != c:
                continue
            pi0 = sum(r["ue"]) / len(r["ue"])
            for k, vs, pn, pe, po in zip(r["k"], r["vs"], r["pn"], r["pe"], r["po"]):
                par = k % 2
                pev = (1 - pi0) if par == 1 else pi0
                pp = pev * pe + (1 - pev) * po
                a = acc[par]
                a[0] += vs; a[1] += pn; a[2] += pn * (1 - pn); a[3] += pp; a[4] += pp * (1 - pp)
        ro, re_ = acc[1][0] / acc[1][1], acc[0][0] / acc[0][1]
        p2a &= ro < 1 < re_
        for par, name in ((1, "k odd "), (0, "k even")):
            O, E, V, Ep, Vp = acc[par]
            okb = 0.85 <= O / Ep <= 1.15
            p2b &= okb
            print(f"  {c:5s} {name}: naive {ratio_line(O, E, V)} | parity model {ratio_line(O, Ep, Vp)} "
                  f"{'ok' if okb else 'FAIL'}")
    print(f"  P2a {'HELD' if p2a else 'FAILED'}; P2b {'HELD' if p2b else 'FAILED'}")
    # P3
    print("-- P3: residue distributions")
    p3a = p3b = p3c = True
    for cell in cells:
        for seed in sorted({r["seed"] for r in rows if r["cell"] == cell}):
            rs = {r["constructor"]: r for r in rows if r["cell"] == cell and r["seed"] == seed}
            if "next" not in rs:
                continue
            counts = rs["next"]["counts"]
            K, h = rs["next"]["K"], rs["next"]["h"]
            dlt = sum(counts) / (K * (h + 1))
            mu = dlt * (h + 1)
            pempty = sum(1 for x in counts if x == 0) / K
            re = pempty / math.exp(-mu)
            ok_a = 1 <= re <= 2.5
            p3a &= ok_a
            Dn = ks_trunc_exp(rs["next"]["v"], dlt, h)
            Dl = ks_trunc_exp([h - v for v in rs["last"]["v"]], dlt, h)
            Dr = ks_uniform([v / h for v in rs["rand"]["v"]])
            Dc = ks_uniform([v / u0 for v, u0 in zip(rs["cofactor"]["v"], rs["cofactor"]["u0"])])
            p3b &= Dn <= 0.10 and Dl <= 0.10
            p3c &= Dr <= 0.05 and Dc <= 0.05
            # mean count drift in k (toy: density falls as kN grows)
            q = K // 4
            drift = [sum(counts[i * q:(i + 1) * q]) / q for i in range(4)]
            print(f"  {cell:4s} s{seed}: mu={mu:.2f} P(empty)={pempty:.3f} exp(-mu)={math.exp(-mu):.3f} ratio={re:.2f}"
                  f" {'ok' if ok_a else 'FAIL'} | KS next={Dn:.3f} last={Dl:.3f} rand={Dr:.3f} cofactor={Dc:.3f}"
                  f" | mean count by k-quarter {[round(x, 2) for x in drift]}")
    print(f"  P3a {'HELD' if p3a else 'FAILED'}; P3b {'HELD' if p3b else 'FAILED'}; P3c {'HELD' if p3c else 'FAILED'}")
    # P4
    print("-- P4: duplicates / primitivity among v-smooth outputs (pooled)")
    dupf = {}
    p4 = True
    for c in cons:
        nv = dup = npr = 0
        for r in rows:
            if r["constructor"] != c:
                continue
            for vs, g, d in zip(r["vs"], r["g"], r["dup"]):
                if vs:
                    nv += 1
                    dup += d
                    npr += int(g > 1)
        dupf[c] = dup / nv
        p4 &= npr / nv >= 0.2
        print(f"  {c:8s} v-smooth={nv} duplicate={dup / nv:.3f} nonprimitive={npr / nv:.3f}")
    p4 &= dupf["next"] > dupf["rand"] and 0.05 <= dupf["next"] <= 0.5
    print(f"  P4 {'HELD' if p4 else 'FAILED'}")
    # descriptive: Dickman
    print("-- descriptive: observed v-smooth rate vs Psi(h)/h-type null and Dickman rho(log h/log y)")
    for cell in cells:
        r0 = [r for r in rows if r["cell"] == cell and r["constructor"] == "rand"]
        if not r0:
            continue
        y, h = r0[0]["y"], r0[0]["h"]
        O = sum(sum(r["vs"]) for r in r0)
        n = sum(len(r["vs"]) for r in r0)
        print(f"  {cell:4s} rand: observed {O / n:.4f}  mean null {sum(sum(r['pn']) for r in r0) / n:.4f}  "
              f"Dickman rho(log h/log y)={dickman_rho(math.log(h) / math.log(y)):.4f}")


# ----------------------------------------------------------------------------------------------------------------------
def score_E3():
    rows = load("E3")
    if not rows:
        print("no E3 data")
        return
    print("== E3 ==")
    honest = ["next", "last", "rand", "cofactor"]
    H = [r for r in rows if r["constructor"] in honest]
    Hd = [r for r in H if r["done"]]
    nz = sum(r["phi_nonzero"] for r in Hd)
    rt = sum(r["random_tried"] for r in Hd)
    rf = sum(r["random_fruitful"] for r in Hd)
    r1 = nz / len(Hd) >= 0.95 and 0.42 <= rf / rt <= 0.58
    print(f"R1: honest runs done {len(Hd)}/{len(H)}; Phi!=0 in {nz}/{len(Hd)}; random deps fruitful {rf}/{rt} = "
          f"{rf / rt:.3f} -> {'HELD' if r1 else 'FAILED'}")
    for c in honest:
        rc = [r for r in Hd if r["constructor"] == c]
        print(f"   {c:8s} Phi!=0 {sum(r['phi_nonzero'] for r in rc)}/{len(rc)}  fruitful "
              f"{sum(r['random_fruitful'] for r in rc)}/{sum(r['random_tried'] for r in rc)}  "
              f"dups {sum(r['dups'] for r in rc)}  mean k_last {sum(r['k_last'] for r in rc) / max(1, len(rc)):.0f}")
    Q = [r for r in rows if r["constructor"] == "qsfb"]
    qb = [r for r in Q if r["blum_diag"] and r["dim_deps"] >= 1]
    qn = [r for r in Q if not r["blum_diag"] and r["dim_deps"] >= 5]
    r2a = all(not r["phi_nonzero"] and r["random_fruitful"] == 0 for r in qb)
    r2b = (sum(r["phi_nonzero"] for r in qn) >= len(qn) / 2) if qn else None
    print(f"R2: J+ base Blum runs with deps {len(qb)} (of {sum(r['blum_diag'] for r in Q)}), all trivial: {r2a}; "
          f"non-Blum runs with >=5 deps {len(qn)} (of {sum(not r['blum_diag'] for r in Q)}), Phi!=0 in "
          f"{sum(r['phi_nonzero'] for r in qn)} -> Blum part: "
          f"{'UNTESTABLE (no Blum J+ run produced a dependency)' if not qb else ('HELD' if r2a else 'FAILED')}; "
          f"non-Blum part: {'UNTESTABLE' if r2b is None else ('HELD' if r2b else 'FAILED')}")
    for r in Q:
        print(f"   qsfb n={r['nbits']} blum={r['blum_diag']} p,q mod4={r['p_mod4_diag']},{r['q_mod4_diag']} "
              f"fb={r['fb']} rels={r['relations']} deps={r['dim_deps']} Phi!=0={r['phi_nonzero']} "
              f"fruit={r['random_fruitful']}/{r['random_tried']} k_last={r['k_last']}")
    B0 = [r for r in rows if r["constructor"] == "beta0"]
    r3 = all(not r["phi_nonzero"] and r["random_fruitful"] == 0 for r in B0)
    print(f"R3: beta0 runs {len(B0)} (done {sum(r['done'] for r in B0)}), all trivial: {r3} -> "
          f"{'HELD' if r3 else 'FAILED'}")
    bl = [r for r in rows if r["blum_diag"]]
    mism = sum(r["beta_mismatch"] for r in bl)
    chk = sum(r["beta_checked"] for r in bl)
    nb = [r for r in rows if not r["blum_diag"]]
    nbm = sum(1 for r in nb if r["beta_mismatch"] > 0)
    r4 = mism == 0 and nbm >= 1
    print(f"R4: Blum: {mism} mismatches in {chk} checked dependencies; non-Blum runs with mismatches {nbm}/{len(nb)} "
          f"(mismatch rate {sum(r['beta_mismatch'] for r in nb)}/{sum(r['beta_checked'] for r in nb)}) -> "
          f"{'HELD' if r4 else 'FAILED'}")
    r5 = True
    for d in (1, 2, 3, 4):
        runs = [r for r in H if len([x for x in r["phi_seq"] if x is not None]) >= d]
        fl = sum(1 for r in runs if all(x == 0 for x in [x for x in r["phi_seq"] if x is not None][:d]))
        if len(runs) >= 40:
            ok = 0.5 * 2 ** -d <= fl / len(runs) <= 2 * 2 ** -d
            r5 &= ok
        else:
            ok = None
        print(f"R5: d={d}: fruitless {fl}/{len(runs)} = {fl / max(1, len(runs)):.3f} vs 2^-d = {2 ** -d:.3f} "
              f"{'' if ok is None else ('ok' if ok else 'FAIL')}")
    print(f"R5 {'HELD' if r5 else 'FAILED'}")
    rk = all(r["rank"] <= r["fb"] - 2 for r in rows)
    eq = sum(1 for r in Hd if r["rank"] == r["fb"] - 2)
    r6 = rk and eq / len(Hd) >= 0.9
    print(f"R6: rank <= fb-2 in all runs: {rk}; rank = fb-2 in {eq}/{len(Hd)} honest done runs -> "
          f"{'HELD' if r6 else 'FAILED'}")
    for c in ("qsfb", "beta0"):
        rc = [r for r in rows if r["constructor"] == c]
        print(f"   {c}: rank-fb distribution {sorted(r['rank'] - r['fb'] for r in rc)}")


# ----------------------------------------------------------------------------------------------------------------------
def score_E4():
    rows = load("E4")
    if not rows:
        print("no E4 data")
        return
    print("== E4 ==")
    ok = sum(r["factor_ok"] for r in rows)
    acc = sum(1 for r in rows if r["kind"] == "ACCIDENTAL_small_factor")
    print(f"F1: factor found {ok}/{len(rows)}; ACCIDENTAL {acc} -> "
          f"{'HELD' if ok >= 0.95 * len(rows) and acc == 0 else 'FAILED'}")
    O = sum(r["v_smooth"] for r in rows)
    E = sum(r["pred_vsmooth_null"] for r in rows)
    print(f"F2: v-smooth outputs {O} vs null {E:.1f}: ratio {O / E:.3f} -> "
          f"{'HELD' if 0.85 <= O / E <= 1.15 else 'FAILED'}")
    print("F3 (descriptive):")
    for r in rows:
        nonempty = r["intervals"] - r["empty"]
        rho = dickman_rho(math.log(r["h"]) / math.log(r["y"]))
        print(f"  n={r['nbits']} s{r['seed']} {r['constructor']:4s} fb={r['fb_size']} intervals={r['intervals']} "
              f"empty={r['empty'] / r['intervals']:.2f} v-smooth={r['v_smooth']} (null {r['pred_vsmooth_null']:.1f}; "
              f"rate {r['v_smooth'] / max(1, nonempty):.3f}, rho(log h/log y)={rho:.3f}) dups={r['duplicates']} "
              f"rels={r['relations']} deps={r['dim_deps']} fruitful {r['random_fruitful']}/{r['random_deps']} "
              f"ok={r['factor_ok']} {r['secs']}s rss={r['rss_mb']}MB")


if __name__ == "__main__":
    which = sys.argv[1:] or ["E1", "E3", "E4"]
    for w in which:
        {"E1": score_E1, "E3": score_E3, "E4": score_E4}[w]()
