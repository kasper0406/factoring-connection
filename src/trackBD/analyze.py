"""Aggregate Track B/D raw results into the tables quoted in notes/trackBD_analysis.md.

Usage: .venv/bin/python -m src.trackBD.analyze > experiments/analysis/trackBD_tables.txt
"""
from __future__ import annotations

import json
import os
from collections import Counter, defaultdict

RAW = "experiments/raw"


def rd(name):
    import gzip
    path = os.path.join(RAW, name)
    if os.path.exists(path):
        return [json.loads(l) for l in open(path)]
    if os.path.exists(path + ".gz"):
        return [json.loads(l) for l in gzip.open(path + ".gz", "rt")]
    return []


def xor(a, b):
    return a + b - 2 * a * b


def tiny():
    rows = rd("trackBD_tiny.jsonl")
    if not rows:
        return
    print("## B1 tiny exhaustive (all composites 4 <= N < 5000)")
    groups = {"odd_semiprime": Counter(), "odd_3plus": Counter(), "prime_power": Counter(),
              "even": Counter(), "carmichael": Counter()}
    carm = {561, 1105, 1729, 2465, 2821}
    viol = 0
    for r in rows:
        fs = r["factors"]
        if r["N"] % 2 == 0:
            g = "even"
        elif len(set(fs)) == 1:
            g = "prime_power"
        elif len(fs) == 2:
            g = "odd_semiprime"
        else:
            g = "odd_3plus"
        groups[g].update(r["agg"])
        if r["N"] in carm:
            groups["carmichael"].update(r["agg"])
        viol += sum(v for k, v in r["agg"].items() if "VIOL" in k)
    print(f"Theorem B1 root-level check violations: {viol}")
    for g, agg in groups.items():
        print(f"\n### group {g}")
        for fam in ("binom", "cyclo", "rand"):
            for reg in ("cmp", "unc"):
                pre = f"{fam}|{reg}|"
                nsk = agg.get(pre + "sketches", 0)
                if not nsk:
                    continue
                print(f"{fam:5s} {reg} sketches={nsk}")
                for ro in ("coef", "content", "euclid"):
                    cls = Counter()
                    for k, v in agg.items():
                        if k.startswith(pre + ro + "|"):
                            lab = k[len(pre + ro) + 1:]
                            # collapse TD reasons
                            cls["TD" if lab.startswith("TD") else lab] += v
                            if lab.startswith("TD") and "res_phi3" in lab:
                                cls["(TD incl. res_phi3)"] += v
                    tot = sum(v for k, v in cls.items() if not k.startswith("("))
                    s = ", ".join(f"{k}:{v}" for k, v in sorted(cls.items()))
                    print(f"   {ro:8s} hits={tot:7d} ({tot / nsk:.4f}/sketch)  {s}")


def scalar():
    for fn in ("trackBD_scalar.jsonl", "trackBD_scalar_diag.jsonl"):
        rows = rd(fn)
        if not rows:
            continue
        print(f"\n## B2/D2 scalar sketches: {fn}")
        print("key  inst  K/inst | F: TD_0 TD_m1 TD_phi3 FROB SPOR all (pred) | G1: ORDER TD_half all (pred) | NULL | "
              "mean p | E[frob+spor on F_p] | E[g] | spor mod 6")
        by = defaultdict(list)
        for r in rows:
            by[(r["bits"], r["extra"])].append(r)
        for k in sorted(by):
            R = by[k]
            h = Counter()
            for r in R:
                h.update(r["hits"])

            def s(pre):
                return sum(v for kk, v in h.items() if kk.startswith(pre))
            pF = pG = 0.0
            fs, gs, m6 = [], [], Counter()
            for r in R:
                if "census" not in r:
                    continue
                p, q = sorted(r["factors"])
                cp, cq = r["census"]["p"], r["census"]["q"]
                pF += r["K"] * xor((cp["triv"] + cp["frob"] + cp["spor"]) / p,
                                   (cq["triv"] + cq["frob"] + cq["spor"]) / q)
                pG += r["K"] * xor(cp["G1"] / p, cq["G1"] / q)
                fs.append(cp["frob"] + cp["spor"])
                gs.append(cp["g"])
                m6[cp["spor"] % 6] += 1
                m6[cq["spor"] % 6] += 1
            mp = sum(sorted(r["factors"])[0] for r in R) / len(R)
            print(f"{k} {len(R)} {R[0]['K']} | F: {s('F|TD_0')} {s('F|TD_m1')} {s('F|TD_phi3')} {s('F|FROB')} "
                  f"{s('F|SPOR')} {h['F|any']} ({pF:.1f}) | G1: {s('G1|ORDER')} {s('G1|TD_half')} {h['G1|any']} "
                  f"({pG:.1f}) | {h['NULL|any']} | {mp:.0f} | {sum(fs) / max(1, len(fs)):.2f} | "
                  f"{sum(gs) / max(1, len(gs)):.2f} | {dict(m6)}")
    rows = rd("trackBD_scalar_carm.jsonl")
    if rows:
        tot = Counter()
        for r in rows:
            tot.update(r["hits"])
        print(f"\nCarmichael (n={len(rows)}, K=2000 each): {dict(sorted(tot.items()))}")


def rankdist():
    rows = rd("trackBD_rankdist.jsonl")
    if not rows:
        return
    print("\n## B3b exact rank distributions, X^r - c, uniform admissible c")
    for r in sorted(rows, key=lambda x: (x["N"], x["r"])):
        if r["tag"].startswith("codex"):
            print(f"{r['N']}={r['p']}*{r['q']} r={r['r']} P_gap={r['P_gap_num']}/{r['P_gap_den']} "
                  f"F_l-roots(p,q)=({r['side_p']['n_Fl_roots']},{r['side_q']['n_Fl_roots']}) "
                  f"attrib p={r['side_p']['attrib']} q={r['side_q']['attrib']}")
    agg = defaultdict(list)
    for r in rows:
        if not r["tag"].startswith("balanced"):
            continue
        sp, sq = r["side_p"], r["side_q"]
        a = sp["attrib"].get("Fl", 0) / sp["adm"]
        b = sq["attrib"].get("Fl", 0) / sq["adm"]
        agg[(r["tag"], r["r"])].append((r["P_gap"], xor(a, b),
                                        sp["attrib"].get("ext", 0) + sq["attrib"].get("ext", 0),
                                        r["P_gap"] > 0, r["p"]))
    print("tag r n mean_Pgap mean_pred_from_F_l_roots frac_inst_with_gap ext_c_total mean(Pgap*p)")
    for k in sorted(agg):
        v = agg[k]
        n = len(v)
        print(k[0], k[1], n, f"{sum(x[0] for x in v) / n:.5f}", f"{sum(x[1] for x in v) / n:.5f}",
              f"{sum(x[3] for x in v) / n:.2f}", sum(x[2] for x in v),
              f"{sum(x[0] * x[4] for x in v) / n:.2f}")


def rank():
    rows = rd("trackBD_rank.jsonl")
    if not rows:
        return
    print("\n## B3 rank-split (column-gcd elimination) on balanced semiprimes")
    n, sp, acc = Counter(), Counter(), Counter()
    bad = Counter()
    for r in rows:
        k = (r["mat"], r["bits"])
        n[k] += 1
        sp[k] += int(r["profile_split"])
        acc[k] += r["acc_nonunits"]
        if r["mat"] == "M1" and (r["rank_p"] != r["pred_rank_p"] or r["rank_q"] != r["pred_rank_q"]):
            bad["M1 rank != Theorem B1 prediction"] += 1
        if bool(r["factor"]) != r["profile_split"]:
            bad["Prop B3 (factor iff profile split)"] += 1
    print("violations:", dict(bad) or "none")
    for k in sorted(n):
        print(f"{k[0]:4s} bits={k[1]} n={n[k]} splits={sp[k]} rate={sp[k] / n[k]:.4f} "
              f"accidental_nonunit_entries/matrix={acc[k] / n[k]:.3f}")


def bm():
    rows = rd("trackBD_bm.jsonl")
    if not rows:
        return
    print("\n## D5 Berlekamp-Massey over Z/N")
    by = defaultdict(list)
    for r in rows:
        by[(r["seq"], r["bits"])].append(r)
    for k in sorted(by):
        R = by[k]
        st = [r for r in R if r.get("class") == "STRUCT"]
        ac = [r for r in R if r.get("class") == "ACC"]
        lc = sum(r["LC_p"] == min(r["e_p"], r["T"] // 2) for r in R) if k[0] == "S1" else "-"
        ratio = (sum(r["step"] / (2 * (r["LC_p"] if r["revealed"] == "p" else r["LC_q"])) for r in st)
                 / max(1, len(st)))
        print(f"{k} n={len(R)} STRUCT={len(st)} ACC={len(ac)} none={len(R) - len(st) - len(ac)} "
              f"LC_p==e_p:{lc} step/(2LC)[STRUCT]={ratio:.2f} meanstep ACC="
              f"{sum(r['step'] for r in ac) / max(1, len(ac)):.0f} mean p={sum(r['p'] for r in R) / len(R):.0f} T={R[0]['T']}")


def sketch_bal():
    rows = rd("trackBD_sketch_bal.jsonl")
    if not rows:
        return
    print("\n## B2b degree-r sketches on balanced semiprimes (per-sketch rates; hits by class)")
    by = defaultdict(list)
    for r in rows:
        by[r["bits"]].append(r)
    for b in sorted(by):
        R = by[b]
        tot = Counter()
        for r in R:
            tot.update(r["counts"])
        mp = sum(r["p"] for r in R) / len(R)
        null = sum(r["null_per_coef"] for r in R) / len(R)
        print(f"\n### bits={b} instances={len(R)} mean p={mp:.0f} null/coef={null:.2e}")
        names = sorted({k.split("|")[0] for k in tot})
        for nm in names:
            nsk = tot[nm + "|sketches"]
            line = []
            for ro in ("coef", "content", "euclid"):
                cls = Counter()
                for k, v in tot.items():
                    if k.startswith(nm + "|" + ro + "|"):
                        lab = k.split("|", 2)[2]
                        cls["TD" if lab.startswith("TD") else lab] += v
                if cls:
                    line.append(f"{ro}: " + ",".join(f"{k}={v}" for k, v in sorted(cls.items())))
            pred = ""
            if nm.startswith("binom"):
                rr = nm[5:]
                pv = [r["pred"][rr] for r in R if "pred" in r]
                if pv:
                    pred = f" pred_split/sketch={sum(pv) / len(pv):.2e}"
                    rcoef = int(rr)
                    pred += f" null_coef/sketch~{min(1.0, rcoef * null):.2e}"
            print(f"  {nm:9s} n={nsk}{pred} | " + " | ".join(line))


if __name__ == "__main__":
    tiny()
    scalar()
    rankdist()
    rank()
    bm()
    sketch_bal()
