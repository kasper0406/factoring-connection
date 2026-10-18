"""Markdown tables from experiments/analysis/eps_battery_summary.json (formatting only, added after the freeze)."""
from __future__ import annotations

import json
import sys

from src.eps.operational import MENU


def fmt(z):
    return "—" if z is None else f"{z:+.2f}"


def main(path="experiments/analysis/eps_battery_summary.json", curve="experiments/analysis/eps_pccurve_summary.json",
         out="experiments/analysis/eps_battery_tables.md"):
    S = json.load(open(path))
    keys = sorted(S, key=lambda k: S[k]["n"])
    L = []
    for tgt in ("eps", "qr2"):
        L.append(f"## Target {tgt}: stratified z per statistic (D = degenerate within strata)\n")
        L.append("| stat | group | role | " + " | ".join(f"n={S[k]['n']} z (raw)" for k in keys) + " |")
        L.append("|---|---|---|" + "---|" * len(keys))
        for j, (name, grp, role) in enumerate(MENU):
            cells = []
            for k in keys:
                r = S[k][tgt]["rows"][j]
                assert r["stat"] == name
                cells.append("D" if r["degenerate"] else f"{r['z']:+.2f} ({r['raw_z']:+.1f})")
            L.append(f"| {name} | {grp} | {role} | " + " | ".join(cells) + " |")
        L.append("")
        L.append("| size | instances | P(y=1) | strata | ρ_MDE | hits | null mean z² | null max|z| | shuffle mean z² | planted z (exp) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for k in keys:
            e = S[k][tgt]
            L.append(f"| {S[k]['n']} | {e['n_inst']} | {e['frac_pos']:.4f} | {e['strata']} | {e['mde_rho']:.4f} | "
                     f"{', '.join(e['hits']) or 'none'} | {e['null_summary']['mean_z2']:.2f} | "
                     f"{e['null_summary']['max_abs_z']:.2f} | {e['shuffle']['mean_z2']:.2f} | "
                     f"{e['planted']['z']:.1f} ({e['planted']['z_exp']:.1f}) |")
        L.append("")
        L.append("| size | PC | fires | freq | violations | z_obs | z_pred | powered | detected | consistent |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for k in keys:
            for pc, v in S[k][tgt]["pc"].items():
                L.append(f"| {S[k]['n']} | {pc} | {v['fires']} | {v['freq']:.5f} | {v['violations']} | {v['z_obs']:.2f} | "
                         f"{v['z_pred']:.2f} | {v['powered']} | {v['detected']} | {v['consistent']} |")
        L.append("")
        L.append("| size | floor k | classes | CV accuracy |")
        L.append("|---|---|---|---|")
        for k in keys:
            for row in S[k][tgt]["floor"]:
                if row["k"] in (3, 4, 6, 8, 10, 12):
                    L.append(f"| {S[k]['n']} | {row['k']} | {row['classes']} | {row['cv_acc']:.4f} ± {row['se']:.4f} |")
            c = S[k][tgt]["floor_chi2"]
            L.append(f"| {S[k]['n']} | χ² vs N mod 2^8 | dof {c['dof']} | χ² = {c['chi2']:.1f} (z_WH {c['z_wh']:+.2f}) |")
        L.append("")
    L.append("## Nonunits and cost\n")
    L.append("| size | gcd checks | events | proper | null expected (Σ 1/p+1/q) | cf nonunit | mean s | total CPU-h |")
    L.append("|---|---|---|---|---|---|---|---|")
    for k in keys:
        nu, c = S[k]["nonunits"], S[k]["cost"]
        L.append(f"| {S[k]['n']} | {nu['gcd_checks']} | {nu['events']} | {nu['proper']} | {nu['null_expected']:.1f} | "
                 f"{nu['cf_nonunit']} | {c['mean_s']:.4f} | {c['total_cpu_h']:.2f} |")
    L.append("")
    try:
        C = json.load(open(curve))
        L.append("## Supplementary: positive-control firing frequency vs n and T\n")
        Ts = ["64", "256", "1024", "4096", "16384"]
        L.append("| n | PC | violations | " + " | ".join(f"T={T}" for T in Ts) + " |")
        L.append("|---|---|---|" + "---|" * len(Ts))
        for c in C:
            for pc in ("pc1", "pc2", "pc3", "pcq2", "pcq3"):
                v = c[pc]
                L.append(f"| {c['n']} | {pc} | {v['violations']} | " + " | ".join(
                    f"{v['freq'][T]:.4f}" for T in Ts) + " |")
    except FileNotFoundError:
        pass
    open(out, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
