"""Part 2 analysis: ambiguity tables and conditional mutual information with the hidden bit.  DIAGNOSTIC.

Input: JSONL from src.tower.run_tower.  Public key K0 = (N mod 128, eps, QR_2).  Hidden bit h (s mod 128 <=> h | K0).

* ambiguity(key): number of distinct key values, how many are ambiguous (two pairs, same key, different s mod 128),
  and how many pairs lie in ambiguous keys -- same convention as src/repair/verify_r5.py.
* cmi(X): plug-in I(h; X | K0) in bits, over all pairs.  Null: permute X within each K0 cell (seeded), 200 times;
  report the null mean and 99th percentile.  "excess" = observed - null mean.  Also the fraction of pairs in
  K0-ambiguous cells whose extended key (K0, X) is h-pure ("decided").

Run:  .venv/bin/python -m src.tower.analyze --inp experiments/raw/tower128_invariants.jsonl.gz \
          --out experiments/analysis/tower128_tables.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
from collections import Counter, defaultdict
from math import log2
from typing import Callable, Dict, List

DISC_NAMES = ["-4N", "-8N", "-12N", "-20N", "-24N", "-28N", "-3N", "-7N"]


def load(path: str) -> List[dict]:
    """Read the JSONL written by run_tower (plain or gzip-compressed, by extension)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [json.loads(line) for line in fh]


def K0(r: dict) -> tuple:
    return (r["N_mod_128"], r["eps"], r["QR2"])


def rk(r: dict, D: str, j: int) -> int:
    ranks = r["cg"][D]["ranks"]
    return ranks[j - 1] if j - 1 < len(ranks) else 0


def invariant_extractors() -> Dict[str, Callable[[dict], object]]:
    ex: Dict[str, Callable[[dict], object]] = {}
    for D in DISC_NAMES:
        for j, nm in ((2, "r4"), (3, "r8"), (4, "r16"), (5, "r32")):
            ex[f"{D}:{nm}"] = (lambda D, j: lambda r: rk(r, D, j))(D, j)
        ex[f"{D}:(r4,r8,r16)"] = (lambda D: lambda r: (rk(r, D, 2), rk(r, D, 3), rk(r, D, 4)))(D)
        ex[f"{D}:type"] = (lambda D: lambda r: tuple(r["cg"][D]["type"]))(D)
        ex[f"{D}:v2(h)"] = (lambda D: lambda r: len(r["cg"][D]["ranks"]) and sum(r["cg"][D]["ranks"]))(D)
        for l in (2, 3, 5, 7):
            ex[f"{D}:depth(amb{l})"] = (lambda D, l: lambda r: (r["cg"][D].get(f"amb{l}") or {}).get("depth"))(D, l)
    return ex


def ambiguity(recs: List[dict], keyf: Callable[[dict], tuple]) -> Dict[str, int]:
    groups = defaultdict(set)
    for r in recs:
        groups[keyf(r)].add(r["s_mod_128"])
    amb = {k for k, v in groups.items() if len(v) > 1}
    return {"keys": len(groups), "ambiguous_keys": len(amb),
            "pairs_in_ambiguous_keys": sum(1 for r in recs if keyf(r) in amb), "pairs": len(recs)}


def _cmi(hs: List[int], xs: list, cells: List[tuple]) -> float:
    n = len(hs)
    cnt_cxh = Counter(zip(cells, xs, hs))
    cnt_cx = Counter(zip(cells, xs))
    cnt_ch = Counter(zip(cells, hs))
    cnt_c = Counter(cells)
    tot = 0.0
    for (c, x, h), k in cnt_cxh.items():
        tot += k * log2(k * cnt_c[c] / (cnt_cx[(c, x)] * cnt_ch[(c, h)]))
    return tot / n


def cmi_with_null(recs: List[dict], f: Callable[[dict], object], nperm: int = 200, seed: str = "tower128-perm"):
    hs = [r["h_bit"] for r in recs]
    cells = [K0(r) for r in recs]
    xs = [f(r) for r in recs]
    obs = _cmi(hs, xs, cells)
    rng = random.Random(seed)
    idx_by_cell = defaultdict(list)
    for i, c in enumerate(cells):
        idx_by_cell[c].append(i)
    null = []
    for _ in range(nperm):
        xp = list(xs)
        for idx in idx_by_cell.values():
            vals = [xs[i] for i in idx]
            rng.shuffle(vals)
            for i, v in zip(idx, vals):
                xp[i] = v
        null.append(_cmi(hs, xp, cells))
    null.sort()
    # fraction of K0-ambiguous pairs decided by (K0, X)
    hk = defaultdict(set)
    for c, x, h in zip(cells, xs, hs):
        hk[(c, x)].add(h)
    amb_cells = {c for c in set(cells) if len({h for cc, h in zip(cells, hs) if cc == c}) > 1}
    in_amb = [i for i, c in enumerate(cells) if c in amb_cells]
    decided = sum(1 for i in in_amb if len(hk[(cells[i], xs[i])]) == 1)
    return {"cmi_bits": round(obs, 6), "null_mean": round(sum(null) / len(null), 6),
            "null_p99": round(null[int(0.99 * (len(null) - 1))], 6), "excess": round(obs - sum(null) / len(null), 6),
            "n_values": len(set(xs)), "decided_frac_of_ambiguous": round(decided / max(1, len(in_amb)), 4)}


def cond_entropy_h(recs: List[dict]) -> float:
    cells = Counter(K0(r) for r in recs)
    ch = Counter((K0(r), r["h_bit"]) for r in recs)
    n = len(recs)
    return -sum(k / n * log2(k / cells[c]) for (c, h), k in ch.items())


def depth2(r: dict, D: str) -> int:
    """2-depth of the ramified prime above 2 (T_4 = (2,2,(N+1)/2) for -4N, T_8 = (2,0,N) for -8N), capped at 2."""
    return min(r["cg"][D]["amb2"]["depth"], 2)


def r_prime(r: dict) -> int:
    """For p = q = 3 (mod 4): the prime r in {p, q} with (other / r) = +1 (exactly one exists)."""
    p, q = r["p"], r["q"]
    return p if pow(q, (p - 1) // 2, p) == 1 else q


def depth_rule(recs: List[dict]) -> Dict[str, object]:
    """Theorems T3/T4 of notes/tower_results.md, checked on every pair with QR_2 true (p = q = +-1 mod 8).

    T3 (eps = -1):  [T_4 in C^4] XOR [T_8 in C^4] = [r = 7 (mod 16)]  (r as in r_prime), hence = h when p = q (16).
    T4 (eps = +1, (p/q) = 1):  with z4 = [T_4 in C^4], z8 = [T_8 in C^4] and p = q (mod 16):
        (z4, z8) = (1, 1) -> h = 0;  (1, 0) or (0, 1) -> h = 1;  (0, 0) -> undetermined.
    T4' (eps = +1, (p/q) = -1):  z4 = z8 XOR [psi(N) = -1]  (public), so the depths carry nothing about h."""
    out = {"T3_checked": 0, "T3_fail": 0, "T3_ambiguous_pairs": 0, "T3_ambiguous_fail": 0,
           "T4_checked": 0, "T4_fail": 0, "T4_decided": 0, "T4_undecided_h": Counter(),
           "T4p_checked": 0, "T4p_fail": 0, "T3_rank_form_fail": 0, "C3_s128_fail": 0,
           "T4_r4_is_2_fail": 0, "C4_s128_fail": 0}
    for r in recs:
        if r["QR2"] != 1:
            continue
        z4, z8 = int(depth2(r, "-4N") >= 2), int(depth2(r, "-8N") >= 2)
        if r["eps"] == -1:
            out["T3_checked"] += 1
            # rank form: r4 = 1 for both, r8 = [T in C^4], and Corollary C3 for s mod 128
            r8a, r8b = rk(r, "-4N", 3), rk(r, "-8N", 3)
            if rk(r, "-4N", 2) != 1 or rk(r, "-8N", 2) != 1 or r8a != z4 or r8b != z8:
                out["T3_rank_form_fail"] += 1
            pred = (-(r["N"] + 1) + 64 * int(r["N"] % 16 == 1) * (r8a ^ r8b)) % 128
            out["C3_s128_fail"] += int(pred != r["s_mod_128"])
            if z4 ^ z8 != int(r_prime(r) % 16 == 7):
                out["T3_fail"] += 1
            if r["psi_p"] == r["psi_q"]:
                out["T3_ambiguous_pairs"] += 1
                out["T3_ambiguous_fail"] += int(z4 ^ z8 != r["h_bit"])
            continue
        p, q = r["p"], r["q"]
        lam = pow(q, (p - 1) // 2, p) == 1
        if not lam:
            out["T4p_checked"] += 1
            out["T4p_fail"] += int(z4 != z8 ^ int(r["N"] % 16 == 9))
            continue
        if r["psi_p"] != r["psi_q"]:
            continue
        out["T4_checked"] += 1
        out["T4_r4_is_2_fail"] += int(rk(r, "-4N", 2) != 2 or rk(r, "-8N", 2) != 2)
        if (z4, z8) == (0, 0):
            out["T4_undecided_h"][r["h_bit"]] += 1
            continue
        out["T4_decided"] += 1
        pred = 0 if (z4, z8) == (1, 1) else 1
        out["T4_fail"] += int(pred != r["h_bit"])
        # Corollary C4: s = (N + 1) + 64 h (mod 128) in this cell (N = 1 mod 16, eps = +1, 2 a square)
        out["C4_s128_fail"] += int((r["N"] + 1 + 64 * pred) % 128 != r["s_mod_128"])
    out["T4_undecided_h"] = dict(out["T4_undecided_h"])
    return out


def ambiguity_by_cell(recs: List[dict], keyf: Callable[[dict], tuple], nperm: int = 20,
                      seed: str = "tower128-ambnull") -> Dict[str, object]:
    """Pairs in K0-ambiguous cells that the extended key leaves ambiguous, per (eps, QR2) cell, with a null:
    the extension part of the key is permuted within K0 cells (destroys any link to h, keeps cardinalities)."""
    ext = [keyf(r) for r in recs]
    cells = [K0(r) for r in recs]
    hs = [r["h_bit"] for r in recs]

    def count(ext_list):
        g = defaultdict(set)
        for c, x, h in zip(cells, ext_list, hs):
            g[(c, x)].add(h)
        res = Counter()
        for r, c, x in zip(recs, cells, ext_list):
            if len(g[(c, x)]) > 1:
                res[(r["eps"], r["QR2"])] += 1
        return res

    obs = count(ext)
    rng = random.Random(seed)
    idx_by_cell = defaultdict(list)
    for i, c in enumerate(cells):
        idx_by_cell[c].append(i)
    nulls = []
    for _ in range(nperm):
        xp = list(ext)
        for idx in idx_by_cell.values():
            vals = [ext[i] for i in idx]
            rng.shuffle(vals)
            for i, v in zip(idx, vals):
                xp[i] = v
        nulls.append(count(xp))
    amb0 = count([()] * len(recs))
    out = {}
    for cell in sorted(set(amb0) | set(obs)):
        nv = sorted(n[cell] for n in nulls)
        out[f"eps={cell[0]},QR2={cell[1]}"] = {
            "K0_ambiguous_pairs": amb0[cell], "still_ambiguous": obs[cell],
            "null_still_ambiguous_mean": round(sum(nv) / len(nv), 1), "null_min": nv[0]}
    return out


def analyze(recs: List[dict], nperm: int = 200) -> Dict[str, object]:
    ex = invariant_extractors()
    out: Dict[str, object] = {"pairs": len(recs), "H(h|K0)_bits": round(cond_entropy_h(recs), 4)}
    amb = {}
    amb["K0"] = ambiguity(recs, K0)
    amb["K0+r4r8(-4N,-8N)"] = ambiguity(recs, lambda r: K0(r) + tuple(rk(r, D, j) for D in ("-4N", "-8N") for j in (2, 3)))
    for D in DISC_NAMES:
        amb[f"K0+{D}:(r4,r8,r16)"] = ambiguity(recs, lambda r, D=D: K0(r) + (rk(r, D, 2), rk(r, D, 3), rk(r, D, 4)))
    amb["K0+all:(r4,r8,r16)"] = ambiguity(recs, lambda r: K0(r) + tuple(rk(r, D, j) for D in DISC_NAMES for j in (2, 3, 4)))
    amb["K0+all:type"] = ambiguity(recs, lambda r: K0(r) + tuple(tuple(r["cg"][D]["type"]) for D in DISC_NAMES))
    amb["K0+all:type+depths"] = ambiguity(recs, lambda r: K0(r) + tuple(
        (tuple(r["cg"][D]["type"]),) + tuple((r["cg"][D].get(f"amb{l}") or {}).get("depth", -1) for l in (2, 3, 5, 7))
        for D in DISC_NAMES))
    out["ambiguity"] = amb
    out["cmi"] = {name: cmi_with_null(recs, f, nperm) for name, f in ex.items()
                  if len({f(r) for r in recs}) > 1}
    # breakdown of the ambiguous population by (eps, QR2)
    by = {}
    for e in (1, -1):
        for b in (0, 1):
            sub = [r for r in recs if r["eps"] == e and r["QR2"] == b]
            by[f"eps={e},QR2={b}"] = {"pairs": len(sub), "H(h|K0)": round(cond_entropy_h(sub), 4) if sub else None}
    out["cells"] = by
    out["depth_rule"] = depth_rule(recs)
    keysets = {
        "-4N:type": lambda r: (tuple(r["cg"]["-4N"]["type"]),),
        "-8N:type": lambda r: (tuple(r["cg"]["-8N"]["type"]),),
        "(-4N,-8N):(r4,r8)": lambda r: tuple(rk(r, D, j) for D in ("-4N", "-8N") for j in (2, 3)),
        "(-4N,-8N):(r4,r8,r16)": lambda r: tuple(rk(r, D, j) for D in ("-4N", "-8N") for j in (2, 3, 4)),
        "(-4N,-8N):depth(T)": lambda r: (depth2(r, "-4N"), depth2(r, "-8N")),
        "(-4N,-8N):(r4,depth(T))": lambda r: (rk(r, "-4N", 2), rk(r, "-8N", 2), depth2(r, "-4N"), depth2(r, "-8N")),
        "(-4N,-8N):type+depth(T)": lambda r: tuple((tuple(r["cg"][D]["type"]), r["cg"][D]["amb2"]["depth"])
                                                   for D in ("-4N", "-8N")),
        "all:(r4,r8,r16)": lambda r: tuple(rk(r, D, j) for D in DISC_NAMES for j in (2, 3, 4)),
        "all:type+depths": lambda r: tuple((tuple(r["cg"][D]["type"]),) + tuple(
            (r["cg"][D].get(f"amb{l}") or {}).get("depth", -1) for l in (2, 3, 5, 7)) for D in DISC_NAMES),
        "all-but-(-4N,-8N):type+depths": lambda r: tuple((tuple(r["cg"][D]["type"]),) + tuple(
            (r["cg"][D].get(f"amb{l}") or {}).get("depth", -1) for l in (2, 3, 5, 7)) for D in DISC_NAMES[2:]),
    }
    out["ambiguity_by_cell"] = {nm: ambiguity_by_cell(recs, f) for nm, f in keysets.items()}
    out["same_2sylow"] = {
        "-12N_vs_-3N": sum(r["cg"]["-12N"]["type"] != r["cg"]["-3N"]["type"] for r in recs),
        "-28N_vs_-7N": sum(r["cg"]["-28N"]["type"] != r["cg"]["-7N"]["type"] for r in recs)}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inp", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--nperm", type=int, default=200)
    ap.add_argument("--bound", type=int, default=None, help="restrict to q < bound")
    args = ap.parse_args()
    recs = load(args.inp)
    if args.bound:
        recs = [r for r in recs if r["q"] < args.bound]
    recs.sort(key=lambda r: (r["p"], r["q"]))
    res = analyze(recs, args.nperm)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, sort_keys=True, indent=1)
            fh.write("\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("cmi", "ambiguity")}, indent=1, sort_keys=True))
    top = sorted(res["cmi"].items(), key=lambda kv: -kv[1]["excess"])[:25]
    for k, v in top:
        print(f"{k:28s} {v}")


if __name__ == "__main__":
    main()
