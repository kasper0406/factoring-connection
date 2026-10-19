"""E1 harness (pre-registered in notes/lv6_prereg.md): the LV6 pipeline on toy non-monic fields.

Usage:
  .venv/bin/python -m src.lv6.e1 experiments/configs/lv6_e1.json experiments/raw/lv6_e1.jsonl
  .venv/bin/python -m src.lv6.e1 --analyze experiments/raw/lv6_e1.jsonl

HARNESS: instance generation (_semiprime, seeded) knows p, q; the operational pipeline (src/lv6/pipeline.py, sqrt.py, chars.py) receives only N, config and a seed.
The factors are used only to classify gcd outcomes after the fact.
"""
from __future__ import annotations

import json
import random
import sys
import time
from math import gcd
from typing import Dict, List

from . import chars, numfield as nf
from .pipeline import atom_columns, collect_relations, left_kernel, toy_field
from .sqrt import choose_prime, hensel_sqrt


def _semiprime(bits: int, seed: int):
    import sympy
    rng = random.Random(10_000 + 97 * bits + seed)
    while True:
        p = int(sympy.randprime(2 ** (bits // 2 - 1), 2 ** (bits // 2)))
        q = int(sympy.randprime(2 ** (bits - bits // 2 - 1), 2 ** (bits - bits // 2)))
        if p != q and (p * q).bit_length() == bits:
            return p * q, p, q


def square_test(atoms_in: List[tuple], fld: dict, g: List[int], p: int, rng, stats):
    """atoms_in: list of (atom, e).  Returns (eta or None, Gamma)."""
    f = fld["f"]
    d = len(f) - 1
    fd = f[d]
    lin = [nf.linear(at["a"], at["b"], fd, d) for at, e in atoms_in]
    gp = nf.gprime_elem(g)
    Gamma = nf.product(lin + [gp, gp], g)
    A = nf.coeff_bound_sqrt(g, [nf.norm1(x) for x in lin])
    eta = hensel_sqrt(Gamma, g, A, p, rng, stats)
    if eta is not None:
        assert max(abs(c) for c in eta) <= A
    return eta, Gamma, A


def nonsquare_witness(Gamma, fld, rng, bits=40, tries=200):
    f = fld["f"]
    fd = f[-1]
    for _ in range(tries):
        r, s, _ = chars.sample_first_degree_prime(f, bits, rng)
        try:
            if chars.chi_elem_value(r, s, fd, Gamma) == 1:
                return True
        except ValueError:
            continue
    return False


def congruence(atoms_in, eta, fld, g, N):
    """X from the rational side (exponent halving), Y from eta; returns (X, Y, ok, notes)."""
    f, g0, g1 = fld["f"], fld["g0"], fld["g1"]
    d = len(f) - 1
    fd = f[d]
    E = sum(e for _, e in atoms_in)
    assert E % 2 == 0
    sign = 1
    v: Dict[int, int] = {}
    for at, e in atoms_in:
        if at["G"] < 0:
            sign = -sign
        for p_, k in at["fG"].items():
            v[p_] = v.get(p_, 0) + e * k
    assert sign == 1 and all(k % 2 == 0 for k in v.values())
    X = pow(g1, -E // 2, N) if E <= 0 else pow(pow(g1, -1, N), E // 2, N)
    for p_, k in v.items():
        X = X * pow(p_, k // 2, N) % N if k >= 0 else X * pow(pow(p_, -1, N), -k // 2, N) % N
    m = g0 * pow(g1, -1, N) % N
    bimg = fd * m % N
    psi_gp = nf.psi(nf.gprime_elem(g), bimg, N)
    gg = gcd(psi_gp, N)
    if gg != 1:
        return None, None, False, "psi_gprime_nonunit"
    psi_minus = 1
    for at, e in atoms_in:
        if e < 0:
            psi_minus = psi_minus * (fd * at["G"] * pow(g1, -1, N)) % N
    Y = nf.psi(eta, bimg, N) * pow(fd, -E // 2, N) % N if E <= 0 else \
        nf.psi(eta, bimg, N) * pow(pow(fd, -1, N), E // 2, N) % N
    Y = Y * pow(psi_minus * psi_gp % N, -1, N) % N
    return X, Y, (X * X - Y * Y) % N == 0, ""


def run_instance(cfg: dict, bits: int, d: int, seed: int) -> dict:
    t0 = time.time()
    N, p_, q_ = _semiprime(bits, seed)
    rng = random.Random(1_000_003 * seed + 7919 * bits + d)          # operational coins
    events: List[str] = []
    par = cfg["per_degree"][str(d)]
    fld = toy_field(N, d, par["hbits"], rng, events)
    f = fld["f"]
    fd = f[d]
    g = nf.monic_g(f)
    D = fld["disc_f"] * fd
    atoms = collect_relations(N, fld, cfg["B"], par["A"], cfg["max_atoms"], events)
    # rows: disjoint signed groups of 1..gmax atoms
    idx = list(range(len(atoms)))
    rng.shuffle(idx)
    rows = []
    i = 0
    while i < len(idx):
        k = rng.randint(1, par["gmax"])
        grp = [(atoms[j], rng.choice([1, -1])) for j in idx[i:i + k]]
        rows.append(grp)
        i += k
    # columns
    colmap: Dict[tuple, int] = {}
    def cidx(c):
        if c not in colmap:
            colmap[c] = len(colmap)
        return colmap[c]
    atom_vec = []
    for at in atoms:
        v = 0
        for c, e in atom_columns(at, D, fd).items():
            if e % 2:
                v ^= 1 << cidx(c)
        atom_vec.append(v)
    aid = {id(at): j for j, at in enumerate(atoms)}
    base_rows = []
    for grp in rows:
        v = 0
        for at, e in grp:
            v ^= atom_vec[aid[id(at)]]
        base_rows.append(v)
    nbase = len(colmap)
    # characters
    t = cfg["n_chars"]
    char_primes = [chars.sample_first_degree_prime(f, cfg["char_bits"], rng)[:2] for _ in range(t)]
    char_rows = []
    for grp in rows:
        v = 0
        for ci, (r, s) in enumerate(char_primes):
            bit = 0
            for at, e in grp:
                bit ^= chars.chi(r, s, at["a"], at["b"])
            v |= bit << ci
        char_rows.append(v)
    full_rows = [br | (cr << nbase) for br, cr in zip(base_rows, char_rows)]
    ker_full = left_kernel(full_rows, nbase + t)
    ker_par = left_kernel(base_rows, nbase)
    hb, hparts = nf.class_group_log2_bound(f)
    p = choose_prime(g, cfg["B"], fd * nf.disc_poly(g))
    out = {"bits": bits, "d": d, "seed": seed, "N": N, "f": f, "g0": fld["g0"], "g1": fld["g1"],
           "atoms": len(atoms), "rows": len(rows), "cols_base": nbase, "n_chars": t,
           "ker_full_dim": len(ker_full), "ker_par_dim": len(ker_par), "events": events,
           "log2H_bound": hb, "log2H_parts": hparts, "sqrt_prime": p,
           "bad_prime_incidences": sum(1 for at in atoms for r in at["fF"] if D % r == 0)}

    def atoms_of(S):
        lst = []
        for ri in range(len(rows)):
            if (S >> ri) & 1:
                lst.extend(rows[ri])
        return lst

    def decide(S, label, rec):
        st: dict = {}
        ats = atoms_of(S)
        eta, Gamma, A = square_test(ats, fld, g, p, rng, st)
        r = {"label": label, "n_atoms": len(ats), "square": eta is not None, "t": st.get("t"),
             "sign_tries": st.get("sign_tries"), "log2A": A.bit_length()}
        if eta is not None:
            r["log2_eta"] = max(abs(c) for c in eta).bit_length()
            X, Y, ok, note = congruence(ats, eta, fld, g, N)
            r["congruence_ok"] = ok
            r["note"] = note
            if ok:
                gg = gcd((X - Y) % N, N)
                r["gcd_class"] = "nontrivial" if 1 < gg < N else "trivial"
        else:
            r["cert"] = "residue_nonsquare" if st.get("residue_nonsquare") else (
                "char_witness" if nonsquare_witness(Gamma, fld, rng) else "UNCERTIFIED")
        rec.append(r)
        return eta is not None

    K = cfg["K_basis"]
    tests: List[dict] = []
    for j, S in enumerate(ker_full[:K]):
        decide(S, f"full{j}", tests)
    for j, S in enumerate(ker_par[:K]):
        decide(S, f"par{j}", tests)
    # character rank on the parity kernel (dimension of the detected part of the image in H)
    crow = []
    for S in ker_par:
        v = 0
        for ri in range(len(rows)):
            if (S >> ri) & 1:
                v ^= char_rows[ri]
        crow.append(v)
    out["char_rank_on_par_kernel"] = len(ker_par) - len(left_kernel(crow, t)) if ker_par else 0
    # Algorithm LV6 combination search over the first mprime full-kernel vectors
    mp = cfg["mprime"]
    found = None
    if len(ker_full) >= 1:
        base = ker_full[:mp]
        for mask in range(1, 2 ** len(base)):
            S = 0
            for j in range(len(base)):
                if (mask >> j) & 1:
                    S ^= base[j]
            if decide(S, f"combo{mask}", tests):
                found = mask
                break
    out["combo_found_mask"] = found
    out["tests"] = tests
    # harness-only classification of gcd outcomes
    out["factors_bits"] = [p_.bit_length(), q_.bit_length()]
    out["seconds"] = round(time.time() - t0, 2)
    return out


def main(argv):
    if argv[0] == "--analyze":
        return analyze(argv[1])
    cfg = json.load(open(argv[0]))
    with open(argv[1], "w") as fh:
        for d, bits in cfg["cells"]:
            for seed in cfg["seeds"]:
                    rec = run_instance(cfg, bits, d, seed)
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                    print(d, bits, seed, rec["atoms"], rec["rows"], rec["ker_full_dim"], rec["seconds"],
                          flush=True)


def analyze(path):
    recs = [json.loads(l) for l in open(path)]
    summ = {"instances": len(recs)}
    alltests = [t for r in recs for t in r["tests"]]
    sq = [t for t in alltests if t["square"]]
    summ["P1_congruence_violations"] = sum(1 for t in sq if not t.get("congruence_ok"))
    summ["P1_squares_checked"] = len(sq)
    summ["P1_notes"] = sorted({t.get("note", "") for t in sq if t.get("note")})
    ns = [t for t in alltests if not t["square"]]
    summ["P2_uncertified_nonsquares"] = sum(1 for t in ns if t.get("cert") == "UNCERTIFIED")
    summ["P2_nonsquares"] = len(ns)
    per = []
    for r in recs:
        full = [t for t in r["tests"] if t["label"].startswith("full")]
        par = [t for t in r["tests"] if t["label"].startswith("par")]
        per.append({"d": r["d"], "bits": r["bits"], "seed": r["seed"], "atoms": r["atoms"], "rows": r["rows"],
                    "ker_full": r["ker_full_dim"], "ker_par": r["ker_par_dim"],
                    "full_sq_frac": (sum(t["square"] for t in full) / len(full)) if full else None,
                    "par_sq_frac": (sum(t["square"] for t in par) / len(par)) if par else None,
                    "char_rank": r["char_rank_on_par_kernel"], "log2H_bound": r["log2H_bound"],
                    "combo_found": r["combo_found_mask"] is not None,
                    "eta_over_A": max([t["log2_eta"] / t["log2A"] for t in r["tests"] if t["square"]] or [None]
                                      , key=lambda z: -1 if z is None else z)})
    summ["per_instance"] = per
    summ["P3_all_full_square"] = sum(1 for x in per if x["full_sq_frac"] == 1.0)
    summ["P3_eligible"] = sum(1 for x in per if x["full_sq_frac"] is not None)
    summ["P4_combo_found"] = sum(1 for x in per if x["combo_found"])
    summ["P5_char_rank_le_bound"] = sum(1 for x in per if x["char_rank"] <= x["log2H_bound"])
    summ["P6_par_frac_le_half"] = sum(1 for x in per if x["par_sq_frac"] is not None and x["par_sq_frac"] <= 0.5)
    gc = [t.get("gcd_class") for t in sq if t.get("gcd_class")]
    summ["gcd_nontrivial"] = gc.count("nontrivial")
    summ["gcd_trivial"] = gc.count("trivial")
    print(json.dumps(summ, indent=1))
    return summ


if __name__ == "__main__":
    main(sys.argv[1:])
