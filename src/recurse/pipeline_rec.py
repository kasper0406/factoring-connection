"""Toy RECURSIVE Theorem-A pipeline (agent recurse; notes/recurse_results.md §2, prereg notes/recurse_prereg.md).

OPERATIONAL: factor(N, cfg, seed, oracle) sees only N, the public config and the seed.  It never calls an
integer factorisation routine on N and never uses factors of N.

It is src/thmA/pipeline.py (imported, not edited) with ONE change, the shape test of DESCEND_ONCE:
  * the complete factorisation of the non-B-smooth cofactor of a split half is delegated to an ORACLE
    with its own coins, drawn from a private stream that is independent of every hidden root r;
  * oracle "pb"  : Pollard-Brent (the thmA toy stand-in for Lenstra-Pomerance), private coins;
  * oracle "rec" : this very pipeline, factor(cof, cfg, subseed, "rec")  -- the recursive shape test.
Both oracles return the unique prime factorisation, so the outer transcript (fields, roots, splits,
atoms, dependencies, gcd outcomes) is a function of (N, seed) only and is IDENTICAL for the two oracles
(oracle-substitution Lemma OS, notes/recurse_results.md §2.1).  Only the running time differs.
The private shape stream is consumed exactly once per shape test, whatever the oracle.
"""
from __future__ import annotations

import random
import time
from math import isqrt
from typing import Dict, List, Optional

from src.thmA import arith
from src.thmA import numfield as nf
from src.thmA import pipeline as tp
from src.thmA.smooth import factor_smooth, primorial, smooth_parts

FactorFound = tp.FactorFound
_gcd_check, _inc = tp._gcd_check, tp._inc


def _oracle_pb(n: int, subseed: int, cfg: Dict, ctx: Dict, depth: int, stats: Dict) -> Dict[int, int]:
    return arith.factor_small(n, random.Random(subseed), ctx["primes_B"])


def _oracle_rec(n: int, subseed: int, cfg: Dict, ctx: Dict, depth: int, stats: Dict) -> Dict[int, int]:
    res = _factor(n, cfg, subseed, "rec", depth + 1, stats)
    out: Dict[int, int] = {}
    for p in res:
        out[p] = out.get(p, 0) + 1
    return out


ORACLES = {"pb": _oracle_pb, "rec": _oracle_rec}


def descend_once(N, A, P, ctx, rng, ctr, shape_rng, oracle, cfg, depth, stats):
    """thmA DESCEND_ONCE with the shape test delegated to `oracle` (private coins from shape_rng)."""
    _inc(ctr, "r_trials")
    r = rng.randrange(1, N)
    if not _gcd_check(r, N, "r", ctr):
        return None
    R = r * r % N
    (u0, v0), _ = arith.gauss_reduce((N, 0), (R, 1))
    if v0 == 0 or u0 == 0:
        _inc(ctr, "split_degenerate"); return None
    assert (u0 - R * v0) % N == 0
    sides = []
    for side, x in ((+1, u0), (-1, v0)):
        if not _gcd_check(x, N, "split_half", ctr):
            return None
        x = abs(x)
        s = smooth_parts([x], ctx["primorial"])[0]
        cof = x // s
        large: List[int] = []
        if cof > 1:
            if cof > P["B1"] ** P["k_large"]:
                _inc(ctr, "shape_reject_size"); return None
            subseed = shape_rng.getrandbits(64)                      # private coins, independent of r
            fac = ORACLES[oracle](cof, subseed, cfg, ctx, depth, stats)   # the shape test
            _inc(ctr, "shape_tests")
            if len(fac) > P["k_large"] or any(e > 1 or q > P["B1"] or q <= P["B"] for q, e in fac.items()):
                _inc(ctr, "shape_reject"); return None
            large = sorted(fac)
        sides.append((side, factor_smooth(s, ctx["primes_B"]) if s > 1 else {}, large))
    if sum(len(l) for _, _, l in sides) < P["min_large"]:
        _inc(ctr, "shape_reject_min_large"); return None
    _inc(ctr, "shape_ok")
    atoms = []
    for side, _, large in sides:
        for ell in large:
            at = tp.sq_search(N, A, ell, P, ctx, ctr)
            if at is None:
                _inc(ctr, "descent_fail"); return None
            atoms.append((ell, side, at))
    t, bimg = A["t"], A["beta_img"]
    q = (1 if u0 * v0 > 0 else N - 1)
    for p, e in sides[0][1].items():
        q = q * pow(p, e, N) % N
    for p, e in sides[1][1].items():
        q = q * pow(pow(p, -1, N), e, N) % N
    rhs = q
    for ell, e, at in atoms:
        n = abs(at["G"]) // ell
        fac = at["sigma"] * A["g1"] * pow(n, -1, N) * pow(t, -1, N) * nf.psi(nf.z_atom(at["a"], at["b"], t, P["d"]), bimg, N)
        rhs = rhs * (fac % N if e > 0 else pow(fac % N, -1, N)) % N
    if rhs != R:
        _inc(ctr, "ROW_IDENTITY_VIOLATION")
    _inc(ctr, "rows")
    return {"r": r, "R": R, "u0": u0, "v0": v0, "sign": 1 if u0 * v0 > 0 else -1,
            "su": sides[0][1], "sv": sides[1][1], "atoms": atoms}


def split(N, cfg, rng, ctr, dep_log, shape_rng, oracle, depth, stats) -> int:
    """thmA SPLIT (verbatim control flow) calling the oracle-parametrised descend_once."""
    P = tp.params_for(N, cfg)
    primes_B = arith.primes_upto(P["B"])
    small_set = set(arith.primes_upto(max(64, isqrt(2 * P["B1"]) + 2 * P["S"] * isqrt(P["B1"]) + 2)))
    ctx = {"primes_B": primes_B, "primorial": primorial(primes_B),
           "isprime_small": lambda b: (b in small_set) if b < max(small_set) else arith.is_prime(b)}
    hp_cache: Dict[int, int] = {}

    def hensel_p(A):
        key = id(A)
        if key not in hp_cache:
            hp_cache[key] = nf.hensel_prime(A["g"], A["D"], P["B"])
        return hp_cache[key]
    ctx["hensel_p"] = hensel_p
    ctr["params"] = {k: P[k] for k in ("d", "M", "B", "B1", "S", "H", "W")}

    def D1(A):
        return descend_once(N, A, P, ctx, rng, ctr, shape_rng, oracle, cfg, depth, stats)
    A, first = None, None
    for _ in range(P["max_fields"]):
        _inc(ctr, "fields_sampled")
        Af = tp.sample_field(N, P, rng, ctr)
        if Af is None:
            continue
        _inc(ctr, "fields_accepted")
        row = D1(Af)
        if row is not None:
            A, first = Af, row
            break
    if A is None:
        raise RuntimeError("field selection cap reached")
    rows = [first]
    for rnd in range(P["max_rounds"]):
        _inc(ctr, "lv6_rounds")
        tr = 0
        colset = set()
        for rv in rows:
            if "_vec" not in rv:
                rv["_vec"] = tp.row_vectors(rv, A, ctx)
            colset.update(k for k, v in rv["_vec"][0].items() if v % 2)
        while True:
            need = len(colset) + P["n_chars"] + P["extra_rows"]
            if len(rows) >= need:
                break
            row = D1(A)
            tr += 1
            if tr > P["max_rtrials"]:
                raise RuntimeError("row cap reached")
            if row is not None:
                row["_vec"] = tp.row_vectors(row, A, ctx)
                colset.update(k for k, v in row["_vec"][0].items() if v % 2)
                rows.append(row)
        ctr["rows_last"] = len(rows)
        g = tp.lv6_round(N, A, rows, P, ctx, rng, ctr, dep_log)
        if g is not None:
            return g
        rows = []
        row = None
        while row is None:
            row = D1(A)
        rows = [row]
    raise RuntimeError("lv6 rounds cap reached")


def _row_digest(row: Dict) -> tuple:
    return (row["r"], row["R"], row["u0"], row["v0"],
            tuple((ell, side, at["a"], at["b"]) for ell, side, at in row["atoms"]))


def _factor(N: int, cfg: Dict, seed: int, oracle: str, depth: int, stats: Dict) -> List[int]:
    """Complete factorisation (thmA FACTOR control flow).  stats collects per-call logs (all depths)."""
    rng = random.Random(f"thmA|{N}|{seed}")
    shape_rng = random.Random(f"recurse-shape|{N}|{seed}")       # private shape-test coin stream
    stats["max_depth"] = max(stats.get("max_depth", 0), depth)
    out: List[int] = []
    stack = [N]
    while stack:
        n = stack.pop()
        if n == 1:
            continue
        if arith.is_prime(n):
            out.append(n); continue
        pp = arith.perfect_power(n)
        if pp:
            stack.extend([pp[0]] * pp[1]); continue
        P = tp.params_for(n, cfg)
        found = False
        for p in arith.primes_upto(P["B"]):
            if n % p == 0:
                while n % p == 0:
                    n //= p
                    out.append(p)
                found = True
        if found:
            stack.append(n); continue
        if n < cfg["N0"]:
            q = 2
            while n % q:
                q += 1
            stack.extend([q, n // q]); continue
        ctr: Dict = {"n_bits": n.bit_length(), "_continue": bool(cfg.get("continue_on_accidental"))}
        deps: List[Dict] = []
        t0 = time.time()
        rows_seen: List[tuple] = []
        try:
            g = split(n, cfg, rng, ctr, deps, shape_rng, oracle, depth, stats)
            kind = "STRUCTURAL_DD"
        except FactorFound as ff:
            g, kind = ff.g, "ACCIDENTAL_" + ff.kind
        ctr.pop("_continue", None)
        stats.setdefault("splits", []).append({
            "depth": depth, "n": n, "factor_kind": kind, "g": g, "deps": deps,
            "seconds": round(time.time() - t0, 3),
            "ctr": {k: v for k, v in ctr.items() if k != "params"}})
        stack.extend([g, n // g])
    return sorted(out)


def factor(N: int, cfg: Dict, seed: int, oracle: str = "rec") -> Dict:
    stats: Dict = {}
    t0 = time.time()
    facs = _factor(N, cfg, seed, oracle, 0, stats)
    return {"factors": facs, "stats": stats, "seconds": round(time.time() - t0, 3)}


def outer_transcript(res: Dict) -> list:
    """Depth-0 transcript used for the oracle-substitution check (P2): per outer split, the input, the
    split outcome, the counters that record every accept/reject decision, and the dependency outcomes."""
    keep = ("r_trials", "shape_tests", "shape_reject", "shape_reject_size", "shape_ok", "rows",
            "descent_fail", "fields_sampled", "fields_accepted", "lv6_rounds", "rows_last", "cols_last",
            "kernel_dim_last", "ROW_IDENTITY_VIOLATION")
    tr = []
    for s in res["stats"].get("splits", []):
        if s["depth"] == 0:
            tr.append((s["n"], s["factor_kind"], s["g"], tuple((k, s["ctr"].get(k)) for k in keep),
                       tuple((d["size"], d["fruitful"], d["order"]) for d in s["deps"])))
    return tr
