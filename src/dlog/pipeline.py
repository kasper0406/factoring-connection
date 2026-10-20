"""Toy Theorem-DL pipeline: discrete logarithms in F_p^* (notes/dlog_results.md §3; prereg notes/dlog_prereg.md).

OPERATIONAL: dlog(p, g, h, cfg, seed) sees only p, g, h, the public config and the seed.  It never uses the
secret exponent, never factors p - 1 and never uses the order of any element.  (Only gcd(B, p - 1) is computed.)

Algorithm (dlog_results §3):
  field selection and special-q descent are those of the Theorem-A toy (src/thmA/pipeline.py, imported):
     F(g0, g1) = p w, f = F_hat + (g1 x - g0) R_c, SQ_SEARCH with batch smoothness;
  rows:   T_i = g^{a_i} h^{b_i} (a_i, b_i uniform mod p - 1), Gauss split T_i = u0/v0 (mod p), shape test,
          descent of every large prime ell = sigma g1 psi(a - b alpha)/n_ell;
  columns (EXACT integer exponents): rational primes <= B, one column for g1 (= the ideal (1, alpha)^{-1},
          entry E_i), first-degree prime ideals (r, s) of the integral ideals (a - b alpha)(1, alpha)^{-1};
  linear algebra: exact integer left kernel Lambda_1 (flint nullspace);
  unit elimination: archimedean logs of gamma_x (arb balls) + LLL integer relations; a relation whose
          certified log vector is below the Blanksby-Montgomery bound is a root of unity;
  answer: x* primitive in Lambda_2, W = prod T_i^{x*_i} = +-psi(zeta) in mu_{2w}, A + X B = log W (mod p - 1);
          Bezout over the relations found; verify g^X = h (Las Vegas).
Toy deviations (prereg): relaxed shape and toy sizes as in the thmA toy; Pollard-Brent shape test;
atoms whose norm has a prime factor dividing D = f_d disc(f) are discarded (the theorem uses local maximal
orders there); Miller-Rabin instead of AKS; dense flint nullspace instead of Dixon lifting.
"""
from __future__ import annotations

import math
import random
import time
from math import gcd
from typing import Dict, List, Optional, Tuple

import flint

from src.thmA import arith
from src.thmA import numfield as nf
from src.thmA import pipeline as tp
from src.thmA.smooth import factor_smooth, primorial, smooth_parts


def _inc(ctr, k, v=1):
    ctr[k] = ctr.get(k, 0) + v


# ----------------------------------------------------------------------------- context
def make_ctx(p: int, cfg: Dict):
    P = tp.params_for(p, cfg)
    primes_B = arith.primes_upto(P["B"])
    small_set = set(arith.primes_upto(max(64, math.isqrt(2 * P["B1"]) + 2 * P["S"] * math.isqrt(P["B1"]) + 2)))
    ctx = {"primes_B": primes_B, "primorial": primorial(primes_B),
           "isprime_small": lambda b: (b in small_set) if b < max(small_set) else arith.is_prime(b)}
    return P, ctx


# ----------------------------------------------------------------------------- atoms and rows
def good_atom(at: Dict, A: Dict) -> bool:
    """Toy restriction: every prime of |F_A(P)| is coprime to D = f_d disc f (first-degree, exact valuations)."""
    return all(A["D"] % r for r in at["fF"])


def sq_atom(p, A, ell, P, ctx, ctr):
    at = tp.sq_search(p, A, ell, P, ctx, ctr)
    if at is None:
        return None
    key = ("_good", ell)
    if key not in A["sq_cache"]:
        A["sq_cache"][key] = good_atom(at, A)
        if not A["sq_cache"][key]:
            _inc(ctr, "atom_bad_prime_discarded")
    return at if A["sq_cache"][key] else None


def descend_target(p: int, A: Dict, P: Dict, ctx: Dict, T: int, rng: random.Random, ctr: Dict) -> Optional[Dict]:
    """Split T = u0/v0 (mod p), shape test, descend every large prime.  Public inputs: T, A, coins."""
    (u0, v0), _ = arith.gauss_reduce((p, 0), (T, 1))
    if u0 == 0 or v0 == 0:
        _inc(ctr, "split_degenerate"); return None
    sides = []
    for side, x in ((+1, u0), (-1, v0)):
        x = abs(x)
        s = smooth_parts([x], ctx["primorial"])[0]
        cof = x // s
        large: List[int] = []
        if cof > 1:
            if cof > P["B1"] ** P["k_large"]:
                _inc(ctr, "shape_reject_size"); return None
            fac = arith.factor_small(cof, rng, ctx["primes_B"])
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
            at = sq_atom(p, A, ell, P, ctx, ctr)
            if at is None:
                _inc(ctr, "descent_fail"); return None
            atoms.append((ell, side, at))
    return {"T": T, "u0": u0, "v0": v0, "su": sides[0][1], "sv": sides[1][1], "atoms": atoms}


def row_columns(row: Dict) -> Tuple[Dict[tuple, int], int]:
    """Exact integer exponents.  T = sgn * prod_p p^{rat} * g1^E * prod psi(a - b alpha)^{e}  (mod p)."""
    cols: Dict[tuple, int] = {}
    for q, e in row["su"].items():
        cols[("rat", q)] = cols.get(("rat", q), 0) + e
    for q, e in row["sv"].items():
        cols[("rat", q)] = cols.get(("rat", q), 0) - e
    E = 0
    for ell, e, at in row["atoms"]:
        E += e
        for q, v in at["fn"].items():
            cols[("rat", q)] = cols.get(("rat", q), 0) - e * v
        a, b = at["a"], at["b"]
        for r, v in at["fF"].items():
            key = ("alg", r, a * pow(b, -1, r) % r)
            cols[key] = cols.get(key, 0) + e * v
    if E:
        cols[("g1",)] = E
    return {k: v for k, v in cols.items() if v}, E


def row_identity_ok(p: int, A: Dict, row: Dict) -> bool:
    """Diagnostic: T == sgn * prod p^rat * g1^E * prod psi(a - b alpha)^e (mod p), up to sign."""
    g1, m = A["g1"], A["g0"] * A["g1inv"] % p
    cols, E = row_columns(row)
    rhs = pow(g1, E, p) if E >= 0 else pow(pow(g1, -1, p), -E, p)
    for k, v in cols.items():
        if k[0] == "rat":
            rhs = rhs * (pow(k[1], v, p) if v >= 0 else pow(pow(k[1], -1, p), -v, p)) % p
    for ell, e, at in row["atoms"]:
        ps = (at["a"] - at["b"] * m) % p
        rhs = rhs * (ps if e > 0 else pow(ps, -1, p)) % p
    return rhs in (row["T"] % p, (-row["T"]) % p)


# ----------------------------------------------------------------------------- archimedean logs
def bm_bound(d: int) -> flint.arb:
    """Blanksby-Montgomery (1971): a non-root-of-unity algebraic integer of degree <= d has Mahler measure
    > 1 + 1/(52 d log(6d)).  Hence ||Log u||_inf >= log(1 + 1/(52 d log 6d))/d for every unit u not in mu_K."""
    x = flint.arb(1) / (52 * d * flint.arb(6 * d).log())
    return (1 + x).log() / d


def embeddings(f: List[int], prec: int):
    flint.ctx.prec = prec
    roots = flint.fmpz_poly(f).complex_roots()
    real = [r for r, m in roots if r.imag.is_zero()]
    cplx = [r for r, m in roots if not r.imag.is_zero() and r.imag > 0]
    assert len(real) + 2 * len(cplx) == len(f) - 1
    return real, cplx


def atom_log(a: int, b: int, real, cplx) -> List[flint.arb]:
    out = [abs(flint.acb(a) - b * r).log() for r in real]
    out += [2 * abs(flint.acb(a) - b * r).log() for r in cplx]
    return out


# ----------------------------------------------------------------------------- post-processing (Prop XK)
def exact_kernel(rows: List[Dict]):
    keys = sorted({k for rw in rows for k in rw["_cols"]})
    idx = {k: j for j, k in enumerate(keys)}
    Et = [[0] * len(rows) for _ in keys]
    for i, rw in enumerate(rows):
        for k, v in rw["_cols"].items():
            Et[idx[k]][i] = v
    K, nullity = flint.fmpz_mat(Et).nullspace()
    basis = [[int(K[i, j]) for i in range(len(rows))] for j in range(nullity)]
    basis = [_prim(v) for v in basis]
    return keys, Et, basis


def _prim(v: List[int]) -> List[int]:
    g = 0
    for x in v:
        g = gcd(g, x)
    return [x // g for x in v] if g > 1 else v


def gamma_exponents(rows: List[Dict], x: List[int]) -> Dict[tuple, int]:
    """y_P = sum_i x_i e_{iP} over atoms P = (a, b)."""
    y: Dict[tuple, int] = {}
    for xi, rw in zip(x, rows):
        if xi == 0:
            continue
        for ell, e, at in rw["atoms"]:
            key = (at["a"], at["b"])
            y[key] = y.get(key, 0) + xi * e
    return {k: v for k, v in y.items() if v}


def log_vector(y: Dict[tuple, int], logs: Dict[tuple, List]) -> List:
    n = len(next(iter(logs.values()))) if logs else 0
    out = [flint.arb(0)] * n
    for key, v in y.items():
        lv = logs[key]
        out = [o + v * l for o, l in zip(out, lv)]
    return out


def unit_relations(basis, rows, A, ctr, extra_bits=64):
    """Integer relations n with sum n_j Log(gamma_{x_j}) = 0, certified via the Blanksby-Montgomery gap."""
    d = len(A["f"]) - 1
    maxbits = max((abs(v).bit_length() for x in basis for v in x), default=1)
    ybits = maxbits + max(4, (len(rows) * 4).bit_length()) + 8
    prec = 2 * ybits + 4 * len(basis) + 128 + extra_bits
    real, cplx = embeddings(A["f"], prec)
    atoms = {(at["a"], at["b"]) for rw in rows for _, _, at in rw["atoms"]}
    logs = {k: atom_log(k[0], k[1], real, cplx) for k in atoms}
    ys = [gamma_exponents(rows, x) for x in basis]
    ctr["basis_nontrivial_gamma"] = sum(1 for y in ys if y)
    vs = [log_vector(y, logs) for y in ys]
    s, nc = len(vs), len(real) + len(cplx)
    ctr["r1r2"] = [len(real), len(cplx)]
    # norm check: every gamma_x (x in Lambda_1) is a unit, so sum of the log vector is 0
    for v in vs:
        tot = sum(v, flint.arb(0))
        if not tot.contains(0):
            _inc(ctr, "NONUNIT_KERNEL_VECTOR")
    kappa = prec // 2
    scale = flint.arb(2) ** kappa
    rowsL = []
    for j, v in enumerate(vs):
        r = [0] * s
        r[j] = 1
        r += [int((c * scale).mid().floor().unique_fmpz()) for c in v[:nc]]
        rowsL.append(r)
    Lred = flint.fmpz_mat(rowsL).lll()
    gap = bm_bound(d) / 2
    rels = []
    for i in range(Lred.nrows()):
        n = [int(Lred[i, j]) for j in range(s)]
        if not any(n):
            continue
        comb = [flint.arb(0)] * nc
        for nj, v in zip(n, vs):
            if nj:
                comb = [c + nj * vv for c, vv in zip(comb, v[:nc])]
        if all(abs(c) < gap for c in comb):
            rels.append(n)
        elif all(abs(c) > -1 for c in comb) and any(abs(c) > gap for c in comb):
            _inc(ctr, "lll_rows_nontorsion")
        elif all(not (abs(c) > gap) for c in comb):
            _inc(ctr, "unit_rel_undecided")
    ctr["kernel_dim"] = s
    ctr["unit_relations"] = len(rels)
    ctr["unit_rank_implied"] = s - len(rels)
    ctr["ybits"] = ybits
    return rels


def solve_from(eqs: List[Tuple[int, int, int]], n: int) -> Optional[int]:
    """Solve A_k + X B_k = lam_k (mod n) jointly; returns X mod n if determined, else None."""
    # Bezout: maintain cur_g = gcd(n, B_1..B_k) == sum_k cur_c[k] B_k (mod n)
    cur_g, cur_c = n, [0] * len(eqs)
    for k, (Ak, Bk, lk) in enumerate(eqs):
        gg = gcd(cur_g, Bk % n)
        if gg == cur_g:
            continue
        # gg = s*cur_g + t*Bk
        s_, t_ = _egcd(cur_g, Bk % n)
        cur_c = [s_ * c for c in cur_c]
        cur_c[k] += t_
        cur_g = gg
    if cur_g != 1:
        return None
    # sum c_k B_k = 1 (mod n)  ->  X = sum c_k (lam_k - A_k)
    X = sum(c * (lk - Ak) for c, (Ak, Bk, lk) in zip(cur_c, eqs)) % n
    return X


def _egcd(a: int, b: int) -> Tuple[int, int]:
    """s, t with s*a + t*b = gcd(a, b), where 'cur_g' is represented implicitly as a combination."""
    # cur_g itself = (previous combination of B's) + multiple of n; we track only B-coefficients, and the
    # n-multiple is irrelevant mod n.  So s multiplies the previous combination, t the new B.
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return x0, y0


def aux_primes(A: Dict, rng: random.Random, count: int = 3, bits: int = 40) -> List[Tuple[int, int]]:
    """Diagnostic: primes q with a simple root s of f mod q (q not dividing D) for exact torsion checks."""
    out = []
    while len(out) < count:
        q = rng.randrange(2 ** (bits - 1), 2 ** bits)
        if not arith.is_prime(q) or A["D"] % q == 0:
            continue
        Z = nf.roots_mod(A["f"], q)
        if Z:
            out.append((q, Z[0]))
    return out


def psi_aux(y: Dict[tuple, int], q: int, s: int) -> int:
    v = 1
    for (a, b), e in y.items():
        z = (a - b * s) % q
        if z == 0:
            return -1
        v = v * (pow(z, e, q) if e >= 0 else pow(pow(z, -1, q), -e, q)) % q
    return v


def pair_rows(p: int, A: Dict, P: Dict, ctx: Dict, count: int, box: int, ctr: Dict) -> List[Dict]:
    """NFS pair relations g1*psi(a - b alpha) = G(a, b) (no hidden exponents): T = 1, a_i = b_i = 0."""
    g0, g1, f = A["g0"], A["g1"], A["f"]
    cand = []
    for b in range(1, box + 1):
        for a in range(-box, box + 1):
            if gcd(a, b) != 1:
                continue
            G = g1 * a - g0 * b
            F = nf.hom_eval(f, a, b)
            if G and F:
                cand.append((a, b, G, F))
    sG = smooth_parts([G for _, _, G, _ in cand], ctx["primorial"])
    sF = smooth_parts([F for _, _, _, F in cand], ctx["primorial"])
    out = []
    for (a, b, G, F), x, y in zip(cand, sG, sF):
        if x == abs(G) and y == abs(F):
            fF = factor_smooth(F, ctx["primes_B"])
            if any(A["D"] % r == 0 for r in fF):
                continue
            at = {"a": a, "b": b, "G": G, "F": F, "sigma": 1, "fn": {}, "fF": fF}
            row = {"T": 1, "u0": 1, "v0": 1, "su": factor_smooth(G, ctx["primes_B"]), "sv": {},
                   "atoms": [(None, -1, at)], "a": 0, "b": 0, "pair": True}
            row["_cols"], _ = row_columns(row)
            if not row_identity_ok(p, A, row):
                _inc(ctr, "ROW_IDENTITY_VIOLATION_pair")
            out.append(row)
            if len(out) >= count:
                break
    ctr["pair_rows"] = len(out)
    return out


def postprocess(p: int, g: int, h: int, A: Dict, rows: List[Dict], ctr: Dict, dep_log: List[Dict],
                aux: Optional[List[Tuple[int, int]]] = None) -> Optional[int]:
    keys, Et, basis = exact_kernel(rows)
    ctr["cols"] = len(keys)
    if len(basis) < 2:
        _inc(ctr, "kernel_too_small"); return None
    rels = unit_relations(basis, rows, A, ctr)
    eqs = []
    xs_used: List[List[int]] = []
    n = p - 1
    # diagnostic: class-group-type obstruction  dim ker(E mod 2) - dim ker_Z(E)
    if 2 * len(rows) * len(keys) < 4_000_000:
        rk2 = flint.nmod_mat(Et, 2).rank()
        ctr["ker_mod2_minus_kerZ"] = (len(rows) - rk2) - len(basis)
    for nvec in rels:
        x = [0] * len(rows)
        for nj, xb in zip(nvec, basis):
            if nj:
                x = [xi + nj * yi for xi, yi in zip(x, xb)]
        x = _prim(x)
        if not any(x):
            continue
        # exactness checks (diagnostic): x^T E = 0
        if any(sum(Et[c][i] * x[i] for i in range(len(rows))) for c in range(len(keys))):
            _inc(ctr, "KERNEL_VIOLATION")
        W = 1
        for xi, rw in zip(x, rows):
            if xi:
                W = W * (pow(rw["T"], xi, p) if xi > 0 else pow(pow(rw["T"], -1, p), -xi, p)) % p
        Aa = sum(xi * rw["a"] for xi, rw in zip(x, rows)) % n
        Bb = sum(xi * rw["b"] for xi, rw in zip(x, rows)) % n
        if W == 1:
            lam = 0
        elif W == p - 1:
            lam = n // 2
        else:
            _inc(ctr, "TORSION_VIOLATION_W_not_pm1")   # d odd: mu_K = {+-1}
            dep_log.append({"W_ok": False})
            continue
        ygam = gamma_exponents(rows, x)
        auxv = [psi_aux(ygam, q, s_) for q, s_ in (aux or [])]
        aux_ok = all(v in (1, q - 1) for v, (q, _) in zip(auxv, aux or []))
        if not aux_ok:
            _inc(ctr, "TORSION_VIOLATION_aux")
        dep_log.append({"W_ok": True, "aux_ok": aux_ok, "gamma_nontrivial": bool(ygam), "gamma_atoms": len(ygam), "W": 1 if W == 1 else -1, "gcdB": gcd(Bb, n),
                        "xbits": max(abs(v).bit_length() for v in x), "Bpar": Bb % 2})
        eqs.append((Aa, Bb, lam))
        xs_used.append(x)
    # diagnostic: basis vectors themselves (generally non-torsion units): how often is W = +-1 ?
    nb = 0
    for xb in basis[:8]:
        W = 1
        for xi, rw in zip(xb, rows):
            if xi:
                W = W * (pow(rw["T"], xi, p) if xi > 0 else pow(pow(rw["T"], -1, p), -xi, p)) % p
        nb += W in (1, p - 1)
    ctr["basis_W_pm1"] = [nb, min(8, len(basis))]
    # diagnostic P_XK3: kernel basis vectors whose gamma is nontrivial and NOT certified torsion
    nt, ntpm = 0, 0
    for xb in basis:
        y = gamma_exponents(rows, xb)
        if not y or not aux:
            continue
        q, s_ = aux[0]
        nt += 1
        ntpm += psi_aux(y, q, s_) in (1, q - 1)
    ctr["basis_gamma_nontrivial_auxpm1"] = [ntpm, nt]
    if not eqs:
        _inc(ctr, "no_torsion_relation"); return None
    # P5 data: for small primes l | p-1, B_k mod l over a subset of relations independent mod l
    law = {}
    for l in (2, 3, 5, 7):
        if n % l:
            continue
        piv: Dict[int, List[int]] = {}
        vals = []
        for (Aa, Bb, lam), x in zip(eqs, xs_used):
            v = [xi % l for xi in x]
            for c in sorted(piv):
                if v[c]:
                    f = v[c]
                    v = [(vi - f * pi) % l for vi, pi in zip(v, piv[c])]
            nz = next((i for i, vi in enumerate(v) if vi), None)
            if nz is None:
                continue
            inv = pow(v[nz], -1, l)
            piv[nz] = [vi * inv % l for vi in v]
            vals.append(Bb % l)
            if len(vals) >= 20:
                break
        law[str(l)] = vals
    dep_log.append({"law": law})
    X = solve_from(eqs, n)
    if X is None:
        _inc(ctr, "gcd_B_nontrivial"); return None
    if pow(g, X, p) != h % p:
        _inc(ctr, "WRONG_ANSWER_BEFORE_VERIFY"); return None
    return X


# ----------------------------------------------------------------------------- top level
def dlog(p: int, g: int, h: int, cfg: Dict, seed: int) -> Dict:
    rng = random.Random(f"dlog|{p}|{g}|{h}|{seed}")
    ctr: Dict = {"p_bits": p.bit_length()}
    deps: List[Dict] = []
    t0 = time.time()
    P, ctx = make_ctx(p, cfg)
    n = p - 1
    ctr["params"] = {k: P[k] for k in ("d", "M", "B", "B1", "S", "H", "W")}

    def new_row(A):
        while True:
            _inc(ctr, "T_trials")
            if ctr["T_trials"] > P["max_rtrials"]:
                raise RuntimeError("row cap")
            a, b = rng.randrange(n), rng.randrange(n)
            T = pow(g, a, p) * pow(h, b, p) % p
            row = descend_target(p, A, P, ctx, T, rng, ctr)
            if row is not None:
                row["a"], row["b"] = a, b
                row["_cols"], _ = row_columns(row)
                if not row_identity_ok(p, A, row):
                    _inc(ctr, "ROW_IDENTITY_VIOLATION")
                _inc(ctr, "rows")
                return row

    A = None
    for _ in range(P["max_fields"]):
        _inc(ctr, "fields_sampled")
        Af = tp.sample_field(p, P, rng, ctr)
        if Af is None:
            continue
        a, b = rng.randrange(n), rng.randrange(n)
        T = pow(g, a, p) * pow(h, b, p) % p
        _inc(ctr, "T_trials")
        row = descend_target(p, Af, P, ctx, T, rng, ctr)
        if row is not None:
            row["a"], row["b"] = a, b
            row["_cols"], _ = row_columns(row)
            A, rows = Af, [row]
            break
    if A is None:
        raise RuntimeError("field cap")
    ctr["f"] = A["f"]
    X = None
    for rnd in range(P["max_rounds"]):
        _inc(ctr, "rounds")
        while True:
            cols = {k for rw in rows for k in rw["_cols"]}
            if len(rows) >= len(cols) + P["extra_rows"]:
                break
            rows.append(new_row(A))
        ctr["rows_last"] = len(rows)
        extra = pair_rows(p, A, P, ctx, cfg.get("pair_rows", 0), cfg.get("pair_box", 0), ctr) \
            if cfg.get("pair_rows", 0) else []
        aux = aux_primes(A, rng)
        X = postprocess(p, g, h, A, rows + extra, ctr, deps, aux)
        if X is not None:
            break
        rows = [new_row(A)]           # fresh hidden exponents, same field
    ctr["seconds"] = round(time.time() - t0, 2)
    return {"X": X, "ctr": ctr, "deps": deps}
