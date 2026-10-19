"""Toy end-to-end Theorem-A pipeline in design (a'''_k; mu)  (notes/thmA_proof.md §3; prereg notes/thmA_prereg.md).

OPERATIONAL: factor(N, cfg, seed) sees only N, the public config and the seed.  It never calls an integer
factorisation routine on N and never uses factors of N.  Steps (thmA_proof §3):
  FACTOR: primality (Miller-Rabin, toy stand-in for AKS), perfect powers, trial division to B, recursion.
  SPLIT:  field selection (first field whose single descent attempt succeeds), Dixon-descent rows from
          uniform hidden roots r (R = r^2 mod N, Gauss split R = u0/v0, shape test, special-q descent of
          every large prime with batch smoothness), LV6 columns + quadratic characters, GF(2) kernel,
          Hensel square root in Z[beta], congruence X^2 = Y^2, gcd(X - Y, N).
Toy deviations (disclosed in the prereg): admissible split = (B-smooth part) x (<= k_large primes in
(B, B1]) instead of "<= k primes of the dyadic J"; shape test by Pollard-Brent instead of Lenstra-Pomerance;
w, H, S, the character count and prime sizes are toy-sized; the h/window/Kern floors are not imposed; the
kernel is computed by Gaussian elimination; dependencies are tested by the lv6 §3.7 refinement (pre-sampled
uniform kernel vectors, each new one outside the span of those already tested) before fresh rows.
Every gcd with N is checked; nonunits are recorded with a classification.
"""
from __future__ import annotations

import random
import time
from fractions import Fraction
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

from . import arith
from . import numfield as nf
from .smooth import factor_smooth, primorial, smooth_parts


class FactorFound(Exception):
    def __init__(self, g: int, kind: str):
        super().__init__(f"{kind}:{g}")
        self.g, self.kind = g, kind


def _gcd_check(x: int, N: int, kind: str, ctr: Dict) -> bool:
    """True if x is a unit mod N.  A proper divisor raises FactorFound; gcd = N returns False."""
    g = gcd(x, N)
    if g == 1:
        return True
    if g == N:
        ctr["gcd_eq_N_" + kind] = ctr.get("gcd_eq_N_" + kind, 0) + 1
        return False
    if ctr.get("_continue"):
        # DIAGNOSTIC mode (E1b only): the factor is recorded, the object is treated as invalid and the
        # pipeline continues on the same N so that the DD congruence statistics can be measured.
        ctr.setdefault("accidental_continued", []).append(kind)
        return False
    raise FactorFound(g, kind)


def _inc(ctr: Dict, k: str, v: int = 1):
    ctr[k] = ctr.get(k, 0) + v


# ----------------------------------------------------------------------------- parameters
def params_for(N: int, cfg: Dict) -> Dict:
    n = N.bit_length()
    tier = None
    for t in cfg["tiers"]:
        if t["bits_max"] >= n:
            tier = t
            break
    tier = tier or cfg["tiers"][-1]
    mu = Fraction(cfg["mu"])
    M = arith.iroot(N ** mu.numerator, mu.denominator)
    B = 2 ** tier["Bbits"]
    return {
        "d": cfg["d"], "M": M, "B": B, "B1": B * tier["B1_mult"], "S": tier["S"],
        "H": 2 ** cfg["hbits"], "W": cfg["W_mult"] * M, "k_large": cfg["k_large"],
        "n_chars": cfg["n_chars"], "char_bits": cfg["char_bits"], "extra_rows": cfg["extra_rows"],
        "n_deps": cfg["n_deps"], "max_rounds": cfg["max_rounds"], "max_fields": cfg["max_fields"],
        "max_rtrials": cfg["max_rtrials"], "irr_primes": cfg["irr_primes"],
        "min_large": cfg.get("min_large", 0), "stats_all_deps": cfg.get("stats_all_deps", True),
    }


# ----------------------------------------------------------------------------- field (SAMPLE_FIELD)
def sample_field(N: int, P: Dict, rng: random.Random, ctr: Dict) -> Optional[Dict]:
    d, M, H, W = P["d"], P["M"], P["H"], P["W"]
    while True:
        g0, g1 = rng.randrange(M // 2, M + 1), rng.randrange(max(2, M // 2), M + 1)
        if gcd(g0, g1) == 1:
            break
    while True:
        w = rng.randrange(W + 1, 2 * W + 1)
        if arith.is_prime(w):
            break
    c = [rng.randrange(-H, H + 1) for _ in range(d)]
    if not _gcd_check(g1, N, "field_g1", ctr) or not _gcd_check(w, N, "field_w", ctr):
        _inc(ctr, "rej_gcdN"); return None
    f = nf.make_f(N, w, d, g0, g1, c)
    t = f[d]
    assert t != 0, "degree drop is impossible (a_d != 0 mod g1)"
    if not _gcd_check(t, N, "field_fd", ctr):
        _inc(ctr, "rej_gcdN"); return None
    e = 0
    for x in f:
        e = gcd(e, x)
    if e != 1:
        if gcd(e, N) not in (1, N):
            if ctr.get("_continue"):
                ctr.setdefault("accidental_continued", []).append("field_content")
                return None
            raise FactorFound(gcd(e, N), "field_content")
        _inc(ctr, "rej_content_w" if e == w else "rej_content_other")
        return None
    if not any(nf.irreducible_mod(f, r) for r in P["irr_primes"] if (t * N) % r):
        _inc(ctr, "rej_irreducibility_uncertified"); return None
    g1inv = pow(g1, -1, N)
    m = g0 * g1inv % N
    fp = [(i + 1) * f[i + 1] for i in range(d)]
    fpm = 0
    for cc in reversed(fp):
        fpm = (fpm * m + cc) % N
    if not _gcd_check(fpm, N, "field_fprime", ctr):
        _inc(ctr, "rej_N_divides_fprime"); return None
    D = t * nf.disc(f)
    g = nf.monic_g(f)
    return {"g0": g0, "g1": g1, "w": w, "c": c, "f": f, "t": t, "D": D, "g": g,
            "beta_img": t * m % N, "g1inv": g1inv, "sq_cache": {}}


# ----------------------------------------------------------------------------- SQ_SEARCH
def sq_search(N: int, A: Dict, ell: int, P: Dict, ctx: Dict, ctr: Dict) -> Optional[Dict]:
    cache = A["sq_cache"]
    if ell in cache:
        return cache[ell]
    _inc(ctr, "sq_searches")
    res = _sq_search(N, A, ell, P, ctx, ctr)
    cache[ell] = res
    if res is not None:
        _inc(ctr, "sq_success")
    return res


def _sq_search(N, A, ell, P, ctx, ctr):
    g0, g1, f, M, S = A["g0"], A["g1"], A["f"], P["M"], P["S"]
    if not _gcd_check(ell, N, "sq_ell", ctr):
        return None
    if g1 % ell == 0:
        _inc(ctr, "sq_ell_divides_g1"); return None
    lam = g0 * pow(g1, -1, ell) % ell
    e1, e2 = arith.gauss_reduce((ell, 0), (lam, 1))
    if 4 * S * S * (e1[0] ** 2 + e1[1] ** 2) < ell:
        _inc(ctr, "sq_irregular"); return None
    V = isqrt(ell * S * S)
    det = abs(e1[0] * e2[1] - e1[1] * e2[0])
    assert det == ell
    X = (abs(e2[0]) + abs(e2[1])) * V // ell + 1
    Yb = (abs(e1[0]) + abs(e1[1])) * V // ell + 1
    cands = []
    for x in range(-X, X + 1):
        for y in range(-Yb, Yb + 1):
            a = x * e1[0] + y * e2[0]
            b = x * e1[1] + y * e2[1]
            if not (2 * b > V and b <= V and abs(a) <= V):
                continue
            if gcd(a, b) != 1 or not ctx["isprime_small"](b):
                continue
            G = g1 * a - g0 * b
            if G == 0 or 32 * abs(G) < M * V or abs(G) > 2 * M * V:
                continue
            F = nf.hom_eval(f, a, b)
            if F == 0:
                continue
            cands.append((b, a, G, F))
    cands.sort()
    _inc(ctr, "sq_candidates", len(cands))
    if not cands:
        return None
    keep = []
    for (b, a, G, F) in cands:          # every nonunit is output (rule 3)
        if _gcd_check(G, N, "sq_G", ctr) and _gcd_check(F, N, "sq_F", ctr):
            keep.append((b, a, G, F))
    cands = keep
    if not cands:
        return None
    Prim = ctx["primorial"]
    sG = smooth_parts([abs(G) // ell for (_, _, G, _) in cands], Prim)
    sF = smooth_parts([abs(F) for (_, _, _, F) in cands], Prim)
    hits = [i for i, (bb, aa, G, F) in enumerate(cands) if sG[i] == abs(G) // ell and sF[i] == abs(F)]
    _inc(ctr, "sq_doubly_smooth", len(hits))
    if not hits:
        return None
    b, a, G, F = cands[hits[0]]
    n = abs(G) // ell
    return {"a": a, "b": b, "G": G, "F": F, "sigma": 1 if G > 0 else -1,
            "fn": factor_smooth(n, ctx["primes_B"]), "fF": factor_smooth(F, ctx["primes_B"])}


# ----------------------------------------------------------------------------- DESCEND_ONCE
def descend_once(N: int, A: Dict, P: Dict, ctx: Dict, rng: random.Random, ctr: Dict) -> Optional[Dict]:
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
            fac = arith.factor_small(cof, rng, ctx["primes_B"])      # shape test stand-in (on a split half)
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
            at = sq_search(N, A, ell, P, ctx, ctr)
            if at is None:
                _inc(ctr, "descent_fail"); return None
            atoms.append((ell, side, at))
    # diagnostic self-check of the row identity r^2 = q t^-E prod psi(z_P)^e (public data + r only)
    t, g1i, bimg = A["t"], A["g1inv"], A["beta_img"]
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


# ----------------------------------------------------------------------------- LV6
def row_vectors(row: Dict, A: Dict, ctx: Dict) -> Tuple[Dict[tuple, int], int, Dict[int, int]]:
    """Integer exponents over non-character columns, E_i, and the rational exponent vector of q_i/g1^E."""
    cols: Dict[tuple, int] = {}
    rat: Dict[int, int] = {}
    neg = 0 if row["sign"] > 0 else 1
    for p, e in row["su"].items():
        rat[p] = rat.get(p, 0) + e
    for p, e in row["sv"].items():
        rat[p] = rat.get(p, 0) - e
    E = 0
    D = A["D"]
    for ell, e, at in row["atoms"]:
        E += e
        if at["sigma"] < 0:
            neg ^= 1
        for p, v in at["fn"].items():
            rat[p] = rat.get(p, 0) - e * v
        a, b = at["a"], at["b"]
        for rr, v in at["fF"].items():
            if D % rr == 0:
                continue                                # bad prime: absorbed into H
            assert b % rr != 0
            key = ("alg", rr, a * pow(b, -1, rr) % rr)
            cols[key] = cols.get(key, 0) + e * v
    if neg:
        cols[("sign",)] = 1
    if E % 2:
        cols[("par",)] = 1
    for p, v in rat.items():
        if v % 2:
            cols[("rat", p)] = v
    return cols, E, rat


def sample_chars(A: Dict, P: Dict, rng: random.Random) -> List[Tuple[int, int]]:
    f, d, bits = A["f"], P["d"], P["char_bits"]
    out = []
    while len(out) < P["n_chars"]:
        rr = rng.randrange(2 ** (bits - 1) + 1, 2 ** bits + 1)
        if not arith.is_prime(rr):
            continue
        Z = nf.roots_mod(f, rr)
        if rng.randrange(d) < len(Z):
            out.append((rr, Z[rng.randrange(len(Z))]))
    return out


def left_kernel(vecs: List[int]) -> List[int]:
    pivots: Dict[int, Tuple[int, int]] = {}
    ker = []
    for i, v in enumerate(vecs):
        h = 1 << i
        while v:
            top = v.bit_length() - 1
            if top in pivots:
                pv, ph = pivots[top]
                v ^= pv
                h ^= ph
            else:
                pivots[top] = (v, h)
                break
        if v == 0:
            ker.append(h)
    return ker


def in_span(v: int, basis: List[int]) -> bool:
    piv: Dict[int, int] = {}
    for b in basis:
        while b:
            top = b.bit_length() - 1
            if top in piv:
                b ^= piv[top]
            else:
                piv[top] = b
                break
    while v:
        top = v.bit_length() - 1
        if top not in piv:
            return False
        v ^= piv[top]
    return True


def try_dependency(N: int, A: Dict, rows: List[Dict], Ssel: List[int], hp: int, rng: random.Random,
                   ctr: Dict) -> Optional[Tuple[int, int]]:
    """Square test + congruence for the row set Ssel.  Returns (X, Y) or None (not a square)."""
    d, g, t, g1 = len(A["f"]) - 1, A["g"], A["t"], A["g1"]
    eP: Dict[tuple, Tuple[int, Dict]] = {}
    rat: Dict[int, int] = {}
    E = 0
    neg = 0
    for i in Ssel:
        rv = rows[i]
        _, Ei, rati = rv["_vec"]
        E += Ei
        for p, v in rati.items():
            rat[p] = rat.get(p, 0) + v
        neg ^= (0 if rv["sign"] > 0 else 1)
        for ell, e, at in rv["atoms"]:
            key = (at["a"], at["b"])
            cur = eP.get(key, (0, at))[0]
            eP[key] = (cur + e, at)
            if at["sigma"] < 0:
                neg ^= 1
    assert E % 2 == 0 and neg == 0 and all(v % 2 == 0 for v in rat.values()), "column parities not killed"
    elems, l1 = [], []
    gp = nf.gprime(g)
    elems += [gp, gp]
    for (a, b), (e, at) in eP.items():
        if e:
            z = nf.z_atom(a, b, t, d)
            elems += [z] * abs(e)
            l1.append((sum(abs(x) for x in z), abs(e)))
    Gamma = nf.zprod(elems, g)
    Abound = nf.sqrt_bound(g, l1)
    st: Dict = {}
    eta = nf.zsqrt(Gamma, g, Abound, hp, rng, st)
    if eta is None:
        _inc(ctr, "dep_nonsquare")
        return None
    _inc(ctr, "dep_square")
    # Y = sqrt(q_S) * t^(-E/2) * psi(eta) * (psi(gamma^-) psi(g'(beta)))^(-1)
    sq = pow(g1, E // 2, N) if E >= 0 else pow(pow(g1, -1, N), -E // 2, N)
    for p, v in rat.items():
        sq = sq * (pow(p, v // 2, N) if v >= 0 else pow(pow(p, -1, N), -v // 2, N)) % N
    tinv = pow(t, -1, N)
    tE = pow(tinv, E // 2, N) if E >= 0 else pow(t, -E // 2, N)
    Y = sq * tE % N * nf.psi(eta, A["beta_img"], N) % N
    den = nf.psi(gp, A["beta_img"], N)
    for (a, b), (e, at) in eP.items():
        if e < 0:
            den = den * pow(t * at["G"] * A["g1inv"] % N, -e, N) % N
    _gcd_check(den, N, "psi_den", ctr)
    Y = Y * pow(den, -1, N) % N
    X = 1
    for i in Ssel:
        X = X * rows[i]["r"] % N
    if (X * X - Y * Y) % N:
        _inc(ctr, "CONGRUENCE_VIOLATION")
        return None
    return X, Y


def lv6_round(N: int, A: Dict, rows: List[Dict], P: Dict, ctx: Dict, rng: random.Random, ctr: Dict,
              dep_log: List[Dict]) -> Optional[int]:
    chars = sample_chars(A, P, rng)
    colidx: Dict[tuple, int] = {}
    vecs = []
    for rv in rows:
        cols = rv["_vec"][0]
        v = 0
        for key, val in cols.items():
            if val % 2:
                if key not in colidx:
                    colidx[key] = len(colidx)
                v |= 1 << colidx[key]
        cbits = 0
        for ell, e, at in rv["atoms"]:
            ab = at.get("_chi")
            if ab is None or at.get("_chi_key") != id(chars):
                ab = 0
                for j, (rr, s) in enumerate(chars):
                    val = (at["a"] - at["b"] * s) % rr
                    assert val != 0
                    if pow(val, (rr - 1) // 2, rr) != 1:
                        ab |= 1 << j
                at["_chi"], at["_chi_key"] = ab, id(chars)
            cbits ^= ab                                  # e = +-1: character values add mod 2
        vecs.append((v, cbits))
    ncol = len(colidx)
    full = [v | (cb << ncol) for v, cb in vecs]
    ctr["cols_last"] = ncol + P["n_chars"]
    ker = left_kernel(full)
    ctr["kernel_dim_last"] = len(ker)
    if not ker:
        _inc(ctr, "kernel_empty"); return None
    tested: List[int] = []
    answer = None
    for _ in range(P["n_deps"]):
        s = 0
        while s == 0:
            for kb in ker:
                if rng.randrange(2):
                    s ^= kb
        # verify the kernel vector
        acc = 0
        for i in range(len(rows)):
            if (s >> i) & 1:
                acc ^= full[i]
        assert acc == 0
        if in_span(s, tested):
            _inc(ctr, "dep_in_span_skipped"); continue
        Ssel = [i for i in range(len(rows)) if (s >> i) & 1]
        out = try_dependency(N, A, rows, Ssel, ctx["hensel_p"](A), rng, ctr)
        if out is None:
            continue
        X, Y = out
        tested.append(s)
        gg = gcd(X - Y, N)
        fr = 1 < gg < N
        dep_log.append({"size": len(Ssel), "fruitful": fr, "order": len(tested)})
        if fr and answer is None:
            answer = gg
        if answer is not None and not P["stats_all_deps"]:
            break
    return answer


# ----------------------------------------------------------------------------- SPLIT / FACTOR
def split(N: int, cfg: Dict, rng: random.Random, ctr: Dict, dep_log: List[Dict]) -> int:
    P = params_for(N, cfg)
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
    # Phase F: field selection
    A, first = None, None
    for _ in range(P["max_fields"]):
        _inc(ctr, "fields_sampled")
        Af = sample_field(N, P, rng, ctr)
        if Af is None:
            continue
        _inc(ctr, "fields_accepted")
        row = descend_once(N, Af, P, ctx, rng, ctr)
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
                rv["_vec"] = row_vectors(rv, A, ctx)
            colset.update(k for k, v in rv["_vec"][0].items() if v % 2)
        while True:
            need = len(colset) + P["n_chars"] + P["extra_rows"]
            if len(rows) >= need:
                break
            row = descend_once(N, A, P, ctx, rng, ctr)
            tr += 1
            if tr > P["max_rtrials"]:
                raise RuntimeError("row cap reached")
            if row is not None:
                row["_vec"] = row_vectors(row, A, ctx)
                colset.update(k for k, v in row["_vec"][0].items() if v % 2)
                rows.append(row)
        ctr["rows_last"] = len(rows)
        g = lv6_round(N, A, rows, P, ctx, rng, ctr, dep_log)
        if g is not None:
            return g
        rows = []                                        # fresh hidden roots, same field
        row = None
        while row is None:
            row = descend_once(N, A, P, ctx, rng, ctr)
        rows = [row]
    raise RuntimeError("lv6 rounds cap reached")


def factor(N: int, cfg: Dict, seed: int) -> Dict:
    """Complete factorisation of N.  Returns {'factors': [...], 'log': [...]} (operational)."""
    rng = random.Random(f"thmA|{N}|{seed}")
    log: List[Dict] = []
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
            log.append({"step": "perfect_power", "n": n, "e": pp[1]})
            stack.extend([pp[0]] * pp[1]); continue
        P = params_for(n, cfg)
        found = False
        for p in arith.primes_upto(P["B"]):
            if n % p == 0:
                e = 0
                while n % p == 0:
                    n //= p
                    e += 1
                out.extend([p] * e)
                found = True
        if found:
            log.append({"step": "trial_division", "B": P["B"]})
            stack.append(n); continue
        if n < cfg["N0"]:
            q = 2
            while n % q:
                q += 1
            log.append({"step": "below_N0", "n": n})
            stack.extend([q, n // q]); continue
        ctr: Dict = {"n_bits": n.bit_length(), "_continue": bool(cfg.get("continue_on_accidental"))}
        deps: List[Dict] = []
        t0 = time.time()
        try:
            g = split(n, cfg, rng, ctr, deps)
            kind = "STRUCTURAL_DD"
        except FactorFound as ff:
            g, kind = ff.g, "ACCIDENTAL_" + ff.kind
        ctr["seconds"] = round(time.time() - t0, 2)
        ctr.pop("_continue", None)
        log.append({"step": "split", "n": n, "factor_kind": kind, "ctr": ctr, "deps": deps})
        stack.extend([g, n // g])
    return {"factors": sorted(out), "log": log}
