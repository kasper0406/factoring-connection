"""zroute core (OPERATIONAL: sees only N and public parameters; never uses the factors).

Deterministic Z route (notes/zroute_results.md):
  * relation sources with NO random seeds:
      - nfs_rows:  base-m monic NFS pairs (a, b) with a - b m and F(a, b) both B-smooth (seedless, d >= 3 odd);
      - orbit_rows: d = 1 "deterministic Dixon-Z": targets lambda_j = base^j mod N, Gauss split lambda = u/v,
        both halves B-smooth (row = x(lambda) - e(u) + e(v));
  * XK-Z: exact elimination of the algebraic side (ideal columns cancelled exactly over Z, units removed with
    archimedean logarithms + LLL + the Blanksby-Montgomery gap, one LLL call of dimension <= 2d per kernel vector);
  * PIN test (the small primes < 2 ln^2 N lie in the Q-span of the relation lattice), the exponent D, Miller.

Everything here is deterministic (no random module is imported).  Every inversion mod N is gcd-checked; nonunits
are raised as FactorExit (rule 3).  Floating point is used only through certified ball arithmetic (flint.arb).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from src.detGRH.core import FactorBase, FactorExit, grh_bound, miller, primes_upto  # reuse by import


# ----------------------------------------------------------------------------- helpers

def iroot(n: int, k: int) -> int:
    """floor(n^(1/k)) exactly."""
    lo, hi = 0, 1 << (n.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def gauss_reduce(N: int, lam: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Lagrange-Gauss reduced basis of {(u, v) : u = lam v (mod N)}; first vector is a shortest one."""
    b1, b2 = (N, 0), (lam % N, 1)
    def n2(x):
        return x[0] * x[0] + x[1] * x[1]
    if n2(b1) < n2(b2):
        b1, b2 = b2, b1
    while True:
        # b2 shorter; reduce b1 by b2
        dot = b1[0] * b2[0] + b1[1] * b2[1]
        nn = n2(b2)
        q = (2 * dot + nn) // (2 * nn)          # round(dot/nn)
        b1 = (b1[0] - q * b2[0], b1[1] - q * b2[1])
        if n2(b1) >= n2(b2):
            return b2, b1
        b1, b2 = b2, b1


# ----------------------------------------------------------------------------- base-m polynomial

@dataclass
class Field:
    N: int
    d: int
    m: int
    f: List[int]                 # coefficients low -> high, monic (f[d] = 1)
    disc: int

    def F(self, a: int, b: int) -> int:
        s, bp = 0, 1
        ap = [1]
        for _ in range(self.d):
            ap.append(ap[-1] * a)
        for i in range(self.d, -1, -1):           # F(a, b) = sum f_i a^i b^(d-i)
            s += self.f[i] * ap[i] * bp
            bp *= b
        return s


def base_m_field(N: int, d: int) -> Field:
    """Deterministic monic base-m polynomial: m = floor(N^(1/d)), f = base-m digits of N (leading digit 1).
    If f is reducible, f = g h gives N = g(m) h(m); a proper split is raised as FactorExit."""
    m = iroot(N, d)
    digits, x = [], N
    for _ in range(d + 1):
        digits.append(x % m)
        x //= m
    if x != 0 or digits[d] != 1:
        raise ValueError("base-m expansion not monic of degree d")
    g = gcd(m, N)
    if 1 < g < N:
        raise FactorExit(g, "m")
    P = flint.fmpz_poly(digits)
    c, facs = P.factor()
    if len(facs) > 1 or facs[0][1] > 1:
        for h, _ in facs:
            hv = int(h(m))
            g = gcd(hv, N)
            if 1 < g < N:
                raise FactorExit(g, "reducible base-m polynomial")
        raise ValueError("reducible f without a split")
    # discriminant of monic f: (-1)^(d(d-1)/2) Res(f, f')
    fp = P.derivative()
    disc = int(flint.fmpz_mat(_sylvester(digits, [int(fp[i]) for i in range(d)])).det())
    if (d * (d - 1) // 2) % 2:
        disc = -disc
    return Field(N=N, d=d, m=m, f=digits, disc=disc)


def _sylvester(f: List[int], g: List[int]) -> List[List[int]]:
    """Sylvester matrix of f (deg p) and g (deg q), coefficient lists low -> high."""
    p, q = len(f) - 1, len(g) - 1
    n = p + q
    rows = []
    fh, gh = f[::-1], g[::-1]
    for i in range(q):
        rows.append([0] * i + fh + [0] * (n - p - 1 - i))
    for i in range(p):
        rows.append([0] * i + gh + [0] * (n - q - 1 - i))
    return rows


# ----------------------------------------------------------------------------- rows

@dataclass
class Row:
    rat: List[int]                      # exponent vector over FB (column 0 = -1): the rational part
    alg: Dict[tuple, int]               # ideal columns (r, s) -> valuation (empty for d = 1 rows)
    atom: Optional[Tuple[int, int]] = None   # (a, b) for NFS rows
    src: str = ""


def small_factor(n: int, primes: Sequence[int]) -> Optional[Dict[int, int]]:
    """Full factorisation of |n| over `primes` (trial division), None if not smooth."""
    n = abs(n)
    if n == 0:
        return None
    out: Dict[int, int] = {}
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            c = 0
            while n % p == 0:
                n //= p
                c += 1
            out[p] = c
    if n > 1:
        if n > primes[-1]:
            return None
        out[n] = out.get(n, 0) + 1
    return out


def nfs_rows(fld: Field, fb: FactorBase, A: int, Bmax: int, cap: int, ctr: Dict) -> List[Row]:
    """Seedless two-sided rows, deterministic order b = 1..Bmax, a = -A..A, gcd(a, b) = 1.
    Row identity: psi(a - b alpha) = a - b m (mod N) exactly; rational part = e(a - b m) (sign in column 0),
    algebraic part = first-degree ideals (r, a/b mod r) with valuation v_r(F(a, b)).  Pairs whose F(a, b) has a
    prime factor dividing disc(f) are skipped (bad primes; disclosed toy restriction, notes §XK0)."""
    N, m = fld.N, fld.m
    primes = fb.primes
    rows: List[Row] = []
    tested = 0
    for b in range(1, Bmax + 1):
        for a in range(-A, A + 1):
            if gcd(a, b) != 1:
                continue
            tested += 1
            G = a - b * m
            if G == 0:
                continue
            eG = fb.vector(G)
            if eG is None:
                continue
            Fv = fld.F(a, b)
            if Fv == 0:
                continue
            fac = small_factor(Fv, primes)
            if fac is None:
                continue
            if any(fld.disc % r == 0 for r in fac):
                ctr["bad_prime_skips"] = ctr.get("bad_prime_skips", 0) + 1
                continue
            g = gcd(G, N)                       # gcd check (G is B-smooth, so g = 1 after preprocessing)
            if 1 < g < N:
                raise FactorExit(g, "a - bm")
            alg = {(r, a * pow(b, -1, r) % r): e for r, e in fac.items()}
            rows.append(Row(rat=eG, alg=alg, atom=(a, b), src="nfs"))
            if len(rows) >= cap:
                ctr["pairs_tested"] = tested
                return rows
    ctr["pairs_tested"] = tested
    return rows


def orbit_target(N: int, base: int, j: int) -> int:
    return pow(base, j, N)


def split_row(N: int, fb: FactorBase, lam: int, x: List[int], src: str) -> Optional[Row]:
    """d = 1 row from target lam = ev(x): Gauss split lam = u/v, both halves B-smooth.
    Row = x - e(u) + e(v) (a kernel relation).  Returns None if a half is not smooth or the split is trivial."""
    (u, v), _ = gauss_reduce(N, lam)
    if v == 0:
        return None
    for h in (u, v):
        g = gcd(h, N)
        if 1 < g < N:
            raise FactorExit(g, "split half")
    eu, ev_ = fb.vector(u), fb.vector(v)
    if eu is None or ev_ is None:
        return None
    r = [a - b + c for a, b, c in zip(x, eu, ev_)]
    if not any(r[1:]) and r[0] % 2 == 0:
        return None
    return Row(rat=r, alg={}, src=src)


def orbit_rows(N: int, fb: FactorBase, base: int, j0: int, J: int, cap: int, ctr: Dict) -> List[Row]:
    """Deterministic targets lambda_j = base^j mod N, j = j0 .. j0 + J - 1 (j0 with base^j0 > N for wrap-around).
    base must be a factor-base prime; x(lambda_j) = j e_base."""
    col = fb.index[base]
    rows: List[Row] = []
    lam = pow(base, j0, N)
    tested = 0
    for j in range(j0, j0 + J):
        tested += 1
        x = [0] * fb.n
        x[col] = j
        r = split_row(N, fb, lam, x, "orbit")
        if r is not None:
            rows.append(r)
            if len(rows) >= cap:
                break
        lam = lam * base % N
    ctr["targets_tested"] = tested
    return rows


def ev_fast(fb: FactorBase, x: Sequence[int]) -> int:
    """prod q^{x_q} mod N for x >= 0 via simultaneous binary exponentiation (bit-sliced)."""
    N = fb.N
    k = max(int(v).bit_length() for v in x) if x else 0
    acc = 1
    qs = [(-1 if i == 0 else fb.primes[i - 1]) % N for i in range(fb.n)]
    for bit in range(k - 1, -1, -1):
        acc = acc * acc % N
        for i, v in enumerate(x):
            if (v >> bit) & 1:
                acc = acc * qs[i] % N
    return acc


# ----------------------------------------------------------------------------- XK-Z: exact elimination

def bm_gap(d: int) -> flint.arb:
    """Blanksby-Montgomery: ||Log u||_inf >= log(1 + 1/(52 d log 6d))/d for units u of degree <= d not in mu_K."""
    x = flint.arb(1) / (52 * d * flint.arb(6 * d).log())
    return (1 + x).log() / d


def embeddings(f: List[int], prec: int):
    flint.ctx.prec = prec
    roots = flint.fmpz_poly(f).complex_roots()
    real = [r for r, _ in roots if r.imag.is_zero()]
    cplx = [r for r, _ in roots if not r.imag.is_zero() and r.imag > 0]
    if len(real) + 2 * len(cplx) != len(f) - 1:
        raise ArithmeticError("root isolation failed: increase precision")
    return real, cplx


def left_kernel(rows: List[Row]) -> Tuple[int, List[List[int]]]:
    """Integer basis of Lambda_1 = {x : x^T E_alg = 0} (exact; flint nullspace + LLL for small entries)."""
    keys = sorted({k for r in rows for k in r.alg})
    M = len(rows)
    if not keys:
        return 0, [[1 if i == j else 0 for i in range(M)] for j in range(M)]
    idx = {k: j for j, k in enumerate(keys)}
    Et = [[0] * M for _ in keys]
    for i, r in enumerate(rows):
        for k, v in r.alg.items():
            Et[idx[k]][i] = v
    K, nul = flint.fmpz_mat(Et).nullspace()
    basis = [[int(K[i, j]) for i in range(M)] for j in range(nul)]
    if basis:
        Lb = flint.fmpz_mat(basis).lll()
        basis = [[int(Lb[i, j]) for j in range(M)] for i in range(Lb.nrows())]
    return len(keys), basis


def _log_vec(rows: List[Row], x: List[int], real, cplx) -> List[flint.arb]:
    out = [flint.arb(0)] * (len(real) + len(cplx))
    for xi, r in zip(x, rows):
        if xi == 0:
            continue
        a, b = r.atom
        lv = [abs(flint.acb(a) - b * t).log() for t in real] + [2 * abs(flint.acb(a) - b * t).log() for t in cplx]
        out = [o + xi * l for o, l in zip(out, lv)]
    return out


def xkz(fld: Field, rows: List[Row], ctr: Dict) -> List[List[int]]:
    """Lemma XK-Z (toy implementation of the proof's algorithm).  Returns a list of integer vectors z in Z^M
    spanning (over Q) Lambda_2 = {x in Lambda_1 : gamma_x in mu_K}, one per kernel vector beyond a log-independent
    pivot set.  Each z is certified by the gap: ||sum n_j Log gamma_{x_j}||_inf < beta_d / 2."""
    d = fld.d
    C_alg, basis = left_kernel(rows)
    ctr["C_alg"] = C_alg
    ctr["dim_Lambda1"] = len(basis)
    if not basis:
        return []
    beta = max((abs(v).bit_length() for x in basis for v in x), default=1)
    prec = 4 * beta + 64 * d + 256
    real, cplx = embeddings(fld.f, prec)
    nc = len(real) + len(cplx)                  # r1 + r2 = r_U + 1
    gap = bm_gap(d)
    vecs = [_log_vec(rows, x, real, cplx) for x in basis]
    # norm check: every gamma_x (x in Lambda_1) is a unit: the logs sum to 0
    for v in vecs:
        if not sum(v, flint.arb(0)).contains(0):
            ctr["NONUNIT"] = ctr.get("NONUNIT", 0) + 1
    vv = [v[:nc - 1] if nc > 1 else [] for v in vecs]   # drop one coordinate (sum is 0): rank r_U lattice
    pivots: List[int] = []
    out: List[List[int]] = []
    undecided = 0
    for j in range(len(basis)):
        S = pivots + [j]
        if nc == 1:                              # r_U = 0: every unit is a root of unity
            out.append(basis[j])
            continue
        n_rel = _lll_relation([vv[i] for i in S], gap, prec)
        if n_rel is None:                         # certified: no relation, so v_j independent of the pivots
            if len(pivots) < nc - 1:
                pivots.append(j)
            else:
                undecided += 1                    # cannot happen: rank of Log image <= r_U
            continue
        if n_rel == "undecided":
            undecided += 1
            continue
        z = [0] * len(rows)
        for c, i in zip(n_rel, S):
            if c:
                for t in range(len(rows)):
                    z[t] += c * basis[i][t]
        g = 0
        for t in z:
            g = gcd(g, t)
        if g > 1:
            z = [t // g for t in z]
        out.append(z)
    ctr["unit_rank_seen"] = len(pivots)
    ctr["torsion_vectors"] = len(out)
    ctr["xk_undecided"] = undecided
    return out


def _lll_relation(vs: List[List[flint.arb]], gap: flint.arb, prec: int):
    """Integer relation among log vectors (XK4): LLL on (e_j, round(2^kappa v_j)); certify with the gap.
    Returns the coefficient list (last coefficient nonzero), None if certified independent, 'undecided'."""
    s = len(vs)
    kappa = prec // 3
    scale = flint.arb(2) ** kappa
    L = []
    for j, v in enumerate(vs):
        r = [0] * s
        r[j] = 1
        r += [int((c * scale).mid().floor().unique_fmpz()) for c in v]
        L.append(r)
    Lr = flint.fmpz_mat(L).lll()
    first = [int(Lr[0, j]) for j in range(s)]
    comb = [flint.arb(0)] * len(vs[0])
    for c, v in zip(first, vs):
        if c:
            comb = [x + c * y for x, y in zip(comb, v)]
    half = gap / 2
    if all(abs(x) < half for x in comb):
        if first[-1] == 0:
            return "undecided"                    # relation among pivots: impossible if pivots independent
        return first
    if any(abs(x) > half for x in comb):
        return None
    return "undecided"


# ----------------------------------------------------------------------------- outputs, PIN, Miller

def rational_relations(fb: FactorBase, rows: List[Row], zs: List[List[int]], ctr: Dict) -> List[List[int]]:
    """o(z) = sum_i z_i rat_i; psi(gamma_z) = +-1 so ev(o) = +-1; fix the sign with e_{-1}; verify exactly."""
    N, n = fb.N, fb.n
    out = []
    bad = 0
    for z in zs:
        o = [0] * n
        for zi, r in zip(z, rows):
            if zi:
                for t in range(n):
                    o[t] += zi * r.rat[t]
        val = ev_signed(fb, o)
        if val == N - 1:
            o[0] += 1
            val = 1
        if val != 1:
            if val * val % N == 1:
                g = gcd(val - 1, N)
                if 1 < g < N:
                    raise FactorExit(g, "torsion value")
            bad += 1
            continue
        if any(o):
            out.append(o)
    ctr["rel_verify_fail"] = bad
    return out


def ev_signed(fb: FactorBase, e: Sequence[int]) -> int:
    """prod q^{e_q} mod N for signed exponents (gcd-checked inverses; FB primes are units)."""
    N = fb.N
    x = N - 1 if e[0] % 2 else 1
    for p, i in fb.index.items():
        k = e[i]
        if k:
            x = x * pow(p, k, N) % N
    return x


def pin_and_exponent(fb: FactorBase, rels: List[List[int]], qmax: int) -> Dict:
    """Rank, PIN (every prime q <= qmax lies in the Q-span of the relation lattice), and an exponent D.

    H := HNF of the relations (plus 2 e_{-1}); rk := rank.  e_q is pinned iff rank(H + e_q) = rk.  With P the rk
    pivot columns of H, D := |det H_P| = product of the pivots; by Cramer, D e_q lies in the lattice for every pinned
    q (and D is a multiple of [Z^n : L] when rk = n).  So q^D = 1 (mod N) for every pinned q (Prop PIN)."""
    n = fb.n
    rels = rels + [[2] + [0] * (n - 1)]
    H = flint.fmpz_mat(rels).hnf()
    Hr = [[int(x) for x in r] for r in H.tolist() if any(int(x) for x in r)]
    rk = len(Hr)
    D = 1
    for r in Hr:
        piv = next(x for x in r if x)
        D *= abs(piv)
    small = [q for q in fb.primes if q <= qmax]
    if rk == n:
        return {"rank": rk, "n": n, "pinned": len(small), "n_small": len(small), "PIN": True, "D": D,
                "full": True, "unpinned": []}
    unpinned = []
    for q in small:
        j = fb.index[q]
        e = [0] * n
        e[j] = 1
        if flint.fmpz_mat(Hr + [e]).rank() > rk:
            unpinned.append(q)
    return {"rank": rk, "n": n, "pinned": len(small) - len(unpinned), "n_small": len(small),
            "PIN": not unpinned, "D": D, "full": False, "unpinned": unpinned}


def miller_split(N: int, D: int) -> Tuple[Optional[int], Optional[int]]:
    """Deterministic Miller test with bases a = 2 .. ceil(2 ln^2 N) (the ERH range)."""
    return miller(N, D, grh_bound(N))
