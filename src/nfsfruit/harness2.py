"""HARNESS for the E1b / E2 / E3 addendum (notes/nfsfruit_prereg.md, Addendum A).

E1b: prefix-nullity test.  Process relations in a random order; the k-th kernel vector found
     spans, with the earlier ones, the prefix dependency space of dimension k.  Record whether
     Phi restricted to it is rank-deficient, for k = 1..12.
E2 : Jacobi-row identities for Blum n (Prop. J): QS  c_i = [J(x_i) = -1] + [(2/n) = -1] v2(z_i)/2;
     degree-1 NFS  c_i = floor(t_i/2), t_i = Omega_T(L1 L2), T = {l : J_n(l) = -1}.  These rows
     are computed OPERATIONALLY (Jacobi symbols mod n only); identity c(S) = Phi(S) is checked on
     every basis vector.
E3 : square-multiplier linear sieve (Theorem L1): hidden-root check for every basis vector and
     every e in mu_2(Z/n); one uniform odd and one uniform even dependency per instance.
"""
from __future__ import annotations
import random
import time
from math import gcd

import gmpy2

from . import families as fam
from .gf2 import left_kernel, coords, rank_of
from .harness import build, phi_class
from .instances import composite


def _phi_of(S, kern, phis):
    c = coords(S, kern)
    ph = 0
    j = 0
    while c:
        if c & 1:
            ph ^= phis[j]
        c >>= 1
        j += 1
    return ph


def _basis_phis(sysm, kern, N, fs):
    out = []
    for S in kern:
        cg = sysm.congruence(S)
        if cg is None:
            return None
        X, Y = cg
        assert (X * X - Y * Y) % N == 0
        out.append(phi_class(X, Y, N, fs))
    return out


def prefix_test(rows, kern, phis, w1, rng, orderings, kmax=12):
    """Returns list over orderings of lists [deficient(k) for k=1..kmax] (None if not reached)."""
    M = len(rows)
    res = []
    for _ in range(orderings):
        perm = list(range(M))
        rng.shuffle(perm)
        _rk, kb = left_kernel([rows[i] for i in perm])
        vals = []
        out = []
        for k, v in enumerate(kb[:kmax], start=1):
            Sg = 0
            t = 0
            while v:
                if v & 1:
                    Sg |= 1 << perm[t]
                v >>= 1
                t += 1
            vals.append(_phi_of(Sg, kern, phis))
            out.append(int(rank_of(vals) < w1))
        out += [None] * (kmax - len(out))
        res.append(out)
    return res


def run_e1b(family, cls, bits, seed, cfg, orderings=5):
    N, fs = composite(cls, bits, seed)
    rec = dict(exp="e1b", family=family, cls=cls, bits=bits, seed=seed)
    try:
        sysm = build(family, N, seed, cfg)
    except fam.FactorFound as e:
        rec.update(status="factor_during_setup", g=int(e.g))
        return rec
    _rk, kern = left_kernel(sysm.rows)
    phis = _basis_phis(sysm, kern, N, fs)
    if phis is None:
        rec["status"] = "nonsquare_in_kernel"
        return rec
    rng = random.Random(f"nfsfruit-e1b|{family}|{cls}|{bits}|{seed}")
    rec["dim_K0"] = len(kern)
    rec["prefix"] = prefix_test(sysm.rows, kern, phis, len(fs) - 1, rng, orderings)
    rec["status"] = "ok"
    return rec


def jacobi_row(sysm, N):
    """OPERATIONAL (uses only N): the explicit factor-detecting row of Prop. J."""
    J = lambda x: int(gmpy2.jacobi(x % N, N))
    if sysm.family == "qs":
        j2 = J(2)
        row = []
        for x, fac in zip(sysm.xs, sysm.facs):
            c = int(J(x) == -1)
            if j2 == -1:
                v2 = fac.get(2, 0)
                assert v2 % 2 == 0
                c ^= (v2 // 2) & 1
            row.append(c)
        return row
    if sysm.family in ("lin_base", "lin_rand"):
        T = {p for p in sysm.fb.primes if J(p) == -1}
        row = []
        for f1, f2 in zip(sysm.f1, sysm.f2):
            t = sum(e for p, e in f1.items() if p in T) + sum(e for p, e in f2.items() if p in T)
            row.append((t // 2) & 1)
        return row
    raise ValueError(sysm.family)


def run_e2(family, bits, seed, cfg):
    N, fs = composite("blum", bits, seed)
    rec = dict(exp="e2", family=family, cls="blum", bits=bits, seed=seed)
    try:
        sysm = build(family, N, seed, cfg)
    except fam.FactorFound as e:
        rec.update(status="factor_during_setup", g=int(e.g))
        return rec
    _rk, kern = left_kernel(sysm.rows)
    phis = _basis_phis(sysm, kern, N, fs)
    row = jacobi_row(sysm, N)
    agree = dis = 0
    for S, ph in zip(kern, phis):
        c = 0
        i = 0
        T = S
        while T:
            if T & 1:
                c ^= row[i]
            T >>= 1
            i += 1
        if c == ph:
            agree += 1
        else:
            dis += 1
    rec.update(dim_K0=len(kern), agree=agree, disagree=dis, phi_rank=rank_of(phis),
               row_weight=sum(row), M=len(row))
    # degeneracy illustration: restrict to relations with c_i = 0, and to c_i = 1
    for lab, val in (("sub_c0", 0), ("sub_c1", 1)):
        idx = [i for i in range(len(row)) if row[i] == val]
        _r, kb = left_kernel([sysm.rows[i] for i in idx])
        vals = []
        for v in kb:
            Sg = 0
            t = 0
            while v:
                if v & 1:
                    Sg |= 1 << idx[t]
                v >>= 1
                t += 1
            vals.append(_phi_of(Sg, kern, phis))
        rec[lab] = [len(kb), rank_of(vals)]
    rec["status"] = "ok"
    return rec


def mu2(fs):
    """HARNESS: all square roots of 1 mod prod(fs) (squarefree), via CRT."""
    N = 1
    for p in fs:
        N *= p
    roots = [1]
    mod = 1
    for p in fs:
        new = []
        for r in roots:
            for s in (1, p - 1):
                # combine r mod `mod` with s mod p
                x = r + mod * ((s - r) * pow(mod, -1, p) % p)
                new.append(x % (mod * p))
        roots = new
        mod *= p
    return roots


def run_e3(cls, bits, seed, cfg):
    N, fs = composite(cls, bits, seed)
    rec = dict(exp="e3", family="lin_sq", cls=cls, bits=bits, seed=seed)
    try:
        sysm = fam.lin_sq(N, seed, cfg["B"], cfg["A"], cfg["Bb"], cfg.get("excess"))
    except fam.FactorFound as e:
        rec.update(status="factor_during_setup", g=int(e.g))
        return rec
    _rk, kern = left_kernel(sysm.rows)
    w1 = len(fs) - 1
    G = mu2(fs)
    r0 = sysm._r
    phis = []
    hidden_ok = hidden_bad = 0
    nodd = 0
    for S in kern:
        X, Y = sysm.congruence(S)
        assert (X * X - Y * Y) % N == 0
        phis.append(phi_class(X, Y, N, fs))
        odd = bin(S).count("1") % 2 == 1
        nodd += odd
        classes = []
        for e in G:
            Xe, Ye = sysm.congruence(S, r_override=r0 * e % N)
            classes.append(phi_class(Xe, Ye, N, fs))
        counts = {}
        for c in classes:
            counts[c] = counts.get(c, 0) + 1
        if odd:
            ok = len(counts) == 2 ** w1 and all(v == 2 for v in counts.values())
        else:
            ok = len(counts) == 1
        hidden_ok += ok
        hidden_bad += (not ok)
    rec.update(dim_K0=len(kern), odd_basis=nodd, hidden_ok=hidden_ok, hidden_bad=hidden_bad,
               phi_rank=rank_of(phis))
    even = [S for S in kern if bin(S).count("1") % 2 == 0]
    rec["dim_even"] = len(kern) - (1 if nodd else 0)
    rng = random.Random(f"nfsfruit-e3|{cls}|{bits}|{seed}")
    if kern:
        v = 0
        for S in kern:
            if rng.random() < 0.5:
                v ^= S
        o = next((S for S in kern if bin(S).count("1") % 2), None)
        if o is not None:
            vo = v if bin(v).count("1") % 2 else v ^ o
            ve = v ^ o if bin(v).count("1") % 2 else v
            for lab, S in (("odd", vo), ("even", ve)):
                if S == 0:
                    rec[lab + "_fruit"] = None
                    continue
                X, Y = sysm.congruence(S)
                rec[lab + "_fruit"] = int(phi_class(X, Y, N, fs) != 0)
                rec[lab + "_class"] = phi_class(X, Y, N, fs)
    rec["status"] = "ok"
    return rec
