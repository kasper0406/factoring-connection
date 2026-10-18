"""HARNESS / DIAGNOSTIC: fruitfulness functional Phi on the true-square dependency space K0.

Uses the known factors only to (i) classify X/Y mod each prime (needed for 3-prime N and
for consistency checks) and (ii) label instances.  Fruitfulness itself for semiprimes is
the operational test 1 < gcd(X - Y, N) < N.
"""
from __future__ import annotations
import json
import random
import time
from math import gcd
from typing import Dict, List

from . import families as fam
from .gf2 import left_kernel, coords, rank_of
from .instances import composite


def phi_class(X, Y, N, fs):
    """Class of X/Y in mu_2(Z/N)/{+-1} as a bitmask in F_2^{w-1}; DIAGNOSTIC (uses fs)."""
    if gcd(Y, N) != 1 or gcd(X, N) != 1:
        return None
    r = X * pow(Y, -1, N) % N
    bits = []
    for p in fs:
        rp = r % p
        if rp == 1:
            bits.append(0)
        elif rp == p - 1:
            bits.append(1)
        else:
            raise AssertionError("X/Y is not a square root of 1 mod p")
    last = bits[-1]
    v = 0
    for i, b in enumerate(bits[:-1]):
        v |= (b ^ last) << i
    return v


def build(family: str, N: int, seed: int, cfg: Dict):
    c = dict(cfg)
    if family == "dixon":
        return fam.dixon(N, seed, c["B"], c["excess"], c["max_trials"])
    if family == "trivial":
        return fam.trivial(N, seed, c["B"], c["excess"], c["max_trials"])
    if family == "qs":
        return fam.qs(N, seed, c["B"], c["M"])
    if family == "lin_base":
        return fam.lin_base(N, seed, c["B"], c["A"], c["Bb"], c.get("excess"))
    if family == "lin_rand":
        return fam.lin_rand(N, seed, c["B"], c["A"], c["Bb"], c.get("excess"))
    if family in ("nfs3_base", "nfs3_rand"):
        Hm = 0 if family == "nfs3_base" else c["Hmul"]
        return fam.nfs3(N, seed, c["B"], c["A"], c["Bb"], Hm, c["nchar"], family, c.get("excess"), c.get("Hshift", 0))
    raise ValueError(family)


def run_instance(family, cls, bits, seed, cfg, subset_trials=40, lin_checks=4):
    t0 = time.time()
    N, fs = composite(cls, bits, seed)
    rec = dict(family=family, cls=cls, bits=bits, seed=seed, N=N)
    try:
        sysm = build(family, N, seed, cfg)
    except fam.FactorFound as e:
        rec.update(status="factor_during_setup", g=int(e.g), accidental=True)
        return rec
    M = len(sysm.rows)
    rank, kern = left_kernel(sysm.rows)
    rec.update(M=M, ncols=sysm.ncols, rank=rank, dim_kerA=len(kern),
               info={k: (v if isinstance(v, (int, str, float)) else str(v)) for k, v in sysm.info.items()})
    if hasattr(sysm, "char_base"):        # NFS: how much do the quadratic characters cut?
        mask = (1 << sysm.char_base) - 1
        rk0, _k0 = left_kernel([r & mask for r in sysm.rows])
        rec["dim_ker_parity_only"] = M - rk0
    rng = random.Random(f"nfsfruit-harness|{family}|{cls}|{bits}|{seed}")
    phis: List[int] = []
    nonsq = 0
    fruitful_basis = 0
    for S in kern:
        cg = sysm.congruence(S)
        if cg is None:
            nonsq += 1
            phis.append(None)
            continue
        X, Y = cg
        assert (X * X - Y * Y) % N == 0, "congruence check failed"
        g = gcd(X - Y, N)
        op_fruit = 1 < g < N
        ph = phi_class(X, Y, N, fs)
        if len(fs) == 2:
            assert ph is None or bool(ph) == op_fruit
        phis.append(ph)
        fruitful_basis += int(op_fruit)
    rec["nonsquare_basis"] = nonsq
    if nonsq:
        rec["status"] = "nonsquare_in_kernel"
        rec["time"] = time.time() - t0
        return rec
    rec["dim_K0"] = len(kern)
    rec["phi_rank"] = rank_of([p for p in phis])
    rec["fruitful_basis"] = fruitful_basis
    rec["phi_basis_nonzero"] = sum(1 for p in phis if p)
    # Lemma H linearity checks (direct recomputation)
    lin_ok = lin_bad = 0
    if len(kern) >= 2:
        for _ in range(lin_checks):
            j, k = rng.sample(range(len(kern)), 2)
            cg = sysm.congruence(kern[j] ^ kern[k])
            ph = phi_class(cg[0], cg[1], N, fs)
            if ph == (phis[j] ^ phis[k]):
                lin_ok += 1
            else:
                lin_bad += 1
    rec["lin_ok"], rec["lin_bad"] = lin_ok, lin_bad
    # subset test: K0(sub) = K0 cap span(sub); Phi via coordinates in full basis
    sub = []
    if M >= 4:
        for _ in range(subset_trials):
            r = rng.randint(max(1, rank - 8), min(M, rank + 8))
            idx = sorted(rng.sample(range(M), r))
            _rk, kb = left_kernel([sysm.rows[i] for i in idx])
            vals = []
            for v in kb:
                Sg = 0
                t = 0
                while v:
                    if v & 1:
                        Sg |= 1 << idx[t]
                    v >>= 1
                    t += 1
                cvec = coords(Sg, kern)
                ph = 0
                j = 0
                while cvec:
                    if cvec & 1:
                        ph ^= phis[j]
                    cvec >>= 1
                    j += 1
                vals.append(ph)
            sub.append([len(kb), rank_of(vals)])
    rec["subsets"] = sub
    rec["status"] = "ok"
    rec["time"] = round(time.time() - t0, 3)
    return rec
