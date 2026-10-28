"""HARNESS for E1 (QS-degenerate base + h Dixon-descended rows).  Factors of N are used only to
label instances (residue class) via nfsfruit.instances; fruitfulness is the operational gcd test.

Usage: python -m src.fruit.harness_e1 <config.json> <out.jsonl>
"""
from __future__ import annotations

import json
import random
import sys
import time

import gmpy2

from ..lemmaM import descent as dsc
from ..nfsfruit.instances import composite
from .linalg import hybrid_rank, kernel_basis, sample_forced, sample_nonzero
from .qs_hybrid import QSHybrid, fruitful


def instances(cls: str, bits: int, count: int):
    want = 1 if cls == "blum1" else 5
    out, seed = [], 0
    while len(out) < count:
        N, fs = composite("blum", bits, seed)
        if N % 8 == want:
            out.append((seed, N, fs))
        seed += 1
    return out


def beta_zero(x: int, N: int) -> bool:
    """nfsfruit Cor J1 (Blum N): beta = [J_N(x) = -1] + [N = 5 mod 8][x odd]."""
    b = (1 if gmpy2.jacobi(x, N) == -1 else 0) ^ (1 if (N % 8 == 5 and x & 1) else 0)
    return b == 0


class QSDegenerate(QSHybrid):
    def add_qs_beta0(self, scan: int, limit: int) -> int:
        N, x0 = self.N, self.x0
        added = 0
        for t in range(scan):
            if added >= limit:
                break
            x = x0 + t
            if gmpy2.jacobi(x, N) == 0:
                continue
            if not beta_zero(x, N):
                continue
            Q = x * x - N
            self.counters["qs_tests"] += 1
            f = self._fbfac(Q)
            if f is None:
                continue
            self.X.append(x % N)
            self.E.append(f)
            self.rows.append(self._row(f))
            self.kind.append("Q")
            self.hidden_r.append(None)
            added += 1
        return added


def run_instance(cls, bits, seed, N, fs, p, cfg):
    t0 = time.time()
    s = QSDegenerate(N, p["B"], p["B1"], p["J"])
    Ct = 1 + sum(1 for q in s.fb.primes if q == 2 or gmpy2.legendre(N % q, q) == 1)
    Mlev = {name: k * Ct + add for name, k, add in cfg["levels"]}
    nq = s.add_qs_beta0(10 ** 9, max(Mlev.values()))
    nbase = len(s.rows)
    recs = []
    rng = random.Random(f"fruit-e1|{N}|{seed}")
    for lev, M in Mlev.items():
        M = min(M, nbase)
        base_idx = list(range(M))
        base_rows = s.rows[:M]
        # h = 0: base only
        kern0 = kernel_basis(base_rows)
        out0 = []
        for _ in range(cfg["rounds"]):
            S = sample_nonzero(kern0, rng)
            if S is None:
                out0.append(None)
                continue
            X, Y = s.sub_congruence(base_idx, S)
            assert (X * X - Y * Y) % N == 0
            out0.append(fruitful(X, Y, N))
        recs.append(dict(cls=cls, bits=bits, seed=seed, N=str(N), omega=len(fs), lev=lev, M=M, Ct=Ct,
                         h=0, round=-1, dimK=len(kern0), base_uniform=out0))
        for rd in range(cfg["rounds"]):
            first = len(s.rows)
            try:
                for _ in range(cfg["hmax"]):
                    if not s.add_descent_exact(rng, cfg["max_trials"]):
                        raise RuntimeError("descent budget exhausted")
            except dsc.FactorFound as e:
                recs.append(dict(cls=cls, bits=bits, seed=seed, lev=lev, round=rd, accidental_factor=e.g))
                continue
            didx = list(range(first, first + cfg["hmax"]))
            for h in cfg["hs"]:
                idx = base_idx + didx[:h]
                rows = [s.rows[i] for i in idx]
                hr = hybrid_rank(base_rows, [s.rows[i] for i in didx[:h]])
                kern = kernel_basis(rows)
                rec = dict(cls=cls, bits=bits, seed=seed, N=str(N), omega=len(fs), lev=lev, M=M, Ct=Ct,
                           round=rd, dimK=len(kern), **hr)
                Sf = sample_forced(kern, M, rng)
                if Sf is not None:
                    X, Y = s.sub_congruence(idx, Sf)
                    assert (X * X - Y * Y) % N == 0
                    rec["forced"] = fruitful(X, Y, N)
                else:
                    rec["forced"] = None
                Su = sample_nonzero(kern, rng)
                if Su is not None:
                    X, Y = s.sub_congruence(idx, Su)
                    assert (X * X - Y * Y) % N == 0
                    rec["uniform"] = fruitful(X, Y, N)
                    rec["uniform_uses_desc"] = bool(Su >> M)
                recs.append(rec)
    for r in recs:
        r["secs_instance"] = round(time.time() - t0, 2)
        r["counters"] = dict(s.counters)
    return recs


def main():
    cfg = json.load(open(sys.argv[1]))
    with open(sys.argv[2], "w") as fo:
        for bits in cfg["bits"]:
            p = cfg[str(bits)]
            for cls in cfg["classes"]:
                for seed, N, fs in instances(cls, bits, cfg["count"]):
                    for r in run_instance(cls, bits, seed, N, fs, p, cfg):
                        fo.write(json.dumps(r) + "\n")
                    fo.flush()


if __name__ == "__main__":
    main()
