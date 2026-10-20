"""POST HOC (disclosed): diagnosis of the E2 P7 failure at (X = 2^13, z = 20).

For every truth output (a, b) of that cell we record d = gcd(G(a, b), P) and compare, divisor by divisor, the
observed number of box points of the 2-D lattice Lambda_d = {(a, b) : d | Y1 a + Y0 b} in the E2 box
[-X, X] x {first 200 admissible b} with the random-model expectation (#rows)(2X+1)/d (times the exact-gcd
factor). The shortest vector of Lambda_d (Gauss reduction) explains excesses: a lattice whose shortest vector
is short in the b-direction puts many points in a box only 200 rows tall."""
from __future__ import annotations

import json
from collections import Counter
from math import gcd, log, prod
from pathlib import Path

from src.pilist.decoder import boneh_threshold, primes_upto

ROOT = Path(__file__).resolve().parents[2]


def gauss_reduce(u, v):
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        m = round((u[0] * v[0] + u[1] * v[1]) / n2(u))
        v = (v[0] - m * u[0], v[1] - m * u[1])
        if n2(v) >= n2(u):
            return u, v
        u, v = v, u


def diag(X=8192, z=20, nrows=200):
    meta = json.loads((ROOT / "experiments/raw/cadobench_testbed_c60/meta.json").read_text())
    Y0, Y1 = int(meta["f0_rational_low_to_high"][0]), int(meta["f0_rational_low_to_high"][1])
    Pd = {q: 1 for q in primes_upto(z) if Y1 % q}
    P = prod(Pd)
    bs = [b for b in range(1, 100000) if gcd(b, Y1) == 1][:nrows]
    thr = boneh_threshold(log(2 * X + 1), Pd)
    obs = Counter()
    for b in bs:
        for a in range(-X, X + 1):
            g = gcd(Y1 * a + Y0 * b, P)
            if log(g) >= thr - 1e-12:
                obs[g] += 1
    out = []
    for d, n in sorted(obs.items(), key=lambda t: -t[1]):
        # exact-gcd expectation: points with d | G and gcd(G/d, P/d) = 1
        exp_ = len(bs) * (2 * X + 1) / d * prod(1 - 1 / q for q in Pd if (P // d) % q == 0)
        root = (-Y0 * pow(Y1, -1, d)) % d          # a == root * b (mod d)
        u, v = gauss_reduce((d, 0), (root, 1))
        out.append({"d": d, "obs": n, "exp": exp_, "short": u, "second": v})
    tot_obs = sum(obs.values())
    return {"thr": thr, "total_obs": tot_obs, "by_d": out}


if __name__ == "__main__":  # pragma: no cover
    r = diag()
    (ROOT / "experiments/raw/pilist_E2_posthoc.json").write_text(json.dumps(r, indent=1))
    print("total", r["total_obs"])
    for x in r["by_d"][:12]:
        print(x)


def coprime_recount():
    """POST HOC: E2 truth counts restricted to primitive points gcd(a, b) = 1, against the model prediction
    times the coprime fraction of the cell's box."""
    cfg = json.loads((ROOT / "experiments/configs/pilist.json").read_text())["E2"]
    meta = json.loads((ROOT / "experiments/raw/cadobench_testbed_c60/meta.json").read_text())
    Y0, Y1 = int(meta["f0_rational_low_to_high"][0]), int(meta["f0_rational_low_to_high"][1])
    from src.pilist.decoder import count_model_fraction, log_P
    rows = []
    for X, nrows in cfg["X_rows"]:
        bs = [b for b in range(1, 100000) if gcd(b, Y1) == 1][:nrows]
        ncop = sum(1 for b in bs for a in range(-X, X + 1) if gcd(a, b) == 1)
        for z in cfg["z"]:
            Pd = {q: 1 for q in primes_upto(z) if Y1 % q}
            P = prod(Pd)
            thr = boneh_threshold(log(2 * X + 1), Pd)
            if thr > log_P(Pd):
                continue
            hi = count_model_fraction(Pd, thr)[1]
            pred = hi * ncop
            if pred < 10:
                continue
            obs = sum(1 for b in bs for a in range(-X, X + 1)
                      if gcd(a, b) == 1 and log(gcd(Y1 * a + Y0 * b, P)) >= thr - 1e-12)
            rows.append({"X": X, "z": z, "obs_coprime": obs, "pred_coprime": pred, "ratio": obs / pred})
            print(rows[-1], flush=True)
    return rows


def primitive_model_check():
    """POST HOC: primitive-pair (projective) densities Pr[q | G(a,b)] = 1/(q+1) for gcd(a,b) = 1, linear G."""
    from src.pilist.decoder import count_model_fraction
    prev = json.loads((ROOT / "experiments/raw/pilist_E2_coprime_posthoc.json").read_text())
    cfg = json.loads((ROOT / "experiments/configs/pilist.json").read_text())["E2"]
    meta = json.loads((ROOT / "experiments/raw/cadobench_testbed_c60/meta.json").read_text())
    Y1 = int(meta["f0_rational_low_to_high"][1])
    out = []
    for r in prev:
        X, z = r["X"], r["z"]
        nrows = dict(map(tuple, cfg["X_rows"]))[X]
        bs = [b for b in range(1, 100000) if gcd(b, Y1) == 1][:nrows]
        ncop = sum(1 for b in bs for a in range(-X, X + 1) if gcd(a, b) == 1)
        Pd = {q: 1 for q in primes_upto(z) if Y1 % q}
        dens = {q: [q / (q + 1), 1 / (q + 1)] for q in Pd}
        thr = boneh_threshold(log(2 * X + 1), Pd)
        lo, hi = count_model_fraction(Pd, thr, dens=dens)
        out.append({**r, "pred_primitive": hi * ncop, "ratio_primitive": r["obs_coprime"] / (hi * ncop)})
        print(out[-1], flush=True)
    return out
