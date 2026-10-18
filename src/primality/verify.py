"""HARNESS: exhaustive verification of the exact leak formulas and torus identities.

Run:  .venv/bin/python -m src.primality.verify [--bound 4000] [--lucas-bound 3000]
      [--frob-bound 400] [--workers 32] [--out experiments/raw/primality_verify.json]

For every N in the verification sets and EVERY base, the operational tests in
``ops`` are run and their leaks tallied; the tallies must equal |group| times the
exact probabilities of ``formulas`` (generic order-pattern engine AND closed forms).
"""
from __future__ import annotations

import argparse
import json
import random
import time
from fractions import Fraction
from math import gcd
from multiprocessing import Pool

import sympy

from src.common.instances import CARMICHAEL_SMALL
from src.primality import formulas as F
from src.primality import ops


def odd_primes_upto(B):
    return [p for p in sympy.primerange(3, B)]


def semiprimes(bound):
    ps = odd_primes_upto(bound // 3 + 1)
    out = []
    for i, p in enumerate(ps):
        for q in ps[i + 1:]:
            if p * q >= bound:
                break
            out.append((p, q))
    return out


def three_prime(bound):
    ps = odd_primes_upto(bound // 15 + 1)
    out = []
    for i, p in enumerate(ps):
        for j in range(i + 1, len(ps)):
            q = ps[j]
            if p * q * q >= bound:
                break
            for r in ps[j + 1:]:
                if p * q * r >= bound:
                    break
                out.append((p, q, r))
    return out


# ------------------------------------------------------------- power tests
def check_power_tests(primes):
    N = 1
    for p in primes:
        N *= p
    tal = {"prefix": 0, "sqrt": 0, "pm": 0, "liar": 0, "fermat": 0, "fermat_chain": 0, "ss": 0,
           "sqrt_minus1": 0, "units": 0, "minus1_final": 0}
    levels = {}
    for a in range(1, N):
        if gcd(a, N) != 1:
            continue
        tal["units"] += 1
        mr = ops.miller_rabin(N, a)
        st = [l for l in mr["leaks"] if l["class"] == "structural"]
        if any(l["kind"] == "minus1" for l in st):
            tal["prefix"] += 1
        if st:
            tal["pm"] += 1
        if mr["sqrt_split"] is not None and 1 < mr["sqrt_split"] < N:
            tal["sqrt"] += 1
        if mr["pattern"] != "NONE":
            tal["liar"] += 1
            levels[mr["pattern"]] = levels.get(mr["pattern"], 0) + 1
            if mr["pattern"] not in ("S0", "L1"):
                # b_(j-1) = -1 with j-1 >= 1: b_(j-2) is a square root of -1 mod N
                j = int(mr["pattern"][1:])
                x = mr["seq"][j - 2]
                assert x * x % N == N - 1
                tal["sqrt_minus1"] += 1
        if mr["seq"][-1] == N - 1:
            tal["minus1_final"] += 1
        if ops.fermat(N, a)["leaks"]:
            tal["fermat"] += 1
        if ops.fermat(N, a, chain=True)["leaks"]:
            tal["fermat_chain"] += 1
        if [l for l in ops.solovay_strassen(N, a)["leaks"] if l["class"] == "structural"]:
            tal["ss"] += 1
    U = tal["units"]
    ex = F.mr_exact(primes)
    pred = {"prefix": ex["prefix"] * U, "sqrt": ex["sqrt"] * U, "pm": ex["pm"] * U,
            "liar": ex["liar"] * U, "fermat": F.fermat_exact(primes) * U,
            "fermat_chain": F.fermat_exact(primes, chain=True) * U,
            "ss": F.ss_exact(primes) * U, "minus1_final": 0}
    lv_pred = {("S0" if j == 0 else f"L{j}"): pr * U for j, pr in ex["liar_level"].items()}
    ok = all(pred[k] == tal[k] for k in pred) and all(
        lv_pred.get(k, 0) == levels.get(k, 0) for k in set(lv_pred) | set(levels))
    if len(primes) == 2:
        p, q = primes
        cl = F.closed_mr_semiprime(p, q)
        ok_closed = (cl["prefix"] * U == tal["prefix"] and cl["sqrt"] * U == tal["sqrt"]
                     and cl["liar"] * U == tal["liar"] and cl["sqrt_minus1"] * U == tal["sqrt_minus1"]
                     and F.closed_ss_semiprime(p, q) * U == tal["ss"]
                     and cl["single"] * U == pred["fermat"])
        gl = F.closed_ladder(p - 1, q - 1, N - 1)
        ok_closed &= gl["prefix"] * U == tal["prefix"] and gl["sqrt"] * U == tal["sqrt"]
    else:
        ok_closed = None
    return {"primes": list(primes), "N": N, "ok": ok, "ok_closed": ok_closed, "tally": tal,
            "pred": {k: str(v) for k, v in pred.items()}, "levels": levels}


# ------------------------------------------------------------------- Lucas
def check_lucas(pq):
    p, q = pq
    N = p * q
    sel = ops.selfridge(N)
    if "D" not in sel or "P" not in sel:
        return {"primes": [p, q], "skip": sel}
    D = sel["D"]
    inv4 = pow(4, -1, N)
    tal = {"prefix": 0, "pm": 0, "liar": 0, "valid": 0, "lucas": 0}
    for P in range(N):
        Q = (P * P - D) * inv4 % N
        if gcd(Q, N) != 1:
            continue
        tal["valid"] += 1
        r = ops.strong_lucas(N, P, Q)
        st = [l for l in r["leaks"] if l["class"] == "structural"]
        if any(l["kind"] == "U" for l in st):
            tal["prefix"] += 1
        if st:
            tal["pm"] += 1
        if r["pattern"] != "NONE":
            tal["liar"] += 1
        if any(l["kind"] == "U" and l["where"] == r["t"] for l in st):
            tal["lucas"] += 1
    ex = F.lucas_exact([p, q], D, exclude_identity=True)
    V = tal["valid"]
    dp, dq = int(sympy.legendre_symbol(D % p, p)), int(sympy.legendre_symbol(D % q, q))
    assert V == (p - dp - 1) * (q - dq - 1)
    lucas_plain = F.leak_prob([F.indicator_dist(pp - int(sympy.legendre_symbol(D % pp, pp)), [(N + 1, 1)], True)
                               for pp in (p, q)])
    pred = {"prefix": ex["prefix"] * V, "pm": ex["pm"] * V, "liar": ex["liar"] * V,
            "lucas": lucas_plain * V}
    ok = all(pred[k] == tal[k] for k in pred)
    # closed form for uniform gamma in T_N (identity included), Thm 3 with n = p - (D/p), e = N+1
    gl = F.closed_ladder(p - dp, q - dq, N + 1)
    ex_u = F.lucas_exact([p, q], D, exclude_identity=False)
    ok_closed = gl["prefix"] == ex_u["prefix"] and gl["sqrt"] == ex_u["sqrt"] and gl["liar"] == ex_u["liar"]
    return {"primes": [p, q], "N": N, "D": D, "legendre": [dp, dq], "ok": ok, "ok_closed_uniform": ok_closed,
            "tally": tal, "pred": {k: str(v) for k, v in pred.items()}}


# --------------------------------------------------------------- Frobenius
def check_frobenius(pq):
    p, q = pq
    N = p * q
    sel = ops.selfridge(N)
    if "P" not in sel:
        return {"primes": [p, q], "skip": True}
    Dl = sel["D"]
    tal = {"structural": 0, "accidental": 0, "units": 0}
    for a0 in range(N):
        for a1 in range(N):
            if gcd((a0 * a0 - Dl * a1 * a1) % N, N) != 1:
                continue
            tal["units"] += 1
            r = ops.frobenius_quadratic(N, Dl, (a0, a1))
            cl = {l["class"] for l in r["leaks"]}
            tal["structural"] += "structural" in cl
            tal["accidental"] += ("accidental" in cl and "structural" not in cl)
    pred = F.frobenius_exact(p, q, Dl) * tal["units"]
    return {"primes": [p, q], "N": N, "Delta": Dl, "ok": pred == tal["structural"], "tally": tal,
            "pred": str(pred)}


# -------------------------------------------------------- torus identities
def qmul(A, B, N, Dl):
    return ((A[0] * B[0] + Dl * A[1] * B[1]) % N, (A[0] * B[1] + A[1] * B[0]) % N)


def qpow(A, e, N, Dl):
    R = (1, 0)
    while e:
        if e & 1:
            R = qmul(R, A, N, Dl)
        A = qmul(A, A, N, Dl)
        e >>= 1
    return R


def conj(A, N):
    return (A[0] % N, (-A[1]) % N)


def check_torus(args):
    """Thm 5/6 on Z/N[w]/(w^2 - D): for random units alpha and gamma = alpha/sigma(alpha):
    gamma^N = gamma^(dq p + dp q - dp dq);  alpha^N = sigma(alpha)^q mod inert p, alpha^q mod split p;
    both inert: alpha^(N+1) = sigma(alpha)^s."""
    (p, q), Dl, trials, seed = args
    N = p * q
    dp, dq = int(sympy.legendre_symbol(Dl % p, p)), int(sympy.legendre_symbol(Dl % q, q))
    rng = random.Random(seed)
    E = F.torus_exponent(p, q, dp, dq)
    bad = 0
    for _ in range(trials):
        while True:
            A = (rng.randrange(N), rng.randrange(N))
            nu = (A[0] * A[0] - Dl * A[1] * A[1]) % N
            if gcd(nu, N) == 1:
                break
        inv_nu = pow(nu, -1, N)
        a2 = qmul(A, A, N, Dl)
        gam = (a2[0] * inv_nu % N, a2[1] * inv_nu % N)  # alpha/sigma(alpha) = alpha^2/nu
        lhs = qpow(gam, N, N, Dl)
        rhs = qpow(gam, E, N, Dl) if E >= 0 else qpow(conj(gam, N), -E, N, Dl)
        bad += lhs != rhs
        AN = qpow(A, N, N, Dl)
        for pr, d, other in ((p, dp, q), (q, dq, p)):
            base = conj(A, N) if d == -1 else A
            r = qpow(base, other, N, Dl)
            bad += (AN[0] - r[0]) % pr != 0 or (AN[1] - r[1]) % pr != 0
        if dp == dq == -1:
            bad += qpow(A, N + 1, N, Dl) != qpow(conj(A, N), p + q, N, Dl)
    return {"primes": [p, q], "D": Dl, "pattern": [dp, dq], "exponent": E, "trials": trials,
            "ok": bad == 0}


def torus_cases(bound, per_pattern=2, seed=0):
    rng = random.Random(seed)
    out = []
    for p, q in semiprimes(bound):
        found = {}
        for Dl in [d for d in range(-60, 61) if d not in (0, 1) and not sympy.ntheory.primetest.is_square(abs(d))]:
            if gcd(Dl, p * q) != 1:
                continue
            pat = (int(sympy.legendre_symbol(Dl % p, p)), int(sympy.legendre_symbol(Dl % q, q)))
            if found.get(pat, 0) < per_pattern:
                found[pat] = found.get(pat, 0) + 1
                out.append(((p, q), Dl, 20, rng.randrange(1 << 30)))
    return out


# ------------------------------------------------------------------- AKS
def check_aks(args):
    """Thm 7: Delta_N(a) = tau_p(Delta_q(a)) mod p, tau_p: X -> X^(p mod r);
    Delta_q(a) mod p depends on q only mod L = lcm(r, ord(X+a)) (checked: q -> q + L);
    Norm(X+a)^N identity via the norm on Z/N[X]/(Phi_r) for prime r."""
    import flint
    p, q, r, a = args
    N = p * q

    def delta(e, mod):
        m = flint.nmod_poly([mod - 1] + [0] * (r - 1) + [1], mod)
        P = flint.nmod_poly([a % mod, 1], mod).pow_mod(e, m)
        c = [int(x) for x in P.coeffs()] + [0] * r
        c = c[:r]
        c[e % r] = (c[e % r] - 1) % mod
        c[0] = (c[0] - a) % mod
        return c

    dN = delta(N, N)
    dN_p = [x % p for x in dN]
    dq = delta(q, p)
    tw = [0] * r
    for i, c in enumerate(dq):
        tw[i * p % r] = (tw[i * p % r] + c) % p
    ok_twist = tw == dN_p
    # order of X + a in (F_p[X]/(X^r-1))^*: divides p^d - 1 with d = ord_r(p) (r prime, r != p)
    d = sympy.n_order(p, r)
    M = p ** d - 1
    m = flint.nmod_poly([p - 1] + [0] * (r - 1) + [1], p)
    b = flint.nmod_poly([a % p, 1], p)
    unit = (pow(-a, r, p) - 1) % p != 0
    if unit and d <= 4:
        o = M
        for pr in sympy.factorint(M):
            while o % pr == 0 and b.pow_mod(o // pr, m) == flint.nmod_poly([1], p):
                o //= pr
        L = sympy.ilcm(o, r)
        ok_period = delta(q, p) == delta(q + L, p)
    else:
        o, L, ok_period = None, None, None
    pass_p = not any(dq)
    return {"p": p, "q": q, "r": r, "a": a, "ok_twist": ok_twist, "ord_Xa": o, "L": L,
            "ok_period": ok_period, "identity_mod_p": pass_p}


def aks_cases(n=300, seed=1):
    rng = random.Random(seed)
    ps = list(sympy.primerange(5, 400))
    rs = [r for r in sympy.primerange(5, 60)]
    out = []
    while len(out) < n:
        p, q = rng.sample(ps, 2)
        r = rng.choice(rs)
        if r in (p, q):
            continue
        if len(out) % 2 == 0 and sympy.n_order(p, r) > 4:
            continue  # half the cases have small ord_r(p) so the period check runs
        out.append((p, q, r, rng.randrange(1, 8)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bound", type=int, default=4000)
    ap.add_argument("--three-bound", type=int, default=4000)
    ap.add_argument("--carm-max", type=int, default=63000)
    ap.add_argument("--lucas-bound", type=int, default=3000)
    ap.add_argument("--frob-bound", type=int, default=400)
    ap.add_argument("--torus-bound", type=int, default=1500)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--out", default="experiments/raw/primality_verify.json")
    args = ap.parse_args()
    t0 = time.time()
    res = {"config": vars(args)}
    with Pool(args.workers) as pool:
        sets = {"semiprime": semiprimes(args.bound),
                "three_prime": three_prime(args.three_bound),
                "carmichael": [tuple(sympy.primefactors(c)) for c in CARMICHAEL_SMALL if c <= args.carm_max]}
        for name, S in sets.items():
            rs = pool.map(check_power_tests, S, chunksize=4)
            res[f"power_{name}"] = {"count": len(rs), "fail": [r for r in rs if not r["ok"] or r["ok_closed"] is False],
                                    "examples": rs[:3] if name != "carmichael" else rs}
            print(name, len(rs), "fail", len(res[f"power_{name}"]["fail"]), f"{time.time()-t0:.1f}s", flush=True)
        rs = [r for r in pool.map(check_lucas, semiprimes(args.lucas_bound), chunksize=4) if "skip" not in r]
        res["lucas"] = {"count": len(rs), "fail": [r for r in rs if not (r["ok"] and r["ok_closed_uniform"])],
                        "examples": rs[:3]}
        print("lucas", len(rs), "fail", len(res["lucas"]["fail"]), f"{time.time()-t0:.1f}s", flush=True)
        rs = [r for r in pool.map(check_frobenius, semiprimes(args.frob_bound)) if "skip" not in r]
        res["frobenius"] = {"count": len(rs), "fail": [r for r in rs if not r["ok"]], "examples": rs[:3],
                            "accidental_total": sum(r["tally"]["accidental"] for r in rs),
                            "structural_total": sum(r["tally"]["structural"] for r in rs),
                            "units_total": sum(r["tally"]["units"] for r in rs)}
        print("frobenius", len(rs), "fail", len(res["frobenius"]["fail"]), f"{time.time()-t0:.1f}s", flush=True)
        rs = pool.map(check_torus, torus_cases(args.torus_bound))
        pats = {}
        for r in rs:
            k = str(r["pattern"])
            pats[k] = pats.get(k, 0) + 1
        res["torus"] = {"count": len(rs), "by_pattern": pats, "fail": [r for r in rs if not r["ok"]]}
        print("torus", len(rs), pats, "fail", len(res["torus"]["fail"]), f"{time.time()-t0:.1f}s", flush=True)
        rs = pool.map(check_aks, aks_cases())
        res["aks"] = {"count": len(rs), "fail_twist": [r for r in rs if not r["ok_twist"]],
                      "fail_period": [r for r in rs if r["ok_period"] is False],
                      "identity_mod_p_true": sum(r["identity_mod_p"] for r in rs),
                      "examples": rs[:5]}
        print("aks", len(rs), "twist fail", len(res["aks"]["fail_twist"]), "period fail",
              len(res["aks"]["fail_period"]), f"{time.time()-t0:.1f}s", flush=True)
    res["wall_s"] = time.time() - t0
    with open(args.out, "w") as f:
        json.dump(res, f, indent=1, default=str)


if __name__ == "__main__":
    main()
