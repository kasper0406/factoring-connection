"""HARNESS: list cost, noise tolerance and derived-format demos (oracle -> decoder).

Sub-experiments (each writes experiments/raw/trackE_<name>.jsonl):
  cost      decode wall time / dimension at k = k_min+delta; wrong-candidate cost
  flips     random bit flips in a low-bits hint: plain vs truncation decoder
  highnoise additive noise on the interval hint, decoder widens by E_tol
  crtnoise  p mod l residues with e wrong ones; random-subset decoder
  sumM      format (ii): s mod M for M = prod of primes of a given size; labels
  sum2k     format (i): s mod 2^k; effective modulus 2^(k-v)
  hsbits    random known bits of p and q, branch-and-prune to depth t, then (a)
  pm1       N mod (p-1) (Track B/D D3) -> Pollard exponent
  legendre  Legendre-only format: survivors in the promise interval (toy n)

Run: .venv/bin/python -m src.trackE.exp_formats <name> [<name> ...]
"""
from __future__ import annotations

import math
import os
import random
import sys
import time
from math import comb, prod
from multiprocessing import Pool

from src.common.instances import balanced_semiprime
from src.common.metrics import Metrics, append_jsonl
from src.trackE import oracles
from src.trackE.coppersmith import (AffineHint, balanced_promise, best_m, decode_affine,
                                    decode_interval, decode_low_bits, decode_residue, choose_params)
from src.trackE.formats import (decode_N_mod_pm1, decode_random_bits, decode_residues,
                                decode_residues_noisy, decode_sum_mod_2k, decode_sum_mod_M,
                                hs_branch, legendre, small_primes, sum_mod_2k_cosets)

RAW = "experiments/raw/trackE_{}.jsonl"
SEEDS = 5


def _out(name):
    path = RAW.format(name)
    if os.path.exists(path):
        os.remove(path)
    return path


def _isfactor(g, N):
    return g is not None and 1 < g < N and N % g == 0


# empirical k_min(d) used to size the noise experiments (from trackE_threshold, low fmt);
# recomputed by analyze.py, these are only starting points.
def kmin_guess(n, d):
    L = 2 ** (n / 2 - 1)
    ly = best_m(n, math.log2(L), d, hermite=1.02)[1]
    return math.ceil(n / 2 - 2 - ly)


# --------------------------------------------------------------------------- cost
def _cost_task(args):
    n, seed, delta = args
    inst = balanced_semiprime(n, seed); N = inst.N
    k = n // 4 + delta
    recs = []
    # empirical-minimal d: smallest d in ladder that succeeds
    for d in (2, 3, 4, 5, 7, 9, 11, 15, 21, 29, 37, 49):
        h = oracles.oracle_low_bits(inst, k)
        mt = Metrics(); t0 = time.perf_counter()
        g = decode_low_bits(N, h["k"], h["r"], d=d, metrics=mt)
        dt = time.perf_counter() - t0
        if _isfactor(g, N):
            # wrong candidate at the same (k, d): fresh random residue
            rng = random.Random(f"wrong|{n}|{seed}|{delta}")
            rw = rng.getrandbits(k) | 1
            if rw == h["r"]:
                rw ^= 2
            mw = Metrics(); t1 = time.perf_counter()
            gw = decode_low_bits(N, k, rw, d=d, metrics=mw)
            dw = time.perf_counter() - t1
            recs.append(dict(n=n, seed=seed, k=k, delta=delta, mode="empirical_min_d", d=d,
                             m=mt.fields["m"], decode_s=round(dt, 5), lll_s=round(mt.fields["time_lll"], 5),
                             max_operand_bits=mt.maxima.get("max_operand_bits"),
                             rigorous=mt.fields["rigorous_condition"],
                             wrong_decode_s=round(dw, 5), wrong_found=_isfactor(gw, N)))
            break
    # rigorous parameters (smallest d meeting the exact HG+LLL inequality), if <= 80
    L, U = balanced_promise(N)
    hint = AffineHint(1 << k, inst.labels["p"] % (1 << k), L, U)
    lo, hi = hint.x_range(); xm = (lo + hi) // 2; Y = max(hi - xm, xm - lo)
    pr = choose_params(N, L, Y, dmax=80)
    if pr is not None:
        d, m = pr
        mt = Metrics(); t0 = time.perf_counter()
        g = decode_affine(N, hint, d=d, m=m, metrics=mt)
        recs.append(dict(n=n, seed=seed, k=k, delta=delta, mode="rigorous_d", d=d, m=m,
                         decode_s=round(time.perf_counter() - t0, 5), success=_isfactor(g, N),
                         max_operand_bits=mt.maxima.get("max_operand_bits")))
    else:
        recs.append(dict(n=n, seed=seed, k=k, delta=delta, mode="rigorous_d", d=None))
    return recs


def exp_cost(pool):
    out = _out("cost")
    jobs = [(n, s, dl) for n in (64, 128, 256, 512) for s in range(3) for dl in (2, 4, 8, 16, 32)]
    for recs in pool.imap_unordered(_cost_task, jobs):
        for r in recs:
            append_jsonl(out, r)


# --------------------------------------------------------------------------- flips
def _flips_task(args):
    n, seed, d, extra, e, trial = args
    inst = balanced_semiprime(n, seed); N = inst.N
    rng = random.Random(f"flips|{n}|{seed}|{d}|{extra}|{e}|{trial}")
    k_tr = kmin_guess(n, d) + 1             # decoder's truncation point (public choice)
    k = k_tr + extra
    pos = sorted(rng.sample(range(k), e))
    h = oracles.oracle_low_bits(inst, k, flip_positions=pos)       # noisy oracle
    t0 = time.perf_counter()
    g_plain = decode_low_bits(N, h["k"], h["r"], d=d)                # plain decoder
    t_plain = time.perf_counter() - t0
    t0 = time.perf_counter()
    g_tr = decode_low_bits(N, k_tr, h["r"] % (1 << k_tr), d=d)      # truncation decoder
    t_tr = time.perf_counter() - t0
    return dict(n=n, seed=seed, d=d, k=k, k_trunc=k_tr, extra=extra, e=e, trial=trial,
                flip_positions=pos, lowest_flip=pos[0],
                plain_ok=_isfactor(g_plain, N), trunc_ok=_isfactor(g_tr, N),
                predicted_trunc_ok=pos[0] >= k_tr, t_plain=round(t_plain, 4), t_trunc=round(t_tr, 4))


def exp_flips(pool):
    out = _out("flips")
    jobs = [(n, s, d, extra, e, t) for (n, d) in ((128, 11), (256, 15)) for s in range(SEEDS)
            for extra in (0, 16, 64) for e in (1, 2, 4) for t in range(4)]
    for r in pool.imap_unordered(_flips_task, jobs):
        append_jsonl(out, r)
    # single flip position sweep (n=256, d=15)
    n, d = 256, 15
    for s in range(3):
        inst = balanced_semiprime(n, s); N = inst.N
        k_tr = kmin_guess(n, d) + 1
        k = k_tr + 16
        for j in list(range(0, k, 4)) + [k - 1]:
            h = oracles.oracle_low_bits(inst, k, flip_positions=[j])
            ok_plain = _isfactor(decode_low_bits(N, k, h["r"], d=d), N)
            ok_tr = _isfactor(decode_low_bits(N, k_tr, h["r"] % (1 << k_tr), d=d), N)
            ok_at_j = _isfactor(decode_low_bits(N, j, h["r"] % (1 << j), d=d), N) if j > 0 else False
            append_jsonl(out, dict(sweep=True, n=n, seed=s, d=d, k=k, k_trunc=k_tr, flip=j,
                                   plain_ok=ok_plain, trunc_ok=ok_tr, trunc_at_flip_ok=ok_at_j))


# --------------------------------------------------------------------------- high-bits noise
def exp_highnoise(pool):
    out = _out("highnoise")
    n, d = 256, 15
    for s in range(SEEDS):
        inst = balanced_semiprime(n, s); N = inst.N
        bits_p = inst.labels["p"].bit_length()
        for slack in (0, 2, 4, 6):          # interval made 2^slack narrower than needed max
            h_bits = bits_p - (n // 4 - 6 - slack)
            for noise_exp in (None, 0, 2, 4, 6, 8):
                rng = random.Random(f"hn|{s}|{slack}|{noise_exp}")
                W = (1 << (bits_p - h_bits)) - 1
                noise = 0 if noise_exp is None else rng.choice((-1, 1)) * (W << noise_exp) // 2
                hb = oracles.oracle_high_bits(inst, h_bits, noise=noise)
                for tol_exp in (None, 0, 2, 4, 6, 8):
                    tol = 0 if tol_exp is None else (W << tol_exp) // 2
                    g = decode_interval(N, hb["A"] - tol, hb["W"] + 2 * tol, d=d)
                    append_jsonl(out, dict(n=n, seed=s, d=d, log2W=W.bit_length(), slack=slack,
                                           noise_over_W=None if noise_exp is None else 2 ** noise_exp / 2,
                                           tol_over_W=None if tol_exp is None else 2 ** tol_exp / 2,
                                           log2_width=(hb["W"] + 2 * tol).bit_length(),
                                           ok=_isfactor(g, N)))


# --------------------------------------------------------------------------- noisy CRT residues
def _crt_task(args):
    n, seed, K, e, d = args
    inst = balanced_semiprime(n, seed); N = inst.N
    primes = [l for l in small_primes(10_000) if l > 2 and N % l][:K]
    rng = random.Random(f"crt|{n}|{seed}|{K}|{e}")
    wrong = set(rng.sample(range(K), e))
    o = oracles.oracle_residues(inst, primes, wrong=wrong, rng=rng)
    need = kmin_guess(n, d) + 1
    mt = Metrics(); t0 = time.perf_counter()
    g = decode_residues_noisy(N, o["primes"], o["residues"], need_bits=need, trials=4000,
                              seed=seed, d=d, metrics=mt)
    total_bits = sum(math.log2(l) for l in primes)
    return dict(n=n, seed=seed, K=K, e=e, d=d, need_bits=need, total_bits=round(total_bits, 1),
                ok=_isfactor(g, N), decodes=mt.counts.get("decodes"), wall=round(time.perf_counter() - t0, 3))


def exp_crtnoise(pool):
    out = _out("crtnoise")
    jobs = [(128, s, K, e, 11) for s in range(SEEDS) for K in (40, 80) for e in (0, 2, 4, 8, 16)]
    for r in pool.imap_unordered(_crt_task, jobs):
        append_jsonl(out, r)


# --------------------------------------------------------------------------- format (ii)
def _sumM_task(args):
    n, seed, lbits, use_labels, d = args
    inst = balanced_semiprime(n, seed); N = inst.N
    need = kmin_guess(n, d) + 2
    rng = random.Random(f"sumM|{n}|{seed}|{lbits}|{use_labels}")
    # public auxiliary primes of ~lbits bits; with labels: only primes with (N/l) = -1
    from src.common.instances import is_prime   # public primality of auxiliary moduli only
    primes, Mb, tries = [], 1, 0
    while Mb.bit_length() < need:
        tries += 1
        if tries > 200000:
            return dict(n=n, seed=seed, lbits=lbits, labels=use_labels, status="too_few_primes")
        l = rng.getrandbits(lbits) | (1 << (lbits - 1)) | 1
        if l < 3 or not is_prime(l) or N % l == 0 or l in primes:
            continue
        if use_labels and legendre(N, l) != -1:
            continue
        primes.append(l); Mb *= l
    o = oracles.oracle_sum_mod(inst, Mb)                 # oracle: s mod M
    s_res = {l: o["s"] % l for l in primes}
    labels = None
    if use_labels:
        lg = oracles.oracle_legendre(inst, primes)       # oracle: (p/l)
        labels = dict(zip(lg["primes"], lg["legendre"]))
    mt = Metrics(); t0 = time.perf_counter()
    g = decode_sum_mod_M(N, s_res, labels=labels, d=d, metrics=mt, order_seed=seed,
                         max_decodes=1 << 12)
    return dict(n=n, seed=seed, lbits=lbits, labels=use_labels, d=d, M_bits=Mb.bit_length(),
                omega=mt.fields.get("omega"), omega_ambiguous=mt.fields.get("omega_ambiguous"),
                list_size=mt.fields.get("list_size"), decodes=mt.counts.get("decodes"),
                ok=_isfactor(g, N), wall=round(time.perf_counter() - t0, 3))


def exp_sumM(pool):
    out = _out("sumM")
    jobs = [(128, s, lb, lab, 7) for s in range(SEEDS) for lb in (48, 24, 16, 12, 10, 9, 8, 7)
            for lab in (False, True)]
    for r in pool.imap_unordered(_sumM_task, jobs):
        append_jsonl(out, r)


# --------------------------------------------------------------------------- format (i)
def exp_sum2k(pool):
    out = _out("sum2k")
    for n, d in ((128, 11), (256, 15)):
        kmin = kmin_guess(n, d)
        for s in range(20):
            inst = balanced_semiprime(n, s); N = inst.N
            p, q = inst.labels["p"], N // inst.labels["p"]
            v = ((q - p) & -(q - p)).bit_length() - 1           # harness label only
            for k in (kmin, kmin + 1, kmin + 2, kmin + 3, kmin + 5, kmin + 8):
                o = oracles.oracle_sum_mod(inst, 1 << k)
                mt = Metrics()
                g = decode_sum_mod_2k(N, o["s"], k, d=d, metrics=mt)
                append_jsonl(out, dict(n=n, seed=s, d=d, k=k, v=v, kmin_guess=kmin,
                                       n_cosets=mt.fields.get("n_cosets"),
                                       coset_bits=mt.fields.get("coset_bits"),
                                       decodes=mt.counts.get("decodes"), ok=_isfactor(g, N)))


# --------------------------------------------------------------------------- random bits
def _hs_task(args):
    n, seed, delta, d = args
    inst = balanced_semiprime(n, seed); N = inst.N
    t = kmin_guess(n, d) + 1
    rng = random.Random(f"hs|{n}|{seed}|{delta}")
    o = oracles.oracle_random_bits(inst, delta, t, rng)
    mt = Metrics(); t0 = time.perf_counter()
    try:
        g = decode_random_bits(N, t, o["known_p"], o["known_q"], depth=t, d=d, metrics=mt)
        status = "ok" if _isfactor(g, N) else "fail"
    except RuntimeError:
        g, status = None, "width_exceeded"
    return dict(n=n, seed=seed, delta=delta, d=d, depth=t,
                known_bits=len(o["known_p"]) + len(o["known_q"]),
                max_width=mt.maxima.get("max_width"), final_width=mt.fields.get("width_final"),
                decodes=mt.counts.get("decodes"), status=status, wall=round(time.perf_counter() - t0, 3))


def exp_hsbits(pool):
    out = _out("hsbits")
    jobs = [(n, s, dl, d) for (n, d) in ((128, 11), (256, 15)) for s in range(SEEDS)
            for dl in (0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7)]
    for r in pool.imap_unordered(_hs_task, jobs):
        append_jsonl(out, r)


# --------------------------------------------------------------------------- N mod (p-1)
def exp_pm1(pool):
    out = _out("pm1")
    for n in (64, 128, 256, 512):
        for s in range(20):
            inst = balanced_semiprime(n, s); N = inst.N
            o = oracles.oracle_N_mod_pm1(inst)
            mt = Metrics()
            g = decode_N_mod_pm1(N, o["e"], metrics=mt)
            append_jsonl(out, dict(n=n, seed=s, ok=_isfactor(g, N), gcds=mt.counts.get("gcd"),
                                   base=mt.fields.get("base")))


# --------------------------------------------------------------------------- Legendre-only
def exp_legendre(pool):
    out = _out("legendre")
    n = 36
    for s in range(5):
        inst = balanced_semiprime(n, s); N = inst.N
        L, U = balanced_promise(N)
        primes = [l for l in small_primes(400) if l > 2 and N % l][:40]
        lg = oracles.oracle_legendre(inst, primes)
        # survivors among ALL integers in [L, U] (this is exhaustive: exponential; toy only)
        alive = list(range(L | 1, U + 1, 2))
        for K, (l, sg) in enumerate(zip(primes, lg["legendre"]), start=1):
            alive = [x for x in alive if legendre(x, l) == sg]
            if K in (1, 2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40):
                append_jsonl(out, dict(n=n, seed=s, K=K, survivors=len(alive),
                                       log2_survivors=round(math.log2(max(1, len(alive))), 3),
                                       divisors_among=sum(1 for x in alive if N % x == 0),
                                       log2_interval=round(math.log2(U - L), 3)))
    # coset-list size for the Coppersmith route: prod (l-1)/2 over smallest primes with M >= 2^kneed
    for n, d in ((128, 11), (256, 15), (512, 21)):
        need = kmin_guess(n, d) + 1
        ps, M = [], 1
        for l in small_primes(100_000)[1:]:
            ps.append(l); M *= l
            if M.bit_length() >= need:
                break
        T = prod((l - 1) // 2 for l in ps)
        append_jsonl(out, dict(coset_route=True, n=n, d=d, need_bits=need, n_primes=len(ps),
                               B=ps[-1], log2_list=round(math.log2(T), 2), info_bits=len(ps)))


EXPS = dict(cost=exp_cost, flips=exp_flips, highnoise=exp_highnoise, crtnoise=exp_crtnoise,
            sumM=exp_sumM, sum2k=exp_sum2k, hsbits=exp_hsbits, pm1=exp_pm1, legendre=exp_legendre)

if __name__ == "__main__":
    with Pool(int(os.environ.get("TRACKE_WORKERS", "8"))) as pool:
        for name in sys.argv[1:]:
            t0 = time.time()
            EXPS[name](pool)
            print(name, "done", round(time.time() - t0, 1), flush=True)
