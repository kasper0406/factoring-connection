"""prepread experiments E1-E5 (pre-registered in notes/prepread_prereg.md).

Usage: .venv/bin/python -m src.prepread.run_experiments {e1,e2,e3,e4,e5}
Raw output: experiments/raw/prepread_<e>.jsonl
"""
from __future__ import annotations
import json, math, os, random, resource, sys, time
from math import lcm
import numpy as np
from sympy import primerange, factorint

from . import slp as S
from . import spectra as SP
from . import compsearch as CSR
from . import dlred

RAW = os.path.join(os.path.dirname(__file__), '..', '..', 'experiments', 'raw')


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


def out(name):
    os.makedirs(RAW, exist_ok=True)
    return open(os.path.join(RAW, f'prepread_{name}.jsonl'), 'a')


def seeded_primes(tag, lo, hi, k):
    rng = random.Random(f'prepread|{tag}')
    return sorted(rng.sample(list(primerange(lo + 1, hi)), k))


# ---------------------------------------------------------------- E1 members
def e1_members():
    mem = []
    for size in (10, 14, 18):
        for i in range(20):
            rng = random.Random(f'prepread|F1|{size}|{i}')
            mem.append((f'F1_s{size}_{i}', 'F1', S.random_slp(size, rng)))
    for c in (1, -1, 3):
        for k in (8, 12, 16, 20):
            mem.append((f'F2_c{c}_k{k}', 'F2', S.composition_slp([c] * k)))
    for k in (10, 16, 20):
        mem.append((f'F3_x2^{k}+x', 'F3', S.power_slp(2 ** k, plus_x=True)))
    for k in (10, 16):
        mem.append((f'F3_dickson2^{k}', 'F3', S.composition_slp([-2] * k)))
    mem.append(('F3_x^lcm30', 'F3', S.power_slp(lcm(*range(1, 31)))))
    return mem


def e1():
    t0 = time.time()
    small = seeded_primes('E1small', 2 ** 13, 2 ** 14, 4)
    large = seeded_primes('E1large', 2 ** 15, 2 ** 16, 2)
    f = out('e1')
    members = e1_members()
    for p in small + large:
        g = SP.primitive_root(p)
        orb = SP.orbit(p, g)
        mask = SP.large_order_mask(p)
        maskall = np.arange(p - 1) != 0
        if p < 2 ** 14:
            cs = range(1, p)
        else:
            cs = random.Random(f'prepread|E1cs|{p}').sample(range(1, p), 256)
        cs = list(cs)
        # baseline: 10 random maps
        base = []
        for b in range(10):
            v = np.random.default_rng([p, b]).integers(0, p, size=p - 1)
            base.append(SP.spectrum_max(v, p, cs, mask)[0])
        Bp = float(np.median(base))
        f.write(json.dumps({'kind': 'baseline', 'p': p, 'g': g, 'B_p': Bp, 'base': base,
                            'ncs': len(cs), 'ell': max(factorint(p - 1))}) + '\n'); f.flush()
        rows = []
        for name, fam, prog in members:
            v = S.eval_mod(prog, orb, p)
            img = int(len(np.unique(v)))
            M, cm, jm = SP.spectrum_max(v, p, cs, mask)
            row = {'kind': 'member', 'name': name, 'family': fam, 'p': p,
                   'deg': S.formal_degree(prog), 'size': len(prog[1]), 'image': img,
                   'nondeg': img >= p / 4, 'M': M, 'ratio': M / Bp, 'c': cm, 'j': jm}
            if p < 2 ** 14 or S.formal_degree(prog) < p - 1:
                w = S.weil_width(S.reduced_poly(prog, p), p)
                row['width'] = w
                if 0 < w < math.sqrt(p):
                    Mall = SP.spectrum_max(v, p, cs, maskall)[0]
                    row['weil_bound'] = w * math.sqrt(p) / (p - 1)
                    row['M_allnonzero'] = Mall
                    row['weil_ok'] = Mall <= row['weil_bound'] * (1 + 1e-9)
            f.write(json.dumps(row) + '\n'); f.flush()
        # planted DL control (HARNESS)
        logs = np.empty(p, dtype=np.int64); logs[orb] = np.arange(p - 1)
        v = (p * np.arange(p - 1)) // (p - 1)
        M, cm, jm = SP.spectrum_max(v, p, cs, mask)
        f.write(json.dumps({'kind': 'member', 'name': 'F4_planted_DL', 'family': 'F4', 'p': p,
                            'image': int(len(np.unique(v))), 'nondeg': True, 'M': M,
                            'ratio': M / Bp, 'c': cm, 'j': jm}) + '\n'); f.flush()
        print(p, 'done', round(time.time() - t0, 1), 's, rss', round(rss_mb()), 'MB', flush=True)
    f.write(json.dumps({'kind': 'meta', 'seconds': time.time() - t0, 'rss_mb': rss_mb()}) + '\n')


# ---------------------------------------------------------------- E2
def v2(n):
    return (n & -n).bit_length() - 1


def e2():
    t0 = time.time()
    f = out('e2')
    for b in (7, 8, 9, 10):
        for p in seeded_primes(f'E2single|{b}', 2 ** b, 2 ** (b + 1), 4):
            k = b + 2
            r = CSR.search(k, [p])
            r.update({'mode': 'single', 'b': b, 'v2m': v2(p - 1), 'v2p': v2(p + 1),
                      'H4a_bound': max(2 * k * k, 4 * 2 ** max(v2(p - 1), v2(p + 1)))})
            f.write(json.dumps(r) + '\n'); f.flush()
            print('single', b, p, r['best'], round(time.time() - t0, 1), flush=True)
    for b in (9, 10):
        ps = seeded_primes(f'E2joint|{b}', 2 ** b, 2 ** (b + 1), 8)
        for k in (b + 2, b + 4):
            r = CSR.search(k, ps)
            r.update({'mode': 'joint', 'b': b})
            f.write(json.dumps(r) + '\n'); f.flush()
            print('joint', b, k, r['best'], round(time.time() - t0, 1), flush=True)
    f.write(json.dumps({'kind': 'meta', 'seconds': time.time() - t0}) + '\n')


# ---------------------------------------------------------------- E3
def tonelli(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, bb * bb % p, t * bb * bb % p, r * bb % p
    return r


def preimages(T, c, p):
    out = []
    for u in T:
        w = (u - c) % p
        if w == 0:
            out.append(0)
        else:
            r = tonelli(w, p)
            if r is not None:
                out.extend((r, p - r))
    return out


def haar_wreath(K):
    pk = [1.0]
    for _ in range(K):
        pk.append((1 - (1 - pk[-1]) ** 2) / 2)
    return pk


def e3():
    t0 = time.time()
    K = 24
    primes = list(primerange(2 ** 14 + 1, 2 ** 16))
    Z = np.zeros((len(primes), K + 1), dtype=np.int64)
    for a, p in enumerate(primes):
        T = [0]
        for k in range(1, K + 1):
            T = preimages(T, 1, p)
            Z[a, k] = len(T)
            if not T:
                break
    pk = haar_wreath(K)
    f = out('e3')
    for k in range(1, K + 1):
        frac = float(np.mean(Z[:, k] > 0))
        tol = 3 * math.sqrt(pk[k] * (1 - pk[k]) / len(primes)) + 0.01
        f.write(json.dumps({'k': k, 'nprimes': len(primes), 'frac_root': frac, 'haar': pk[k],
                            'tol': tol, 'H6_frac_ok': abs(frac - pk[k]) <= tol,
                            'mean': float(Z[:, k].mean()), 'max': int(Z[:, k].max()),
                            'var': float(Z[:, k].var())}) + '\n')
    f.write(json.dumps({'kind': 'meta', 'seconds': time.time() - t0}) + '\n')
    f.close()
    print('e3 done', round(time.time() - t0, 1))


# ---------------------------------------------------------------- E4
def e4():
    t0 = time.time()
    rows = [json.loads(l) for l in open(os.path.join(RAW, 'prepread_e1.jsonl'))]
    base = {r['p']: r for r in rows if r.get('kind') == 'baseline'}
    members = {name: prog for name, fam, prog in e1_members()}
    f = out('e4')
    for p, br in sorted(base.items()):
        g = br['g']
        orb = SP.orbit(p, g)
        rng = np.random.default_rng([p, 4])
        targets = [int(x) for x in rng.integers(0, p - 1, size=20)]
        # (a) planted control: f(y) = e(floor(p*log y/(p-1))/p). HARNESS-built oracle.
        logs = np.empty(p, dtype=np.int64); logs[orb] = np.arange(p - 1)
        vplant = (p * logs) // (p - 1)
        fa = lambda ys, vp=vplant, pp=p: np.exp(2j * np.pi * vp[ys] / pp)
        Fj = SP.A_value((p * np.arange(p - 1)) // (p - 1), p, 1, 1)
        ok, qs = 0, []
        for x in targets:
            xr, q = dlred.recover_dlog(p, g, int(orb[x]), fa, 1, Fj, rng, gpow=orb)
            ok += (xr == x); qs.append(q)
        f.write(json.dumps({'part': 'a', 'p': p, 'ok': ok, 'n': 20, 'Fj': abs(Fj),
                            'queries_per_recovery': int(np.mean(qs))}) + '\n'); f.flush()
        # (b) best non-control member at this prime
        cand = [r for r in rows if r.get('kind') == 'member' and r['p'] == p and r['family'] != 'F4'
                and r['nondeg']]
        best = max(cand, key=lambda r: r['ratio'])
        prog = members[best['name']]
        v = S.eval_mod(prog, orb, p)
        c, j = best['c'], best['j']
        Fj = SP.A_value(v, p, c, j)
        fb = lambda ys, pr=prog, pp=p, cc=c: np.exp(2j * np.pi * ((cc * S.eval_mod(pr, ys, pp)) % pp) / pp)
        ok, qs = 0, []
        for x in targets:
            xr, q = dlred.recover_dlog(p, g, int(orb[x]), fb, j, Fj, rng, gpow=orb)
            ok += (xr == x); qs.append(q)
        f.write(json.dumps({'part': 'b', 'p': p, 'member': best['name'], 'ratio': best['ratio'],
                            'c': c, 'j': j, 'Fj': abs(Fj), 'ok': ok, 'n': 20,
                            'queries_per_recovery': int(np.mean(qs))}) + '\n'); f.flush()
        print('e4', p, round(time.time() - t0, 1), flush=True)
    f.write(json.dumps({'kind': 'meta', 'seconds': time.time() - t0}) + '\n')


# ---------------------------------------------------------------- E5
def e5():
    t0 = time.time()
    primes = seeded_primes('E5', 2 ** 13, 2 ** 14, 12)
    xs = {p: np.arange(p, dtype=np.int64) for p in primes}
    f = out('e5')
    nhit = 0
    for i in range(2000):
        prog = S.random_slp(20, random.Random(f'prepread|E5|{i}'))
        dens, zero = [], []
        for p in primes:
            v = S.eval_mod(prog, xs[p], p)
            nz = int(np.count_nonzero(v == 0))
            dens.append(nz / p); zero.append(nz == p)
        row = {'i': i, 'deg': S.formal_degree(prog), 'min': min(dens), 'median': float(np.median(dens)),
               'max': max(dens), 'n_identically_zero': sum(zero)}
        if min(dens) >= 1 / 64:
            nhit += 1
            row['hit'] = True
            row['dens'] = dens
            row['prog'] = prog
        f.write(json.dumps(row) + '\n')
    f.write(json.dumps({'kind': 'meta', 'seconds': time.time() - t0, 'hits': nhit,
                        'primes': primes}) + '\n')
    print('e5 hits', nhit, round(time.time() - t0, 1))


if __name__ == '__main__':
    {'e1': e1, 'e2': e2, 'e3': e3, 'e4': e4, 'e5': e5}[sys.argv[1]]()
