"""Run with the repository Python -B. All classification decisions are certified."""
import argparse
import itertools
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

ROOT = Path("/Volumes/git/factoring-connection")
sys.path[:0] = [str(ROOT / "experiments/masullo_ext/mext_c")]
from cover import ex, pair_units
from flint import acb, arb
from sympy import integer_nthroot, primerange


def orbit(z, units, E):
    reps = []
    for u, v in units:
        w = tuple(t for a, b in zip(z[::2], z[1::2])
                  for t in (u*a+2*v*b, u*b+v*a))
        if max(map(abs, w)) <= E:
            reps.append(w)
    return min(reps)


def row_prefix(row, points, m, E):
    # Precompute half sums as balls. Sort by the certified lower endpoints.
    h = len(row)//2
    def table(vals):
        return {z: sum((a*b for a,b in zip(vals,z)), acb(0))
                for z in itertools.product(range(-E,E+1), repeat=len(vals))}
    left, right = table(row[:h]), table(row[h:])
    def size(z):
        return abs(left[z[:h]] + right[z[h:]])
    ordered = []
    for z in points:
        v = size(z)
        ordered.append((v.lower(), v.upper(), z))
    ordered.sort(key=lambda it: it[0])
    trial = []
    for low, high, z in ordered:
        if ex.kranks(trial+[z]) > len(trial):
            trial.append(z)
            if len(trial) == m:
                lo = low
                hi = max(size(w).upper() for w in trial)
                break
    prefix = []
    for low, high, z in ordered:
        if low >= lo:
            break
        if high < lo and ex.kranks(prefix+[z]) > len(prefix):
            prefix.append(z)
    assert len(prefix) == m-1
    # True strict span is proper and contains this hyperplane, hence equals it.
    return prefix, lo, hi, size


def run(bits, seed, n, m, E):
    ex.n, ex.d = n, 2*n
    N = ex.fixture(bits,seed)
    _, b, P, Q = ex.select(N)
    pols = [X*X-2 for X in (P,Q)]
    gens = [ex.expand(X, f, ex.generators(X,"divided",m))
            for X,f in zip((P,Q),pols)]
    residues = [int(g(b)) % N for g in gens[0]]
    assert residues == [int(g(b)) % N for g in gens[1]]
    points = [z for z in itertools.product(range(-E,E+1),repeat=2*m)
              if next((v for v in z if v),0) > 0 and
              math.gcd(sum(a*v for a,v in zip(z,residues)),N)==1]
    rows,_ = ex.rows_at(ex.coeffs(pols[0]),list(map(ex.coeffs,gens[0])),
                        ex.coeffs(pols[1]),list(map(ex.coeffs,gens[1])),192)
    spaces = []
    w0 = (-2,1,-1,1,-1,0)
    witness = bits==72 and seed==19 and (n,m,E)==(8,3,2)
    ratios, local = [], []
    tl, th = arb(1), arb(1)
    for row in rows:
        basis, lo, hi, size = row_prefix(row,points,m,E)
        if witness:
            v = size(w0)
            ratios.append((v.lower()/hi,v.upper()/lo))
            local.append([str(lo),str(hi),str(v),basis,ex.kranks(basis+[w0])])
            tl,th = tl*lo, th*hi
        if all(ex.kranks(old+basis) > m-1 for old in spaces):
            spaces.append(basis)
    units = pair_units()
    cap = 3**(4*n) * 2**((bits+2*n-1)//(2*n))
    root,exact = integer_nthroot(bits*(bits-1).bit_length()**3,4)
    B = 2**(int(root)+int(not exact))
    primes = list(primerange(B+1))
    primeprod = math.prod(primes)
    fields = [ex.Ideals(f,primes) for f in pols]
    covered, escape, columns = [], [], {}
    count, escaping = 0, 0
    for z in sorted({orbit(z,units,E) for z in points}):
        elts = [ex.element(G,z) for G in gens]
        ns = [abs(int(f.resultant(a))) for f,a in zip(pols,elts)]
        F = math.prod(ns)
        if not 0 < F <= cap:
            continue
        count += 1
        inside = any(ex.kranks(s+[z]) == m-1 for s in spaces)
        escaping += not inside
        if not ex.smooth(F,primeprod):
            continue
        mask = 0
        for side in range(2):
            for prime in fields[side].parity(elts[side],ns[side],primes):
                mask |= 1 << columns.setdefault((side,prime),len(columns))
        (covered if inside else escape).append(mask)
    rec = dict(bits=bits,seed=seed,n=n,m=m,E=E,N=str(N),B=B,cap=str(cap),
               planes=len(spaces),valid=2*len(points),capped=count,esc=escaping,
               smooth=len(covered)+len(escape),smooth_esc=len(escape),
               covered_rank=ex.rank(covered),fresh=ex.rank(covered+escape)-ex.rank(covered))
    if witness:
        witnessdata = dict(local=local,threshold_product=[str(tl),str(th)],
                           min_ratio=str(min(p[0] for p in ratios)),
                           max_ratio=str(max(p[1] for p in ratios)),
                           alternate_ternary=[p for p in points if max(map(abs,p))==1 and
                                              orbit(p,units,E)==orbit(w0,units,E)])
        Path("archived_details.json").write_text(json.dumps(witnessdata,indent=1)+"\n")
    return rec


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("seed",type=int)
    args = p.parse_args()
    for bits, (n,m) in zip((60,72,84,96,112),((7,3),(8,3),(9,4),(10,4),(12,4))):
        for degree,rank,E in dict.fromkeys(((8,3,2),(8,3,1),(n,m,1))):
            print(json.dumps(run(bits,args.seed,degree,rank,E)),flush=True)
