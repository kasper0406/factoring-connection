"""Run with the repository .venv/bin/python -B experiment.py > runs.jsonl."""
import json
import math
import random
import sys
from collections import Counter
from pathlib import Path
from time import perf_counter

from flint import fmpz, fmpz_mat, fmpz_poly
from sympy import primerange, sqrt_mod

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "experiments/subthird/codex"))
import round2_crt as fixtures

V = 32768
K = 64
TILE = 64
LAYERS = {2: [4, 4, 4, 2, 2, 2, 2], 3: [4, 4, 2, 2], 5: [2, 2], 7: [2]}
DIM = 6


def cap(p, maximum):
    power, e = 1, 0
    while power * p <= maximum:
        power *= p
        e += 1
    return e


def primorial(primes):
    return fmpz(math.prod(primes))


def smooth(n, product):
    d = math.gcd(n, int(product % n))
    while d > 1:
        n //= d
        d = math.gcd(n, d)
    return n == 1


def curves(pairs, primes, maximum, collision, choices, messages, projective=False):
    weight = sum(cap(p, maximum) * math.log2(p) for p in primes)
    alphabet_bits = weight + (sum(math.log2(1 + 1/p) for p in primes) if projective else 0)
    johnson = math.sqrt(choices * collision * weight)
    fractions, multiples = [0] * 4, [0] * 3
    best = 0
    for n, part in pairs:
        for i in range(4):
            fractions[i] += part**4 >= n**(i+1)
        score = math.log2(part)
        for i, multiple in enumerate((0.25, 0.5, 1.0)):
            multiples[i] += score >= multiple * johnson
        best = max(best, score)
    return dict(W=weight, R=math.log2(messages)/alphabet_bits, J=johnson,
                max_score=best, fraction_counts=fractions, johnson_counts=multiples)


def interval_sieve(start, primes):
    residuals = list(range(start, start+V))
    hits = divisions = 0
    for p in primes:
        for i in range(-start % p, V, p):
            hits += 1
            value = residuals[i]
            while value % p == 0:
                value //= p
                divisions += 1
            residuals[i] = value
    witnesses = {i for i, rem in enumerate(residuals) if rem == 1}
    return residuals, witnesses, dict(incidences=hits, divisions=divisions, tested=V)


def lattice_setup():
    # For valuation j the i-th Taylor term must supply the missing divisibility.
    scales = [math.prod(p**sum(max(h-i, 0) for h in profile)
                        for p, profile in LAYERS.items()) for i in range(DIM)]
    # Translating the Taylor center by a multiple of every maximal p^e
    # preserves the conditions (the caps of the scoring profile suffice).
    period = math.prod(p**len(profile) for p, profile in LAYERS.items())
    return scales, period


def forced_divisor(n):
    result = 1
    for p, profile in LAYERS.items():
        for increment in profile:
            if n % p:
                break
            result *= p**increment
            n //= p
    return result


def interval_decoder(start, primes):
    product = primorial(primes)
    scales, period = lattice_setup()
    forced = extra = smooth_forced = 0
    witnesses = set()
    records = []
    for tile in range(0, V, TILE):
        center = start + tile + TILE//2
        c = center % period
        rows = [[scales[i]*math.comb(i,j)*c**(i-j)*(TILE//2)**j
                 if j <= i else 0 for j in range(DIM)] for i in range(DIM)]
        reduced = fmpz_mat(rows).lll()
        short = [reduced[0,j] for j in range(DIM)]
        bound = sum(abs(x) for x in short)
        poly = fmpz_poly([x // (TILE//2)**j for j, x in enumerate(short)])
        roots = set()
        for val, _ in poly.roots():
            z = int(val)
            if -TILE//2 <= z < TILE//2:
                roots.add(z)
                guaranteed = forced_divisor(center+z) > bound
                forced += guaranteed
                extra += not guaranteed
                if smooth(center+z, product):
                    witnesses.add(tile+TILE//2+z)
                    smooth_forced += guaranteed
        records.append((center, int(bound), roots))
    return witnesses, records, dict(lll=V//TILE, forced=forced, extra=extra,
                                    smooth_forced=smooth_forced, hits=len(witnesses))


def run_interval(start, primes):
    then = perf_counter()
    rem, truth, counts = interval_sieve(start, primes)
    counts.update(seconds=perf_counter()-then, hits=len(truth))
    then = perf_counter()
    hits, records, decoded = interval_decoder(start, primes)
    decoded['seconds'] = perf_counter()-then
    assert hits <= truth
    product = primorial(primes)
    # Every residual is coprime to the entire base, tested in blocks.
    for i in range(0,V,256):
        assert product.gcd(fmpz(math.prod(rem[i:i+256]))) == 1
    for center, bound, roots in records:
        for z in range(-TILE//2,TILE//2):
            if forced_divisor(center+z) > bound:
                assert z in roots
    curve = curves(((start+i, (start+i)//r) for i,r in enumerate(rem)),
                   primes, start+V-1, math.log2(V), 1, V)
    return dict(X=start, B=max(primes), baseline=counts, decoder=decoded, curve=curve)


def norm(form, u, v):
    a,b,c = form
    return a*u*u+b*u*v+c*v*v


def line_roots(form, p):
    a,b,c = form
    if p == 2:
        return [r for r in range(2) if norm(form,r,1)%p == 0]
    if a % p == 0:
        return [(-c*pow(b,-1,p))%p] if b%p else []
    delta = b*b-4*a*c
    return [int((z-b)*pow(2*a,-1,p)%p) for z in sqrt_mod(delta%p,p,all_roots=True)]


def full_lattice_sieve(form, primes):
    a,b,c = form
    residuals = {(u,v): norm(form,u,v) for v in range(1,K+1)
                 for u in range(1,K+1) if math.gcd(u,v)==1}
    removed_bad = Counter()
    incidences = 0
    for p in primes:
        roots = line_roots(form,p)
        assert any(t%p for t in form)
        bad = (2*a*(b*b-4*a*c))%p == 0

        def strip(u,v):
            nonlocal incidences
            key = (u,v)
            if key not in residuals:
                return
            n = residuals[key]
            assert n%p == 0
            incidences += 1
            original = n
            while n%p == 0:
                n //= p
            residuals[key] = n
            if bad:
                removed_bad[key] = (removed_bad[key] or 1)*(original//n)

        for v in range(1,K+1):
            if v%p == 0:
                if a%p == 0:
                    for u in range(1,K+1):
                        strip(u,v)
                continue
            for r in roots:
                for u in range((r*v-1)%p+1,K+1,p):
                    strip(u,v)
    truth = {key for key,v in residuals.items() if v==1}
    return truth, residuals, removed_bad, dict(incidences=incidences, tested=len(residuals))


def composite_roots(form, primes):
    a,b,c = form
    limit = 4*K*K
    # Only good prime powers have two simple roots and an affine chart.
    pieces = []
    for p in primes:
        if (2*a*(b*b-4*a*c))%p == 0:
            continue
        roots = line_roots(form,p)
        if not roots:
            continue
        assert len(roots)==2
        powers = []
        m = p
        while m<=limit:
            powers.append((m,roots))
            if m*p > limit:
                break
            roots = [r + m*((-(norm(form,r,1)//m)*pow(2*a*r+b,-1,p))%p)
                     for r in roots]
            m *= p
        pieces.append((p,powers))

    def visit(i, modulus, residues):
        if modulus > 2*K*K:
            yield modulus, residues
        for j in range(i,len(pieces)):
            p, powers = pieces[j]
            if p*modulus > limit:
                break
            for m, roots in powers:
                if m*modulus > limit:
                    break
                inv = pow(modulus,-1,m)
                combined = [r+modulus*(((z-r)*inv)%m) for z in roots for r in residues]
                yield from visit(j+1, modulus*m, combined)
    yield from visit(0,1,[0])


def reconstruct(r,d):
    old, new = d, r
    old_v, new_v = 0, 1
    while new > K:
        quotient, rem = divmod(old,new)
        old,new = new,rem
        old_v,new_v = new_v,old_v-quotient*new_v
    if new > 0 and 1 <= new_v <= K and math.gcd(new,new_v)==1:
        return new,new_v
    return None


def forced_search(form, primes, use_reconstruction):
    product = primorial(primes)
    seen, hits = set(), set()
    branches = moduli = before_dedup = 0
    for d, residues in composite_roots(form,primes):
        moduli += 1
        branches += len(residues)
        for r in residues:
            if use_reconstruction:
                point = reconstruct(r,d)
                candidates = [point] if point else []
            else:
                candidates = []
                for v in range(1,K+1):
                    u = r*v%d
                    if 1 <= u <= K and math.gcd(u,v)==1:
                        candidates.append((u,v))
            before_dedup += len(candidates)
            for point in candidates:
                if point not in seen:
                    seen.add(point)
                    if smooth(norm(form,*point),product):
                        hits.add(point)
    return seen,hits,dict(moduli=moduli, branches=branches, matched=before_dedup,
                          tested=len(seen), hits=len(hits))


def run_form(form, primes):
    then = perf_counter()
    truth, residuals, bad, base = full_lattice_sieve(form,primes)
    base.update(seconds=perf_counter()-then, hits=len(truth))
    outputs = {}
    for key, method in [('recovery',True), ('matched',False)]:
        then = perf_counter()
        seen, hits, work = forced_search(form,primes,method)
        work['seconds'] = perf_counter()-then
        outputs[key] = work
        if method:
            expected = seen,hits
        else:
            assert expected == (seen,hits)
        assert hits == truth & seen
    product = primorial(primes)
    assert truth == {z for z in residuals if smooth(norm(form,*z),product)}
    a,b,c = form
    good = [p for p in primes if (2*a*(b*b-4*a*c))%p and len(line_roots(form,p))==2]
    # The two linear components have trivial common kernel modulo every good p.
    for p in good:
        r,s = line_roots(form,p)
        assert math.gcd(r-s,p)==1
    pairs = ((n,n//rem//(bad.get(z,1))) for z,rem in residuals.items()
             for n in [norm(form,*z)])
    curve = curves(pairs,good,(a+abs(b)+c)*K*K,math.log2(K*K),2,len(residuals),True)
    return dict(form=form,B=max(primes),base=base,curve=curve,**outputs)


def main():
    for bits,b1,b2 in [(40,256,256),(64,2048,1024),(80,16384,2048),(96,65536,4096)]:
        for rep in range(6):
            rng = random.Random(617030+1000*bits+rep)
            start = rng.randrange(1<<(bits-1),(1<<bits)-V+1)
            primes = list(primerange(b2+1))
            n = fixtures.fixture(bits,rng)
            form = fixtures.form(n,primes,rng)
            interval = run_interval(start,list(primerange(b1+1)))
            m2 = run_form(form,primes)
            print(json.dumps(dict(bits=bits,rep=rep,N=n,bounds=[b1,b2],
                                  interval=interval,m2=m2)),flush=True)


if __name__=='__main__':
    main()
