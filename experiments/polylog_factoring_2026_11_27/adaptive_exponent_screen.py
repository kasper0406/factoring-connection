"""Screen public scale-dependent power predictors on fresh balanced moduli."""

import json
from math import gcd, isqrt, log, sqrt
from random import Random

from sympy import nextprime

from check_quarter_fourier import quarter


TIERS = (24, 36, 48, 64)
MODULI_PER_TIER = 8
SAMPLES_PER_MODULUS = 8192
SEED = 261130
NAMES = ("scale", "inverse", "quotient", "remainder", "euclid",
         "bit_count", "previous_convergent", "fixed_65537", "planted")


def public_exponents(n, scale):
    quotient, remainder = divmod(n, scale)
    x, y = n, scale
    previous, denominator = 0, 1
    cutoff = isqrt(n)
    euclidean_steps = 0
    current = 1
    while y:
        step, residue = divmod(x, y)
        x, y = y, residue
        previous, denominator = denominator, step*denominator+previous
        if denominator <= cutoff:
            current = denominator
        euclidean_steps += 1
    return (
        scale, pow(scale,-1,n), quotient, remainder,
        (1 << euclidean_steps) + 1, (1 << scale.bit_count()) + 1, current, 65_537,
    )


def fixtures(bits, rng, count=MODULI_PER_TIER):
    high = 1 << (bits//2)
    low = 3*high//4
    chosen = set()
    while len(chosen) < count:
        p,q = sorted(int(nextprime(rng.randrange(low,high))) for _ in range(2))
        n = p*q
        if p != q and q <= 2*p and n.bit_length() == bits and n not in chosen:
            chosen.add(n)
            yield p,q


def trial(p, q, rng):
    n = p*q
    hidden = 1+p*(-2*pow(p,-1,q) % q)
    correct = [0]*len(NAMES)
    for _ in range(SAMPLES_PER_MODULUS):
        c = rng.randrange(1,n)
        while gcd(c,n) != 1:
            c = rng.randrange(1,n)
        actual = quarter(n,c*hidden)
        for i,exponent in enumerate(public_exponents(n,c)):
            correct[i] += (quarter(n,pow(c,exponent,n)) == actual)
        planted = actual if rng.random() < .55 else -actual
        correct[-1] += planted == actual
    return {"N": n, "correct": dict(zip(NAMES,correct))}


def main():
    fixture_rng = Random(SEED)
    sample_rng = Random(SEED+1)
    records = {bits: [trial(p,q,sample_rng) for p,q in fixtures(bits,fixture_rng)]
               for bits in TIERS}
    confidence_cells = len(TIERS)*len(NAMES)
    sample_count = MODULI_PER_TIER*SAMPLES_PER_MODULUS
    radius = sqrt(log(2*confidence_cells/.05)/(2*sample_count))
    output = {
        "fixture_seed": SEED, "sample_seed": SEED+1,
        "samples_per_modulus": SAMPLES_PER_MODULUS,
        "simultaneous_accuracy_radius": radius, "tiers": records,
    }
    print(json.dumps(output,indent=2))


def confirm_scale():
    fixture_rng = Random(261230)
    sample_rng = Random(261231)
    samples = 16_384
    records = []
    for bits in (48, 64):
        for p, q in fixtures(bits, fixture_rng, count=16):
            n = p*q
            hidden = 1+p*(-2*pow(p,-1,q) % q)
            correct = 0
            for _ in range(samples):
                c = sample_rng.randrange(1,n)
                while gcd(c,n) != 1:
                    c = sample_rng.randrange(1,n)
                correct += quarter(n,pow(c,c,n)) == quarter(n,c*hidden)
            records.append({"N": n, "bits": bits, "correct": correct})
    radius = sqrt(log(2/.05)/(2*samples*len(records)))
    print(json.dumps({"fixture_seed": 261230, "sample_seed": 261231,
                      "samples_per_modulus": samples, "pooled_accuracy_radius": radius,
                      "records": records},indent=2))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--confirm-scale", action="store_true")
    if parser.parse_args().confirm_scale:
        confirm_scale()
    else:
        main()
