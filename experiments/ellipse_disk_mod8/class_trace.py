"""Character-trace experiment for a cyclic class group of prime discriminant."""

import json
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt
from time import perf_counter

from class_prefix import branch_from_n, reduced_forms_at
from flint import fmpz_mod_poly_ctx, fmpz_poly, fq_default_ctx
from sympy import factorint, isprime

from genusleak.forms import compose, identity, inverse, power


def band_label(form, cutoff):
    """One representative per sign; below the cap these are the only primitive directions."""
    a, b, d = form
    return (
        sum(
            cutoff < 4 * t and t <= cutoff and isqrt(t) ** 2 != t
            for t in (a, d, a + b + d, a - b + d)
        )
        % 2
    )


def field_bits(value):
    return hex(sum(int(bit) << i for i, bit in enumerate(value.to_list())))


def field_trace_bit(x, degree):
    total = x
    for _ in range(1, degree):
        x = x * x
        total += x
    return int(total)


class UnsupportedClassGroup(ValueError):
    """The public discriminant has a trivial or noncyclic class group."""


class ClassTraces:
    """Fourier data indexed by (character exponent, character order, trace degree)."""

    def __init__(self, delta):
        if delta <= 3 or delta % 8 != 3 or not isprime(delta):
            raise ValueError(
                "requires a prime discriminant magnitude greater than 3 and 3 mod 8"
            )
        self.delta = delta
        self.work = Counter()
        proper = set()
        for a in range(1, isqrt(delta // 3) + 1):
            self.work["leading_coefficients"] += 1
            for f in reduced_forms_at(delta, a, self.work):
                proper.update((f, inverse(f)))
        h = self.h = len(proper)
        if h == 1:
            raise UnsupportedClassGroup("trivial_class_group")
        unit = identity(-delta)
        prime_divisors = factorint(h)
        for g in sorted(proper):
            self.work["generator_candidates"] += 1
            if all(power(g, h // l) != unit for l in prime_divisors):
                break
        else:
            raise UnsupportedClassGroup("noncyclic_class_group")
        self.forms = [unit]
        for _ in range(1, h):
            self.forms.append(compose(self.forms[-1], g))
        if set(self.forms) != proper:
            raise ArithmeticError("incomplete cyclic class enumeration")
        self.work["classes"] = h

        coeffs = fmpz_poly.cyclotomic(h)
        polynomial = fmpz_mod_poly_ctx(2)([int(c) % 2 for c in coeffs])
        _, factors = polynomial.factor()
        irr = min((f for f, _ in factors), key=lambda f: tuple(map(int, f)))
        self.field = fq_default_ctx(modulus=irr)
        self.work["cyclotomic_degree"] = polynomial.degree()
        self.work["field_degree"] = irr.degree()
        zeta = self.field.gen()
        self.powers = [self.field.one()]
        for _ in range(1, h):
            self.powers.append(self.powers[-1] * zeta)
        self.orbits = []
        unused = set(range(1, h))
        while unused:
            k = min(unused)
            j, d = k, 0
            while j in unused:
                unused.remove(j)
                unused.remove((-j) % h)
                j = (2 * j) % h
                d += 1
            self.orbits.append((k, h // gcd(k, h), d))

    def unordered(self, k, e):
        return self.powers[(k * e) % self.h] + self.powers[(-k * e) % self.h]

    def transform(self, values):
        """Fourier representatives from one value per inverse pair, including the identity."""
        coefficients = []
        for k, _, _ in self.orbits:
            value = self.field(values[0] & 1)
            for e in range(1, len(values)):
                if values[e] & 1:
                    value += self.unordered(k, e)
            coefficients.append(value)
        return coefficients

    def calibration(self, cutoff):
        if not 0 <= cutoff <= isqrt(self.delta // 3):
            raise ValueError("cutoff outside the decoder cap")
        labels = [band_label(f, cutoff) for f in self.forms[: (self.h + 1) // 2]]
        return labels, self.transform(labels)

    def norm_profile(self, n):
        """Norm-N ideal multiplicities, via positive representations up to sign; no factors."""
        work = Counter()
        counts = []
        delta = self.delta
        for a, b, _ in self.forms[: (self.h + 1) // 2]:
            target = 4 * a * n
            r = pow(target, (delta + 1) // 4, delta)
            work["public_root_tests"] += 1
            count = 0
            if (r * r - target) % delta == 0:
                high = isqrt(target)
                for residue in {r, -r % delta}:
                    first = -high + (residue + high) % delta
                    for w in range(first, high + 1, delta):
                        v2 = (target - w * w) // delta
                        v = isqrt(v2)
                        work["integer_square_tests"] += 1
                        if v * v == v2 and (v or w > 0) and (w - b * v) % (2 * a) == 0:
                            count += 1
            counts.append(count)
        return counts, work

    def contract(self, constant, weights, observation):
        answer = constant
        for (_, _, degree), w, u in zip(self.orbits, weights, observation, strict=True):
            if w and u:
                answer ^= field_trace_bit(w * u, degree)
        return answer


def evaluate(n, cutoffs):
    br, shortcut = branch_from_n(n)
    if shortcut:
        return {"N": n, "branch_divisor": shortcut}
    return evaluate_branch(br, cutoffs)


def evaluate_branch(br, cutoffs):
    """Evaluate a promised public branch with c² + ell² delta = 4N and c a unit."""
    n = br.n
    if br.c**2 + br.ell**2 * br.delta != 4 * n or gcd(br.c * br.ell, n) != 1:
        raise ValueError("invalid or nonunit public branch")
    start = perf_counter()
    model = ClassTraces(br.delta)
    middle = perf_counter()
    counts, work = model.norm_profile(n)
    observation = model.transform(counts)
    read = perf_counter()
    results = []
    for b in cutoffs:
        labels, weights = model.calibration(b)
        a = model.contract(labels[0], weights, observation)
        results.append(
            {
                "B": b,
                "band_parity": a,
                "positive_inverse_pairs": sum(labels),
                "weights": [field_bits(v) for v in weights],
            }
        )
    return {
        "N": n,
        "delta": br.delta,
        "h": model.h,
        "generator": model.forms[1],
        "orbits": model.orbits,
        "observation": [field_bits(u) for u in observation],
        "setup": dict(model.work),
        "readout": dict(work),
        "setup_seconds": middle - start,
        "readout_seconds": read - middle,
        "cutoffs": results,
        "calibration_and_query_seconds": perf_counter() - read,
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("cutoffs", nargs="+", type=int)
    args = parser.parse_args()
    print(json.dumps(evaluate(args.N, args.cutoffs), indent=2))
