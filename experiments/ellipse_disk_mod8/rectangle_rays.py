"""Count a public rectangle of norm levels and positive rows modulo 2 or 4."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from itertools import product
from math import gcd, isqrt, prod

from bisector import branch_from_n, public_factorization
from row_batches import first_visible


def signed_divisors(primes):
    terms = [(1, 1)]
    for p in primes:
        terms += [(p * d, -sign) for d, sign in terms]
    return terms


def multiplicity(low, high, terms):
    if high < low:
        return 0
    return sum(sign * (high // d - (low - 1) // d) for d, sign in terms)


class Rectangles:
    def __init__(self, br):
        self.br = br
        if br.delta % 4 == 0:
            self.level, self.d = br.n, br.delta // 4
        else:
            self.level, self.d = 4 * br.n, br.delta
        self.work = Counter()
        factors = public_factorization(self.d, self.work)
        self.powers = sorted(p**e for p, e in factors.items())
        self.bad = {p for p, e in factors.items() if e >= 2}
        self.terms = signed_divisors(self.bad)
        self.work["inclusion_terms_built"] = len(self.terms)
        divisors = [1]
        for power in self.powers:
            divisors += [power * d for d in divisors]
        self.work["unitary_divisors_built"] = len(divisors)
        self.seed_moduli = sorted(divisors)
        self.level_rules = {}
        self.level_factor_work = Counter()

    def primitive_rule(self, t):
        if t not in self.level_rules:
            factors = public_factorization(t, self.level_factor_work)
            stride = prod(p for p, e in factors.items() if e >= 2 and p not in self.bad)
            inert = False
            for p in factors:
                if p > 2 and self.level % p and self.d % p:
                    self.work["legendre_tests"] += 1
                    if pow(-self.d, (p - 1) // 2, p) == p - 1:
                        inert = True
                        break
            self.work["base_levels_factored"] += 1
            self.level_rules[t] = (stride, inert)
        return self.level_rules[t]

    def metrics(self):
        return {**self.work, "base_level_factorization": dict(self.level_factor_work)}

    def count(self, low, high, ymin=1, ymax=None, *, modulus=4, rays=True, screen=True):
        """Positive counts. modulus=None gives integers; rays=False tests raw candidates."""
        if modulus not in (None, 2, 4):
            raise ValueError("expected exact, parity or mod 4")
        if not (1 <= low <= self.br.cap + 1 and 0 <= high <= self.br.cap and ymin >= 1):
            raise ValueError("rectangle outside the supported range")
        self.work["queries"] += 1
        if low > high or (ymax is not None and ymax < ymin):
            return 0
        # A scaled point in a narrow level interval can come from below its lower end.
        first = 1 if rays else low
        width = high - first + 1
        seed = next(d for d in self.seed_moduli if d >= width)
        other = [q for q in self.powers if seed % q]
        levels = defaultdict(list)
        inv = pow(self.level, -1, seed)
        for r in range(seed):
            t = first + (inv * r * r - first) % seed
            if t <= high:
                levels[t].append(r)
        self.work["seed_residues"] += seed
        self.work["distinct_seed_levels"] += len(levels)

        # Every rejection here applies to all x residue lifts and all positive rows.
        active = {}
        clipped_rows = ymin != 1 or ymax is not None
        for t, seeds in levels.items():
            hmin = 1 + isqrt((low - 1) // t) if rays else 1
            hmax = isqrt(high // t) if rays else 1
            weight = multiplicity(hmin, hmax, self.terms)
            self.work["inclusion_terms_visited"] += len(self.terms)
            if not weight or (not clipped_rows and modulus and weight % modulus == 0):
                self.work["levels_dropped_by_weight"] += 1
                continue
            active[t] = (seeds, hmin, hmax, weight)
        if not active:
            return 0

        tables = []
        for q in other:
            targets = {(self.level * t) % q for t in active}
            table = defaultdict(list)
            for r in range(q):
                key = r * r % q
                if key in targets:
                    table[key].append(r)
            tables.append(table)
            self.work["local_residues"] += q
            self.work["local_target_requests"] += len(active)
        components = [seed, *other]
        crt = [(self.d // q) * pow(self.d // q, -1, q) for q in components]
        total = 0
        for t, (seeds, hmin, hmax, weight) in active.items():
            target = self.level * t
            roots = [seeds] + [
                table.get(target % q, ())
                for q, table in zip(other, tables, strict=True)
            ]
            if not all(roots):
                continue
            primitive = 1
            if rays:
                primitive, inert = self.primitive_rule(t)
                if screen and inert:
                    self.work["inert_levels_dropped"] += 1
                    continue
            before = self.work["candidate_square_tests"]
            xmax = isqrt(target - self.d) if target >= self.d else 0
            for residues in product(*roots):
                r = sum(a * b for a, b in zip(crt, residues, strict=True)) % self.d
                self.work["terminal_residues"] += 1
                for x in range(r or self.d, xmax + 1, self.d):
                    self.work["positions"] += 1
                    if rays and gcd(x, primitive) != 1:
                        self.work["scaled_positions"] += 1
                        continue
                    v = (target - x * x) // self.d
                    copies = weight
                    if clipped_rows:
                        left = max(hmin, 1 + isqrt((ymin * ymin - 1) // v))
                        right = (
                            hmax if ymax is None else min(hmax, isqrt(ymax * ymax // v))
                        )
                        copies = multiplicity(left, right, self.terms)
                        self.work["inclusion_terms_visited"] += len(self.terms)
                        if not copies or (modulus and copies % modulus == 0):
                            self.work["positions_dropped_by_clip"] += 1
                            continue
                    y = isqrt(v)
                    self.work["candidate_square_tests"] += 1
                    if y * y == v:
                        total += copies
            if 4 * t > high and isqrt(t) ** 2 != t:
                self.work["upper_nonsquare_tests"] += (
                    self.work["candidate_square_tests"] - before
                )
        return total if modulus is None else total % modulus


def recover(br):
    """Use positive mod-4 disk queries, then mod-4 prefixes of the first shell."""
    rectangles = Rectangles(br)

    def visible(cutoff):
        return (rectangles.count(1, cutoff) - isqrt(cutoff)) % 4

    a, _, _ = first_visible(br, visible)
    k = isqrt(a)
    lo, hi = 0, isqrt(4 * br.n * a // br.delta)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        pub = k * k == a and k * br.ell <= mid
        if (rectangles.count(a, a, ymax=mid) - pub) % 4:
            hi = mid
        else:
            lo = mid
    y = hi
    w2 = 4 * br.n * a - br.delta * y * y
    w = isqrt(w2)
    numerator = br.c * y - br.ell * w
    d = gcd(numerator // 2, br.n)
    if w * w != w2 or numerator % 2 or not (1 < d < br.n and br.n % d == 0):
        raise ArithmeticError("first-shell divisor failed verification")
    return a, d, rectangles.metrics()


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("--cutoff", type=int)
    args = parser.parse_args()
    br, direct = branch_from_n(args.N)
    if direct is not None:
        print("branch divisor", direct)
    elif args.cutoff is None:
        a, divisor, work = recover(br)
        print("first mixed level", a, "factor", divisor, "work", dict(work))
    else:
        rectangles = Rectangles(br)
        print("signed residue mod 16", 4 * rectangles.count(1, args.cutoff))
        print("work", rectangles.metrics())
