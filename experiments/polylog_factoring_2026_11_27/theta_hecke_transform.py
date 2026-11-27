"""Check the root-lattice Hecke decomposition behind the fixed-cap disk."""

from itertools import combinations
from math import gcd, isqrt


def primes(limit):
    return [
        p
        for p in range(3, limit + 1, 2)
        if all(p % d for d in range(3, isqrt(p) + 1, 2))
    ]


def branch(n):
    c = isqrt(4 * n)
    delta = 4 * n - c * c
    if delta <= 4:
        c -= 1
        delta = 4 * n - c * c
    b = delta % 2
    return c, delta, b, (delta + b) // 4


def root_form(n, delta, root):
    """Return F with Q0(N*u+s*v,v)=N*F(u,v)."""
    b = delta % 2
    d = (delta + b) // 4
    s = (root - b) * pow(2, -1, n) % n
    value, rem = divmod(s * s + b * s + d, n)
    assert rem == 0
    return n, 2 * s + b, value, s


def representations(form, value):
    a, b, d = form
    delta = 4 * a * d - b * b
    ymax = isqrt(4 * a * value // delta)
    answer = []
    for y in range(-ymax, ymax + 1):
        radius = 4 * a * value - delta * y * y
        width = isqrt(radius)
        lo = (-b * y - width) // (2 * a) - 1
        hi = (-b * y + width) // (2 * a) + 1
        for x in range(lo, hi + 1):
            if a * x * x + b * x * y + d * y * y == value:
                answer.append((x, y))
    return answer


def roots_from_factors(n, p, q, c):
    roots = set()
    for rp in (c % p, -c % p):
        for rq in (c % q, -c % q):
            roots.add(rp + p * ((rq - rp) * pow(p, -1, q) % q))
    assert len(roots) == 4
    return sorted(roots)


def check_case(p, q):
    n = p * q
    c, delta, b, d = branch(n)
    if gcd(c, n) != 1:
        return 0, 0
    principal = (1, b, d)
    roots = roots_from_factors(n, p, q, c)
    forms = [root_form(n, delta, root) for root in roots]
    cap = isqrt(delta // 3)
    levels = points = 0
    for j in range(1, cap + 1):
        source = representations(principal, n * j)
        partition = [[] for _ in roots]
        for x, y in source:
            assert gcd(y, n) == 1
            root = (2 * x + b * y) * pow(y, -1, n) % n
            partition[roots.index(root)].append((x, y))
        assert sum(map(len, partition)) == len(source)
        lifted_count = 0
        lifted_abs = 0
        counts = []
        for index, (a, middle, last, s) in enumerate(forms):
            target = representations((a, middle, last), j)
            counts.append(len(target))
            image = [(n * u + s * v, v) for u, v in target]
            assert sorted(image) == sorted(partition[index])
            lifted_count += len(target)
            lifted_abs += sum(abs(2 * n * u + middle * v) for u, v in target)
        assert lifted_count == len(source)
        assert lifted_abs == sum(abs(2 * x + b * y) for x, y in source)
        known = [counts[i] for i, root in enumerate(roots) if root in (c, n - c)]
        hidden = [counts[i] for i, root in enumerate(roots) if root not in (c, n - c)]
        assert known == [len(representations(principal, j))] * 2
        assert len(set(hidden)) == 1
        assert len(source) == 2 * known[0] + 2 * hidden[0]
        levels += 1
        points += len(source)
    # The four complete lines are Fourier self-dual in their slope variable.
    dual_values = range(n) if n <= 200 else range(0, n, max(1, n // 31))
    for a in dual_values:
        if gcd(a, n) != 1:
            continue
        for dual_b in dual_values:
            count = sum((a * root + dual_b) % n == 0 for root in roots)
            assert count == int((-dual_b * pow(a, -1, n) % n) in roots)
    return levels, points


def wrong_prime_form_counterexample():
    p, q, j = 127, 227, 7
    n = p * q
    c, delta, b, d = branch(n)
    actual = len(representations((1, b, d), n * j))
    wrong = 2 * len(representations((1, b, d), j)) + 2 * len(
        representations((p, c, q), j)
    )
    mixed = next(root for root in roots_from_factors(n, p, q, c) if root not in (c, n - c))
    hidden = root_form(n, delta, mixed)[:3]
    correct = 2 * len(representations((1, b, d), j)) + 2 * len(
        representations(hidden, j)
    )
    assert (actual, wrong, correct) == (0, 4, 0)
    return {"N": n, "j": j, "actual": actual, "prime_form_prediction": wrong}


def main():
    cases = levels = points = 0
    for p, q in combinations(primes(79), 2):
        if q <= 2 * p:
            got_levels, got_points = check_case(p, q)
            cases += got_levels > 0
            levels += got_levels
            points += got_points
    print(
        {
            "balanced_branches": cases,
            "levels": levels,
            "source_representations": points,
            "counterexample": wrong_prime_form_counterexample(),
        }
    )


if __name__ == "__main__":
    main()
