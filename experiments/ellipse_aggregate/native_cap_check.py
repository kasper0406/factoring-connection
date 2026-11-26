"""Check public-disk geometry, row-telescoping witnesses and shifted ideal powers."""

from math import comb, isqrt
from random import Random

from sympy import Poly, symbols


def check_polynomials():
    w, y = symbols("w y")
    random = Random(24817)
    checked = 0
    for n in (15, 77, 143):
        modulus = 4 * n
        for m in range(1, 5):
            for _ in range(4):
                h, k, delta = (random.randint(1, 12) for _ in range(3))
                form = Poly((w + h) ** 2 + delta * (y + k) ** 2, w, y)
                auxiliary = Poly(0, w, y)
                for power in range(m + 2):
                    shift = Poly(
                        sum(
                            random.randint(-2, 2) * w**i * y**j
                            for i in range(3)
                            for j in range(3)
                        ),
                        w,
                        y,
                    )
                    auxiliary += modulus ** max(0, m - power) * shift * form**power
                # Subtract the two leading w-coefficients via their monic F-adic
                # basis terms. This tests cancellations/degree drops in ordinary shifts.
                for _ in range(2):
                    degree = auxiliary.degree(w)
                    exponent, remainder = divmod(degree, 2)
                    content = modulus ** max(0, m - exponent)
                    leading_terms = {
                        (0, j): c for (i, j), c in auxiliary.terms() if i == degree
                    }
                    assert all(c % content == 0 for c in leading_terms.values())
                    leading = Poly.from_dict(leading_terms, (w, y))
                    auxiliary -= leading * Poly(w**remainder, w, y) * form**exponent
                    checked += 1
    return checked


def check_row():
    n, p = 38_942_974_393, 153_247
    c = isqrt(4 * n)
    delta = 4 * n - c * c
    cap = isqrt(delta // 3)
    w, y = 2_038_973, 8_499

    def mixed_increment(x):
        assert 1 <= x and x * x + delta * y * y <= 4 * n * cap
        assert (x - c * y) % n and (x + c * y) % n
        return x if (x * x + delta * y * y) % (4 * n) == 0 else 0

    w2 = w - 2 * p
    assert n % p == 0 and w % p == w2 % p
    hit, miss = mixed_increment(w), mixed_increment(w2)
    assert hit % p == 46_762 and miss == 0
    # Exact cancellation in a complete period holds, but does not give the partial row.
    for v in range(1, 16):
        assert sum(x for x in range(77) if (x * x - v * v) % 77 == 0) % 77 == 0
    return {
        "N": n,
        "cap": cap,
        "y": y,
        "row_end": isqrt(4 * n * cap - delta * y * y),
        "hit": hit,
        "miss_coordinate": w2,
        "residue_p": hit % p,
    }


def check_differences():
    # For a degree-i univariate polynomial, its i-th difference kills all
    # lower terms. This certificate also tests integer-grid rather than complex norms.
    for i in range(1, 13):
        for h in range(1, 7):
            vals = [(j * h + 5) ** i + (j * h - 3) ** (i - 1) for j in range(i + 1)]
            diff = sum((-1) ** (i - j) * comb(i, j) * v for j, v in enumerate(vals))
            factorial = 1
            for j in range(1, i + 1):
                factorial *= j
            assert diff == factorial * h**i


if __name__ == "__main__":
    check_differences()
    print({"leading_coefficients": check_polynomials(), "row": check_row()})
