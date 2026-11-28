"""Native collection controls and an exhaustive two-dependency label check."""

from itertools import product
from random import Random

from shanks_smooth_batch import (
    coefficients_mod, factor_smooth_batch, gamma, kernel_basis, mixed_factor,
    primary_generator, shortest_rational_lift, smooth_row, split_seed,
)


def two_dependencies():
    n, prime, t = 713, 37, -7
    generator = primary_generator(prime)
    _, pair = shortest_rational_lift(n, t)
    row = smooth_row(pair, [(prime, generator)])
    assert row is not None and row[0] % 3 != 0
    deps = kernel_basis([[row[0]] * 3], 3)
    assert len(deps) == 2
    roots = [coefficients_mod(gamma((x, -1)), n)
             for x in range(n) if split_seed(n, x)[1] == t]
    group = [coefficients_mod(gamma(generator), n)]
    outcomes = dict.fromkeys(product((False, True), repeat=2), 0)
    for witnesses in product(roots, repeat=3):
        result = [mixed_factor(n, [row] * 3, witnesses, group, dep) for dep in deps]
        assert all(d is None or (1 < d < n and n % d == 0) for d in result)
        outcomes[tuple(d is not None for d in result)] += 1
    assert outcomes == {(False, False):81, (False, True):162,
                        (True, False):162, (True, True):324}
    print('Two independent dependencies: ', outcomes)


if __name__ == '__main__':
    two_dependencies()
    # Seeds chosen to exercise the matrix branch before an incidental seed gcd.
    for n, bound, seed in [(10002200057, 241, 7), (1000036000099, 499, 8)]:
        factor, stats = factor_smooth_batch(n, bound, rng=Random(seed))
        assert 1 < factor < n and n % factor == 0
        assert stats['source'] == 'mixed dependency' and stats['nullity'] >= 2
        print(f'N={n}, B={bound}, factor={factor}, {stats}')
