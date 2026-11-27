"""Recover a modular root with a variable zero-prefix existence oracle."""

from math import gcd

from bounded_quadratic_roots import interval_limit, quadratic_roots


def recover_root(n, y, exists):
    """Assume n>1 coprime to 6 and y a unit square. Return a root mod n."""
    bits = n.bit_length()
    bound = 1 << bits
    period = 3
    residue = 1
    limit = interval_limit(n)
    while (bound - 1) // period > limit:
        next_period = 3 * period
        modulus = n * next_period
        zeros = modulus.bit_length() - bits
        inverse = pow(n, -1, next_period)
        for digit in range(2):
            candidate = residue + digit * period
            square = y + n * (((candidate*candidate-y)*inverse) % next_period)
            if exists(modulus, square, zeros):
                residue = candidate
                break
        else:
            # The three signed children partition the occupied parent.
            residue += 2 * period
        period = next_period
    return finish_progression(n, y, bound, period, residue)


def recover_root_binary(n, y, exists):
    """Recover a square root for odd n > 1 using even oracle moduli."""
    bound = 1 << n.bit_length()
    period, residue = 4, 1
    limit = interval_limit(n)
    while (bound - 1) // period > limit:
        # Squaring modulo 4D selects one signed class modulo 2D.
        scale = 4 * period
        square = y + n * (((residue**2 - y) * pow(n, -1, scale)) % scale)
        if not exists(n * scale, square, scale.bit_length() - 1):
            residue += period
        period *= 2
    return finish_progression(n, y, bound, period, residue)


def finish_progression(n, y, bound, period, residue):
    inverse = pow(period, -1, n)
    for origin in (residue, period - residue):
        if origin >= bound:
            continue
        maximum = (bound - 1 - origin) // period
        offsets = quadratic_roots(n, 2 * origin * inverse,
                                  (origin**2 - y) * inverse**2, maximum)
        if offsets:
            return (origin + period * offsets[0]) % n
    raise ValueError("Selected progression contains no root")


if __name__ == "__main__":
    calls = 0

    def direct_oracle(modulus, square, zeros):
        global calls
        calls += 1
        bound = 1 << (modulus.bit_length() - zeros)
        return any(x*x % modulus == square for x in range(bound))

    cases = 0
    for n in (35, 55, 77, 175, 385, 1001):
        for y in {h*h % n for h in range(n) if gcd(h, n) == 1}:
            root = recover_root(n, y, direct_oracle)
            assert root*root % n == y
            cases += 1
    print(f"Verified {cases} unit squares with {calls} prefix queries.")

    cases = calls = 0
    for n in range(3, 256, 2):
        for y in {h*h % n for h in range(n)}:
            root = recover_root_binary(n, y, direct_oracle)
            assert root*root % n == y
            cases += 1
    print(f"Binary decoder verified {cases} squares with {calls} prefix queries.")
