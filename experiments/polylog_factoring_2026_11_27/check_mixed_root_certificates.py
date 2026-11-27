"""Compare mixed-lift vanishing with the quadratic shift ideal."""

from random import Random


def evaluate(coefficients, x, modulus):
    total = 0
    for value in reversed(coefficients):
        total = (total*x+value) % modulus
    return total


def shift_expansion(coefficients, center, scale):
    # Divide by (x-center)^2-scale^2 repeatedly, keeping linear remainders.
    remaining = coefficients[:]
    constant, linear = center*center-scale*scale, -2*center
    while len(remaining) > 2:
        quotient = [0]*(len(remaining)-2)
        for degree in range(len(remaining)-1, 1, -1):
            top = remaining[degree]
            quotient[degree-2] = top
            remaining[degree-1] -= linear*top
            remaining[degree-2] -= constant*top
        yield tuple(remaining[:2])
        remaining = quotient
    yield tuple(remaining)


def check(n, smallest_prime, rng):
    cases = 0
    max_degree = smallest_prime-1
    for m in (1, 2):
        modulus = n**m
        for center, scale in ((0, 1), (2, 2), (4, 4)):
            inverse = pow(scale, -1, n)
            roots = [
                x for x in range(n)
                if (x-center)*inverse % n not in (1, n-1) and
                ((x-center)**2-scale*scale) % n == 0
            ]
            for _ in range(60):
                coefficients = [0]*(max_degree+1)
                basis = [1]
                for i in range(max_degree//2+1):
                    for shift in range(min(2, max_degree-2*i+1)):
                        weight = rng.randint(-3, 3) * n**max(0, m-i)
                        for j, entry in enumerate(basis):
                            coefficients[j+shift] += weight*entry
                    product = [0]*(len(basis)+2)
                    for j, entry in enumerate(basis):
                        for k, f in enumerate((center*center-scale*scale, -2*center, 1)):
                            product[j+k] += entry*f
                    basis = product
                if rng.randrange(3) == 0:
                    coefficients[rng.randrange(max_degree+1)] += rng.randrange(1, modulus)
                member = all(
                    all(v % n**max(0,m-i) == 0 for v in row)
                    for i, row in enumerate(shift_expansion(coefficients, center, scale))
                )
                vanishes = all(
                    evaluate(coefficients, x+lift*n, modulus) == 0
                    for x in roots
                    for lift in range(n**(m-1))
                )
                assert member == vanishes
                cases += 1
    return cases


if __name__ == "__main__":
    rng = Random(261127)
    cases = sum(check(n, p, rng) for n, p in [(15, 3), (35, 5), (45, 3), (77, 7), (105, 3)])
    print(f"{cases} polynomial certificate checks")
