"""Factor from one predicted hitting exponent per selected auxiliary prime."""

from flint import fmpz, fmpz_mat
from sympy import Poly, nextprime, symbols

if __package__:
    from .quartic_collision_powers import divisors_from_roots
else:
    from quartic_collision_powers import divisors_from_roots


def split_from_hits(n, predict, inverse_agreement=1):
    """Return a proper divisor, or n for a prime, under the weighted agreement promise."""
    if n < 3 or n % 2 == 0 or inverse_agreement < 1 or not isinstance(inverse_agreement, int):
        raise ValueError("Expected odd n >= 3 and positive integral inverse agreement")
    bits = n.bit_length()
    maximum_factors, minimum_product = 0, 1
    while minimum_product * (2 * maximum_factors + 3) <= n:
        minimum_product *= 2 * maximum_factors + 3
        maximum_factors += 1
    degree = 2 * maximum_factors * inverse_agreement
    target_bits = 2 * inverse_agreement * degree * (bits + 5 * maximum_factors + 2)
    prime = max(bits, 5)
    residue, period = 0, 1
    while period.bit_length() <= target_bits:
        prime = int(nextprime(prime))
        if n % prime in (0, 1):
            continue
        answer = predict(n, prime) % prime
        correction = (answer - residue) * pow(period, -1, prime) % prime
        residue += period * correction
        period *= prime

    # This row basis spans all polynomials G with G(residue) = 0 modulo period.
    rows = [[period] + [0] * degree]
    for power in range(1, degree + 1):
        row = [0] * (degree + 1)
        row[0] = -pow(residue, power, period)
        row[power] = 1
        rows.append(row)
    vector = [int(value) for value in fmpz_mat(rows).lll().tolist()[0]]
    norm_squared = sum(value * value for value in vector)
    dimension = degree + 1
    if norm_squared**dimension > (1 << (degree * dimension)) * period**2:
        raise ArithmeticError("Short-vector certificate exceeds determinant bound")
    polynomial = Poly.from_list(vector[::-1], symbols("j"))
    for divisor in divisors_from_roots(n, polynomial.ground_roots()):
        if divisor < n:
            return divisor
    if fmpz(n).is_prime():
        return n
    raise ValueError("Predictions did not yield a verified split")
