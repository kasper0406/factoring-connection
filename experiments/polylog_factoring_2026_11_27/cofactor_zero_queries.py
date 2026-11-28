"""Decode prime divisors using only whole-modulus zero tests on dyadic factorial ratios."""


if __package__:
    from .factored_gaussian_packing import certified_prime
else:
    from factored_gaussian_packing import certified_prime


def _odd_candidate(n, signature):
    lower = (n+signature-1)//signature
    return lower if lower % 2 else lower+1


def zero_from_complement(n, mask, *, prime_test=certified_prime, stats=None):
    """Answer a whole-modulus query by enumerating only unselected mask positions.

    A call has at most 2**k candidates for k unselected positions. For the
    bit-time estimate use a certified polynomial-time primality test.
    """
    if n <= 1 or n % 2 == 0:
        raise ValueError("expected an odd input above one")
    full = (1 << (n.bit_length()-1))-1
    allowed = full ^ 1
    if mask < 0 or mask & ~allowed:
        raise ValueError("expected positive dyadic bit positions with top at least two")
    free = allowed ^ mask
    subset = free
    examined = 0
    while True:
        examined += 1
        candidate_mask = subset | 1
        p = _odd_candidate(n, candidate_mask)
        if p >= 3 and n % p == 0 and prime_test(p):
            quotient = n // p
            while quotient and not (quotient & mask):
                quotient //= p
            if not quotient:
                answer = False
                break
        if subset == 0:
            answer = True
            break
        subset = (subset-1) & free
    if stats is not None:
        stats.update(unselected=free.bit_count(), candidates=examined)
    return answer


def largest_prime(n, zero):
    """Return the largest prime of an odd n>1.

    zero(n, mask) tests divisibility by n of (product F_(2^j))^bitlength(n)
    over the set bit positions j of mask. The decoder never sees the product.
    """
    if n <= 1 or n % 2 == 0:
        raise ValueError("expected an odd input above one")
    # Higher positions have top <=1 and give the constant ratio 1.
    last = n.bit_length()-2
    mask = 0
    for j in range(last, 0, -1):
        trial = mask | (1 << j)
        if not zero(n, trial):
            mask = trial
    full = (1 << (last+1))-1
    signature = full ^ mask
    # The signature is < N/s + 2N/s² for its prime s.
    factor = _odd_candidate(n, signature)
    if factor < 3 or n % factor:
        raise ValueError("the answers do not decode to an odd divisor")
    return factor


def factorization(n, zero):
    if n < 1:
        raise ValueError("expected a positive input")
    factors = {}
    while n % 2 == 0:
        factors[2] = factors.get(2, 0)+1
        n //= 2
    while n > 1:
        factor = largest_prime(n, zero)
        power = 0
        while n % factor == 0:
            n //= factor
            power += 1
        factors[factor] = power
    return factors


def semiprime_from_pairs(n, zero):
    """Use only singleton and two-ratio zero tests on a product of distinct odd primes.

    The list of possible queries is fixed by n. Reconstruction may short-circuit
    after the singleton tests or the first cross-edge.
    """
    if n < 15 or n % 2 == 0:
        raise ValueError("expected a product of two distinct odd primes")
    bits = [1 << j for j in range(1, n.bit_length()-1)]
    shared = 1
    other = []
    for bit in bits:
        if zero(n, bit):
            shared |= bit
        else:
            other.append(bit)

    # If the smaller mask is contained in the larger it is N/q=p exactly.
    if 1 < shared < n and n % shared == 0:
        return shared
    # On balanced inputs, shared = p AND q. The least bit of N-shared²
    # belongs to exactly one of p or q, so it is a guaranteed cross-edge pivot.
    difference = n - shared*shared
    first = difference & -difference
    if first in other:
        other.remove(first)
        other.insert(0, first)
    pair_answers = {}
    for pivot in other:
        opposite = shared
        for bit in other:
            if bit == pivot:
                continue
            pair = pivot | bit
            if pair not in pair_answers:
                pair_answers[pair] = zero(n, pair)
            if pair_answers[pair]:
                opposite |= bit
        if opposite != shared:
            # With higher-power contributions this is a signature, not
            # necessarily an unaltered cofactor. Its reciprocal isolates a prime.
            factor = _odd_candidate(n, opposite)
            if 1 < factor < n and n % factor == 0:
                return factor
            break
    raise ValueError("the answers do not decode to a semiprime factor")
