"""Check gated band probabilities exactly and compare bounded-search backends."""

from collections import Counter
from fractions import Fraction
from importlib.util import find_spec
from itertools import combinations, product
from math import gcd, prod
from random import Random

if __package__:
    from .sqrt_prime_marks import (
        _DecimalRing, _confirmed_band, _dyadic_convolution, _least_nonunit, _reveal,
        sqrt_mark_probability, sqrt_prime_marks,
    )
else:
    from sqrt_prime_marks import (
        _DecimalRing, _confirmed_band, _dyadic_convolution, _least_nonunit, _reveal,
        sqrt_mark_probability, sqrt_prime_marks,
    )


def exact_band_law(primes, constant, gate):
    """Replay both possible results of each rational coin, including the outer gate."""
    class NeedCoin(Exception):
        pass

    law = Counter({(): 1-gate})
    pending = [((), gate)]
    while pending:
        path, mass = pending.pop()
        decisions = iter(path)

        def coin(p):
            assert 0 <= p <= 1
            if p == 0 or p == 1:
                return bool(p)
            decision = next(decisions, None)
            if decision is None:
                raise NeedCoin(p)
            return decision

        try:
            marks = tuple(_confirmed_band(primes, constant, gate, coin))
        except NeedCoin as request:
            probability, = request.args
            pending.append((path+(True,), mass*probability))
            pending.append((path+(False,), mass*(1-probability)))
        else:
            law[marks] += mass
    return +law


def run():
    backends = ["py", "decimal", "dyadic"]
    if find_spec("flint") is not None:
        backends.append("flint")
    checked = 0
    for n in (2, 4, 25, 30, 77, 121, 1009, 41*47, 97*101, 1009**2, 5003*5011):
        for upper in (1, 2, 7, 18, 40, 101, min(n, 2050)):
            expected = next((u for u in range(2,upper+1) if gcd(u,n) > 1), None)
            assert all(_least_nonunit(n,upper,b) == expected for b in backends)
            checked += 1
    for backend in backends:
        assert _reveal(3*5**2*7*11, 5, 7, backend) == (11,[5,7])
        assert _reveal(97, 96, 97, backend) == (1,[97])
        assert _reveal(2**11*43**3, 40, 46, backend) == (1,[43])

    rng = Random(5681)
    for n, a_len, b_len in ((91,31,36), ((1<<127)-1,53,27), ((1<<1999)+235,28,25)):
        a = [rng.randrange(n) for _ in range(a_len)]
        b = [rng.randrange(n) for _ in range(b_len)]
        expected = [sum(a[i]*b[j-i] for i in range(max(0,j-b_len+1),min(a_len,j+1))) % n
                    for j in range(a_len+b_len-1)]
        for ring in (_DecimalRing(n), _DecimalRing(n, native_digit_limit=0)):
            assert ring.mul(a,b) == ring.make(expected)
            assert ring.mullow(a,b,25) == ring.make(expected[:25])

    # Direct convolution avoids masking an off-by-one error by the modulus.
    small = [tuple(v) for n in range(4) for v in product((-1,0,1), repeat=n)]
    for a,b in product(small,repeat=2):
        reference = [0] * (len(a)+len(b)-1) if a and b else []
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                reference[i+j] += x*y
        assert _dyadic_convolution(a,b,101) == [x % 101 for x in reference]
    # A long, constant-coefficient triangle exercises twiddle drift at 2048
    # points. Use modulus greater than every exact coefficient.
    height, length = (1 << 77)-1, 521
    big = [height]*length
    triangle = [height**2 * min(i,2*length-i) for i in range(1,2*length)]
    assert _dyadic_convolution(big,big,max(triangle)+1) == triangle

    cases = 0
    for primes,constant,gate in (
        ((), 1, Fraction(1,3)), ((37,), 1, Fraction(1,5)),
        ((37,41), 1, Fraction(1,3)), ((67,71,73), 1, Fraction(3,8)),
        ((2,3), 8, Fraction(1)), ((89,97,101), 4, Fraction(1)),
    ):
        law = exact_band_law(primes,constant,gate)
        for k in range(len(primes)+1):
            for subset in combinations(primes,k):
                target = prod(sqrt_mark_probability(p,constant) if p in subset
                              else 1-sqrt_mark_probability(p,constant) for p in primes)
                assert law[subset] == target
        cases += 1

    # Identical independent draws must give identical marks for each backend.
    n = 2**4*13*37**2*41*101
    for seed,constant in product(range(20), (0,1,8)):
        marks = [list(sqrt_prime_marks(n,Random(seed),constant=constant,backend=b)) for b in backends]
        assert all(value == marks[0] for value in marks)
        assert sorted(set(marks[0])) == marks[0]
        assert all(n % p == 0 for p in marks[0])
    assert list(sqrt_prime_marks(1,Random(8))) == []
    print(f"{checked} least-nonunit bounds using {', '.join(backends)}; three native/crossover comparisons")
    print(f"{len(small)**2} signed dyadic convolutions and one 2048-point triangle")
    print(f"{cases} complete band laws and 60 matched native draws")


if __name__ == "__main__":
    run()
