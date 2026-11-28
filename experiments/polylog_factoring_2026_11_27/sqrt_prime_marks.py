"""Independent exact marks of probability min(1,C/ceil(sqrt(p))) for primes p|N.

Opening a dyadic band costs soft-O(sqrt(P)) times a polynomial in log N,
using native products at bounded sizes and a certified dyadic FFT beyond
the native crossover. The probability of
opening it is O(C log(N)/sqrt(P)). The 'py' Karatsuba backend is an exact
cross-check and does not carry the same expected polynomial bit bound.
"""

from fractions import Fraction
from importlib.util import find_spec
from math import gcd, isqrt
from pathlib import Path
from random import SystemRandom
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.trackC.baselines import ps_interval
from src.trackC.polyback import PolyRing, multipoint_eval, product_tree_linear


# A fixed upper bound on the native fast path avoids extrapolating a library's
# transform-size thresholds to unbounded inputs.
_NATIVE_DIGIT_LIMIT = 1 << 26


def sqrt_mark_probability(prime, constant=8):
    """The independent per-prime mark probability (prime is a supplied prime base)."""
    if prime < 2 or not isinstance(constant, int) or constant < 0:
        raise ValueError("expected a prime base and nonnegative integer constant")
    return min(Fraction(1), Fraction(constant, isqrt(prime-1) + 1))


def _complex_product(a, b, bits):
    x, y = a
    u, v = b
    return ((x*u-y*v) >> bits, (x*v+y*u) >> bits)


def _twiddles(length, bits):
    """Positive half-circle of Lth roots, scaled by S=2^bits, L a power of two."""
    scale = 1 << bits
    cosine, sine = 0, scale  # Fourth root i.
    for _ in range(length.bit_length()-3):
        # Halve the angle. Use sin(theta)/(2*cos(theta/2)) to avoid
        # subtracting nearly equal cosines to find a very small sine.
        cosine = isqrt(scale*(scale+cosine)//2)
        sine = sine*scale//(2*cosine)
    roots = [(scale, 0)]
    for _ in range(1, length//2):
        roots.append(_complex_product(roots[-1], (cosine, sine), bits))
    return roots


def _fft(data, roots, bits, inverse=False):
    """In-place, unnormalized radix-two FFT in dyadic fixed-point arithmetic."""
    length = len(data)
    j = 0
    for i in range(1, length):
        bit = length//2
        while j & bit:
            j ^= bit
            bit //= 2
        j ^= bit
        if i < j:
            data[i], data[j] = data[j], data[i]
    width = 2
    while width <= length:
        half = width//2
        for offset in range(half):
            u, v = roots[offset*(length//width)]
            root = (u, -v if inverse else v)
            for i in range(offset, length, width):
                x, y = data[i]
                z, w = _complex_product(root, data[i+half], bits)
                data[i], data[i+half] = (x+z, y+w), (x-z, y-w)
        width *= 2


def _dyadic_convolution(a, b, modulus):
    """Return exact coefficients modulo modulus with no native transform ceiling.

    Let L be the padded length, M>=1 the coefficient height, S the scale.
    Angle halving has maximum coordinate error e' <= 2e+3/S. The cosine
    error is <=e/2+1/S by the square root's Lipschitz constant; its rounded
    value stays >=1/2. Sine division then costs <=e+2(e/2+1/S)+1/S.
    After log2(L)-2 halvings the primitive root has complex error <=3L/S.
    All complex errors below are Euclidean. Rounding a product by flooring
    its two coordinates costs <2/S. Twiddle powers have error recurrence
    E' <= (1+3L/S)E+(3L+2)/S, starting at E=0. Since (1+3L/S)^L <= exp(1)
    for S>=3L^2, all twiddle errors are <=16L^2/S.

    For true inputs bounded by B>=1 and initial complex error E, a butterfly
    stage maps error e to at most 3e+(16L^3 B+2)/S (true intermediate
    magnitude <=LB). Since 3^log2(L)<=L^2, a transform has final error
    <=L^2 E+18L^5 B/S. Each forward transform starts exact. The pointwise
    products then have error <=56L^6 M^2/S and true size <=L^2 M^2.
    Applying the inverse bound and dividing by L gives <=74L^7 M^2/S.
    Our scale exceeds 1024L^8 M^2, so nearest-integer rounding is exact.

    The O(L log L) butterflies use operands of O(log L + log M) bits;
    ordinary integer multiplication and square root suffice for a uniform
    O(L log L poly(log L + log M)) bit bound.
    """
    if not a or not b:
        return []
    size = len(a)+len(b)-1
    if size == 1:
        return [a[0]*b[0] % modulus]
    depth = (size-1).bit_length()
    length = 1 << depth
    height = max(1, max(map(abs,a)), max(map(abs,b)))
    bits = 10 + 8*depth + 2*height.bit_length()
    roots = _twiddles(length, bits)
    left = [(x << bits, 0) for x in a] + [(0,0)]*(length-len(a))
    right = [(x << bits, 0) for x in b] + [(0,0)]*(length-len(b))
    _fft(left, roots, bits)
    _fft(right, roots, bits)
    products = [_complex_product(x,y,bits) for x,y in zip(left,right)]
    _fft(products, roots, bits, inverse=True)
    # Combine inverse normalization and nearest-integer rounding in one shift.
    shift = bits + depth
    half = 1 << (shift-1)
    return [((x+half) >> shift) % modulus for x,_ in products[:size]]


class _DecimalRing(PolyRing):
    """List polynomials with bounded-size decimal packing and a uniform FFT fallback."""

    def __init__(self, n, native_digit_limit=_NATIVE_DIGIT_LIMIT):
        super().__init__(n, "py")
        self.native_digit_limit = native_digit_limit
        try:
            import _decimal
        except ImportError:
            _decimal = None
        self.decimal = _decimal

    def _kara(self, a, b):
        if min(len(a), len(b)) <= 20:
            return self._school(a, b)
        decimal = self.decimal
        # 10^width exceeds the greatest possible integer convolution coefficient.
        width = (min(len(a), len(b)) * (self.N-1)**2).bit_length()
        digit_count = width * (len(a) + len(b))
        if decimal is None or digit_count > self.native_digit_limit:
            return _dyadic_convolution(a, b, self.N)
        context = decimal.Context(prec=digit_count, Emin=decimal.MIN_EMIN,
                                  Emax=decimal.MAX_EMAX, traps=[decimal.Inexact, decimal.Overflow])

        def pack(coefficients):
            # Convert each small coefficient separately. Packing and unpacking
            # the full long operand through Python binary integers would lose
            # the quasi-linear bound in the length of the polynomial.
            return decimal.Decimal("".join(format(decimal.Decimal(x), "f").zfill(width)
                                           for x in reversed(coefficients)))

        blocks = len(a) + len(b) - 1
        digits = format(context.multiply(pack(a), pack(b)), "f").zfill(blocks*width)
        return [int(decimal.Decimal(digits[i:i+width])) % self.N
                for i in range((blocks-1)*width, -1, -width)]


def _least_nonunit(n, upper, backend):
    """Least nonunit integer in [2, upper], hence a prime, or None."""
    if upper < 2:
        return None
    length = upper-1
    block = isqrt(length)
    packed_size_bound = (4*block+16) * (2*n.bit_length() + block.bit_length() + 4)
    if backend == "py" or (backend == "flint" and packed_size_bound <= _NATIVE_DIGIT_LIMIT):
        return ps_interval(n, 2, upper+1, backend=backend)["witness"]
    count = length//block
    ring = _DecimalRing(n, native_digit_limit=0 if backend == "dyadic" else _NATIVE_DIGIT_LIMIT)
    polynomial = product_tree_linear(ring, range(2, 2+block))[-1][0]
    products = multipoint_eval(ring, polynomial, [j*block for j in range(count)])
    for j, value in enumerate(products):
        if gcd(n, value) != 1:
            start = 2 + j*block
            return next(x for x in range(start, start+block) if gcd(n, x) != 1)
    return next((x for x in range(2+count*block, upper+1) if gcd(n, x) != 1), None)


def _reveal(remainder, low, high, backend):
    """Remove all prime powers with bases <= high, retaining band bases only."""
    band = []
    while remainder > 1:
        cutoff = min(high, isqrt(remainder))
        p = _least_nonunit(remainder, cutoff, backend)
        if p is None:
            if remainder > high:
                break
            # The full square-root prefix was coprime, certifying primality.
            p = remainder
        if p >= low:
            band.append(p)
        while remainder % p == 0:
            remainder //= p
    return remainder, band


def _confirmed_band(primes, constant, gate, coin):
    """The outer gate has opened. Return an independently marked subset in law
    after accounting for the alternative that the outer gate stayed closed.
    """
    probabilities = [sqrt_mark_probability(p, constant) for p in primes]
    empty = [Fraction(1)] * (len(primes)+1)
    for i in reversed(range(len(primes))):
        empty[i] = empty[i+1] * (1-probabilities[i])
    nonempty = 1-empty[0]
    if not nonempty or not coin(nonempty/gate):
        return []
    marked = []
    for i, (p, h) in enumerate(zip(primes, probabilities)):
        # Until a first mark occurs, condition the remaining suffix on a hit.
        if coin(h if marked else h/(1-empty[i])):
            marked.append(p)
    return marked


def sqrt_prime_marks(modulus, rng=None, *, constant=8, backend="auto"):
    """Yield each marked distinct prime divisor of the public modulus, in order.

    For fixed C this takes expected polynomial bit time with 'auto', 'flint',
    'decimal' or 'dyadic'. Native FLINT or libmpdec accelerate bounded sizes.
    Larger convolutions use the fixed-point transform in this module. An
    individual opened high band can cost exponentially more than the mean.
    """
    if modulus < 1 or not isinstance(constant, int) or constant < 0:
        raise ValueError("expected positive modulus and nonnegative integer constant")
    if backend == "auto":
        backend = "flint" if find_spec("flint") is not None else "decimal"
    if backend not in ("flint", "decimal", "dyadic", "py"):
        raise ValueError("expected auto, flint, decimal, dyadic, or py backend")
    if not constant:
        return
    rng = SystemRandom() if rng is None else rng

    def coin(probability):
        if not 0 <= probability <= 1:
            raise ArithmeticError("invalid band or subset probability")
        return rng.randrange(probability.denominator) < probability.numerator

    remainder = modulus
    bits = modulus.bit_length()-1
    for j in range(1, bits+1):
        low = 1 << j
        # Since N >= low^count, there are at most floor(log2(N)/j) band primes.
        gate = min(Fraction(1), Fraction(constant*(bits//j), isqrt(low-1)+1))
        if coin(gate):
            remainder, primes = _reveal(remainder, low, 2*low-1, backend)
            yield from _confirmed_band(primes, constant, gate, coin)
            if remainder == 1:
                return
