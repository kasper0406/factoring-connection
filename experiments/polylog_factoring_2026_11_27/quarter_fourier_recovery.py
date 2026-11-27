"""Sample-heavy-frequency recovery from a scaled quarter-window predictor."""

from fractions import Fraction
from math import ceil, cos, gcd, isqrt, pi, sin

from mpmath import mp


def shells(n):
    half = n//2
    return [(0, 0)] + [(1 << j, min(half, (2 << j)-1)) for j in range(half.bit_length())]


def kernel_weight(n, low, high, signed_t, importance, precise=False):
    size = high-low+1
    if signed_t == 0:
        return size*importance
    phase = ((high+low)*signed_t) % (2*n)
    oscillation = (size*signed_t) % (2*n)
    if precise:
        return mp.cospi(mp.mpf(phase)/n)*mp.sinpi(mp.mpf(oscillation)/n) / mp.sinpi(mp.mpf(signed_t)/n)*importance
    return cos(pi*phase/n)*sin(pi*oscillation/n)/sin(pi*signed_t/n)*importance


def estimator(n, predict, rng, samples, precise=False):
    bins = shells(n)
    shell_count = len(bins)

    products = {}
    for _ in range(samples):
        left, right = rng.choice(bins)
        distance = rng.randint(left, right)
        t = -distance if rng.getrandbits(1) else distance
        count = 2*(right-left+1) if distance else 1
        x = rng.randrange(n)
        first = predict(x)
        if not first:
            continue
        second = predict((x+t) % n)
        if second:
            products[t] = products.get(t, 0) + first*second*shell_count*count

    weights = [
        (t, (mp.mpf(value) if precise else value)/n)
        for t, value in products.items() if value
    ]

    def energy(low, high):
        return sum(
            kernel_weight(n, low, high, t, weight, precise)
            for t, weight in weights
        )/samples

    return energy


def heavy_frequencies(n, threshold, energy, precise=False):
    threshold = Fraction(threshold)
    cap = ceil(4/threshold**2)
    cutoff = threshold**2/2
    cutoff = mp.mpf(cutoff.numerator)/cutoff.denominator if precise else float(cutoff)
    intervals, frequencies = [(0, n-1)], []
    while intervals:
        children = []
        retained = 0
        for low, high in intervals:
            if energy(low, high) < cutoff:
                continue
            retained += 1
            if retained > cap:
                return []
            if low == high:
                frequencies.append(low)
                continue
            mid = (low+high)//2
            children.extend(((low, mid), (mid+1, high)))
        intervals = children
    return frequencies


def extract_divisor(n, frequencies, harmonics):
    for h in range(1, harmonics+1):
        divisor = gcd(h, n)
        if 1 < divisor < n:
            return divisor
        for frequency in frequencies:
            for h_signed in (h, -h):
                divisor = gcd(frequency-h_signed, n)
                if 1 < divisor < n:
                    return divisor
    return None


class FoundDivisor(Exception):
    def __init__(self, divisor):
        self.divisor = divisor


def recover(n, predictor, rng, threshold, harmonics, samples, probe_samples=True, precise=False):
    def predict(x):
        d = gcd(x, n)
        if d == 1:
            return predictor(x)
        if probe_samples and 1 < d < n:
            raise FoundDivisor(d)
        return 0

    try:
        energy = estimator(n, predict, rng, samples, precise)
        frequencies = heavy_frequencies(n, threshold, energy, precise)
    except FoundDivisor as answer:
        return answer.divisor, []
    return extract_divisor(n, frequencies, harmonics), frequencies


def factor_from_prediction(n, predictor, rng, advantage, failure=Fraction(1, 10)):
    """Use the proven bounds, which favor a conservative sample size."""
    advantage = Fraction(advantage)
    failure = Fraction(failure)
    if not 0 < advantage <= Fraction(1, 2) or not 0 < failure < 1:
        raise ValueError("Expected a positive advantage and a failure probability in (0,1)")
    harmonics = ceil(32/advantage**2)
    if n <= 2*harmonics+1:
        return next((d for d in range(2, isqrt(n)+1) if n % d == 0), None)
    threshold = advantage / (4*(1+harmonics.bit_length()))
    energy_error = threshold**2/8
    budget = 2*n-1  # All nodes of the fixed interval tree.
    confidence_bits = ceil(2*budget/failure).bit_length()
    samples = ceil(2*len(shells(n))**2*confidence_bits/energy_error**2)
    precision = n.bit_length() + (1/energy_error).numerator.bit_length() + samples.bit_length() + 64
    with mp.workprec(precision):
        divisor, _ = recover(n, predictor, rng, threshold, harmonics, samples, precise=True)
    return divisor


if __name__ == "__main__":
    from random import Random
    rng = Random(261128)
    for p, q in [(11, 13), (101, 109), (503, 937), (101**2, 109), (153247, 254119)]:
        n = p*q
        hidden = 1+p*(-2*pow(p, -1, q) % q)
        calls = [0]

        def noisy_predictor(scale):
            calls[0] += 1
            value = scale*hidden % n
            result = 1 if min(value, n-value) <= n//4 else -1
            return -result if rng.random() < .12 else result

        divisor, frequencies = recover(n, noisy_predictor, rng, Fraction(1, 4), 1, 16_384, probe_samples=False)
        assert divisor in (p, q)
        assert any(frequency in (hidden, n-hidden) for frequency in frequencies)
        print(n, divisor, calls[0], len(frequencies))
