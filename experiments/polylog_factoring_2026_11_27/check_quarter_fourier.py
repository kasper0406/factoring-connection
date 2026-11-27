"""Check Fourier localization for a planted quarter-window predictor."""

from cmath import exp, pi
from math import gcd, sin
from random import Random


def quarter(n, value):
    value %= n
    return 1 if min(value, n-value) <= n//4 else -1


def shells(n):
    half = n//2
    return [(0, 0)] + [(1 << j, min(half, (2 << j)-1)) for j in range(half.bit_length())]


def kernel(n, low, high, signed_t):
    size = high-low+1
    if signed_t == 0:
        return size
    phase = ((high+low)*signed_t) % (2*n)
    oscillation = (size*signed_t) % (2*n)
    return exp(pi*1j*phase/n) * sin(pi*oscillation/n) / sin(pi*signed_t/n)


def energy_estimate(n, low, high, predict, rng, trials):
    bins = shells(n)
    estimate = 0.0
    for _ in range(trials):
        left, right = rng.choice(bins)
        distance = rng.randint(left, right)
        signed_t = rng.choice((-distance, distance))
        size = 2*(right-left+1) if distance else 1
        x = rng.randrange(n)
        value = predict(x) * predict((x+signed_t) % n)
        weight = kernel(n, low, high, signed_t).real * len(bins)*size/n
        assert abs(weight) <= len(bins)+1e-8
        estimate += value*weight
    return estimate/trials


def locate(n, threshold, energy):
    intervals = [(0, n-1)]
    result = []
    while intervals:
        low, high = intervals.pop()
        if energy(low, high) < threshold**2/2:
            continue
        if low == high:
            result.append(low)
        else:
            mid = (low+high)//2
            intervals.extend([(low, mid), (mid+1, high)])
    return result


def exact_spectrum(values):
    n = len(values)
    phases = [exp(-2j*pi*x/n) for x in range(n)]
    return [
        sum(value*phases[k*x % n] for x, value in enumerate(values) if value)/n
        for k in range(n)
    ]


def check(p, q, rng):
    n = p*q
    hidden = 1 + p*(-2*pow(p, -1, q) % q)
    values = [
        quarter(n, u*hidden) * (-1 if rng.random() < .25 else 1)
        if gcd(u, n) == 1 else 0
        for u in range(n)
    ]
    spectrum = exact_spectrum(values)
    prefix = [0.0]
    for value in spectrum:
        prefix.append(prefix[-1]+abs(value)**2)

    def energy(low, high):
        return prefix[high+1]-prefix[low]

    found = locate(n, .15, energy)
    assert any(1 < gcd(value-1, n) < n for value in found)
    for low, high in [(0, n//3), (n//3, n-1), (hidden, hidden)]:
        measured = energy_estimate(n, low, high, values.__getitem__, rng, 65_536)
        assert abs(measured-energy(low, high)) < .055
    print(n, len(found), sorted(gcd(value-1, n) for value in found if gcd(value-1,n) > 1))


if __name__ == "__main__":
    rng = Random(261127)
    for p, q in [(11, 13), (13, 19), (17, 23), (29, 37)]:
        check(p, q, rng)
