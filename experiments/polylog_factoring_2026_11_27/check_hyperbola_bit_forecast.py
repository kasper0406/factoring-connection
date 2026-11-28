"""Independent grid calibration, including an unfitted comparison on distinct prime pairs."""

from itertools import combinations
from math import pi
from random import Random

import numpy as np
from sympy import primerange

if __package__:
    from .hyperbola_bit_forecast import upper_bit_forecast
else:
    from hyperbola_bit_forecast import upper_bit_forecast


def grid_forecast(z, r, d, h):
    """Integrate the log-prime density on each exact product, splitting at all jumps."""
    lam = h*r//d
    out = []
    indices = np.arange(2*lam,4*lam+1)
    cells = indices*d/(2*h)
    for n in z:
        low, high = max(r,n/(2*r)), min(2*r,n/r)
        knots = np.unique(np.clip(np.r_[low, high, cells, n/cells],low,high))
        mids = (knots[1:]+knots[:-1])/2
        log_x = np.log(knots)
        primitive = np.log(log_x/(np.log(n)-log_x))
        inside = ((2*h*mids/d).astype(int) & (2*h*n/mids/d).astype(int) & 1).astype(bool)
        weights = np.diff(primitive)
        out.append(weights[inside].sum()/weights.sum())
    return np.asarray(out)


def small():
    cases = 0
    for r in (16,32,64):
        primes = list(primerange(r,2*r))
        pairs = list(combinations(primes,2))
        ns = np.asarray([p*q for p,q in pairs])
        for h in (1,3):
            for d in (r,r//2,r//4):
                expected = grid_forecast(ns,r,d,h)
                actual = [float(upper_bit_forecast(int(n),r,d,h)) for n in ns]
                assert np.max(np.abs(expected-actual)) < 1e-11
                for (p,q),v in zip(pairs,actual):
                    if v == 0 or v == 1:
                        assert bool((h*p & h*q) & (d//2)) == v
                cases += len(ns)
    return cases


def cohort():
    r = 1 << 19
    primes = list(primerange(r,2*r))
    rng = Random(55741)
    pairs = [rng.sample(primes,2) for _ in range(20000)]
    p,q = np.array(pairs).T
    n = p*q
    for h in (1,3):
        for ratio in (4,16,64):
            d = r//ratio
            y = ((h*p & h*q & (d//2)) != 0).astype(float)
            mu = grid_forecast(n,r,d,h)
            score = -np.cos(4*pi*h*np.sqrt(n)/d+pi/4)
            forecast = 0.25+0.25*np.sqrt(d/h)*n.astype(float)**(-0.25)*score
            truth_const = np.mean((y-y.mean())**2)
            def gain(s):
                return truth_const-np.mean((y-s)**2)
            residual = (y-mu)*score
            print(h,ratio,*(f"{x:.6f}" for x in
                            (gain(forecast),gain(mu),np.mean(residual),
                             np.std(residual,ddof=1)/np.sqrt(len(y)))))


if __name__ == "__main__":
    print("decimal vs independent interval enumerator:",small())
    print("h R/D square-root gain, hyperbola gain, residual covariance term, sample SE")
    cohort()
