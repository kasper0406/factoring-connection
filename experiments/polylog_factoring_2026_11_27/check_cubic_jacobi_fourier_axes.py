"""Check additive Fourier axes of public cubic Jacobi signatures."""

from itertools import product
from math import gcd

import numpy as np

from eisenstein_torus import cubic_jacobi


def symbol(n, x):
    denominator = (x*x + x + 1) % n
    if gcd(denominator, n) != 1:
        return None
    scale = pow(denominator, -1, n)
    return cubic_jacobi(n, ((x*x - 1)*scale % n, (-2*x - 1)*scale % n))


def signature(local, x, shifts):
    values = tuple(local[(x + a) % len(local)] for a in shifts)
    return None if None in values else values


def label_law(local, shifts):
    j = len(shifts)
    law = np.zeros((3,)*j)
    for x in range(len(local)):
        key = signature(local, x, shifts)
        if key is not None:
            law[key] += 1 / len(local)
    valid = law.sum()
    assert valid > 0
    law /= valid
    bias = np.abs(np.fft.fftn(law)).ravel()[1:].max()
    return law, valid, bias


def verify_axes(p, q, shifts):
    n = p*q
    at_p = [symbol(p, x) for x in range(p)]
    at_q = [symbol(q, x) for x in range(q)]
    public = [symbol(n, x) for x in range(n)]
    for x in range(n):
        a, b = at_p[x % p], at_q[x % q]
        assert public[x] == (None if None in (a, b) else (a+b) % 3)

    j = len(shifts)
    plaw, tp, bp = label_law(at_p, shifts)
    qlaw, tq, bq = label_law(at_q, shifts)
    keys = list(product(range(3), repeat=j))
    k = len(keys)
    actual_keys = [signature(public, x, shifts) for x in range(n)]
    freq = np.arange(n)
    axis_p = (freq % q == 0) & (freq != 0)
    axis_q = (freq % p == 0) & (freq != 0)
    tables = []
    for key in keys:
        tables.append({a: int(a == key) for a in keys})
    root = np.exp(2j*np.pi/3)
    tables.append({a: root**(sum(a) % 3) for a in keys})
    tables.append({a: 1 if (sum(v*v for v in a) % 3 == 1) else -1 for a in keys})
    verified = nonvacuous = 0
    for table in tables:
        state = np.array([0 if key is None else table[key] for key in actual_keys])
        total = np.mean(np.abs(state)**2)
        if total == 0:
            continue
        energies = np.abs(np.fft.fft(state)/n)**2
        from_p, from_q = energies[axis_p].sum()/total, energies[axis_q].sum()/total
        grid = np.zeros((p, q), dtype=complex)
        for x in range(n):
            grid[x % p, x % q] = state[x]
        assert abs(np.var(grid.mean(axis=1))/total - from_p) < 1e-11
        assert abs(np.var(grid.mean(axis=0))/total - from_q) < 1e-11
        eta = (k-1)*bp*bq
        if eta < 1:
            bound = (1-tp) + (1+(k-1)*bp)*bq*bq/(1-eta)
            swapped = (1-tq) + (1+(k-1)*bq)*bp*bp/(1-eta)
            assert from_p <= bound + 1e-11 and from_q <= swapped + 1e-11
            nonvacuous += bound + swapped < 1
        if j == 1 and all(abs(table[(v,)] - root**v) < 1e-11 for v in range(3)):
            assert abs(from_p - (tp - 1/p**2)/(q**2*tp*tq)) < 1e-11
            assert abs(from_q - (tq - 1/q**2)/(p**2*tp*tq)) < 1e-11
        verified += 1

    combined = np.fft.ifftn(np.fft.fftn(plaw)*np.fft.fftn(qlaw)).real
    predicted = np.array([actual_keys.count(a)/n for a in keys]).reshape((3,)*j)
    assert np.max(np.abs(combined*tp*tq - predicted)) < 1e-11
    return verified, nonvacuous


def local_word_spectrum(keys, exponents):
    roots = np.exp(2j*np.pi*np.arange(3)/3)
    word = [0 if key is None else roots[sum(a*b for a, b in zip(exponents, key)) % 3]
            for key in keys]
    return np.fft.fft(word)/len(keys)


def verify_quarter(p, q, shifts):
    """Check local twists, CRT, and the optimal label-table prediction on unit scales."""
    n, j = p*q, len(shifts)
    loc = []
    for r in (p, q):
        values = [symbol(r, x) for x in range(r)]
        loc.append([signature(values, x, shifts) for x in range(r)])
    public = [signature_key(loc[0][x % p], loc[1][x % q]) for x in range(n)]
    positions = np.arange(n)
    valid = np.array([key is not None for key in public])
    units = (positions % p != 0) & (positions % q != 0)
    freq_p = positions * pow(q, -1, p) % p
    freq_q = positions * pow(p, -1, q) % q
    max_word = np.zeros(n)
    for exponents in list(product(range(3), repeat=j))[1:]:
        p_hat, q_hat = [local_word_spectrum(keys, exponents) for keys in loc]
        for r, spectrum in ((p, p_hat), (q, q_hat)):
            assert np.abs(spectrum).max() <= 4*j/r**0.5 + 6*j/r + 1e-12
        crt = p_hat[freq_p] * q_hat[freq_q]
        assert np.max(np.abs(local_word_spectrum(public, exponents) - crt)) < 1e-11
        max_word = np.maximum(max_word, np.abs(crt))

    mixed = (1 + p*((-2*pow(p, -1, q)) % q)) % n
    assert mixed % p == 1 and mixed % q == q-1
    results = []
    for u in (1, 2, mixed):
        argument = positions*u % n
        target = 2*(np.minimum(argument, n-argument) <= n//4).astype(int)-1
        bins = {}
        for x in positions[valid & units]:
            bins[public[x]] = bins.get(public[x], 0) + target[x]
        # Every valid signature may choose its own majority sign. On singular
        # unit scales give the predictor the correct target for free.
        optimum = (sum(abs(v) for v in bins.values()) + np.sum(units & ~valid))/units.sum()
        g_hat = np.abs(np.fft.fft(target)/n)
        capacity = ((3**j-1)**0.5 * np.dot(max_word, g_hat)
                    + abs(np.mean(target*valid)) + np.mean(valid & ~units)
                    + np.mean(units & ~valid))/np.mean(units)
        assert optimum <= capacity + 1e-11
        results.append((optimum, capacity))
    return results


def signature_key(a, b):
    if a is None or b is None:
        return None
    return tuple((x+y) % 3 for x, y in zip(a, b))


def verify_divisors(r, shifts):
    """For split and inert r use pairs of F_r coordinates for a geometric root."""
    split_root = next((a for a in range(r) if (a*a+a+1) % r == 0), None)
    if split_root is None:
        root, conjugate = (0, 1), (-1, -1)
    else:
        root, conjugate = (split_root, 0), (-1-split_root, 0)
    checked = 0
    for e in product(range(3), repeat=len(shifts)):
        divisor = {}
        local_weights = {}
        for a, coefficient in zip(shifts, e):
            local_weights[a % r] = (local_weights.get(a % r, 0) + coefficient) % 3
            for point, sign in ((root, 1), (conjugate, -1)):
                key = ((point[0]-a) % r, point[1] % r)
                divisor[key] = (divisor.get(key, 0) + sign*coefficient) % 3
        assert any(divisor.values()) == any(local_weights.values())
        checked += 1
    return checked


def check_injective(p, q, count):
    n = p*q
    labels = [symbol(n, x) for x in range(n)]
    assert None not in labels
    for j in range(1, count + 1):
        signatures = {signature(labels, x, range(j)) for x in range(n)}
        assert (len(signatures) == n) == (j == count)


if __name__ == "__main__":
    tested = nonvacuous = 0
    quarter_cases = []
    for p, q in ((7, 11), (13, 17), (37, 53), (101, 103)):
        for shifts in ((0,), (0, 1), (0, 1, 3)):
            cases, informative = verify_axes(p, q, shifts)
            tested += cases
            nonvacuous += informative
            quarter_cases.extend(verify_quarter(p, q, shifts))
    divisor_checks = sum(verify_divisors(r, shifts) for r in (5, 7, 11, 13)
                         for shifts in ((0, 1, 2), (0, 2, r, 2*r+1)))
    check_injective(17, 23, 11)
    check_injective(23, 29, 13)
    print(f"{tested} nonlinear Fourier/variance checks, {nonvacuous} with combined bound <1.")
    print(f"{divisor_checks} divisor-profile checks and two exact injectivity thresholds.")
    print(f"{len(quarter_cases)} optimal quarter-table checks, "
          f"{sum(bound < 1 for _, bound in quarter_cases)} with bound <1.")
    print(f"101 x 103, mixed root, three shifts: optimal correlation including singular gift "
          f"{quarter_cases[-1][0]:.6f}, finite Fourier bound {quarter_cases[-1][1]:.6f}.")
