"""Check a single phase interval and constant-bin approximation."""

from math import gcd, sqrt

import numpy as np
from sympy import jacobi_symbol

from check_salie_phase_weights import table


def interval_values(n, u):
    residues = np.arange(n)
    units = [c for c in range(n) if gcd(c,n)==1]
    character = [int(jacobi_symbol(c,n)) for c in range(n)]
    chi = np.array(character)
    if n%4==1:
        indicator = np.minimum(residues, n-residues) <= n//4
    else:
        indicator = residues < (n-1)//2
    all_shifted = np.fft.ifft(np.fft.fft(chi)*np.conj(np.fft.fft(indicator))).real
    f = sum(all_shifted[(residues*r)%n] for r in (1,-1,u,-u))/sqrt(n)
    return f, character, units, indicator


def check(p, q):
    n = p*q
    u = 1+p*(-2*pow(p,-1,q) % q)
    f, character, units, indicator = interval_values(n,u)
    residues = np.arange(n)
    spectrum = np.fft.fft(f)/n
    assert abs(spectrum[u]) > .28
    uf = np.zeros(n)
    uf[units] = f[units]
    assert abs(np.fft.fft(uf)[u]/n) > .18
    errors = []
    for count in (16,64,256):
        weight = np.exp(2j*np.pi*((count*residues)//n)/count)
        exact = np.exp(2j*np.pi*residues/n)
        squared_errors = []
        for c in units[:24]:
            coefficients = table(n,c,units,character)
            assert abs(sum(np.array(coefficients)[indicator])/sqrt(n)-f[c]) < 1e-9
            squared_errors.append(abs(np.dot(coefficients,weight-exact))**2/n)
        errors.append(sqrt(sum(squared_errors)/len(squared_errors)))
    print(n, round(abs(spectrum[u]),6), [round(e,6) for e in errors])


if __name__ == "__main__":
    for p,q in [(11,13),(13,17),(31,37),(101,109)]:
        check(p,q)
