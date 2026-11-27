"""Check the quadratic-scale Salié spectrum and its autocorrelation."""

from math import gcd, sqrt

import numpy as np
from sympy import jacobi_symbol


def check(p, q, squarefree=True):
    n = p*q
    residues = np.arange(n, dtype=np.int64)
    units = np.array([gcd(int(c),n)==1 for c in residues])
    signal = np.zeros(n, dtype=complex)
    for m in residues[units]:
        value = int(m)
        inverse = pow(value,-1,n)
        signal[inverse] = int(jacobi_symbol(value,n))*np.exp(2j*np.pi*value/n)
    gauss = sqrt(n)*(1 if n % 4==1 else 1j)
    salie = np.fft.ifft(signal)*n/gauss
    values = salie[(residues**2*pow(4,-1,n))%n]
    spectrum = np.fft.fft(values)/n
    ramanujan = np.fft.ifft(units)*n
    assert np.max(np.abs(spectrum-ramanujan[(1-residues**2)%n]/n)) < 1e-12
    autocorrelation = np.fft.ifft(abs(np.fft.fft(values))**2)/n
    probability = sum(gcd(int(m)+1,n)==1 for m in residues[units])/n
    assert np.max(np.abs(autocorrelation[units]-probability*values[units])) < 1e-12
    if squarefree:
        assert abs(autocorrelation[0]-(4-(6*(p+q)-9)/n)) < 1e-12
    print(n, int(units.sum()), round(float(autocorrelation[0].real),8))


if __name__ == "__main__":
    for p,q in [(11,13),(31,37),(101,109)]:
        check(p,q)
    for p,q in [(9,25),(27,49)]:
        check(p,q,squarefree=False)
