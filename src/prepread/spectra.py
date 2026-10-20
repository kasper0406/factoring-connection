"""Mixed multiplicative/additive spectra A_p(j, c) = (1/(p-1)) sum_t e(c*v_t/p) e(-j t/(p-1)),
v_t = Phi(g^t) mod p, computed by float64 FFT (disclosed; relative error ~1e-12).

HARNESS/diagnostic code: uses the factorisation of p - 1 (the prime's own data) for the
primitive root and for the large-prime-order frequency set J_p.
"""
from __future__ import annotations
import numpy as np
from sympy import factorint


def primitive_root(p: int) -> int:
    fac = factorint(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // l, p) != 1 for l in fac):
            return g
    raise ValueError


def orbit(p: int, g: int) -> np.ndarray:
    """y_t = g^t mod p for t = 0..p-2 (int64)."""
    out = np.empty(p - 1, dtype=np.int64)
    y = 1
    for t in range(p - 1):
        out[t] = y
        y = y * g % p
    return out


def large_order_mask(p: int) -> np.ndarray:
    """J_p: frequencies j whose character order (p-1)/gcd(j,p-1) is divisible by l = P^+(p-1)."""
    fac = factorint(p - 1)
    l = max(fac)
    lv = l ** fac[l]
    j = np.arange(p - 1)
    return (j % lv) != 0


def spectrum_max(v: np.ndarray, p: int, cs, mask: np.ndarray, batch: int = 256):
    """max over c in cs and j in mask of |A_p(j, c)|; also return the argmax (c, j) and the
    max over all nonzero j (for reference)."""
    n = p - 1
    table = np.exp(2j * np.pi * np.arange(p) / p)
    cs = np.asarray(list(cs), dtype=np.int64)
    best = (-1.0, None, None)
    for s in range(0, len(cs), batch):
        cb = cs[s:s + batch]
        idx = (cb[:, None] * v[None, :]) % p
        A = np.abs(np.fft.fft(table[idx], axis=1)) / n
        Am = np.where(mask[None, :], A, 0.0)
        k = int(np.argmax(Am))
        r, jj = divmod(k, n)
        if Am[r, jj] > best[0]:
            best = (float(Am[r, jj]), int(cb[r]), int(jj))
        del idx, A, Am
    return best


def A_value(v: np.ndarray, p: int, c: int, j: int) -> complex:
    n = p - 1
    t = np.arange(n)
    ph = np.exp(2j * np.pi * ((c * v) % p) / p) * np.exp(-2j * np.pi * (j * t % n) / n)
    return complex(ph.sum() / n)
