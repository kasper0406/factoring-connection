"""HARNESS: orbits, Fourier spectra on Z/r, metrics, and the transfer identity (uses factors).

Floating point: spectra are computed with numpy's float64 FFT. Its error is <= ~ c*log2(r)*2^-53
relative to ||F||_2 (Higham), i.e. < 1e-12 relative mass at r <= 2^26. Every reported comparison has
orders of magnitude of margin, which is the only use made of these numbers (AGENT_RULES rule 1:
floats only for measurement, never inside an operational decision).
"""
from __future__ import annotations

import math
from typing import Dict

import gmpy2
import numpy as np
import sympy


def order_mod_prime(a: int, p: int) -> int:
    o = p - 1
    for ell, e in sympy.factorint(p - 1).items():
        for _ in range(e):
            if pow(a, o // ell, p) == 1:
                o //= ell
            else:
                break
    return o


def order_info(a: int, p: int, q: int) -> Dict[str, int]:
    rp, rq = order_mod_prime(a % p, p), order_mod_prime(a % q, q)
    return {"r_p": rp, "r_q": rq, "r": rp * rq // math.gcd(rp, rq), "g": math.gcd(rp, rq)}


def orbit(a: int, N: int, r: int) -> np.ndarray:
    """[a^0, ..., a^{r-1}] mod N by vectorised doubling (needs N < 2^31)."""
    assert N < 2 ** 31
    out = np.empty(r, dtype=np.int64)
    out[0] = 1
    k = 1
    while k < r:
        step = pow(a, k, N)
        m = min(k, r - k)
        out[k:k + m] = (out[:m] * step) % N
        k += m
    return out


def spectrum(F: np.ndarray):
    """Relative masses |F^(j)|^2 / ||F||^2 on Z/r, F^(j) = E_t F(t) e(-jt/r).

    Real F: returns (P_half, True) for j = 0..r//2 (mass at j and r-j are equal); complex F: (P, False).
    """
    r = len(F)
    norm = float(np.mean(np.abs(F) ** 2))
    if np.isrealobj(F):
        P = np.abs(np.fft.rfft(F) / r) ** 2
        half = True
    else:
        P = np.abs(np.fft.fft(F) / r) ** 2
        half = False
    if norm > 0:
        P /= norm
    return P, half, norm


def metrics(spec, r: int, jac_a: int, chunk: int = 1 << 22) -> Dict[str, float]:
    """M_peak over j with eta_j outside {1, Jacobi|_H}; H_mass on orders >= sqrt(r). Chunked."""
    P, half, norm = spec
    s = math.isqrt(r)
    best, bj, hmass, exc = -1.0, 0, 0.0, 0.0
    L = len(P)
    for lo in range(0, L, chunk):
        hi_ = min(L, lo + chunk)
        j = np.arange(lo, hi_, dtype=np.int64)
        Pc = P[lo:hi_]
        w = np.ones(len(j))
        if half:
            w = np.where((j == 0) | (2 * j == r), 1.0, 2.0)
        excl = (j == 0) | ((jac_a == -1) & (2 * j == r))
        exc += float((Pc * w)[excl].sum())
        orders = r // np.gcd(j, r)
        hmass += float((Pc * w)[orders >= s].sum())
        Pm = np.where(excl, -1.0, Pc)
        k = int(np.argmax(Pm))
        if Pm[k] > best:
            best, bj = float(Pm[k]), int(j[k])
        del j, Pc, w, excl, orders, Pm
    return {"M_peak": best, "M_peak_abs": best * norm, "argmax_j": bj,
            "argmax_order": r // math.gcd(bj, r), "H_mass": hmass, "exceptional_mass": exc,
            "rand_baseline": math.log(r) / r, "norm2": norm}


# ---------------------------------------------------------------- transfer identity (tiny N)

def primitive_root(p: int) -> int:
    return int(sympy.primitive_root(p))


def transfer_check(N: int, p: int, q: int, a: int, gfun) -> float:
    """Max |F^(j) - sum_{chi(a)=e(j/r)} g^(chi)| over j. Characters chi_{u,v} of (Z/p)^* x (Z/q)^*."""
    gp, gq = primitive_root(p), primitive_root(q)
    P1, Q1 = p - 1, q - 1
    pw_p = [pow(gp, i, p) for i in range(P1)]
    pw_q = [pow(gq, k, q) for k in range(Q1)]
    # CRT grid y(i,k)
    ip = pow(p, -1, q)
    Y = np.empty((P1, Q1), dtype=np.int64)
    for i in range(P1):
        yp = pw_p[i]
        for k in range(Q1):
            Y[i, k] = yp + p * (((pw_q[k] - yp) * ip) % q)
    G = gfun(Y.ravel(), N).reshape(P1, Q1)
    phi = P1 * Q1
    # ghat(u,v) = (1/phi) sum g(y) conj(chi_{u,v}(y)), chi_{u,v}(y(i,k)) = e(ui/P1 + vk/Q1)
    ghat = np.fft.fft2(G) / phi
    log_p = {v: i for i, v in enumerate(pw_p)}
    log_q = {v: k for k, v in enumerate(pw_q)}
    ia, ka = log_p[a % p], log_q[a % q]
    info = order_info(a, p, q)
    r = info["r"]
    acc = np.zeros(r, dtype=np.complex128)
    for u in range(P1):
        for v in range(Q1):
            num = r * (u * ia * Q1 + v * ka * P1)
            assert num % (P1 * Q1) == 0
            acc[(num // (P1 * Q1)) % r] += ghat[u, v]
    F = gfun(orbit(a, N, r), N)
    Fh = np.fft.fft(F) / r
    return float(np.max(np.abs(Fh - acc)))
