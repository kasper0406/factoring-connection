"""Hankel window bound for coefficient extraction (OPERATIONAL math helper; no N-specific data).

Prop HW (notes/nativeana_results.md): if m_n = sum_{i<T} c_i z_i^n (|z_i| <= 1, any c_i) then
    sum_{|n-N| <= W} |m_n - delta_{nN}| >= 1   whenever T <= W  (and N >= W).
Tight: T = W+1 nodes z_i = e(i/(W+1)), c_i = e(-iN/(W+1))/(W+1) give m_n = [n == N mod (W+1)].
"""
from __future__ import annotations

import numpy as np


def moments(nodes: np.ndarray, coeffs: np.ndarray, ns: np.ndarray) -> np.ndarray:
    return (coeffs[None, :] * nodes[None, :] ** ns[:, None]).sum(axis=1)


def window_l1_error(nodes, coeffs, N: int, W: int) -> float:
    ns = np.arange(N - W, N + W + 1)
    m = moments(np.asarray(nodes, dtype=complex), np.asarray(coeffs, dtype=complex), ns)
    m[W] -= 1.0
    return float(np.abs(m).sum())


def window_hankel_opnorm(nodes, coeffs, N: int, W: int) -> float:
    """Operator norm of the (W+1)x(W+1) block [m_{j+k} - delta_{j+k,N}]_{0<=j<=W, N-W<=k<=N}; Prop HW says >= 1
    if T <= W, and it is <= the window l1 error."""
    nodes = np.asarray(nodes, dtype=complex)
    coeffs = np.asarray(coeffs, dtype=complex)
    j = np.arange(0, W + 1)
    k = np.arange(N - W, N + 1)
    idx = j[:, None] + k[None, :]
    flat = np.unique(idx)
    m = dict(zip(flat.tolist(), moments(nodes, coeffs, flat)))
    H = np.vectorize(lambda t: m[t])(idx).astype(complex)
    H = H - (idx == N)
    return float(np.linalg.norm(H, 2))


def least_squares_coeffs(nodes, N: int, W: int) -> np.ndarray:
    """Best l2 fit of delta_N on the window by sum c_i z_i^n with the given nodes."""
    nodes = np.asarray(nodes, dtype=complex)
    ns = np.arange(N - W, N + W + 1)
    A = nodes[None, :] ** ns[:, None]
    b = (ns == N).astype(complex)
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    return c


def tight_example(N: int, W: int):
    T = W + 1
    i = np.arange(T)
    nodes = np.exp(2j * np.pi * i / T)
    coeffs = np.exp(-2j * np.pi * i * N / T) / T
    return nodes, coeffs
