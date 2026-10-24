"""DIAGNOSTIC oracles: geometric quantities computed FROM THE FACTORS p, q (never imported by decoders.py)."""
from math import gcd

from src.geometry import modcurve as mc


def genus_oracle(p: int, q: int) -> int:
    return mc.genus({p: 1, q: 1})


def index_oracle(p: int, q: int) -> int:
    """psi(N) = |P^1(Z/N)| = [SL_2(Z):Gamma_0(N)] = (3/pi) * hyperbolic area of X_0(N)."""
    return (p + 1) * (q + 1)


def pillai_oracle(p: int, q: int) -> int:
    """sum_{g=1}^{N} gcd(g, N) = sum_{d|N} d phi(N/d) = total number of components of the star polygons {N/g}."""
    return (2 * p - 1) * (2 * q - 1)


def sigma_oracle(p: int, q: int) -> int:
    """sigma(N) = number of sublattices of Z^2 of index N (= connected N-sheeted covers of the torus, Hecke)."""
    return (1 + p) * (1 + q)


def necklace_oracle(p: int, q: int, k: int, r: int) -> int:
    """sum_{d|N} phi(d) k^{N/d} mod r  (Burnside: N times the number of k-colour necklaces of length N)."""
    N = p * q
    return (pow(k, N, r) + (p - 1) * pow(k, q, r) + (q - 1) * pow(k, p, r) + (p - 1) * (q - 1) * k) % r


def lucas_V(t: int, n: int, r: int) -> int:
    """V_n(t) = lam^n + lam^{-n}, lam + 1/lam = t, mod r (V_0 = 2, V_1 = t)."""
    # binary ladder on (V_m, V_{m+1})
    a, b = 2 % r, t % r
    for bit in bin(n)[2:]:
        if bit == "1":
            a, b = (a * b - t) % r, (b * b - 2) % r
        else:
            a, b = (a * a - 2) % r, (a * b - t) % r
    return a


def catmap_oracle(p: int, q: int, t: int, r: int) -> int:
    """Primitive periodic points of period N of a hyperbolic toral automorphism with trace t, mod r:
    sum_{d|N} mu(N/d) |det(A^d - I)| = sum mu(N/d) (V_d(t) - 2)   (t >= 3)."""
    N = p * q
    F = lambda d: (lucas_V(t, d, r) - 2) % r
    return (F(N) - F(p) - F(q) + F(1)) % r


def x0_count_oracle(p: int, q: int, l: int) -> int:
    """#X_0(N)(F_l) = l + 1 - Tr(T_l | S_2(Gamma_0(N)))  (Eichler-Shimura; trace formula with the factors)."""
    return l + 1 - mc.trace_formula({p: 1, q: 1}, l, brute=False)


def x0_count_sq_oracle(p: int, q: int, l: int) -> int:
    """#X_0(N)(F_{l^2}) from the factors (amendment A1)."""
    return mc.x0_count_sq({p: 1, q: 1}, l, brute=False)
