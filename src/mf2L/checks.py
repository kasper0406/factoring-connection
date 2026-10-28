"""mf2L: exact checks of the finite statements in notes/mf2L_results.md (diagnostic; no factoring, no N).

- homogeneity(): Prop INV — for EVERY integral binary form H of degree e and primitive non-proportional P1, P2 with
  D = |det(P1, P2)|, H(P1) == lam^e H(P2) (mod D), lam the unit with P1 == lam P2 (mod D).  So no randomisation of the
  coefficients (g, w, c, extra multipliers, shifts) changes the pair lattice {x1 == lam^e x2 (mod D)}.
- iota_average(): Prop AVG — averaging the lattice pair count over ALL units iota gives exactly A* B* (no character
  term); characters enter only through the non-uniformity of iota over the pair family.
- dth_power_kernel(): the characters trivial on d-th powers (order | d) never average out on the algebraic side
  (iota = lam^d), while for a quadratic character and odd d, chi(lam^d) = chi(lam).
- nid_example(): Prop NID — two weightings with identical interval sums at every scale >= k but different pair counts
  on the lattice {x1 == x2 (mod k)}: interval densities do not determine pair covariances.
- gauss_transfer(): Lemma GS — for primitive chi mod k: |sum_n a_n chi(n)|^2 = |sum_b conj(chi)(b) S(b/k)|^2 / k and
  sum_{chi mod k} |sum_b conj(chi)(b) S(b/k)|^2 = phi(k) sum_{(b,k)=1} |S(b/k)|^2 — checked exactly in Q(zeta_m).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import gcd
from typing import Dict, List, Sequence, Tuple

import sympy


# ---------------------------------------------------------------- Prop INV
def _form_value(coeffs: Sequence[int], P: Tuple[int, int]) -> int:
    a, b = P
    e = len(coeffs) - 1
    return sum(c * a ** i * b ** (e - i) for i, c in enumerate(coeffs))


def lattice_unit(P1: Tuple[int, int], P2: Tuple[int, int]) -> Tuple[int, int]:
    """(D, lam) with D = |det(P1, P2)| and P1 == lam P2 (mod D) (pc3 Lemma LAT(a))."""
    D = abs(P1[0] * P2[1] - P1[1] * P2[0])
    if D == 1:
        return 1, 0
    for i in (0, 1):
        if gcd(P2[i], D) == 1:
            lam = P1[i] * pow(P2[i], -1, D) % D
            if (P1[0] - lam * P2[0]) % D == 0 and (P1[1] - lam * P2[1]) % D == 0:
                return D, lam
    # general primitive P2: use a combination coprime to D
    for s in range(1, D + 1):
        t = (s * P2[0] + P2[1]) % D
        if gcd(t, D) == 1:
            lam = (s * P1[0] + P1[1]) * pow(t, -1, D) % D
            assert (P1[0] - lam * P2[0]) % D == 0 and (P1[1] - lam * P2[1]) % D == 0
            return D, lam
    raise ValueError("no unit")


def homogeneity(coeffs: Sequence[int], P1, P2) -> bool:
    D, lam = lattice_unit(P1, P2)
    e = len(coeffs) - 1
    return (_form_value(coeffs, P1) - pow(lam, e, D) * _form_value(coeffs, P2)) % D == 0


# ---------------------------------------------------------------- Prop AVG
def iota_average(a: Dict[int, int], b: Dict[int, int], k: int) -> Tuple[int, int]:
    """Return (sum over units iota of the lattice pair count, A* B*) for finitely supported a, b >= 0."""
    units = [i for i in range(1, k) if gcd(i, k) == 1] if k > 1 else [0]
    total = 0
    for iota in units:
        for x1, w1 in a.items():
            for x2, w2 in b.items():
                if (x1 - iota * x2) % k == 0 and gcd(x2, k) == 1:
                    total += w1 * w2
    Astar = sum(w for x, w in a.items() if gcd(x, k) == 1)
    Bstar = sum(w for x, w in b.items() if gcd(x, k) == 1)
    return total, Astar * Bstar


# ---------------------------------------------------------------- d-th powers
def dth_power_kernel(k: int, d: int) -> Dict[str, int]:
    units = [i for i in range(1, k) if gcd(i, k) == 1]
    image = {pow(u, d, k) for u in units}
    # number of characters trivial on the d-th powers = index of the image = #{chi : chi^d = chi_0}
    return {"phi": len(units), "image": len(image), "chars_order_dividing_d": len(units) // len(image)}


def quadratic_odd_power(p: int, d: int) -> bool:
    """Legendre symbol (lam^d / p) == (lam / p) for all units lam, d odd."""
    return all(sympy.legendre_symbol(pow(l, d, p), p) == sympy.legendre_symbol(l, p) for l in range(1, p))


# ---------------------------------------------------------------- Prop NID
def nid_example(k: int, L: int, eps: Fr) -> Dict[str, Fr]:
    """Weights on [1, L]: s(n) = 1 (n coprime to k) versus t(n) = 1 + eps*chi(n) with chi the real character
    n -> +1 if n == 1 (mod k), -1 if n == -1 (mod k), 0 otherwise (needs k = 3, 4 or 6 so that this is a character).
    Returns max over intervals of length >= k of |sum t - sum s| / length, and the pair counts on {x1 == x2 (mod k)}."""
    assert k in (3, 4, 6)
    chi = {n: (1 if n % k == 1 else (-1 if n % k == k - 1 else 0)) for n in range(1, L + 1)}
    s = {n: Fr(1) if gcd(n, k) == 1 else Fr(0) for n in range(1, L + 1)}
    t = {n: s[n] * (1 + eps * chi[n]) for n in range(1, L + 1)}
    worst = Fr(0)
    for lo in range(1, L + 1):
        for hi in range(lo + k - 1, L + 1):
            dev = abs(sum(t[n] - s[n] for n in range(lo, hi + 1))) / (hi - lo + 1)
            worst = max(worst, dev)
    pair_s = sum(s[x] * s[y] for x in s for y in s if (x - y) % k == 0)
    pair_t = sum(t[x] * t[y] for x in t for y in t if (x - y) % k == 0)
    return {"max_interval_density_gap": worst, "pair_s": pair_s, "pair_t": pair_t,
            "sum_s": sum(s.values()), "sum_t": sum(t.values())}


# ---------------------------------------------------------------- Lemma GS (exact cyclotomic arithmetic)
class Cyc:
    """Elements of Q(zeta_m) as coefficient lists modulo the m-th cyclotomic polynomial."""

    def __init__(self, m: int):
        self.m = m
        x = sympy.Symbol("x")
        self.phi = [Fr(int(c)) for c in reversed(sympy.Poly(sympy.cyclotomic_poly(m, x), x).all_coeffs())]
        self.n = len(self.phi) - 1

    def reduce(self, c: List[Fr]) -> List[Fr]:
        c = list(c) + [Fr(0)] * max(0, self.n - len(c))
        for i in range(len(c) - 1, self.n - 1, -1):
            q = c[i]
            if q:
                for j in range(self.n + 1):
                    c[i - self.n + j] -= q * self.phi[j]
        return c[: self.n]

    def zeta(self, j: int) -> List[Fr]:
        c = [Fr(0)] * self.m
        c[j % self.m] = Fr(1)
        return self.reduce(c)

    def add(self, u, v):
        return [a + b for a, b in zip(u, v)]

    def scal(self, s, u):
        return [s * a for a in u]

    def mul(self, u, v):
        c = [Fr(0)] * (2 * self.n)
        for i, a in enumerate(u):
            if a:
                for j, b in enumerate(v):
                    c[i + j] += a * b
        return self.reduce(c)

    def conj(self, u):
        c = [Fr(0)] * self.m
        for i, a in enumerate(u):
            c[(-i) % self.m] += a
        return self.reduce(c)

    def zero(self):
        return [Fr(0)] * self.n


def _characters(k: int):
    """All Dirichlet characters mod an odd prime power k (cyclic group), as exponent maps; and primitivity."""
    p = sympy.primefactors(k)[0]
    g = sympy.primitive_root(k)
    ph = sympy.totient(k)
    log = {}
    x = 1
    for t in range(ph):
        log[x] = t
        x = x * g % k
    chars = []
    for j in range(ph):
        primitive = (k == p and j != 0) or (k != p and j % p != 0)
        chars.append((j, primitive))
    return int(ph), log, chars


def gauss_transfer(k: int, a: Dict[int, int]) -> Dict[str, bool]:
    ph, log, chars = _characters(k)
    m = k * ph
    K = Cyc(m)

    def chi_val(j, n, inv=False):  # chi_j(n) = zeta_ph^{j log n}; values as zeta_m powers
        if gcd(n, k) != 1:
            return K.zero()
        e = j * log[n % k] * (m // ph)
        return K.zeta(-e if inv else e)

    def S(b):  # sum_n a_n e(b n / k)
        out = K.zero()
        for n, w in a.items():
            out = K.add(out, K.scal(Fr(w), K.zeta(b * n * (m // k))))
        return out

    Sb = {b: S(b) for b in range(k)}
    ok_prim, lhs_total = True, K.zero()
    for j, primitive in chars:
        T = K.zero()
        for b in range(1, k):
            if gcd(b, k) == 1:
                T = K.add(T, K.mul(chi_val(j, b, inv=True), Sb[b]))
        T2 = K.mul(T, K.conj(T))
        lhs_total = K.add(lhs_total, T2)
        if primitive:
            A = K.zero()
            for n, w in a.items():
                A = K.add(A, K.scal(Fr(w), chi_val(j, n)))
            ok_prim &= K.scal(Fr(k), K.mul(A, K.conj(A))) == T2
    rhs = K.zero()
    for b in range(1, k):
        if gcd(b, k) == 1:
            rhs = K.add(rhs, K.mul(Sb[b], K.conj(Sb[b])))
    return {"primitive_identity": ok_prim, "orthogonality": lhs_total == K.scal(Fr(ph), rhs)}
