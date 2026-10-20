"""Pre-registered factor classification (notes/slpcirc_prereg.md): LIN, CYC, KUM, CHEB, LOW, OTHER."""
from __future__ import annotations
from fractions import Fraction
import flint

TYPES = ("LIN", "CYC", "KUM", "CHEB", "LOW", "OTHER")


def _compose_affine(g: flint.fmpz_poly, s: Fraction, c: int) -> flint.fmpz_poly:
    """Primitive integer polynomial proportional to g(s*x + c)."""
    coeffs = [int(a) for a in g.coeffs()]
    d = len(coeffs) - 1
    # g(s x + c) * den^d where s = num/den
    num, den = s.numerator, s.denominator
    lin = flint.fmpz_poly([c * den, num])  # den*(s x + c)
    res = flint.fmpz_poly(0)
    powl = flint.fmpz_poly(1)
    for k in range(d + 1):
        res += coeffs[k] * powl * (den ** (d - k))
        powl *= lin
    if res == 0:
        return res
    ct = res.content()
    return res / ct if ct != 1 else res


def _is_cyc(h: flint.fmpz_poly) -> bool:
    if h.degree() < 1:
        return False
    if h.leading_coefficient() < 0:
        h = -h
    return h.is_cyclotomic() != 0


def is_cyc(g):
    for s in (Fraction(1), Fraction(-1), Fraction(2), Fraction(-2), Fraction(1, 2), Fraction(-1, 2)):
        for c in range(-3, 4):
            if _is_cyc(_compose_affine(g, s, c)):
                return True
    return False


def is_kum(g, kmax=128):
    for c in range(-3, 4):
        h = _compose_affine(g, Fraction(1), c)
        d = h.degree()
        if d < 1:
            continue
        # work over Q: x^k mod h via fmpq_poly
        hq = flint.fmpq_poly([int(a) for a in h.coeffs()])
        xq = flint.fmpq_poly([0, 1])
        r = flint.fmpq_poly(1)
        for k in range(1, kmax + 1):
            r = (r * xq) % hq
            if r.degree() <= 0 and r != 0:
                return True
    return False


def is_cheb(g):
    coeffs = [int(a) for a in g.coeffs()]
    d = len(coeffs) - 1
    u = flint.fmpz_poly([0, 1])
    for a in (1, -1):
        for kap2 in (1, 2):  # kappa = 1/kap2
            for c in range(-2, 3):
                # x = (u^2 + kap2*c*u + a)/(kap2*u); H = sum g_k (u^2+kap2 c u + a)^k (kap2 u)^(d-k)
                q = flint.fmpz_poly([a, kap2 * c, 1])
                ku = flint.fmpz_poly([0, kap2])
                H = flint.fmpz_poly(0)
                qp = flint.fmpz_poly(1)
                for k in range(d + 1):
                    H += coeffs[k] * qp * ku ** (d - k)
                    qp *= q
                if H == 0 or H.degree() < 1:
                    continue
                _, facs = H.factor()
                ok = True
                for f, _e in facs:
                    if f.degree() == 1 and f.coeffs()[0] == 0:  # factor u (u = 0) cannot occur unless a=0
                        ok = False
                        break
                    if not _is_cyc(f):
                        ok = False
                        break
                if ok:
                    return True
    return False


def classify_factor(g: flint.fmpz_poly) -> str:
    d = g.degree()
    if d == 1:
        return "LIN"
    if is_cyc(g):
        return "CYC"
    if is_kum(g):
        return "KUM"
    if is_cheb(g):
        return "CHEB"
    if d <= 4:
        return "LOW"
    return "OTHER"
