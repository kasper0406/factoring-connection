"""Relative-resultant height penalty for frame-less type-D (subfield) modules.  (grdeg, Proposition RR.)

DIAGNOSTIC code. No factoring; no secret inputs. The invariant is computed from the module and N only.

Setting. The frame-less degeneracy is a module over a common subfield k = Q(w), [k:Q] = e, realised as relative
modules on K = k(s), K' = k(s') with RELATIVE minimal polynomials f, g in O_k[x] whose reductions mod N share the
linking root v. Let R = Res_k(f, g) in O_k be the relative resultant. Then (Lemma RES over k):

    * the maps give a surjection psi: O_k -> Z/N, and psi(R) = Res_{Z/N}(psi f, psi g) = 0 (common root v);
    * so R in ker psi, an ideal of index N, hence  N | Norm_{k/Q}(R).

Consequently, by the product formula / Hadamard on the resultant,
    |Norm_{k/Q}(R)| <= exp( n' h1 + n h2 + O(d) ),   h_i = (relative) archimedean generator heights,
so  h1 + h2 >= x/n - O(d) = e*x/d - O(d).  For a rank r = 2e module this is the N-saving penalty ~ r*x/D that
restores alpha >= 1/3 (ledger.py).  If instead R = 0, the two sides descend to the SAME relative form (repeated
norms): products are squares already in k, giving trivial congruences of squares (useless; rankdegen Lemma CC).

SCOPE (Opus critic). RR is rigorous for the HEIGHT/divisibility statement (RR1, RR2) for a relative BINARY module
(relative rank 2, r = 2e) or a relative power basis. The passage from a height bound to a box small-value count (RR3 =>
alpha >= 1/3) is model/heuristic (G-height does NOT imply G-count), and only the relative-binary / power-basis case is
covered: a type-D module over a FIXED small k with relative rank m >= 3 and non-power generators gets only
kappa_gen ~ (r/m)*x/d = c=0 (L[1/4]) from RR, and closing it needs an unproved "relative R1" over k. See the note.

rr_quadratic handles f = x^2 - A, g = x^2 - B with A, B in Q (codex's N=65 and the N-dependent escape), where
Res_k(f,g) = (A-B)^2 in Q and Norm_{k/Q}((A-B)^2) = (A-B)^{2e}. rr_over_Qi handles a GENUINE k=Q(i) example (R not
rational, O_k does real work). psi_resultant_is_zero verifies psi(R)=0 mod N via the integer resultant of the reduced
polynomials (the adjugate identity R = Af + Bg maps under psi and is evaluated at v, so it holds for any relative
polynomials; the reduced resultant being 0 mod N is the computable witness).
"""
from __future__ import annotations

import math
from typing import Dict, List

from sympy import Poly, resultant, symbols


def rr_quadratic(A: int, B: int, e: int, N: int) -> Dict[str, object]:
    """f = x^2 - A, g = x^2 - B over k with A,B in Q, [k:Q]=e. R = (A-B)^2 in Q; Norm_{k/Q}(R) = (A-B)^{2e}.
    Returns R, Norm, and whether N | Norm."""
    R = (A - B) ** 2
    norm = R ** e  # Norm_{k/Q}(rational) = rational^e
    return {
        "R_relative_resultant": R,
        "Norm_k/Q(R)": norm,
        "N_divides_Norm": (norm % N == 0),
        "Norm/N": norm // N if norm % N == 0 else None,
    }


def psi_resultant_is_zero(f_coeffs: List[int], g_coeffs: List[int], v: int, N: int) -> Dict[str, object]:
    """Verify psi(R) = 0 mod N directly: reduce f,g mod N, check they share root v, and that Res mod N == 0.
    f_coeffs, g_coeffs are integer coefficient lists low->high (the reduced relative polys, coeffs already in Z via
    the linking)."""
    x = symbols("x")
    def ev(c, s):
        return sum((ci % N) * pow(s, i, N) for i, ci in enumerate(c)) % N
    fp = Poly(list(reversed(f_coeffs)), x)
    gp = Poly(list(reversed(g_coeffs)), x)
    res = int(resultant(fp, gp))
    return {
        "f(v) mod N": ev(f_coeffs, v),
        "g(v) mod N": ev(g_coeffs, v),
        "Res(f,g) mod N": res % N,
        "psi(R)=0": (res % N == 0),
    }


def rr_over_Qi(f_gauss, g_gauss, psi_i: int, v: int, N: int) -> Dict[str, object]:
    """Genuine O_k example with k=Q(i) (O_k does real work: R is not rational). f_gauss, g_gauss are sympy Polys over
    Z[i]; psi: i -> psi_i mod N; v the shared linking root mod N. Verifies ψ(R)=0 and N | Norm_{Q(i)/Q}(R)."""
    import sympy as sp
    R = sp.expand(sp.resultant(f_gauss, g_gauss))
    a, b = int(sp.re(R)), int(sp.im(R))
    norm = a * a + b * b
    # reduce f,g mod (N, i - psi_i) and check common root v
    def red_coeffs(poly):
        cs = poly.all_coeffs()[::-1]  # low->high
        return [(int(sp.re(sp.expand(c))) + int(sp.im(sp.expand(c))) * psi_i) % N for c in cs]
    def ev(c):
        return sum(ci * pow(v, i, N) for i, ci in enumerate(c)) % N
    fr, gr = red_coeffs(f_gauss), red_coeffs(g_gauss)
    return {"R": (a, b), "Norm_k/Q(R)": norm, "N_divides_Norm": norm % N == 0,
            "psi_f(v)": ev(fr), "psi_g(v)": ev(gr), "psi(R)=0_consistent": (ev(fr) == 0 and ev(gr) == 0)}


def height_penalty_report(name: str, A: int, B: int, e: int, N: int) -> Dict[str, object]:
    out = rr_quadratic(A, B, e, N)
    x = math.log(N)
    # generator s has s^2 = A so height log|A|/2; the N-saving exponent realised:
    out.update({
        "name": name, "e": e, "N_bits": int(N).bit_length(),
        "relative_gen_height_log|A|": math.log(abs(A)) if A else 0.0,
        "note": "N | Norm(R) forces the generator height; small generators => R~0 => repeated norm (useless)",
    })
    return out


if __name__ == "__main__":
    # codex's N=65 Weil-restriction example: k=Q(i), K=Q(i,sqrt14), K'=Q(i,sqrt79), f=x^2-14, g=x^2-79, e=2.
    print("codex N=65 example (k=Q(i), f=x^2-14, g=x^2-79):")
    print("  RR:", rr_quadratic(14, 79, 2, 65))
    print("  psi(R)=0 check (linking root v=53 mod 65):", psi_resultant_is_zero([-14, 0, 1], [-79, 0, 1], 53, 65))

    # codex's 'escape' family: k=Q(t) t^2=-N, s^2=1-N, s'^2=1-2N, e=2.
    N = 101 * 103
    print(f"\ncodex escape family (N={N}, k=Q(sqrt-N), f=x^2-(1-N), g=x^2-(1-2N)):")
    print("  RR:", rr_quadratic(1 - N, 1 - 2 * N, 2, N))
    print("  => R = N^2, Norm = N^4, N | Norm trivially; generators s~sqrt(N): height x/2, norms ~ N^4.")
