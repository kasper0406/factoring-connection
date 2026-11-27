# coeffam: uniform coefficient families at N = pq (2026-11-27)

**Verdict: no poly(log N) evaluator for any factor-revealing family. Every member that carries a hidden bit
costs N^{1/2-o(1)} or more with known methods, or needs a local datum at p.** L[1/3] not broken.

Pre-registration `experiments/coeffam/prereg.md` (f1dea58). Code and raw output in `experiments/coeffam/`.
Codex gpt-6.2-astra: `codex/coeffam_1` (mod-2 theory, session 01a21e3b-fcc4-7dd3-b98f-d24e79bb996a),
`codex/coeffam_2` (non-Frobenius routes, 01a21e3b-fcc4-7dd0-ae4e-e40d12f42e19),
`codex/coeffam_3` (Hecke-trace open step, 01a21e49-b09e-77e0-a65c-6b99d61377c3).

## Yield ledger

| family | hidden content at pq | bits/query (measured) | best factor-blind evaluator | status |
|---|---|---|---|---|
| Delta^k mod 2, h(k) <= 1 | none (a = 0) | 0 | trivial | proved (NS) + exp |
| Delta^7, Delta^9, Delta^17 (h = 2) | {p,q} mod 8 = {3,5}, {3,3}, {5,5}: QR bits (-1/p), (2/p) | <= 0.63 (H=2 total) | Õ(N^{1/2}): 2.6/5.1/10.6 ms at 36/38/40 bits | exp + proved (codex 1) |
| Delta^11, Delta^13 (h = 3) | quartic character of 2 / p = a^2 + 32b^2 (D4 Frobenius, non-abelian) | 1.52 (H=3 total) | Õ(N^{1/2}) with factoring of N - 8z^2 (never N): 0.11/0.18/0.53 s | exp, 36/36 at 36-40 bits |
| all Delta^k with h(k) <= H, i.e. (H+1)(H+2)/2 forms, k <~ 4H^2 | T_p T_q mod m^{H+1} | 0.63, 1.52, 2.23, 3.44, 4.60, 5.64, 6.96 for H = 2..8 | Õ(N^{1/2}) per member at best | exp |
| binary members x^2 + 2^m y^2 (Delta^{1+2^m}) | class of p vs q in Cl(Z[sqrt(-2^m)]) | ≈ H(Pr[a=1]), Pr[a=1] ≈ (1-2)·cost/sqrt N | enumeration, cost sqrt(N/2^m)/2 | exp (3,000 N at 40 bits) + proved average bound |
| P_N(1) = sum_k a(N,Delta^k) | Tr H_{-7N^2} mod 2 (CM trace, conductor N, Q(sqrt-7)) | 1 (fair coin, not a function of N or (p,q) mod 2^10) | O(N log N) series | identity proved (codex 3), 102/102 |
| h(-4Nm) mod 8 (= r_3 / 12) | Rédei: 4((a+b)s + eb), s = (p/q) mixed symbol | 1 per m | class group, L[1/2] | exp 16/16 cells |
| HW_N = [x^{N-1}] f^{(N-1)/2} | C_p C_q mod N, C_r integer (not a_r) | 0 useful | Õ(N) | proved, 435/435 |

Conditioning on N mod 2^12 instead of 2^8 leaves the H <= 6 entropies unchanged (4.54 vs 4.60 at H = 6), so these are
bits about the factors, not public 2-adic digits of N. Given (p mod 2^12, q mod 2^12), >= 1.07 bits remain at H = 6
(non-abelian part). The trend is about +1 bit per height step up to the sample cap. Uniform recovery of p mod 2^r
from this family needs K = Θ(4^r) (codex 1, proved), so polynomial k gives only O(log n) abelian bits; the rest is
non-abelian 2-tower Frobenius data. In the Chebotarev limit the vector v_H has Θ(H^2) bits (codex 3, proved); our
linear growth is sample-capped. No decoder from v_H, even at H = poly(n), to factor bits is known (open).

## Direction 1: mod-2 forms (kill)
- a(N, Delta^k) = parity of #{N = sum_{i in bits k} 2^i x_i^2, x_i odd}; T_p mod m^4 table in codex 1.
- Minimal non-abelian member: k = 11 (proved by codex 1; found independently by E1c). On N = 3 mod 16:
  a = 1 iff {p,q} = {1|9, 3|11} mod 16 and [2 is a 4th power mod p] XOR [p = 9 mod 16] = 0.
- The lowest hidden layer of every member is a fixed-D quadratic-residuosity bit (Blum-type). Small-D QR oracles are
  not known to factor (codex 1 §4).
- Lacunary cutoff idea, refuted: an evaluator that scans y <= Y sees a nonzero coefficient with average probability
  <= 0.41·Y/sqrt(N) + Y/N (lattice count, proved), and measured Pr[a=1] ≈ (1-2)·Y/sqrt(N). So polylog members carry
  poly(n)·N^{-1/2} bits.
- No automaton on all indices (codex 1, via Christol and Ritchie). The semiprime promise is open.

## Direction 2: non-Frobenius routes (kill; codex 2)
- HW_N = C_N = C_p C_q (mod N) is proved, but the residue mixes a_p with the integer C_q mod p, so it gives nothing.
  The direct cost is Õ(N); BGS does not apply mod N because of division by (N-1)!.
- Dwork truncation: F_{<N}(t) = F_{<p}(t) F_{<q}(t^p) mod p (90/90). The unit-root ratio fails at N
  (u_35 vs u_1 mod 7).
- Res(P, g(x^N)) = Res(P, g(x^q)) mod p. A straight-line program's output mod p depends only on its inputs mod p and
  the program. A Monsky–Washnitzer lift of x -> x^N is impossible: the x^p-coefficient of (1+x)^N is q, not 0, mod p.
- The integer sum over x mod N of (f(x)/N) equals a_p a_q, but costs Õ(N). Codex's S_3-bit test runs in Õ(N^{1/4}) via
  a Pollard–Strassen root sieve, which is factoring-strength.

## Direction 3 (own idea, post-hoc E5, kill pending one open step)
The whole mod-2 level-1 family at N is the single polynomial P_N(X) = J_N(X) mod 2 (Faber polynomial). Its value
at an ordinary j0 is the Hecke correspondence sum over C ⊂ E[N] of j(E/C). At j0 = 1, i.e. y^2 + xy = x^3 + 1 with CM
by Q(sqrt-7) and 2 split, it collapses to Tr H_{-7N^2} mod 2 (proved by codex 3 via isogeny volcanoes; 102/102). One evaluation point gives a
fair hidden bit. A point u in F_{2^e} gives a value in F_{2^e}, up to e bits, but the plain CM-trace form fails
there (P_15 = X^7 + X^15).
- Every known evaluator is at least N: the q-series (N log N), Cl(O_N) of order (p - chi(p))(q - chi(q)), the psi(N)
  correspondence terms, and the exact Rademacher/Salié series, whose output Tr(7N^2) has Θ(N) digits.
- Cl(O_N) = (O_K/N)^x/(Z/N)^x has explicit elements but no explicit j-map.

## Sharpest open steps
1. Decide in poly(log N) time the parity of #{x, y odd > 0 : x^2 + 8y^2 = N}, i.e. [p = q = 3 mod 8]. This is a fixed-D
   Blum-type bit, O(1) per form.
2. Family version: evaluate the Hecke trace J_N(1) mod 2 = Tr H_{-7N^2} mod 2 at one ordinary point without the
   order-N class group. Uniform in the point j0 in F_{2^e}, this would give up to e bits per query. No N^{o(1)}
   method or lower bound is known (codex 3).
3. Decoder: turn v_H (Θ(H^2) non-abelian bits) into ~n/4 leading bits of p (Coppersmith). Even an oracle
   would not help until such a decoder exists.
