# lemmaM round 2 (adversarial review): Dixon descent replaces NFS fruitfulness

You are a hostile but fair referee in analytic and algorithmic number theory. Find errors, gaps, overclaims and prior art.
Label each verdict CORRECT / CORRECT WITH FIXES / WRONG / UNVERIFIABLE, and give fixes. Use web search for prior art.
Scratch files go in the current directory only.

## Context
Lee–Venkatesan (LV, arXiv 1805.08873) rigorously produce NFS congruences of squares in L[1/3,(64/9)^{1/3}]. Non-triviality
(fruitfulness) is conditional on their Conj. 7.1 (Blum n). Lenstra–Pomerance (JAMS 1992) rigorous L[1/2, 1+o(1)] is the
record. Earlier rounds (codex/nfsfruit_1, nfsfruit_2, lemmaM_1) showed:
- a uniform Dixon relation r² mod n that is smooth over an L[1/3] factor base costs L[2/3];
- a uniform hidden root λ = r² has small ratio representatives c/d (|c|,|d| ≤ H) with probability ≤ 4H²2^ω/φ(n);
- re-randomizing the root of a known relation is already factoring.

## Claim 1 (Theorem DD, "Dixon descent")
Setting: n odd with ω ≥ 2 distinct prime factors. Relations are pairs (X_i ∈ (Z/n)^×, c_i ∈ G), where G is an abelian group of
"certificates" (e.g. Q^× × K^× for NFS) with a homomorphism ev: G → (Z/n)^× (defined on the relevant subgroup) and X_i² = ev(c_i).
For S with ∏_S c_i = g² and g computable, put X_S = ∏_S X_i and Y = ev(g), so X_S² = Y².
A *descent relation* is produced as follows. Draw r uniform in (Z/n)^× (if gcd(r,n) > 1, output the factor). Put R = r² mod n.
Run a randomized procedure D(R; ξ) (fresh coins ξ, plus coins shared across the run, e.g. the polynomial (m,f), all
independent of r) that either fails or returns (W, c) with W² · R = ev(c), i.e. (rW)² = ev(c). Set X = rW.
**Claim:** if the dependency S and the square root g are computed from (n, all R_i, all coins, all W_i, c_i), and from any
other relations, but not from the r_i beyond R_i, and S contains at least one descent relation, then
Pr[1 < gcd(X_S − Y, n) < n] = 1 − 2^{1−ω} exactly (conditional on everything except the hidden fibre elements). No Blum
condition, no Conj. 7.1, and nothing about the other relations is needed.
*Proof idea:* conditional on R_i and everything else, r_i is uniform on its fibre r_i μ₂(Z/n); X_S/Y is multiplied by ε when
r_{i0} ↦ ε r_{i0} for one i0 ∈ S ∩ D, so X_S/Y is uniform on μ₂(Z/n).

## Claim 2 (descent costs, heuristic smoothness model; the implication from a hypothesis is the claim)
Descent of R: (u,v) := Gauss-reduced shortest vector of {u ≡ R v (mod n)} (|u|,|v| ≤ √(2n)), accept iff u, v are B1-smooth.
For each prime ℓ ∈ (B, B1] dividing uv, a special-q relation: (a,b) in the lattice Λ_ℓ = {a ≡ bm (mod ℓ)} (reduced basis
~√ℓ), with (a − bm)/ℓ and F(a,b) both B-smooth. Then ℓ ≡ φ(a − bα)/c (mod n) with c = (a − bm)/ℓ, and r² ≡ u/v becomes a
relation over the factor base (rational primes ≤ B, degree-1 prime ideals of norm ≤ B, LV's characters).
Use C + O(1) descent relations only (C = number of columns ≈ L[·, a]). Then:
(a) Fixed odd degree d, B = L[1/2, a], B1 = L[1/2, b1] with b1 → ∞ slowly (or B1 = exp((log n)^{3/4}), where HECM detects
B1-smoothness rigorously in L[1/2]^{o(1)}). Special-q values have size n^{θ+o(1)}, θ = 2/(d+1) (the ℓ^{d/2} and box factors are
negligible at this scale). Cost L[1/2, max(a + θ/(2a), 2a)] → **L[1/2, 2/√(d+1)]**, so c < 1 for every odd d ≥ 5
(d = 5: 0.816; d = 7: 0.707).
(b) Growing degree, one-level descent: B = L[2/5], B1 = L[3/5], d ≍ (log n/log log n)^{1/5}. The split costs L[2/5]. Special-q
values have log-size ≍ (log n)^{4/5}, so B-smoothness costs L[2/5]. HECM tests B1-smoothness in L_{B1}[2/3] = L[2/5]. Total **L[2/5, c]**.
Generally, one level of descent gives α ≥ 2/5: the split forces log B1 ≥ (log n)^{1−α}, and the special-q values force
(1 + β1)/2 − α ≤ α, where β1 is the exponent of B1. L[1/3] needs multi-level descent.
(c) QS analogue: special-q x² − n = ℓ c (x ≡ ±√n mod ℓ) gives a heuristic L[1/2,1] QS whose dependencies are provably
fruitful, removing the "H2" (nontrivial-square-root) hypothesis from deterministic-QS analyses, at the price of randomness.

## Claim 3 (what remains for a rigorous theorem: "Hypothesis D" and "Lemma SQ")
Hypothesis D: for the LV-randomized NFS (random m, f = f̂ + (x − m)R), with probability ≥ L^{−o(1)} over (m,f), a uniform
square R descends with probability ≥ L^{−o(1)} in time T_D. Then factoring of every n in time C·T_D·L^{o(1)} + (LV's linear
algebra and square-root costs), rigorously. The core of Hyp. D is a *per-target* special-q statement (Lemma SQ):
Pr_{(m,f), ℓ∼ν}[the special-q region of ℓ contains a doubly-smooth pair] ≥ η, with η^k ≥ L^{−o(1)} (k = number of large
primes per split: k ≍ (log n)^{1/4} for B1 = exp((log n)^{3/4}) in regime (a); k ≍ (log n)^{2/5} in regime (b)). The ℓ arising in
a split depend only on R, so they are independent of (m,f). We claim:
- the first moment E_{m,f,ℓ}[density] is of LV type: LV Lemmas 5.9/5.20/5.21, averaged also over ℓ with weights ∝ 1/ℓ, which is
  comparable to a uniform average over moduli a − bm in dyadic boxes;
- the missing ingredient is anti-concentration (a second moment / pair-correlation upper bound for smoothness of F(a,b),
  F(a',b') under random R and of the rational cofactors under random m), which LV explicitly avoid via stochastic deepening.
  Stochastic deepening cannot replace it, because one field must serve many prescribed targets.
- With two randomized algebraic sides f_j = f̂_j + (x − m)R_j and only completely split large primes descended (density
  (1/d!)^k = L^{−o(1)} for fixed d), every value is randomized. The pair correlations then reduce to (i) a zonotope density bound
  (joint density of (⟨c,v⟩,⟨c,v'⟩) ≤ C_d/(H²‖v∧v'‖) for the Veronese vectors v = v(a,b)) and (ii) Bombieri–Vinogradov-type
  estimates for smooth numbers in progressions mod a − bm over thin sectors of (a,b).
Assess: is Hyp. D genuinely more tractable than Conj. F (first-order smoothness rather than hidden-character equidistribution)?
Is the reduction correct? Are (a) and (b) correct as heuristic costs, and as implications from Hyp. D?

## Claim 4 (barriers)
(B1) Degree-1 two-form sieve, Blum n: if every prime ℓ ≤ B has (ℓ/n) = +1, every congruence is trivial (Φ = Σ⌊Ω_T/2⌋ with T
empty). So "Conj. F for degree-1 sieves for all Blum n" implies that the least prime ℓ with (ℓ/n) = −1 is ≤ B = L[1/2]. That is far
beyond Burgess (n^{1/(4√e)+ε}) unconditionally. For almost all n it is harmless (Erdős–Linnik).
(B2) Unit twist: in the CF/QS relation set of √n, k ↦ k + P (one period) is multiplication by the fundamental unit
ε = a + b√n. It preserves Q_k and multiplies the root by a, a² = (−1)^P. Exact runs: a is a nontrivial √1 in 80/80 Blum
semiprimes, 58/80 mixed (p ≡ 1, q ≡ 3 mod 4) and 31/80 with p ≡ q ≡ 1 (mod 4) (49 of the latter have norm −1), at 28–40 bits.
log₂P has slope 0.511 in bits. So the canonical x-preserving, β-flipping involution exists but displaces relations by ~√n.
By the square-twist identity, any *explicit* β-flipping pair factors n directly. So involutions can enter only
non-constructively (as symmetries of the search law), which is what Claim 1 realises through the coins.

## Tasks
1. Referee Claims 1–4 line by line. In particular: hidden dependence of S or g on r_i; LV's character / square-root machinery
   applied to special-q elements outside their box; whether the split's acceptance or the HECM test leaks fibre information;
   the regime (a)/(b) arithmetic.
2. Prior art (web): is "Dixon relations via descent of random squares in an NFS/QS factor base" known? Compare Stange,
   *Factoring using multiplicative relations modulo n* (arXiv 2211.06821; index calculus mod n with random exponents g^x).
   Also Pomerance's surveys, Schirokauer, Commeine–Semaev, Joux–Lercier, Barbulescu, "rigorous descent". Is a one-level-descent
   L[2/5] heuristic noted anywhere?
3. Is any of the heuristic claims (a) L[1/2, 2/√(d+1)] with proven fruitfulness, or (b) L[2/5], already implied by a
   *rigorous* existing result plus Theorem DD? If Hyp. D can be proved in regime (a) with known tools (LV + Pascadi /
   Harper-type BV + zonotope bounds), sketch how. If not, state precisely the missing lemma.
Verdict first.
