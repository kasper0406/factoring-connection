# rkrank: Hypothesis RK (GF(2) rank of NFS relations) and rigorous small-square samplers

Agent `rkrank`, phase 8, 2026-10-19.
- Pre-registration: `notes/rkrank_prereg.md`, committed in `7a2d0fe` before the E1 run.
- Code: `src/rkrank/`. Tests: `tests/test_rkrank.py`. Raw data: `experiments/raw/rkrank_e1.jsonl`.
- GPT-6.1-sol rounds: `codex/rkrank_1.{md,out}` (blind, lever A), `codex/rkrank_2.{md,out}` (lever B) and
  `codex/rkrank_3.{md,out}` (adversarial review of this note).
- Reviews: two independent read-only Opus critics, one on proofs and one on literature and novelty (§7).

## Verdict

**Neither lever succeeded. We did not prove Hypothesis RK. We found no rigorous small-square sampler
with θ < 2/3, either in the literature or by our own attempts. Lenstra–Pomerance's L[1/2, 1+o(1)]
(1992) stands, and no published rigorous c < 1 was found.** Design V remains conditional at c = √(2/3).

What is new is a sharper map of the problem: what RK reduces to, which weaker form already suffices,
which obstructions are forced, and why the natural proof templates fail. The toy experiments support a
very strong form of RK.

**Lever A (RK).** Labels in brackets.
- **Theorem V_κ [proved implication].** RK may be weakened to codimension D^κ, where D is the number of
  factor-base columns. With a sampler exponent θ and LV rows of size exponent η = 3/(d+1), the cost
  constant is c = √(2·max{η, θ/(2−κ)}).
  - For Vallée (θ = 2/3) and d ≥ 9, c < 1 holds iff **κ < 2/3**.
  - So deficiency D^{2/3−ε} already beats Lenstra–Pomerance, conditionally on the other Design V
    inputs. RK as originally stated is κ = 0.
- **Theorem U [proved].** Uniform rank growth implies RK. Suppose every surviving dual space V (meeting
  the forced space trivially) is escaped by a random relation with probability ≥ dim V/(D·L^{o(1)}).
  Then D·L^{o(1)} relations reach deficiency L^{o(1)}.
- **Prop G′ [proved, existential; GPT-6.1-sol].** Raw matrix membership and genuine certificates differ by
  at most h = dim H_K = O(log N) dimensions, where H_K is the "square ideal" Selmer group. So r + h + 1
  sampler rows force a *genuine* dependency with a nonzero sampler projection. Design V's
  "genuine-certificate" subtlety is therefore harmless existentially; LV §6 is still needed for efficient
  extraction.
- **Prop C [proved; elementary plus Weil].** This classifies the *algebraic* obstructions: quadratic
  character identities χ(a − bm)·ψ(a − bα) ≡ 1 that hold on all coprime pairs.
  - Away from primes ℓ ≤ C_d and ℓ | disc(f)·f_d, they are exactly the ω characters at p | N.
  - Archimedean and |S|-parity obstructions come only from the box geometry, and number at most r₁ + 1.
  - Correction to hypD §5 (found by calibration): a real root of f "near m", hence an archimedean sign
    obstruction, exists only when the R-height H exceeds m. It is not automatic.
  - Every other obstruction is an accident of the smooth set X.
- **Lemma S [proved; GPT-6.1-sol blind].** For uniform B-smooth integers, a switching argument gives
  exactly the escape bound that Theorem U needs. This is the rational-side model of RK. It fails for the
  joint law, because multiplying a − bm by p does not map relations to relations.
- **Barriers [proved].**
  - (B-F) Averaging over fields cannot prove RK via per-functional bounds. There are at most N^{O(1)}
    fields against 2^D functionals, and an explicit spread construction (GPT-6.1-sol) has every fixed
    functional failing in at most one field while every field has deficiency D/2. Any proof must be a
    per-field structural statement about the smooth set X_{m,f}.
  - (B-M) Fixed low moments (pair correlations) of the relation law cannot give near-full span. One needs
    moments of order D/L^{o(1)}.
  - (B-C) Coverage does not imply span, even with arbitrarily many occurrences per prime.
- **Experiment E1 [experimentally supported; pre-registered; toy cubic NFS at 64/80 bits].**
  - In all 72 records, W ⊆ Λ/2, and W contains all of Λ/2 on primes occurring ≥ 3 times (c₃ = 0).
  - Small-square vectors on that support lie in W: 1141/1141.
  - Dense and Liouville-type functionals have δ ∈ [0.45, 0.55].
  - Two sharper predictions failed narrowly (P1: peeling explains every rare-prime obstruction; P2b: c₁ = 0
    at ≥ 2·cols). The cause was ≤ 5 rare-prime obstructions, never dense or fake ones.
  - 16 of 20 instances exhausted their registered box before 3·cols (disclosed).

**Lever B (samplers).**
- **Prop Q [proved implications].** Take a fibre-symmetric sampler of exponent θ, with per-call cost
  L^{s+o(1)}.
  - Dixon-style: c = s + √(s² + 2θ), which is < 1 iff θ + s < 1/2. No rank hypothesis is needed.
  - Design V, under RK and its other inputs: c = (s + √(s² + 4θ))/2, which is < 1 iff θ + s < 1.
  - So "θ_V < 2/3 gives √θ_V" is correct (s = 0, d + 1 ≥ 6/θ). Note also that θ < 1/2 alone would beat
    LP unconditionally.
- **Theorem L (linearisation barrier) [proved, for a stated model].** Take any affine cell
  x₁ + Σ t_i w_i whose quadratic part is uniformly small as an integer form. Such a cell collapses to a
  rank-one cell x₁ + vτ with v² ≡ σ small, so it holds ≤ 2√h + 1 distinct x.
  - If x₁ is chosen obliviously (density ≤ λ/N), the expected yield per cell is ≤ λ·2^ω·(2√h+1)(2h+1)/N.
  - For θ < 2/3 this forces N^{Ω(1)} cells per output.
  - The natural Coppersmith shift lattice works only when T²h < N. That is incompatible with the needed
    Th ≥ N^{1−o(1)}; derived independently by us and GPT-6.1-sol.
- **Literature.**
  - The post-1991 progress on Vallée's sampler is in *uniformity*, not exponent: Coron (CRYPTO 2002)
    widened the interval with statistical closeness, and Gentry (CRYPTO 2004) gave exactly uniform sampling
    at width 8N^{2/3}. With Gentry, κ = 1 in Design V.
  - We found no θ < 2/3, and no rigorous c < 1, through 2026-10-19.

---

## 1. Setting and notation

These follow hypD §5.
- A field is A = (m, f), with f = f̂_{N,m} + (x − m)R, f(m) = N, fixed odd d and B = L[1/2, a].
- A relation is a coprime (a, b) in a box of size L^{O(1)} with a − bm and F(a, b) both B-smooth.
- Its row consists of:
  - rational parities q(a, b) ∈ F₂^{D_rat}: the sign and v_p(a − bm), p ≤ B;
  - the "high" block z(a, b): degree-one ideal parities, real-embedding signs, LV characters, and the
    |S|-parity column when f is non-monic.
- R is the row space and W = {u : (u, 0) ∈ R}.
- Λ/2 = {u : ∏p^{u_p} ∈ (Z/N)^{×2}}. Gen = (Λ/2)^⊥, spanned by the ω functionals w_i(p) = [(p/p_i) = −1].
- U is the target space: Λ/2 restricted to a set of well-covered primes, or the span of the sampler
  supports.
- "Deficiency" means codim_U(W ∩ U).

## 2. Lever A: Hypothesis RK

### 2.1 Duality [proved; standard, stated in codex rkrank_1]
- W^⊥ = proj_rat(R^⊥).
- codim_U(W ∩ U) = dim res_U(W^⊥).
- Proof: a rational functional w kills W iff it extends to a functional (w, y) on R. That is, w∘q lies in
  the span of the coordinate functions of z on the relation set.

So a **dual obstruction** is a pair (w, y) with w(q(a, b)) = y(z(a, b)) for every relation. RK bounds the
dimension of the rational projections of the non-forced obstructions.

Linear algebra also gives dim W = rank R − rank(z-block). Hence
codim W = D_rat − rank R + rank(z-block). This identity is asserted in the code (`harness.level_record`).

### 2.2 Raw versus genuine certificates [proved, existential; codex rkrank_1, checked by us]
- Let E₀ be the combinations S whose algebraic parities and chosen character columns vanish.
- Then ∏_S(a − bα) has even valuation at every prime, so it defines a class
  ρ(S) ∈ H_K := {x ∈ K^× : all v_𝔭(x) even}/K^{×2}. H_K sits in the exact sequence
  0 → O_K^×/O_K^{×2} → H_K → Cl(K)[2] → 0.
- So h := dim H_K ≤ r₁ + r₂ + log₂ h_K = O_d(log N), since |disc f| ≤ N^{O(d)}.
- Order and index bookkeeping adds O(log N) more.
- Genuine certificates form ker ρ, which has codim ≤ h in E₀, so dim W_raw − dim W_true ≤ h.
- **Consequence.** If the raw deficiency on U is ≤ r, then any t > r + h sampler rows supported in U give a
  genuine dependency with nonzero sampler projection. That dependency is fruitful with probability
  ≥ 1/(2κ) (hypD §5).
- Efficient extraction of the genuine combination still needs LV §6 (characters). This is an
  algorithmic, not an existential, requirement.

### 2.3 Theorem U: uniform rank growth implies RK [proved]
**Setup.** Z₁, Z₂, … are i.i.d. with law μ on F₂^s. G is the space of functionals vanishing μ-a.s. (the
forced space). Δ_M = span(Z₁, …, Z_M)^⊥.

**Hypothesis (Exp_T).** For every subspace V of functionals with V ∩ G = 0 and dim V = k > r,
Pr_μ[Z ∉ V^⊥] ≥ k/T.

**Theorem U.** Let τ_r = min{M : dim Δ_M − dim G ≤ r}. Then E τ_r ≤ T·Σ_{k=r+1}^{s} 1/k ≤ T(1 + ln s).
- Hence Pr[dim Δ_M − dim G > r] ≤ T(1 + ln s)/M.
- And codim_U(W ∩ U) ≤ dim res_U(proj_rat G) + r.

*Proof.*
1. G ⊆ Δ_M always. Write Δ_M = G ⊕ V_M with dim V_M = k_M.
2. If Z_{M+1} ∉ V_M^⊥ then Z_{M+1} ∉ Δ_M^⊥, so dim Δ drops by exactly one.
3. Given the past, this happens with probability ≥ k_M/T. So the time spent at level k is stochastically
   dominated by a geometric variable of mean T/k.
4. Sum over k, and apply Markov's inequality.
5. Finally, W^⊥ = proj_rat(Δ). ∎

*Remarks.*
- For one field, the relation law is the kept law on the finite smooth set X = X_{m,f} (nfsfruit Prop R).
  So (Exp_T) is a **deterministic expansion property of X**: for every V with V ∩ G = 0,
  #{x ∈ X : x ∉ V^⊥} ≥ |X|·dim V/T.
- The natural scale is T = D·L^{o(1)}. A prime p ≈ B occurs in a relation with probability ≈ u/B, so
  rare-prime spaces have escape rate ≈ k·u/B. That is the coupon-collector regime, and it is why the
  rate must be allowed to scale with k rather than be a constant p₀.
- (Exp_T) is the GF(2) analogue of de Boer–Ducas–Fehr's "p-evenly distributed" sampling condition, used
  by de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588) for lattice generation under ERH. Their condition
  asks for a *constant* escape probability from every proper sublattice. NFS rows cannot satisfy that,
  because of rare primes. The k-dependent rate is what NFS rows can at best satisfy.

### 2.4 Fourier form and the moment barrier [proved; codex rkrank_1, checked]
- E 2^{s − rank R_M} = Σ_w ((1 + μ̂(w))/2)^M, where μ̂(w) = E(−1)^{w·Z}.
- Also Pr[S_k = S′_k] = 2^{−s}Σ_w μ̂(w)^{2k}, where S_k, S′_k are independent sums of k rows.
- A single row has support size J ≤ L^{O(1)} (the box size), so Pr[S_k = S′_k] ≥ J^{−k}.
- Certifying a collision probability near 2^{−s+r} therefore needs k ≥ (s − r)log 2/log J = D/L^{o(1)}.
- **(B-M):** pair-correlation or other fixed-order moment information about relations cannot by itself
  certify near-full span. This is the RK analogue of hypD's Prop P.

### 2.5 Prop C: the forced obstructions [proved]
Call a pair (ψ₁, ψ₂) an **algebraic identity** if:
- ψ₁ is a quadratic character on integers coprime to some q₁ (a Dirichlet character, of any conductor);
- ψ₂ is a quadratic character of (O_K/𝔮)^× for some ideal 𝔮;
- ψ₁(a − bm)ψ₂(a − bα) = 1 for **all** coprime (a, b) with both arguments coprime to the moduli.

By Kummer theory and Chebotarev, every dual functional on the smooth groups is the restriction of such a
pair, of enormous conductor. So "algebraic identity" is a property of the realisation, not of the
functional. Prop C says which realisations hold identically.

**Prop C.** Let C_d = 4(d+1)². In an algebraic identity, the local component at every prime
ℓ ∤ 2·f_d·disc(f) with ℓ > C_d is trivial unless ℓ | N. At ℓ | N it is a power of the pair
(χ_ℓ, χ_{𝔓_ℓ}), 𝔓_ℓ = (ℓ, α − m).
- Consequently the identities span a space of dimension ≤ ω + O_d(#{ℓ ≤ C_d} + ω(f_d·disc f)) =
  ω + O_d(log N).
- There are no archimedean identities on all of Z²: as (a, b) ranges over Z², a/b crosses m and every
  real root.

*Proof.*
1. **CRT.** The moduli factor over primes, and (a, b) mod ∏ℓ^{e_ℓ} ranges over a product of local sets.
   Fixing all other local coordinates at (1, 0) (value 1) shows each local component is identically 1.
2. **Local reduction.** For odd ℓ ∤ f_d·disc(f), quadratic characters of (Z/ℓ^e)^× and (O_K/𝔩^e)^×
   factor through residue fields. The local component is (a, b) ↦ χ_ℓ(G(a, b)) with
   G = (a − bm)^{ε₀}·∏_{𝔩|ℓ} N_𝔩(a − bα)^{ε_𝔩}. This uses (x/𝔩) = χ_ℓ(N_{k(𝔩)/F_ℓ}(x)). G is a binary form
   of degree ≤ d + 1 over F_ℓ.
3. **The irreducible factors.**
   - N_𝔩(a − bα) is the homogenised irreducible factor of f mod ℓ belonging to 𝔩.
   - These are pairwise coprime, since ℓ ∤ disc(f).
   - They are coprime to a − bm unless some 𝔩 has degree one with root ≡ m (mod ℓ), i.e. ℓ | f(m) = N.
4. **ℓ ∤ N.** G is then a product of distinct irreducible forms, not a constant times a square, unless all
   exponents vanish. By Weil's bound, |Σ_{(a:b)∈P¹(F_ℓ)} χ_ℓ(G(a, b))| ≤ (deg G − 1)√ℓ + deg G. For
   ℓ > C_d this sum is < ℓ + 1 − deg G. Hence χ_ℓ(G) = −1 at some point, contradicting the identity.
5. **ℓ | N.** The same argument applied to G/(a − bm)² leaves only ε₀ = ε_{𝔓_ℓ} with all other ε = 0.
6. The excluded primes (ℓ ≤ C_d, ℓ | 2f_d·disc f) contribute at most d + 2 local dimensions each. ∎

**Box-geometry obstructions [proved].**
- If a real embedding σ_j has sign(a − bσ_j(α)) constant on the box, its column is constant.
- If, in addition, sign(a − bm) is constant (it is −1 when |a| < m and b ≥ 1), then e_sign plus that
  column, or plus the |S|-parity column for non-monic f, lies in R^⊥. In that case W ⊆ {u_sign = 0}.
- There are at most r₁ + 1 such obstructions.
- They are harmless to Design V: use positive targets, or pair negative ones.
- **Correction to hypD §5 ("archimedean root near m").** That root exists only if H ≳ m. Diagnostically,
  for degree-9 LV-style fields with R-height m^η, a real root within m/4 of m occurred in 0/20 trials for
  η ≤ 1 and 20/20 for η ≥ 1.5. This is consistent with the root displacement N/f′(m) ≈ m²/H.
- In the monic toy no such obstruction arises unless a real-sign column is constant 1 (P4 in §5).

**Reading.** Every non-forced obstruction is an *accident of the smooth set X*. Its realisations hold on a
positive density of all pairs, yet fail on none of the density-L^{−c} smooth pairs. RK asserts that such
accidents span only L^{o(1)} (or, by V_κ, ≤ D^{2/3−ε}) dimensions.

### 2.6 Lemma S: the rational model satisfies uniform rank growth [proved; codex rkrank_1 blind, checked]
**Lemma S.** Let n be uniform among the B-smooth integers in [1, x], and let V be a k-dimensional space of
prime-parity functionals. Then
  Pr[n ∉ V^⊥] ≥ k·Ψ(x/B, B)/((1 + log₂x)·Ψ(x, B)).
For x = N^{c+o(1)} and B = L[1/2, a], this is ≥ k/(D·L^{o(1)}).

*Proof.*
1. Choose primes p₁, …, p_k ≤ B whose evaluation vectors ev(p_i) ∈ F₂^k form a basis. This is possible
   because V is k-dimensional on the prime coordinates.
2. For B-smooth n ≤ x/B, the k + 1 smooth integers n, p₁n, …, p_kn ≤ x have distinct evaluation vectors
   ev(n), ev(n) + ev(p_i). At most one of them is 0.
3. Each integer y ≤ x arises from at most 1 + ω(y) ≤ 1 + log₂x pairs (n, i).
4. So #{y : ev(y) ≠ 0} ≥ k·Ψ(x/B, B)/(1 + log₂x).
5. By Canfield–Erdős–Pomerance, Ψ(x/B, B)/Ψ(x, B) = B^{−1}·exp(o(u log u)) = B^{−1}L^{−o(1)}, and
   D ≍ B/log B. ∎

**Where it breaks for NFS [assessment, agreed with codex].** The switching n ↦ p_i·n is not an operation
on relations of one field. There is no (a′, b′) in the box with a′ − b′m = p_i(a − bm), and F(a′, b′) would
change uncontrollably; F(a, b) ≡ Nb^d (mod a − bm) couples the two sides. A *bounded-congestion switching
on the joint kept law* would prove (Exp_T) via Theorem U. We know no such switching. Varying m changes the
field, and all rows of one dependency must share the field.

### 2.7 Barriers for proofs of RK [proved]
- **(B-F) Field averaging is useless for union bounds.**
  - There are ≤ N^{O(1)} fields, so for a fixed rational functional w, Pr_A[w ∈ W_A^⊥] is 0 or
    ≥ N^{−O(1)}. A union bound over the 2^{D_rat} functionals (D_rat = L^{a}) therefore needs all but
    2^r·N^{O(1)} of them to be obstructions for *no* field. That is a pointwise statement.
  - Abstract example (codex rkrank_1): take a spread {E_A} of pairwise trivially intersecting
    D/2-dimensional subspaces of F₂^D and draw field-A rows uniformly from E_A^⊥. Every fixed nonzero
    functional is an obstruction in at most one field, yet every field has deficiency D/2.
  - So RK must be proved per field, as an expansion property of X_{m,f} (§2.3).
  - The same holds for "many independent fields": W = Σ_i W_i gives (ΣW_i)^⊥ = ∩W_i^⊥. With L^{o(1)}
    fields this helps only under a uniform multiplicative contraction of the surviving dual per field.
    Fixed-functional survival bounds would need Θ(D) fields.
- **(B-M)** is in §2.4.
- **(B-C) Coverage does not imply span.**
  - Rows (e_i, e_i), with rational prime i paired with algebraic ideal i, cover every coordinate any
    number of times, yet W = 0 (codex).
  - So "well-covered support" in RK must mean a span statement, not an occurrence count.
  - E1's P2 tests exactly this for the toy.
- **(B-L) Random-exponent generation is not available at NFS size [assessment].**
  - Hafner–McCurley, de Boer–Pellet-Mary–Wesolowski and Dixon-with-random-exponents prove generation
    because, given the observable, the relation vector is (quasi-)uniform on a coset of the relation
    lattice.
  - A single NFS row has support ≤ L^{O(1)}, far below 2^{Θ(D)} = |Λ/2|.
  - We checked the obvious substitute, the twisted sampler c_e·x² mod N small with c_e a product of k
    factor-base primes. Its relation vectors e + vec(z) are quasi-uniform on (weight-k vectors) ∩ coset.
    They provably satisfy (Exp_T) with T = D·κ·2^ω, by the argument of Lemma S with p uniform in its square
    class.
  - But every such row costs a full sampler smoothness test. Spanning with them is Vallée's own algorithm,
    c = √(2θ).

### 2.8 Theorem V_κ: weaker RK suffices [proved implication; formula confirmed by codex rkrank_1]
**Assumptions.** Design V's inputs: LV supply at fixed odd d, LV §6, and a fibre-symmetric sampler of
exponent θ with cost L^{o(1)} per call. RK is weakened to: raw deficiency on U ≤ D^κ = L^{κa+o(1)}, with
probability ≥ L^{−o(1)} over the field.

**Algorithm.** Collect D^κ + h + 1 sampler rows; by §2.2 this suffices. With η = 3/(d+1), the costs are:
- linear algebra: L^{2a};
- LV rows: L^{a+η/(2a)};
- sampler rows: L^{κa+θ/(2a)}.

**Optimum.**
- c(κ) = √(2·max{η, θ/(2−κ)}), at a = √(θ/(2(2−κ))) when the LV term does not bind.
- κ = 0 recovers √θ (Theorem V). κ = 1 recovers Vallée/Dixon's √(2θ) with no rank hypothesis, because
  D + 1 rows are always dependent.
- **c < 1 iff κ < 2 − 2θ**, given η ≤ θ/(2−κ). For Vallée, θ = 2/3, this means κ < 2/3, with d ≥ 9 for
  every κ ∈ [0, 1].

## 3. Lever B: small-square samplers

### 3.1 Prop Q: consequences of a better sampler [proved implications]
**Definition.** A sampler of exponent θ:
- outputs x with |x² mod N| ≤ N^{θ+o(1)};
- has Pr[εx] ≥ Pr[x]/κ for all ε ∈ μ₂(Z/N), with κ = L^{o(1)};
- outputs z = x² mod N that is B-smooth with probability ≥ ρ(log h/log B)·L^{−o(1)}.

**Inputs.**
- Fruitfulness ≥ 1/(2κ) per dependency containing a sampler row (hypD §5, pointwise bound).
- Smoothness counting in square classes. Take products of two smooth numbers ≤ √h of the same class, with
  Ω ≤ K = o(log L) so that representation multiplicities are 2^{2K} = L^{o(1)} (codex rkrank_2). The loss
  is only L^{o(1)}.
- Preprocessing by trial division to exp(√log N), which gives 2^ω = L^{o(1)}.
- Smoothness of batches by Bernstein's batch algorithm, at L^{o(1)} amortised cost (codex rkrank_2). So no
  HECM assumption is needed.

**Costs with per-call cost L^{s}.**
- Dixon-style: max(2a, a + s + θ/(2a)), so **c_D = s + √(s² + 2θ)**, which is < 1 iff **θ + s < 1/2**.
  - Codex rkrank_2 wrote s + √(2θ). That is an upper bound only; corrected here.
- Design V, conditional: max(2a, a + η/(2a), s + θ/(2a)), so **c_V = (s + √(s² + 4θ))/2**, which is < 1 iff
  **θ + s < 1**. This needs η ≤ c_V²/2, i.e. d + 1 ≥ 6/(c_V²) at s = 0.
- At s = 0 these give **√(2θ)** and **√θ**. So the hypD claim "θ_V < 2/3 gives √θ_V" is correct, under
  RK, the supply hypothesis and d + 1 ≥ 6/θ.

### 3.2 Why 2/3, and Theorem L (linearisation barrier) [proved for the stated model]
**Mechanism (Vallée; reconstructed by us and codex rkrank_2).**
- In the cell x = x₁ + t, (x₁ + t)² ≡ x₁² + 2x₁t + t².
- Absorbing t² ≤ h into the error leaves a 2-dimensional affine lattice of determinant N in a T × h box
  with T ≤ √h.
- A cell is productive only if Th ≳ N, i.e. h ≳ N^{2/3}.
- Vallée's Farey dissection handles the skewed cells. The lattice point counts are uniform exactly down
  to this threshold (Vallée §5.1).

**Lemma (rank collapse) [proved].** Setup:
- w ∈ (Z/N)^r, and the box 𝔅 = ∏[−T_i, T_i] with T_i ≥ 1;
- s_ij = w_iw_j mod N, centred;
- q(t) = Σ s_ij t_it_j is the integer quadratic form;
- hypothesis: |q(t)| ≤ h on 𝔅 and 3h² < N.

Then, unless a gcd with N reveals a factor, there are v ∈ Z/N, σ ∈ Z \ {0} and a primitive u ∈ Z^r with
w ≡ v·u, v² ≡ σ and q(t) = σ(u·t)². The distinct points x₁ + w·t form a subset of
{x₁ + vτ : |τ| ≤ τ_max}, with |σ|τ_max² ≤ h. Hence there are at most 2√h + 1 of them.

*Proof.*
1. q(e_i) = s_ii and q(e_i + e_j) − s_ii − s_jj = 2s_ij, so |s_ij| ≤ 3h/2.
2. The 2×2 minors s_iis_jj − s_ij² are ≡ 0 (mod N), because S ≡ ww^T. In absolute value they are < N, so
   they vanish. Thus S = σuu^T over Q, with σ ∈ Z because u is primitive.
3. Pick i with u_i, w_i units (otherwise gcd gives a factor, or w_i ≡ 0). Put v = w_i/u_i. Then
   w_j ≡ σu_iu_j/w_i ≡ vu_j.
4. τ_max = Σ|u_i|T_i is attained at t_i = sgn(u_i)T_i. ∎

**Theorem L (oblivious-cell barrier) [proved].** Setup:
- a sampler draws a cell centre x₁ from a law with density ≤ λ/N on Z/N, independent of N's arithmetic
  beyond that bound;
- its cell is an affine family with uniformly small integer quadratic part, as in the Lemma.

Then the expected number of x with |x² mod N| ≤ h in the cell is ≤ λ·2^ω(2h+1)(2√h+1)/N. For h = N^θ,
θ < 2/3, a sampler with L^{o(1)} distortion therefore needs N^{2/3−θ−o(1)} cells per output, whatever its
per-cell enumeration cost.

*Proof.*
- The set {x : |x² mod N| ≤ h} has ≤ 2^ω(2h + 1) elements, because each residue has ≤ 2^ω square roots
  for squarefree N.
- Combine this with the Lemma and linearity of expectation over x₁. ∎
- Codex rkrank_2 gives the special case of interval cells.

*Scope.* The barrier does **not** cover:
- arithmetic (non-oblivious) choice of productive cells;
- curved regions solved exactly;
- arbitrary lattices.
No unconditional lower bound on θ can exist without a factoring lower bound: polynomial-time factoring
gives θ = 0.

**Coppersmith [proved failure of the natural construction].**
- For F(t, z) = (x₁ + t)² − z mod N with |t| ≤ T and |z| ≤ h, take the triangular shift lattice
  N^{m−j}t^iF^j, 0 ≤ i ≤ 2(m − j).
- It has dimension (m+1)² and det = (NT)^A h^C, with A = m(m+1)(4m+5)/6 and C = m(m+1)(2m+1)/6.
- The short-vector condition tends to **T²h < N** (our Jochemsz–May count and codex's exact count agree).
- Productive cells need Th ≥ N^{1−o(1)}, so the two never overlap.
- With z known, univariate Coppersmith applies, but z ranges over N^θ values.

### 3.3 Other routes [assessment unless marked]
- **Ratio forms** x² ≡ u/v with |u|, |v| ≤ N^{θ/2}: these are equivalent to (vx)² ≡ uv, so there is no gain.
  [proved]
- **Class groups (LP, Seysen, Hafner–McCurley).**
  - They behave like θ = 1/2 objects: reduced forms have leading coefficient ≤ √(|Δ|/3). That explains
    c = 1.
  - [proved, codex rkrank_2] Coverage obstruction: ideals of norm ≤ |Δ|^θ cover only O(A log A) of the
    |Δ|^{1/2+o(1)} classes. So uniform class sampling cannot have θ < 1/2.
- **Elkies (ANTS 2000)** gives points near curves via a compact cover. Tangent subdivision reproduces the
  curvature threshold, so it supplies no sampler.
- **Schnorr's SVP factoring (ePrint 2021/933)**: no rigorous smoothness, rank or fibre symmetry.
  **Boneh's CRT decoding**: no square structure. **Hittmeir/Harvey**: exponential.
- **Square identities.** A polynomial family x = P(t) with P(t)² ≡ Q(t)² identically would, if
  fibre-symmetric, factor N directly (codex). So an easy identity cannot hide the hard part.

### 3.4 Literature (codex rkrank_2, our literature subagent, and the critic in §7)
- **Vallée.**
  - Math. Comp. 56 (1991) 823–849; STOC 1989. The target interval is [−4N^{2/3}, 4N^{2/3}], with a Farey
    order of about N^{1/3}/4.
  - We did not read the primary PDF; AMS blocked it. Codex rkrank_2 and hypD_2 read it via ResearchGate.
  - Global quasi-uniformity constants: unverified by us.
- **Coron**, CRYPTO 2002 (LNCS 2442, 613–626): wider intervals, N^{2/3+ε}, statistically close to uniform.
- **Gentry**, CRYPTO 2004 (LNCS 3152, 179–200, §4): exactly uniform sampling at width 8N^{2/3} in expected
  polynomial time.
- **No θ < 2/3 found** (16 OpenAlex citers of the 1991 paper checked by the literature subagent).
- **Stange 2022/23 (arXiv 2211.06821v2).**
  - Hyp. 3.1 is saturation of one coordinate line of an *integer* relation lattice; it is heuristic.
  - Thm 3.2 is proved relative to a relation oracle; the output is always a multiple of the order.
  - NFS is mentioned only as a remark. There is no GF(2) statement.
  - Her "1 − 1/ζ(c+1)" should read 1/ζ(c+1) (codex; a typo).
- **Rigorous rank results.**
  - Balister–Bollobás–Morris (Annals 2018) give the first-dependency threshold for uniform integers only.
  - Coja-Oghlan et al. give rank results for configuration models, not NFS rows.
  - We found no saturation theorem for sieve relations.
- **Rigorous c < 1:** none found. We checked the 115 OpenAlex citers of LP 1992 and the citers of LV.

## 4. What would close the gap
1. **(Exp_T) for one field's smooth set**, with T = D·L^{o(1)}. Equivalently, a *joint* switching lemma for
   the NFS kept law with bounded congestion.
   - By Theorem V_κ, it suffices that the non-forced dual has dimension ≤ D^{2/3−ε}, i.e. an expansion
     statement only for V of dimension > D^{2/3−ε}. Large V are the "most generic" ones.
   - A possible lever, not worked out: rows from L^{o(1)} fields *sharing m* but with independent R.
     These have common rational parts, which allows switching (a, b) between fields. But dependencies
     must be per field (§2.7).
2. **A non-oblivious sampler below 2/3.** It would have to choose productive cells arithmetically while
   keeping fibre symmetry. Theorem L shows that oblivious choices cannot work.

## 5. Experiment E1 (pre-registered; exact; toy)
- Commands are in `notes/rkrank_prereg.md`. Analysis:
  `.venv/bin/python -m src.rkrank.analyze experiments/raw/rkrank_e1.jsonl`.
- 20 instances (blum/tri × 64/80 bits × seeds 0–4), 72 level-records. Runtime: 15–40 s at 64 bits,
  5–13 min at 80 bits.
- 0 accidental gcd events.

**Disclosed limitation (not a deviation).** The registered boxes were exhausted before 3.0·cols in 16 of
the 20 instances; the analysis skips unreached levels, as registered.
- Levels reached, per instance: 3.0 in 4, 2.0 in 8, 1.5 in 4, 1.25 in 4.
- So the "levels ≥ 2.0" predictions rest on 16 records rather than the intended ~40.
- 9 of the 10 80-bit instances stopped at b = 12000.

| prediction | registered | result |
|---|---|---|
| P0 W ⊆ Λ/2 | all records | 72/72 — **holds** |
| P1a excess_peel ≤ 4 (levels ≥ 1.5) | all records | 30/32 — **fails** (two records with 5: tri 64 s1 and tri 80 s1, both at level 1.5) |
| P1b excess_peel = 0 (levels ≥ 2.0) | ≥ 80% | 12/16 = 75% — **fails narrowly** (the others have 1, 1, 2 at 2.0 and 1 at 3.0) |
| P2a c₃ = 0 (levels ≥ 1.5); c₃ ≤ 2 | ≥ 90%; all | 32/32 — **holds** (c₃ = 0 in every record at every level except four at 1.0) |
| P2b c₁ = 0 (levels ≥ 2.0) | ≥ 80% | 12/16 = 75% — **fails narrowly** |
| P3 small squares supported on P₃ in W | ≥ 0.95 | 1141/1141 — **holds**; negative z: 514/514 |
| P4 sign obstruction ⟺ constant-1 real column | all records | 72/72 — **holds** (no sign obstruction arose; every constant real column was constant 0) |
| P5 δ ∈ [0.45, 0.55] (3 Liouville-type + 10 dense) | all instances | 20/20 — **holds** |

**Descriptive, not registered.**
- c₂ = 0 in 12/16 records at level 1.5, 11/12 at 2.0 and 4/4 at 3.0. Every residual obstruction therefore
  involves primes occurring at most twice (c₃ = 0 everywhere at ≥ 1.25), and most involve primes occurring
  at most once.
- The excess over the peeling predictor shrinks with M:
  - level 1.0: 4–23;
  - level 1.5: 0–5;
  - level 2.0: 0–2;
  - level 3.0: 0–1.
- c₁ equals dim W^⊥ − ω − #{occ = 0} in every record.
- Small-square vectors *not* supported on P₃ lie in W at rate 416/473 = 0.88. The failures contain a
  rarely covered prime.
- Single-prime functionals have δ between 0.0002 and 0.10; dense functionals stay at 0.5.

**Reading [experimentally supported, toy only].**
- In these monic cubic toys the dual obstructions are exactly:
  - the ω genuine characters;
  - functionals supported on primes that occur at most twice, almost all accounted for by primes missing
    from the algebraic 2-core.
- There are no dense near-obstructions (P5), no archimedean obstruction (P4) and no fake obstruction on
  well-covered primes (P2a, P3).
- The two sharp predictions that failed (P1, P2b) failed because a handful (≤ 5) of rare-prime
  obstructions escape the simple peeling predictor and the occ ≥ 1 cut. They did not fail through any
  dense or "fake" obstruction.
- This is the strong, sparse form of RK, "W = Λ/2 on the well-covered support", at toy scale. It says
  nothing about the asymptotic regime: d ≥ 9, L[1/2]-size factor bases, the LV coefficient law.

## 6. Candidate cards
```text
ID: RK-U  uniform rank growth (Theorem U) for one LV field's smooth set.
Proved: Theorem U (Exp_T => deficiency r after T(1+ln s) relations); Lemma S (rational model satisfies
  Exp_T with T = D L^{o(1)}); Prop G' (raw vs genuine <= O(log N)); Prop C (algebraic identities = Gen +
  O(log N) exceptional + <= r1+1 box-archimedean); Theorem V_kappa.
New lemma required: Exp_T for the joint kept law (joint switching), or deficiency <= D^{2/3-eps}.
Barriers: field averaging (B-F), fixed moments (B-M), coverage != span (B-C), random exponents too costly (B-L).
Status: OPEN; E1 toy supports "W = Lambda/2 on the algebraic 2-core" (see §5).

ID: VS-L  small-square samplers below 2/3.
Proved: Prop Q (Dixon s+sqrt(s^2+2theta); Design V (s+sqrt(s^2+4theta))/2); rank-collapse lemma;
  Theorem L (oblivious linearisable cells need h >= N^{2/3-o(1)}); Coppersmith shift lattice needs T^2 h < N.
Status: refuted for this model (oblivious linearisable cells; natural Coppersmith); open in general.
```

## 7. Reviews
REVIEWS_PLACEHOLDER
