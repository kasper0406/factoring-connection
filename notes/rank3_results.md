# rank3: a rank-r ledger for linked modules with two nonlinear sides (resbarrier §2.6, ptesplit O3)

Agent `rank3`, phase 17, 2026-10-20.
- Code: `src/rank3/forms.py` (exact Δ_F via certified balls, box counts, Liouville check), `src/rank3/ledger.py`
  (constants), `src/rank3/e1.py`, `src/rank3/e2.py`. Tests: `tests/test_rank3.py` (13 tests).
- Pre-registration: `notes/rank3_prereg.md` (commit 865d205, before E1/E2 ran). Raw data:
  `experiments/raw/rank3_E{1,2}.jsonl`.
- GPT-6.1-sol: blind round `codex/rank3_1.{md,out}` (§9); adversarial round `codex/rank3_2.{md,out}` (§10).
- Reviews: two read-only Opus critics, proofs and literature (§10).

Labels: **proved**; **proved (model)** = exact consequence of the stated model (heuristic H, linear algebra B², single
N, RAM); **heuristic** (named); **KNOWN** (cited); **measured**; **open**. Novelty is never claimed.

---

## 0. Verdict

**Plain language.**
- The open case is settled in one direction, for a large and natural class of families: **higher rank does not
  beat rank 2.** A classical volume inequality for products of linear forms (Bean–Thunder 1996) generalises
  Bean's area theorem to r variables. Combined with one elementary observation about how the linking to Z/N forces
  divisibility, it gives a rank-r ledger that charges both norms, for every N, special or not, with no counting
  heuristic about polynomials or fields. This turns resbarrier's heuristic G_r into a theorem wherever it is
  applied inside the scope below.
- The cost of rank r is a norm floor N^{(r−1)/(D−1)} in place of rank 2's N^{1/(D−1)}. It grows with r, while the
  entropy gain (D/r per candidate) does not compensate: **r = 2 is optimal, and α ≥ 1/3**, rigorously in the model.
- The §4.2 "aggregate" L[1/4] family is impossible as specified: its premise (generators of family height
  N^{1/d}) contradicts a rank-r Liouville inequality that holds for every N.
- What is not closed: (a) degenerate configurations where some r of the D linear factors are dependent
  (Δ_F = 0; this includes codex's counterexamples, where the linking is explained in characteristic 0); (b) boxes so
  large, or ranks so high with tiny boxes, that the lattice-point error term (2E+1)^{r−1} exceeds the factor base.
  There a purely geometric lattice-point heuristic replaces G_r. (c) ptesplit's O3 split side is charged only
  through its *size*: the theorem closes the O3 instances whose split-side norm is below an explicit threshold,
  not the others.
- Construction side: nothing beats rank 2. In a pre-registered toy comparison the best rank-3 design's smooth
  yield fell from 0.95× to 0.30× of the best rank-2 design's as N grew from 40 to 60 bits.

**Precisely.**
1. **Theorem R1 [proved; a corollary of Bean–Thunder 1996 + an elementary divisibility; found blind by codex too,
   in the unramified form].** For a linked rank-r module (both sides of any degree ≥ 1, D = d + d′ > r) with
   Δ_F ≠ 0, Δ_F is a nonzero integer, N^{r(r−1)(D−2)!/(D−r)!} divides it for squarefree N (ramified primes
   included), and V_F ≤ C_r·N^{−r(r−1)/(D(D−1))}, with C_r ≤ (2/r)(12(r+1))^{r−1}.
2. **Lemma R2 [proved].** #{a ∈ Z^r, |a_j| ≤ E : 0 < |F(a)| ≤ X} ≤ V_F X^{r/D} + r(D+1)(2E+1)^{r−1}.
3. **Corollary R3 [proved (model)].** Under H5_r (r(D+1)(2E+1)^{r−1} ≤ B/4), a yield ≥ B forces
   log X ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r). If also log C_r ≤ λ/2 (H7_r; e.g. r log(12(r+1)) ≤ ¼ log B), then
   α ≥ 1/3. For fixed r the floor is (32/9)^{1/3} = 1.5263 at r = 2 and 2(2(r−1)/3)^{2/3} ≥ 2.4228 for r ≥ 3. If the
   lattice-point error could be ignored, it would be (64(r−1)/(9r))^{1/3} ≥ 1.6803 for r ≥ 3. Rank 2 is optimal.
4. **Proposition R4 (rank-r Liouville) [proved; box-free].** For every basis of the module,
   ∏_{i≤D} max_j |L_i(e_j)| ≥ (2^r/(r!C_r))^{D/r}·N^{(r−1)/(D−1)}. In height units: h_total ≥ (r−1)/(D−1) − O(D log r)/x,
   for all N. This is the rigorous counterpart of G_r's height floor (r−1)/(d+1), weaker by about a factor 2 for
   two degree-d sides.
5. **Split side (ptesplit O3) [proved (model), partial].** #{|W| ≤ X₁, |W′| ≤ X₂} ≤ C_r N^{−κ}(X₁X₂)^{r/D} + error.
   The split side must therefore carry log X₁ ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r) − log X₂. Below that threshold O3
   is closed. Above it the N-cost can sit on the split side, and no two-norm ledger can say more (codex example, §6).

## 1. Setting

- K, K′ number fields of degrees d, d′ (D = d + d′), orders O ⊂ K, O′ ⊂ K′, ring maps φ: O → Z/N, φ′: O′ → Z/N,
  N squarefree.
- A **linked rank-r module**: generators w_j = (u_j, u′_j) ∈ O × O′ (j ≤ r) with φ(u_j) = φ′(u′_j) =: v_j.
  Equivalently, a rank-r sublattice of Λ = O ×_{Z/N} O′ whose two projections are injective.
- Candidates a ∈ Z^r: θ(a) = Σ a_j u_j, θ′(a) = Σ a_j u′_j. Norm forms W(a) = N_K(θ(a)) = ∏_σ ℓ_σ(a) with
  ℓ_σ = (σ(u_j))_j, and W′ = ∏_τ ℓ′_τ. F = W·W′ is a decomposable form of degree D in r variables with integer
  coefficients. Its D linear factors L_1..L_D are the ℓ_σ and ℓ′_τ.
- **Bean–Thunder invariant** Δ_F = ∏ det(L_{i_1}, …, L_{i_r}) over ordered r-tuples of distinct indices.
  - It is unchanged by rescalings L_i ↦ c_iL_i with ∏c_i = 1, since each index lies in the same number of tuples.
  - Δ_F ≠ 0 iff every r of the D factors are linearly independent ("general position", GP).
  - Under F ↦ F∘T: Δ ↦ det(T)^{D!/(D−r)!}Δ and V_F ↦ V_F/|det T|, so |Δ_F|^{(D−r)!/D!}V_F is GL_r(R)-invariant.
- V_F = vol{x ∈ R^r : |F(x)| ≤ 1}.
- rank 2 is resbarrier's C2 case: there Δ_F = ±disc F and N | Res gives N² | disc.

## 2. Theorem R1 (volume bound for linked modules)

**KNOWN input (Bean–Thunder, J. London Math. Soc. (2) 54 (1996) 39–49).** Let F be a decomposable form of degree
D > r in r variables, a product of complex linear forms, with Δ_F ≠ 0. Then

  V_F ≤ C_r·|Δ_F|^{−(D−r)!/D!},
  C_r = (2/r)∏_{k=1}^{r−1}[B(1/(r+1), k/(r+1)) + B((r−k)/(r+1), k/(r+1)) + B((r−k)/(r+1), 1/(r+1))].

- C_2 = 3B(1/3,1/3) = 15.8997 is Bean's constant; C_3 = 213.71; C_r ~ (2/π)e^{1−γ}(2r)^r (as quoted by Stewart,
  arXiv 2410.01097, and codex). The constant is sharp and is the volume for X₁⋯X_r(X₁+⋯+X_r).
- Elementary bound used here: since B(x, y) ≤ 2/x + 2/y on (0, 1], each bracket is ≤ 12(r+1), so
  log C_r ≤ (r−1)log(12(r+1)) + log(2/r) (`log_C_upper_elementary`, tested).
- Statement and constant were read from the citation in Stewart's 2024 preprint and by codex. The literature critic
  is asked to verify them against the original (§10).

**Theorem R1 [proved].** Let the module be linked, D > r, and Δ_F ≠ 0. Then:
- (i) Δ_F ∈ Z \ {0};
- (ii) for every prime p | N, v_p(Δ_F) ≥ P(r, D) := r(r−1)(D−2)!/(D−r)! = r!·C(D−2, r−2);
- (iii) hence, for squarefree N, |Δ_F| ≥ N^{P(r,D)} and

  V_F ≤ C_r·N^{−κ_r},  κ_r = P(r,D)·(D−r)!/D! = r(r−1)/(D(D−1)).

*Proof.*
- (i) The u_j, u′_j are algebraic integers, so every det is an algebraic integer in a Galois closure L of KK′.
  Gal(L/Q) permutes the D forms (σ ↦ gσ on each side), hence permutes ordered tuples. det(g·rows) = g(det rows),
  so g(Δ_F) = Δ_F. Thus Δ_F ∈ Q ∩ Z̄ = Z. It is nonzero by GP.
- (ii) *Linking gives congruent rows.* Fix p | N and a prime 𝔓 of O_L above p; normalise v(p) = 1.
  - Let 𝔭 = ker(φ mod p), a prime of O with O/𝔭 = F_p. Fix an embedding ι: K → L. By lying-over there is a prime
    𝔓₁ ⊃ ι(𝔭) of O_L. All primes of O_L above p are conjugate, so 𝔓₁ = g𝔓 for some g.
  - Put 𝔭_K := ι^{−1}(𝔓₁) ∩ O_K, a prime of O_K containing 𝔭. Let T_K be the set of embeddings σ: K → L with
    σ^{−1}(𝔓) ∩ O_K = 𝔭_K. It contains g^{−1}ι, and #T_K = e(𝔭_K|p)·f(𝔭_K|p) ≥ e := e(𝔭_K|p) (local degree;
    L_𝔓/Q_p is Galois).
  - For σ ∈ T_K and y ∈ O: y − φ̃(y) ∈ 𝔭 ⊂ 𝔭_K (φ̃ an integer lift), so
    v(σ(y) − φ̃(y)) ≥ e(𝔓|σ𝔭_K)/e(𝔓|p) = 1/e. Apply this to y = u_j: every row ℓ_σ, σ ∈ T_K, is congruent to the
    vector ṽ = (φ̃(u_j))_j with valuation ≥ 1/e.
  - The same holds on the K′ side: a set T_{K′} of size ≥ e′, rows congruent to the *same* ṽ (linking:
    φ(u_j) = φ′(u′_j)) with valuation ≥ 1/e′.
  - Let T = T_K ∪ T_{K′}, t = #T ≥ e + e′ ≥ 2. Any two rows in T differ by a vector of valuation ≥ 1/max(e, e′).
  - *Determinants.* If an r-subset I contains k ≥ 2 rows of T, subtract one of them from the others. This leaves
    k − 1 rows of valuation ≥ 1/max(e, e′), with all entries algebraic integers, so
    v(det_I) ≥ (k − 1)/max(e, e′).
  - *Counting.* Let S(t) := Σ_{|I|=r}(|I ∩ T| − 1)_+. Adding one element to T raises S by
    C(D−1, r−1) − C(D−t, r−1) ≥ C(D−2, r−2) when t ≥ 2. With S(2) = C(D−2, r−2) this gives
    S(t) ≥ (t−1)C(D−2, r−2) (identity S(t) = tC(D−1,r−1) − C(D,r) + C(D−t,r); `test_inclusion_exclusion_monotone`).
  - Hence v_p(Δ_F) = r!·Σ_I v(det_I) ≥ r!(e + e′ − 1)C(D−2, r−2)/max(e, e′) ≥ r!C(D−2, r−2) = P(r, D),
    since e + e′ − 1 ≥ max(e, e′).
- (iii) Multiply over p | N and apply Bean–Thunder. ∎

**Remarks.**
- *Unramified case.* e = e′ = 1 and t = 2: the single congruent cross pair (σ₀, τ₀) already gives P(r, D). This is
  codex's (7). Codex also gave a ramified bound with exponent divided by ⌊r²/4⌋, via 2×2 minors and the Gauss lemma.
  The inclusion–exclusion above removes that loss: ramification only adds congruent rows.
- *Small primes.* No condition p > D is needed.
- *Power basis (t-dimensional sieving).* With u_j = α^{j−1}, u′_j = β^{j−1}, every det is a Vandermonde, so
  |Δ_F| = |disc(fg)|^{P(r,D)/2} (`test_power_basis_delta_is_disc_power`). GP then holds automatically when f, g are
  squarefree and coprime. R1 becomes V_F ≤ C_r·|disc f·disc g·Res(f,g)²|^{−r(r−1)/(2D(D−1))}, the r-variable form of
  resbarrier C2.
- *What linking contributes.* It contributes exactly one coincidence modulo 𝔓 between a K-row and a K′-row. The
  r-wise determinants containing that pair carry it r(r−1)/2 times more often (normalised) than the rank-2
  resultant does. That is why the N-exponent grows like r², not r. A pairwise (2-plane) argument alone gives only the
  r-independent floor x/d of resbarrier §2.6.
- *GP is necessary for any N-saving.* Codex's counterexample: K = Q(ζ₅, √2), K′ = Q(ζ₅, √3), module
  (1, ζ₅, ζ₅²) on both sides. Then F = G⁴ with G = N_{Q(ζ₅)/Q}(a + bζ₅ + cζ₅²), and V_F does not depend on N. Here
  cross rows coincide in characteristic 0, so Δ_F = 0 and every binary cross-resultant vanishes identically. This
  is the rank-r analogue of resbarrier C2's hypothesis disc F ≠ 0.
  - Δ_F = 0 has three sources: proportional rows on one side (module inside a subfield direction), proportional
    cross rows (linking explained in characteristic 0), and r pairwise-independent but dependent rows (r ≥ 3).
  - Δ_F is computable from the module, so GP can be checked algorithmically (`delta_F`).

## 3. Lemma R2 (box count)

**Lemma R2 [proved].** For any decomposable F of degree D in r variables, any X > 0 and E ≥ 0:
#{a ∈ Z^r ∩ [−E,E]^r : |F(a)| ≤ X} ≤ vol({|F| ≤ X} ∩ [−E,E]^r) + r(D+1)(2E+1)^{r−1} ≤ V_F X^{r/D} + r(D+1)(2E+1)^{r−1}.

*Proof.*
- Let S = {x ∈ [−E,E]^r : |F(x)| ≤ X}. Define I_k = Σ_{a_1..a_k ∈ Z} ∫ 1_S(a_1..a_k, t_{k+1}..t_r) dt. Then
  I_r = the count and I_0 = vol S.
- With all coordinates but the k-th fixed, s ↦ F(…, s, …) is a polynomial of degree ≤ D. So
  {s ∈ [−E,E] : |F| ≤ X} is a union of ≤ D + 1 closed intervals (or all of [−E,E] if constant), and it holds at most
  (its length) + (D + 1) integers.
- By Tonelli, I_k ≤ I_{k−1} + (D+1)(2E+1)^{k−1}(2E)^{r−k}. Sum over k. ∎

The same slicing gives Thunder 2001's Lemma 14 (codex quotes error rD(2E+1)^{r−1}) and resbarrier's C2 count for r = 2.

## 4. Corollary R3 (the rank-r ledger) [proved (model)]

**Model.**
- resbarrier §1: heuristic H for examined candidates, cost ≥ #examined + B², single N.
- Candidates are an examined subset (size e^C) of a box of side 2E in a linked rank-r module with GP, D > r.
- **H5_r:** r(D+1)(2E+1)^{r−1} ≤ B/4. This is the rank-r form of resbarrier's H5. In L-units, (r−1)e ≤ β, i.e.
  C ≤ rβ/(r−1).

**Yield step** (as resbarrier C2(e)).
- count(Y) ≤ main(Y) + B/4, with main(Y) = C_r N^{−κ_r}Y^{r/D}.
- Yield ≥ B needs a dyadic Y with min(e^C, main(Y))·ρ(u_Y) ≥ B/(4 log₂ Y_max).
- So, writing λ = log min(e^C, main(Y)):

  **log Y ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r)**,  with λ ≥ log B − log(4 log₂Y_max) + u_Y log u_Y(1+o(1)) and λ ≤ C.

**(a) Exponent [proved (model)].** Add **H7_r:** log C_r ≤ λ/2. Sufficient: r log(12(r+1)) ≤ ¼ log B, since
λ ≥ ½ log B for large x. Then α ≥ 1/3.
- Proof. Then log Y ≥ (r−1)x/(D−1) + (D/(2r))λ ≥ 2√((r−1)xλ/(2r)) ≥ √(xλ), using (r−1)/r ≥ ½.
- If the cost were L[α, c] with α < 1/3, then λ, log B ≤ c x^α(log x)^{1−α}, and u_Y = log Y/log B ≥ √(x/(2 log B)).
- That makes u_Y log u_Y ≳ x^{(1−α)/2}(log x)^{O(1)} ≫ λ, contradicting λ ≥ u_Y log u_Y. ∎
- This holds for every N and every r ≥ 2 with r ≤ D − 1, in scope. Growing rank is allowed up to H7_r.

**(b) Constants for fixed r [proved (model)].**
- Set t = x^{1/3}(log x)^{2/3}, D = δ(x/log x)^{1/3}, log B = βt, λ = γt.
- Norm coefficient η ≥ (r−1)/δ + δγ/r ≥ 2√((r−1)γ/r). H requires γ ≥ β + η/(3β). Cost ≥ max(γ, 2β).
- H5_r: γ ≤ rβ/(r−1). For r ≥ 3 this binds, giving β^{3/2} ≥ 2(r−1)/3, i.e. cost ≥ 2(2(r−1)/3)^{2/3}.
- For r = 2 H5 does not bind, and the floor is resbarrier's 1.5263.
- Closed forms are in `ledger.scoped_constant`/`unscoped_constant` and cross-checked by brute force
  (`numeric_floor`).

| r | floor in scope (H5_r) | formal floor if the lattice error were negligible | size-model value of power-basis sieving |
|---|---|---|---|
| 2 | **1.5263** (= RB2; SNFS attains) | 1.5263 | 1.5263 (SNFS) / 1.9230 (GNFS) |
| 3 | 2.4228 | 1.6803 = (128/27)^{1/3} | 1.6803 (SNFS-type pair) / 2.1165 (GNFS, heurnfs B2) |
| 4 | 3.1748 | 1.7472 | 1.7472 / 2.2013 |
| ∞ | → ∞ | (64/9)^{1/3} = 1.9230 | 1.9230 / 2.4229 |

- **Reading.** Every rank-r ≥ 3 family is strictly worse than the rank-2 floor, both in scope and formally. The
  formal rank-3 floor 1.6803 is what SNFS-type 3-dimensional sieving attains in the size model, so the theorem's
  r-dependence is the right one.
- The comparison with heuristic G_r: G_r's height floor for two degree-d sides is (r−1)/(d+1). R1/R4 give
  (r−1)/(D−1) = (r−1)/(2d−1). That is weaker by about 2, exactly as Bean's N² is weaker than GNFS's κ = 2 at rank 2.
  It is the same factor that separates SNFS from GNFS. In both cases the exponent conclusion is unaffected.

**(c) Special-q sublattices.** As in resbarrier C2(c).
- An index-q sublattice M with q | F has V_{F∘M} = V_F/q.
- The q-lattices on which q | F correspond to degree-1 primes, so there are ≤ D of them per q.
- Summing over q ≤ Q gives a gain ≤ (D²/r)Q^{r/D} = L^{o(1)} for Q ≤ B.

## 5. Proposition R4 (rank-r Liouville inequality) and the aggregate family

**Proposition R4 [proved].**
- Let a^{(1)}, …, a^{(r)} ∈ Z^r be linearly independent, and put 𝓜 := ∏_{i≤D} max_k |L_i(a^{(k)})|. Then

  𝓜 ≥ (2^r/(r!·C_r))^{D/r}·|Δ_F|^{(D−r)!/(r(D−1)!)} ≥ (2^r/(r!C_r))^{D/r}·N^{(r−1)/(D−1)}.

- *Proof.* The cross-polytope conv(±a^{(k)}) has volume 2^r|det(a^{(k)})|/r! ≥ 2^r/r!. On it,
  |L_i(Σt_ka^{(k)})| ≤ max_k|L_i(a^{(k)})| when Σ|t_k| ≤ 1, so it lies in {|F| ≤ 𝓜}. Its volume is therefore
  ≤ V_F𝓜^{r/D}. Now apply R1. ∎
- No box and no counting are used. With a^{(k)} = e_k: for **every** basis (u_j, u′_j) of the module,
  Σ_σ log max_j|σ(u_j)| + Σ_τ log max_j|τ(u′_j)| ≥ ((r−1)/(D−1))·x − D·log(r·12(r+1)/2). Checked in E1 (P3).
- If the r candidates have comparable embedding profiles (|L_i(a^{(k)})| ≤ T·m_i for all i, k, with ∏m_i ≤ X), then
  log(T^D X) ≥ (r−1)x/(D−1) − O(D log r).

**Consequence for resbarrier §4.2 (the aggregate L[1/4] family) [proved].**
- The mechanism asked for r ≫ 3 linked generators with balanced embedding profiles and family height
  H₁H₂ ≈ N^{1/d}, i.e. generator max-heights with log 𝓜 ≈ x/d = 2x/D.
- R4 forces log 𝓜 ≥ (r−1)x/(D−1) − O(D log r). For the §4.2 parameters (r ≍ x^{1/4}(log x)^{3/4},
  D ≍ √(x/log x)) that is ≍ x^{3/4}(log x)^{5/4}. The family claimed log|F| ≍ √(x log x).
- So **for every N** the aggregate family does not exist as specified. G_r had excluded it only for typical N.
- A variant could survive only if its candidates' norms were far below their generators' size
  (|F(a)| ≪ (r^{1/2})^D𝓜). That requires massive simultaneous cancellation in many conjugates. By R1 such candidates
  occupy a region of volume ≤ C_rN^{−κ_r}X^{r/D} = e^{−Ω(x^{1/2})} at these parameters. They could be numerous only
  as lattice points far in excess of volume (the Davenport error regime, §7 (b)).

## 6. The split side (ptesplit O3)

**Corollary R5 [proved (model)].** In scope (GP, H5_r):
#{a in box : 0 < |W(a)| ≤ X₁, 0 < |W′(a)| ≤ X₂} ≤ C_r N^{−κ_r}(X₁X₂)^{r/D} + r(D+1)(2E+1)^{r−1}, because the set lies in
{|F| ≤ X₁X₂}.
- This is the "rank-r ledger that also charges the split side" that ptesplit §2.4 asked for. The split side is
  charged through its *size* X₁, not its smoothness.
- A split-by-construction family (W free for smoothness) with yield ≥ B must satisfy

  log X₁ + log X₂ ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r).

**What it closes.**
- If the split side's norm is bounded by X₁ (in ptesplit's View B, X₁ ≤ ∏_{j≤M}|F̃(a_j, b_j)|), then the charged
  side carries log X₂ ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r) − log X₁.
- So the O3 instances whose split-side norm satisfies log X₁ = o((r−1)x/(D−1)) inherit the full rank-r ledger:
  α ≥ 1/3 and r = 2 optimal.

**What it does not close.**
- When log X₁ ≳ (r−1)x/(D−1), the N-cost can sit on the split side.
- Codex's example (rank3_1): K′ = Q(ζ₅), K = Q(ζ₅, √2), θ′ = a + bζ₅ + cζ₅², θ = a + b(ζ₅ + N√2) + c(ζ₅² + Nζ₅√2).
  This is linked, with coprime norm forms and a fixed charged form W′. It shows that linking alone bounds nothing on
  W′ when all height is on the uncharged side.
- Any closure of that regime must use the split *structure* (how many split members land in the box: ptebarrier
  Theorem W-type rigidity, construction entropy), not a norm-volume ledger. ptesplit's verdict, "G_r alone does not
  suffice", is therefore sharpened: **no** two-norm ledger suffices there, and R5 is the best possible of that kind
  up to constants.
- Codex also noted that the natural rank-3 split parametrisation q(T) = (s + tT)(u + vT) has W′(q) = W′(s+tT)W′(u+vT).
  It produces no relation rank beyond its binary constituents.

## 7. What remains open (precise)

- **(a) Non-general position (Δ_F = 0).**
  - Proportional cross rows (char-0 linking, codex counterexample): no N-saving is possible, and the linking is not
    an "N-dependent" coincidence.
  - Coprime forms whose maximal minors vanish (r ≥ 3, codex's unresolved target).
  - Power bases always satisfy GP. A design would have to be deliberately degenerate.
- **(b) The lattice-point error (H5_r fails).**
  - Boxes with (2E+1)^{r−1} > B/(4r(D+1)), in L-units (r−1)e > β. This includes power-basis 3-dimensional sieving
    at its natural size-model optimum (e = 2β/3), so even that classical family is covered only under a heuristic.
  - It also includes high-rank tiny-box designs (a_j ∈ {0, ±1}, r ≈ log B).
  - The theorem then still bounds the *volume* of the small-|F| region; what is missing is that box lattice points
    do not exceed that volume by more than L^{o(1)} (or that excess lies in lower-rank sub-modules, which R1 again
    covers).
  - Call this the lattice-point heuristic LP. It is geometric and N-free, unlike G_r, which posited the scarcity of
    linked number-field configurations.
  - Thunder's 2005 bound N_F(X) ≪_{r,D}(X/m(F))^{r/D} + X^{(r−1)/D} (with m(F) ≥ |Δ_F|^{(D−r)!/(r(D−1)!)} by Hadamard)
    has no box error. Its constants are not explicit in r and D, so it does not cover growing degree.
- **(c) Very high rank inside H5_r without H7_r:** b/(4 log 12b) ≲ r ≤ b/log 3, b = log B. Here C_r can cancel the
  N-term only if D ≳ ((r−1)x/log(12r))^{1/2} ≳ x^{(1+α)/2}. Such families are excluded under a typical-size
  heuristic (log|F| ≳ D), not rigorously.
- **(d) Non-box sampling** (selection from boxes larger than B^{1/(r−1)}, as resbarrier H5) and the O3 split
  regime above (§6).

**Status of G_r after this note.**
- Inside scope, no counting heuristic about fields or polynomials is needed: R1–R4 hold for every N.
- Outside scope (b)–(c), G_r's role is taken over by LP or typical-size heuristics. They are geometric and say
  nothing about N.

## 8. Literature

- **Bean–Thunder 1996** (JLMS 54, 39–49): V_F|Δ_F|^{(D−r)!/D!} ≤ C_r for Δ_F ≠ 0. Bean 1994 (Compositio 92) is
  r = 2.
- **Thunder 2001** (Ann. Math. 153): finite type ⇔ N_F(m) < ∞ for all m; N_F(m) ≪_{r,D} m^{r/D}; Lemma 14 (box
  discrepancy). **Thunder 2003/2005** (J. Number Theory; Compositio; arXiv math/0305078): V_F ≪ M(F)^{−1/D}(log)^{r−1};
  asymptotics with height-dependent error.
- **Stewart 2024** (arXiv 2410.01097): essentially finite type; quotes Bean–Thunder (19) with C_r ~ c(2r)^r.
- **Evertse 2000** (norm-form counts with explicit constants (16d)^{(r+11)³/3}); **Schmidt 1972** (Ann. Math. 96;
  norm-form finiteness); **Evertse–Győry**, *Unit Equations* Ch. 9 and *Discriminant Equations* §14.9
  (semi-discriminants); **J. Liu** (thesis, Leiden/Bordeaux 2015; p-adic decomposable form inequalities). These are
  from codex rank3_1, to be verified by the literature critic.
- None of these states a linked-module or factoring application. No prior statement of R1 (the N-divisibility of
  Δ_F) or R4 was found by us or by codex. R1 is an immediate corollary of Bean–Thunder plus an elementary
  congruence; no novelty is claimed.

## 9. Experiments (pre-registered, commit 865d205)

### E1: theorem checks (exact; `src/rank3/e1.py`, 1.1 s)

- 30 instances:
  - (A) GNFS-like, N = pq with p, q of 10–11 bits, monic f, g with a common root mod N;
  - (B) SNFS-like, N = the squarefree part of Res(f, g) for small f, g;
  - (C) non-power-basis modules, P_j ≠ Q_j, linked only mod N.
- Cells (d, d′, r) ∈ {(3,3,3), (4,3,3), (4,4,3), (4,4,4)}. Δ_F was certified as an exact integer with arb balls
  (prec up to 2^17 bits).

| Prediction | Result |
|---|---|
| P1: v_p(Δ_F) ≥ P(r,D) for all p \| N | **held**, 0 violations. Equality v_p = P (24/30/36/360) in 29/30 instances; one (C) instance has 48/42 > 36 (extra coincidences) |
| P2: count ≤ C_r\|Δ_F\|^{−(D−r)!/D!}X^{r/D} + r(D+1)(2E+1)^{r−1} at every grid X | **held**, 0 violations |
| P3: R4 Liouville | **held**, 0 violations; margins 10–54 nats |
| P4: max saturation ≥ 10^{−3} in ≥ half of (B) | **held**: 12/12 (0.027–0.18) |
| P5: median saturation (B) > (A) per cell | **held** in all 4 cells: 0.143/0.023, 0.132/0.037, 0.114/0.023, 0.029/0.003 |

- Post-hoc (not pre-registered; `tests/test_rank3.py::test_T1_ramified_prime` and a scratch run): six instances with
  f Eisenstein at p, so the linking prime is totally ramified (e = 3). They gave v_p = 48 ≥ 24, against the proof's
  lower bound r!·S(4)/3 = 40.
- Reading: the volume bound with the exact Δ_F is saturated to ≈ 0.1–0.2 by SNFS-like rank-3 families, about the
  same as resbarrier's rank-2 E1 (≈ 0.2). R1 is sharp up to a shape constant, not only in exponent.
- Caveat: eligible X (main term > error) exist only at large X in these tiny boxes. P4/P5 are therefore weak
  sharpness tests.

### E2: construction comparison (toy sizes; measured, not asymptotic; `src/rank3/e2.py`, 16 s)

- Sizes n = log₂N ∈ {40, 50, 60}, 3 seeds each; B = 2^12. Yield = candidates with both norms B-smooth.
- Equal boxes: 2^17 points. Primitive candidates: 79 791 at rank 2 and 112 753 at rank 3, which favours rank 3.

| n (median over seeds) | R2-rat | R2-nn | R3-rat | R3-nn | R3-lll | best3/best2 |
|---|---|---|---|---|---|---|
| 40: yield / median log₂\|F\| | 3818 / 49 | 695 / 64 | 3711 / 57 | 1219 / 66 | 80 / 83 | 0.92–0.97 |
| 50 | 2795 / 54 | 536 / 69 | 1451 / 67 | 691 / 77 | 11 / 99 | 0.44–0.55 |
| 60 | 1330 / 59 | 206 / 74 | 444 / 77 | 117 / 86 | 0 / 113 | 0.29–0.33 |

- **P6** (best rank-3 yield ≤ 2× best rank-2): **held** 9/9. **P7** (rank-3 medians ≥ best rank-2 median − 2 bits):
  **held** 9/9. The kill/lead signal was not triggered.
- The ratio falls with n, as the ledger predicts.
- R3-lll, the attempt to build a non-power rank-3 module from the shortest vectors of the linked lattice Λ, was the
  worst design throughout. Its short vectors include (1, 1), and the others pay the full lattice covolume.
- "Best rank-2 design" was operationalised as the larger rank-2 yield (always R2-rat). The prereg did not fix this.
- Disclosed design choices:
  - R3-lll's weight search uses rounded real embeddings.
  - R3-lll uses norms of algebraic integers (Res/lc^{deg}). The other designs use the standard integers Res(f, A),
    which include leading-coefficient powers.

## 10. Reviews and adversarial round

(See the end of this note; filled in after the reviews.)

## 11. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_rank3.py      # 13 tests, < 1 s
.venv/bin/python -m src.rank3.e1                        # E1, ~1 s
.venv/bin/python -m src.rank3.e2                        # E2, ~16 s
.venv/bin/python -c "from src.rank3.ledger import *; print([float(scoped_constant(r)) for r in (2,3,4)])"
```
