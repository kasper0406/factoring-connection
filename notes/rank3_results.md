# rank3: a rank-r ledger for linked modules with two nonlinear sides (resbarrier §2.6, ptesplit O3)

Agent `rank3`, phase 17, 2026-10-20. Revised after codex rank3_2 (MAJOR REVISION) and two Opus critics (proofs:
ACCEPT WITH FIXES; literature: verified with citation fixes). §10 lists every applied fix.

**Files.**
- Code:
  - `src/rank3/forms.py`: exact Δ_F via certified arb balls; box counts; Liouville checks.
  - `src/rank3/ledger.py`: constants.
  - `src/rank3/e1.py`, `src/rank3/e2.py`: the pre-registered experiments.
  - `src/rank3/e1_posthoc.py`: post hoc.
- Tests: `tests/test_rank3.py`, 18 tests. The 5 marked post-hoc were added after review.
- Pre-registration: `notes/rank3_prereg.md` (commit 865d205, before E1/E2 ran).
- Raw data: `experiments/raw/rank3_E{1,2}.jsonl`; post hoc `rank3_E1_posthoc.jsonl`.
- GPT-6.1-sol: blind round `codex/rank3_1.{md,out}`; adversarial round `codex/rank3_2.{md,out}`.

**Labels.**
- **proved**: proved here, or a cited published theorem.
- **proved (model)**: an exact consequence of the stated model. The model is heuristic H, linear algebra B², a single
  N and the RAM model, plus the named scope conditions.
- **heuristic**: depends on a named assumption.
- **KNOWN**: cited.
- **measured**: an experimental result.
- **open**: not settled.

No novelty is claimed.

---

## 0. Verdict

**Plain language.**
- **In a large class of families, higher rank provably does not beat rank 2.**
  - The class: boxes of side up to about B^{1/(r−1)}, and modules whose linear factors are in general position.
  - The tool: a classical volume inequality for products of linear forms (Bean–Thunder 1996), which extends Bean's
    area theorem to r variables.
  - The step added here: an elementary observation that the shared reduction to Z/N forces a high power of N to
    divide the invariant that inequality uses.
  - Together they give a rank-r ledger that charges both norms. It holds for every N, special or not, and needs no
    counting heuristic about polynomials or fields.
- **What it costs.** Rank r pays a norm floor of N^{(r−1)/(D−1)}, where rank 2 pays N^{1/(D−1)} (D is the total
  degree). The entropy gain does not compensate. The derived lower-bound constants are minimised at r = 2, and the
  exponent is α ≥ 1/3.
- **Where resbarrier's heuristic G_r stood, there is now a theorem**, with a weaker height floor ((r−1)/(D−1) instead
  of G_r's (r−1)/(d+1)). The weaker floor still suffices for these conclusions.
- **The §4.2 "aggregate" L[1/4] family**, as specified (generator heights ≈ N^{1/d}), is impossible for every N
  whenever its linear factors are in general position.
- **What is not closed.**
  - (a) Degenerate modules, in which some r of the D linear factors are dependent (Δ_F = 0).
  - (b) Boxes larger than about B^{1/(r−1)}. There the lattice-point error term (2E+1)^{r−1} exceeds the factor base.
    This includes classical power-basis 3-dimensional sieving at its size-model optimum. There a lattice-point
    assumption, not G_r, is what is missing.
  - (c) A narrow very-high-rank window.
  - (d) ptesplit's O3 split side, above an explicit size threshold. There linking alone cannot charge the other side.
- **Construction side.** Nothing beats rank 2. In a pre-registered toy comparison, the best rank-3 design's smooth
  yield fell from 0.95× to 0.30× of the best rank-2 design's as N grew from 40 to 60 bits.

**Precisely** (setting in §1; D = d + d′ > r, both sides of any degree ≥ 1).
1. **Theorem R1 [proved].** It is a corollary of Bean–Thunder 1996 plus an elementary divisibility. Codex found the
   unramified form blind.
   - If Δ_F ≠ 0 (general position, GP), then Δ_F ∈ Z and v_p(Δ_F) ≥ r(r−1)(D−2)!/(D−r)! for every p | N.
   - This holds at ramified primes, at index divisors, and at p ≤ D.
   - Hence V_F ≤ C_r N^{−r(r−1)/(D(D−1))} for squarefree N, with C_r ≤ (2/r)(12(r+1))^{r−1}.
2. **Lemma R2 [proved].** #{a ∈ Z^r : |a_j| ≤ E, 0 < |F(a)| ≤ X} ≤ V_F X^{r/D} + rD(2E+1)^{r−1}.
3. **Corollary R3 [proved (model)].** Assume GP and H5_r: r(D+1)(2E+1)^{r−1} ≤ B/4.
   - A yield ≥ B forces log X ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r).
   - With H7_r (log C_r ≤ λ/2), α ≥ 1/3.
   - For fixed r, the in-scope floor is 1.5263 at r = 2 and 2(2(r−1)/3)^{2/3} ≥ 2.4228 at r ≥ 3.
   - If the lattice error were negligible, the floor would be (64(r−1)/(9r))^{1/3}, which is 1.6799 at r = 3.
   - In both versions the constants increase with r. No rank-r ≥ 3 family in scope beats the universal rank-2 (SNFS)
     floor 1.5263. This does not say that a rank-2 construction attaining 1.5263 exists for every N.
4. **Proposition R4 (rank-r Liouville) [proved; elementary; no box].**
   - For any r independent candidates a^{(k)}: ∏_{i≤D} max_k |L_i(a^{(k)})| ≥ r^{−D/2}|Δ_F|^{(D−r)!/(r(D−1)!)}.
   - With R1(ii) the right side is ≥ r^{−D/2}N^{(r−1)/(D−1)}.
   - The bound is sharp in exponent for SNFS with a rational side.
5. **Split side [proved (model), partial].**
   - #{|W| ≤ X₁, |W′| ≤ X₂} ≤ C_r N^{−κ_r}(X₁X₂)^{r/D} + error.
   - A split family must have log X₁ + log X₂ ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r). This closes the O3 instances whose
     split-side norm is below this threshold, under the charged-side analogue of H.
   - Above the threshold no bound on the charged norm alone is possible.

## 1. Setting

- K, K′ are number fields of degrees d, d′, with D = d + d′.
- O ⊂ K and O′ ⊂ K′ are orders, with unital ring maps φ: O → Z/N and φ′: O′ → Z/N. N is squarefree; this is used only
  in R1(iii).
- **Linked rank-r module:** r generator pairs w_j = (u_j, u′_j) ∈ O × O′ such that φ(u_j) = φ′(u′_j) =: v_j.
  - Projections need not be injective: r may exceed d or d′, and a rational side (d′ = 1) is allowed.
  - Candidates with W(a) = 0 or W′(a) = 0 are discarded.
  - The D linear factors form a disjoint, tagged union of the two sides, even when K = K′.
- **Candidates:** a ∈ Z^r, with θ(a) = Σ a_j u_j and θ′(a) = Σ a_j u′_j.
- **Norm forms:** W(a) = N_K(θ(a)) = ∏_σ ℓ_σ(a) with ℓ_σ = (σ(u_j))_j, and likewise W′ = ∏_τ ℓ′_τ.
  - F = W·W′ is a decomposable form of degree D in r variables with integer coefficients.
  - Its linear factors L_1..L_D are the ℓ_σ and the ℓ′_τ.
- **Bean–Thunder invariant:** Δ_F = ∏ det(L_{i_1},…,L_{i_r}), over ordered r-tuples of distinct indices.
  - It is unchanged under rescaling L_i ↦ c_iL_i with ∏c_i = 1, because each index lies in r(D−1)!/(D−r)! tuples.
  - Δ_F ≠ 0 ⇔ every r of the D factors are independent. This is **GP**.
  - Under F ↦ F∘T, Δ ↦ det(T)^{D!/(D−r)!}Δ and V_F ↦ V_F/|det T|.
- **Volume:** V_F = vol{x ∈ R^r : |F(x)| ≤ 1}.
  - GP and D > r imply that at most r − k factors vanish on any k-dimensional subspace. That is Thunder's a(F) < D/r,
    so V_F < ∞.
- **Rank 2** is resbarrier's C2 case: Δ_F = ±disc F.
- **Power bases** (u_j = α^{j−1}, u′_j = β^{j−1}; t-dimensional sieving, ptesplit's residue view): every minor is a
  Vandermonde, and |Δ_F| = |disc(fg)|^{P(r,D)/2}. So GP holds whenever fg is squarefree and D > r
  (`test_power_basis_delta_is_disc_power`).

## 2. Theorem R1

**KNOWN input (Bean–Thunder, J. London Math. Soc. (2) 54 (1996) 39–49, doi:10.1112/jlms/54.1.39).**
- Statement: let F = ∏_{i≤D} L_i, where the L_i are complex linear forms in r < D variables and Δ_F ≠ 0. Then

  V_F ≤ C_r·|Δ_F|^{−(D−r)!/D!}.

- C_r is the value of |Δ|^{(D−r)!/D!}V at X₁⋯X_r(X₁+⋯+X_r), which has degree r + 1, so the bound is sharp.
- C_r ~ (2/π)e^{1−γ}(2r)^r. The source is the Bean–Thunder abstract and zbMATH Zbl 0854.11019 (a review by Bean).
- The Beta-product formula is taken from Stewart, arXiv 2410.01097, eq. (19):
  C_r = (2/r)∏_{k=1}^{r−1}[B(1/(r+1), k/(r+1)) + B((r−k)/(r+1), k/(r+1)) + B((r−k)/(r+1), 1/(r+1))].
- Checks of the formula:
  - C_2 = 3B(1/3,1/3) = 15.8997, Bean's constant (Compositio 92 (1994) Thm 1).
  - The proof critic computed vol{|xyz(x+y+z)| ≤ 1} = 213.7144 by quadrature, against C_3 = 213.7148.
  - The asymptotic ratio agrees.
- Verification status: checked against the abstract and the review; the full text was not consulted.
- Elementary bound used here: B(x, y) ≤ 2/x + 2/y on (0, 1], so log C_r ≤ (r−1)log(12(r+1)) + log(2/r).
- The literature critic measured C_r < (2r)^r for all r ≤ 5000 (measured, unproved; not used).

**Theorem R1 [proved].** Let the module be linked, with D > r and Δ_F ≠ 0. Then:
- (i) Δ_F ∈ Z \ {0}.
- (ii) For every prime p | N, v_p(Δ_F) ≥ P(r, D) := r(r−1)(D−2)!/(D−r)! = r!·C(D−2, r−2).
- (iii) For squarefree N, |Δ_F| ≥ N^{P(r,D)}, and

  **V_F ≤ C_r·N^{−κ_r}, κ_r = P(r,D)(D−r)!/D! = r(r−1)/(D(D−1)).**

*Proof of (i).* Every minor is an algebraic integer in a Galois closure L of KK′. Gal(L/Q) permutes the factors
within each side (σ ↦ gσ), so it permutes the ordered tuples, and g(det rows) = det(g·rows). Hence Δ_F is Galois
invariant, so it lies in Q ∩ Z̄ = Z. It is nonzero by GP.

*Proof of (ii).* Fix p | N and a prime 𝔓 of O_L above p. Normalise v(p) = 1.
- **The linking primes.** 𝔭 := ker(φ mod p) is a prime of O with O/𝔭 = F_p.
  - Fix an embedding ι: K → L. By lying-over there is a prime 𝔓₁ ⊃ ι(𝔭) of O_L, and 𝔓₁ = g𝔓 for some g.
  - Put 𝔭_K := ι^{−1}(𝔓₁) ∩ O_K. It contains 𝔭.
- **The congruent K-rows.** Let T_K be the set of embeddings σ with σ^{−1}(𝔓) ∩ O_K = 𝔭_K.
  - T_K contains g^{−1}ι, and #T_K = e(𝔭_K|p)f(𝔭_K|p) (the local degree, since L_𝔓/Q_p is Galois).
  - For σ ∈ T_K and y ∈ O, y − φ̃(y) ∈ 𝔭 ⊂ 𝔭_K, where φ̃ is an integer lift. Hence
    v(σ(y) − φ̃(y)) ≥ e(𝔓|σ𝔭_K)/e(𝔓|p) = 1/e, where e := e(𝔭_K|p).
  - So every row ℓ_σ with σ ∈ T_K is congruent to ṽ = (φ̃(u_j))_j, with valuation ≥ 1/e.
- **The congruent K′-rows.** On the K′ side there is likewise a set T_{K′} of size ≥ e′. Its rows are congruent, at the
  *same* 𝔓, to the same ṽ, because linking gives φ(u_j) = φ′(u′_j).
  - Let T := T_K ⊔ T_{K′}. Then t := #T ≥ e + e′ ≥ 2.
  - Any two rows of T differ by a vector of valuation ≥ 1/max(e, e′).
- **Minors.** Suppose an r-subset I contains k ≥ 2 rows of T. Subtract one of these rows from the other k − 1. Then
  v(det_I) ≥ (k−1)/max(e, e′), since all entries are integral.
- **Counting.** Let S(t) := Σ_{|I|=r}(|I∩T| − 1)_+.
  - Enlarging T from s to s+1 elements raises S by C(D−1, r−1) − C(D−1−s, r−1) ≥ C(D−2, r−2) when s ≥ 1.
  - S(2) = C(D−2, r−2), so S(t) ≥ (t−1)C(D−2, r−2).
  - Closed form: S(t) = tC(D−1,r−1) − C(D,r) + C(D−t,r) (`test_inclusion_exclusion_monotone`).
- **Conclusion.** v_p(Δ_F) = r!Σ_I v(det_I) ≥ r!(e+e′−1)C(D−2,r−2)/max(e,e′) ≥ P(r, D), since
  e + e′ − 1 ≥ max(e, e′). ∎

*(iii).* Multiply over p | N and apply Bean–Thunder.

**Remarks.**
- **Ramification and index divisors.** In the unramified case (e = e′ = 1) a single congruent cross pair already
  gives P(r, D); this is codex's (7). Codex's ramified bound lost a factor ⌊r²/4⌋ and needed p > D. The
  inclusion–exclusion argument removes both: ramification only adds congruent rows.
  - Adversarial checks found no counterexample:
    - p ∈ {2,3,5,7} ≤ D (60 instances);
    - index divisors with residue degree 3 above F_p (10);
    - Eisenstein on both sides, e = 3 and e′ = 4 (12);
    - codex's N = 2 index-8 example (v₂ = 168 ≥ 24).
  - Post-hoc tests: `test_T1_ramified_prime`, `test_T1_small_prime_and_index_divisor`.
- **What linking contributes.** Exactly one coincidence modulo 𝔓 between a K-row and a K′-row, up to ramification.
  - The r-wise minors carry it r(r−1)/2 times more often (after normalisation) than the rank-2 resultant does. That
    is why κ_r grows like r².
  - A pairwise (2-plane) argument alone gives only the r-independent floor x/d of resbarrier §2.6.
- **GP cannot be dropped.**
  - Codex's example: K = Q(ζ₅, √2), K′ = Q(ζ₅, √3), module (1, ζ₅, ζ₅²) on both sides. Then F = G⁴ with V_F
    independent of N, and cross rows coincide in characteristic 0.
  - This is the analogue of resbarrier C2's hypothesis disc F ≠ 0.
  - It shows that linking alone gives no uniform saving without GP. It does not show that non-GP families are useless.
  - Δ_F is computable from the module (`delta_F`).

## 3. Lemma R2 (box count) [proved]

For any decomposable F of degree D in r variables:

#{a ∈ Z^r ∩ [−E,E]^r : |F(a)| ≤ X} ≤ vol({|F| ≤ X} ∩ [−E,E]^r) + rD(2E+1)^{r−1} ≤ V_F X^{r/D} + rD(2E+1)^{r−1}.

*Proof.*
- Fix all coordinates but one. The function s ↦ F is then a polynomial P of degree ≤ D, or a constant.
- {|P| ≤ X} = {P² − X² ≤ 0} is a union of at most D intervals, because P² − X² has degree 2D and positive leading
  coefficient. Each interval contains at most (length + 1) integers.
- Apply this one coordinate at a time and use Tonelli: I_k ≤ I_{k−1} + D(2E+1)^{k−1}(2E)^{r−k}. ∎

This is Thunder 2001 (Ann. Math. 153) Lemma 14, which gives the same error dr(2B₀+1)^{r−1}, two-sided. The
pre-registered code uses the cruder r(D+1)(2E+1)^{r−1}.

## 4. Corollary R3 (the rank-r ledger) [proved (model)]

**Model.** As in resbarrier §1, plus:
- the module satisfies GP and D > r;
- the examined set has e^C points, inside a box of side 2E;
- **H5_r:** r(D+1)(2E+1)^{r−1} ≤ B/4. In L-units, (r−1)e ≤ β, i.e. C ≤ rβ/(r−1).

**Yield step** (repairing resbarrier C2(e)'s bookkeeping, codex #5).
1. Discard the B/4 examined candidates with smallest |F| once and for all. They yield at most B/4 relations.
2. For the remaining candidates, the cumulative counts satisfy n(Y) ≤ main(Y) := C_r N^{−κ_r}Y^{r/D}, by R1, R2 and
   H5_r.
3. Truncate the large-norm tail. Choose Y* with e^C·ρ(u_{Y*}) ≤ B/4. Then log Y* = O(u*·log B) with u* log u* ≈ C,
   which is polynomial in x, so the number of dyadic bins is poly(x) = L^{o(1)}.
4. Bins (Y/2, Y] use ρ at the lower endpoint.
5. Hence a yield ≥ B needs some bin with min(e^C, main(Y))ρ(u_{Y/2}) ≥ B/L^{o(1)}. Writing
   λ := log min(e^C, main(Y)):

  **log Y ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r)**,  λ ≤ C,  λ ≥ log B + u log u (1+o(1)) − o(log B).

**(a) Exponent.** Add **H7_r:** log C_r ≤ λ/2. A sufficient condition is r log(12(r+1)) ≤ ¼ log B. Then α ≥ 1/3.
- *Proof.* We have log Y ≥ (r−1)x/(D−1) + (D/(2r))λ ≥ 2√((r−1)xλ/(2r)) ≥ √(xλ), using (r−1)/r ≥ ½.
- Suppose the cost were L[α, c] with α < 1/3. Then λ, log B ≤ c x^α(log x)^{1−α}, and λ ≥ ½ log B, so
  u = log Y/log B ≥ √(x/(2 log B)).
- That makes u log u ≳ x^{(1−α)/2} ≫ λ, a contradiction. ∎
- Growing r is allowed, subject to H5_r and H7_r.

**(b) Constants for fixed r.**
- Parametrise t = x^{1/3}(log x)^{2/3}, D = δ(x/log x)^{1/3}, log B = βt, λ = γt.
- The norm coefficient is η ≥ (r−1)/δ + δγ/r ≥ 2√((r−1)γ/r).
- H requires γ ≥ β + η/(3β), and the cost is ≥ max(γ, 2β).
- H5_r means γ ≤ rβ/(r−1). It binds for r ≥ 3, giving β^{3/2} ≥ 2(r−1)/3.
- Code: `ledger.scoped_constant` and `unscoped_constant`. Brute-force checks are `numeric_floor` and the proof critic's
  independent optimisation.

| r | floor in scope (H5_r) | formal floor (lattice error neglected) | size-model value, power-basis sieving: special pair / GNFS |
|---|---|---|---|
| 2 | **1.52629** (= RB2) | 1.52629 | 1.52629 (SNFS) / 1.92300 |
| 3 | 2.42283 | 1.67989 = (128/27)^{1/3} | 1.67989 / 2.11653 (heurnfs B2) |
| 4 | 3.17480 | 1.74716 | 1.74716 / 2.20128 |
| r → ∞ (fixed-r limit) | → ∞ | (64/9)^{1/3} = 1.92300 | 1.92300 / (128/9)^{1/3} = 2.42283 |

**Reading.**
- The lower-bound constants increase with r, in scope and formally.
- The formal rank-3 value 1.67989 is what a special-N two-algebraic power-basis pair attains in the size model. E1(B)
  saturates R1 at about 0.1–0.2, so the r-dependence is right.
- The formal rank-3 floor lies below the GNFS constant 1.923. For a given general N, nothing here shows that a rank-2
  construction beats every rank-3 one. It shows only that no rank-3 family goes below the universal rank-2 (SNFS)
  floor, or below 1/3.
- **Comparison with G_r.**
  - G_r's height floor for two degree-d sides is (r−1)/(d+1). R1 and R4 give (r−1)/(D−1) = (r−1)/(2d−1).
  - That is weaker by about a factor of 2, the same factor that separates SNFS from GNFS at rank 2.
  - The exponent conclusions are unaffected.
- **Folklore context.** That 2D sieving is asymptotically optimal is known heuristically in the analogous
  large-characteristic DLP setting: JLSV, CRYPTO 2006 §3.1, §5; Barbulescu–Pierrot, LMS JCM 17A (2014) §5.2. We know
  of no factoring-specific proof. In medium characteristic, higher-dimensional sieving does help (Gaudry–Grémy–Videau
  2016; Grémy 2019).

**(c) Special-q (main term only).**
- For standard special-q lattices M defined by a degree-1 prime ideal with a nonzero reduction functional:
  V_{F∘M} = V_F/q, and there are ≤ D such lattices per q. The summed main-term gain is ≤ (D²/r)Q^{r/D}.
- This is L^{o(1)} when (r/D)log Q = o(log L).
- The lattice errors of the ≈ Dπ(Q) sublattices must be controlled separately: Σ_M rD(2E_M+1)^{r−1} ≤ B/4.
- Modules scaled by q, where q | F identically, are excluded (codex #9). The same caveat applies to resbarrier C2(c).

## 5. Proposition R4 (rank-r Liouville inequality) and the aggregate family

**Proposition R4 [proved, elementary].**
- Statement: let a^{(1)},…,a^{(r)} ∈ Z^r be independent and 𝓜 := ∏_{i≤D} max_k |L_i(a^{(k)})|. Then

  𝓜 ≥ r^{−D/2}·|Δ_F|^{(D−r)!/(r(D−1)!)} ≥ r^{−D/2}·N^{(r−1)/(D−1)} (the last step uses R1(ii)).

- *Proof.*
  - Let A be the matrix of the a^{(k)}, so |det A| ≥ 1. Put G = F∘A.
  - Then Δ_G = det(A)^{D!/(D−r)!}Δ_F.
  - Hadamard bounds each r-minor of G by ∏‖L_iA‖₂, and ‖L_iA‖₂ ≤ √r·max_k|L_i(a^{(k)})|.
  - Each index lies in r(D−1)!/(D−r)! ordered tuples. ∎
- The proof critic noted that this form is sharper than the cross-polytope proof via Bean–Thunder: it loses
  D·log(r)/2 instead of about D·log(r²).
- *Check.* `liouville_hadamard_check` holds on all 30 E1 instances (post hoc), with margins 3.5–41.6 nats. The
  pre-registered P3 used the cross-polytope form.
- **Special case a^{(k)} = e_k.** For every basis of the module, Σ_σ log max_j|σ(u_j)| + Σ_τ log max_j|τ(u′_j)| ≥
  ((r−1)/(D−1))x − (D/2)log r.
- **Sharpness.** SNFS with a rational side (f small of degree d, g = x − m, m ≈ N^{1/d}, power basis) has
  log 𝓜 = (r−1)log m = (r−1)x/(D−1).

**Consequence for resbarrier §4.2 [proved, for GP families].**
- *What §4.2 needs.* r ≫ 3 linked generators with balanced profiles and family height H₁H₂ ≈ N^{1/d}.
  - Read as generator max-heights, that means log 𝓜 ≈ x/d.
  - Read as coefficient height, it means the same up to O(D + r): each linear factor's Mahler measure is at least its
    largest coefficient, so 𝓜_W ≤ M(W) ≤ ‖W‖₂ ≤ C(d+r−1, r−1)^{1/2}H(W), and likewise for W′.
- *What R4 forces.* log 𝓜 ≥ (r−1)x/(D−1) − (D/2)log r.
  - At §4.2's parameters (r ≍ x^{1/4}(log x)^{3/4}, D ≍ √(x/log x)), this is ≍ x^{3/4}(log x)^{5/4}.
  - The family needs ≍ √(x log x).
- **So for every N, no GP family with these parameters and this height bound exists.** G_r had excluded it only for
  typical N.
- *What remains open.*
  - Non-GP aggregate families (§7(a)).
  - Families described only through selected small-norm candidates rather than their generators. Their norms would
    have to lie far below their generators' size, which needs simultaneous cancellation in many conjugates. By R1
    such candidates occupy volume e^{−Ω(√x)} at these parameters, so they could be numerous only as lattice points in
    excess of volume (§7(b)).

## 6. The split side (ptesplit O3)

**Corollary R5 [proved].** Assume GP. Then

#{a ∈ box : 0 < |W(a)| ≤ X₁, 0 < |W′(a)| ≤ X₂} ≤ C_r N^{−κ_r}(X₁X₂)^{r/D} + rD(2E+1)^{r−1},

because the set lies inside {|F| ≤ X₁X₂}.

**What it gives** [proved (model)].
- Hypotheses: H′ (the charged side W′ is B-smooth as often as a random integer of its size), with λ defined from
  W′-smoothness, plus GP, H5_r and H7_r.
- Conclusion: a family whose W is free for smoothness but bounded in size by X₁ must satisfy

  log X₁ + log X₂ ≥ (r−1)x/(D−1) + (D/r)(λ − log C_r).

- If log X₁ = o((r−1)x/(D−1)), the charged side inherits the full rank-r ledger: α ≥ 1/3, with constants as in §4.
- This is a rank-r ledger that charges the split side through its *size*, which ptesplit §2.4 asked for.
- Applicability: ptesplit's residue-view representation is a power basis, so GP holds whenever fg is squarefree.
- The geometric inequality alone does not supply the transfer of H to box-linear element dictionaries. ptesplit
  left that transfer unproved, and it remains a hypothesis here.

**What it does not give.**
- When log X₁ ≳ (r−1)x/(D−1), linking alone cannot force the charged norm to grow.
- GP example (codex rank3_2; it replaces rank3_1's example, which has Δ_F = 0):
  - K = K′ = Q(ζ₅), with generators (1, ζ₅+N, (ζ₅+N)²) on the K side and (1, ζ₅, ζ₅²) on the K′ side.
  - N is squarefree and supported on primes that split.
  - All rows are Vandermonde rows at distinct parameters, so GP holds, and the charged form W′ stays fixed.
  - This does not make W smooth by construction.
- So in that regime, any closure must use the split *structure*, such as ptebarrier Theorem W-type rigidity or
  construction entropy, rather than a norm-volume ledger.
- Also from codex: the natural rank-3 split parametrisation (s+tT)(u+vT) has
  W′((s+tT)(u+vT)) = W′(s+tT)W′(u+vT). It yields no relation rank beyond its binary constituents.

## 7. What remains open

- **(a) Non-general position, Δ_F = 0.**
  - Proportional cross rows, where the linking holds in characteristic 0 (codex's example).
  - Coprime forms whose maximal minors vanish (r ≥ 3).
  - Natural symmetries or subfield constructions could produce these, so they are not only "deliberate" degeneracies.
  - One unexamined route: a congruent cross pair still divides every *nonzero* minor containing it (a
    semi-discriminant, as in Evertse–Győry, *Discriminant Equations* §14.9). Thunder 2005 Thm 4 gives
    V_F ≪_{r,D} m(F)^{−r/D}(1+|log m(F)|)^{r−1}. Whether the nonzero minors bound m(F) from below is open, and the
    constants are not stated explicitly.
- **(b) Lattice-point error, where H5_r fails.**
  - These are boxes with (r−1)e > β. They include power-basis 3-dimensional sieving at its size-model optimum,
    e = 2β/3, and in general all selection from boxes larger than B^{1/(r−1)}.
  - R1 still bounds the *volume* of the small-|F| region. What is missing is control of the lattice points in excess
    of volume.
  - **LP (an alternative assumption, not shown weaker than G_r).** For thresholds Y at which a family's bulk relations
    are produced, the box points with |F| ≤ Y number at most L^{o(1)}·V_F Y^{r/D}. Points lying in proper rational
    subspaces are excepted, and those are again linked modules of lower rank.
    - LP cannot hold uniformly at every threshold: ±e₁ has F = 1 when u₁ = u′₁ = 1.
    - Restrictions to subspaces need not keep GP.
    - LP is N-free and geometric, whereas G_r posited that linked field configurations are scarce.
  - Thunder 2005 (Compositio 141, Thm 3) gives N_F(X) ≪_{r,D} (X/m(F))^{r/D} + X^{(r−1)/D}, with no box error, for F of
    finite type with a′(F) < D/r (for example gcd(r, D) = 1).
    - With m(F) ≥ |Δ_F|^{(D−r)!/(r(D−1)!)} (Hadamard), this would give a box-free ledger.
    - But a′ < D/r can fail (for example r = 2 with D even), and its constants are not stated with a growth in r and D
      usable at D ≍ x^{1/3}.
- **(c) A very-high-rank window inside H5_r where H7_r fails.**
  - Under H5_r a yield ≥ B forces 2E+1 ≥ 4r(D+1), so in-scope rank is r ≤ 1 + log B/log(4r(D+1)).
  - Tiny-box designs (a_j ∈ {0, ±1}) are therefore entirely case (b).
  - The window is log B/(4 log(12 log B)) ≲ r ≤ log B/log(4r(D+1)).
  - There C_r can cancel the N-term only if D ≳ ((r−1)x/log(12r))^{1/2}. Such families are excluded under a
    typical-size heuristic (log|F| ≳ D), not rigorously.
- **(d) The split regime above the threshold (§6)**, and selection rules that are not boxes.

**Status of G_r after this note.**
- Inside the scope of R1–R4, the barrier conclusions resbarrier drew from G_r (α ≥ 1/3; no rank-r family below the
  rank-2 floor; no aggregate family) are theorems for every N, with a weaker height floor.
- Outside the scope, (a)–(d) remain. There the missing input is LP or a typical-size heuristic, not a
  number-field counting heuristic.

## 8. Literature (checked by the literature critic)

- **Bean 1994.** Compositio Math. 92, 115–131. Thm 1 covers complex binary forms, D ≥ 3, D_F ≠ 0; Cor. 1 is the
  integer case.
- **Bean–Thunder 1996.** See §2.
- **Thunder 2001.** Ann. Math. 153, 767–804, arXiv math/0105267.
  - Finite type ⇔ N_F(m) < ∞ for all m, and then N_F(m) ≪ m^{n/d}, with constants stated to be explicit.
  - Thm 3 gives the asymptotic; Lemma 14 the box discrepancy.
- **Thunder 2003.** J. Number Theory 101, 294–309, Thm 2: M(F)^{−n/d} ≪ V ≪ M(F)^{−1/d}(1+(log M)^{n−1}).
- **Thunder 2005.** Compositio 141, 271–292 = arXiv math/0305078.
  - Thm 1: V ≥ (2/n)^n m(F)^{−n/d}.
  - Thm 3: the count, under finite type and a′ < D/r.
  - Thm 4.
- **Stewart 2024.** arXiv 2410.01097 (essentially finite type).
  - Eq. (19) quotes Bean–Thunder.
  - Its Cor. 3 misprints the exponent as (d−n)!/n!, and its eq. (7) misquotes Evertse's constant.
- **Evertse 2000.** Publ. Math. Debrecen 56, 337–374, Thm 1, with constant (16r)^{(n+11)³/3}.
- **Schmidt.** Math. Ann. 191 (1971) 1–20 (finiteness criterion); Ann. Math. 96 (1972) 526–551.
- **Evertse–Győry.**
  - *Unit Equations*, CUP 2015, Ch. 9, pp. 231–283.
  - *Discriminant Equations*, CUP 2017 (NMM 32), §14.9.
- **J. Liu.** *On p-adic decomposable form inequalities*, PhD thesis, Leiden/Bordeaux I, 2015.
- **Prior art.** None of these mentions linked modules, NFS or factoring. Neither the literature critic nor codex
  found a prior statement of R1's N-divisibility or of R4. R1 is an immediate corollary of Bean–Thunder plus an
  elementary congruence; no novelty is claimed.

## 9. Experiments (pre-registered, commit 865d205)

### E1: theorem checks (exact; `src/rank3/e1.py`, 1 s)

**Instances.** 30 in total:
- (A) GNFS-like: N = pq with 10–11-bit primes, monic f, g sharing a root mod N.
- (B) SNFS-like: N = the squarefree part of Res(f, g) for small f, g.
- (C) non-power modules, with P_j ≠ Q_j linked only mod N.

Cells (d, d′, r) ∈ {(3,3,3), (4,3,3), (4,4,3), (4,4,4)}. Δ_F was certified exactly (arb, up to 2^17 bits).

| Prediction | Result (re-scored independently by codex and the proof critic) |
|---|---|
| P1: v_p(Δ_F) ≥ P(r,D) | **held**, 30/30. Equality in 29/30; one (C) instance has 48/42 > 36 |
| P2: count bound at every grid X | **held**, 3448/3448, but **vacuous**: the error term exceeds the box population (e.g. 3549 > 2197) |
| P3: R4 (cross-polytope form) | **held**, 30/30; margins 10.4–53.5 nats |
| P4: max saturation ≥ 10^{−3} in ≥ half of (B) | **held**, 12/12 (0.027–0.176). Finite-box diagnostic only |
| P5: median saturation (B) > (A) per cell | **held**, 4/4 (0.143/0.023, 0.132/0.037, 0.114/0.023, 0.029/0.003) |

**Post hoc, not pre-registered** (`src/rank3/e1_posthoc.py`, `experiments/raw/rank3_E1_posthoc.jsonl`).
- The 9 rank-3 (3,3) instances were rerun with E = 40 (531 440 box points, error 137 781 < box). That gives 760 checks
  and **0 violations**.
- 117 of the checks are informative, meaning count > error.
- (count − error)/main reaches 0.054–0.061 for the SNFS-like (B) instances, against ≤ 3·10^{−4} for (A) and (C).
- Also post hoc: the Hadamard form of R4 holds on all 30 instances, and the ramified, small-prime, index-divisor,
  rational-side and r > d tests pass.

**Reading.** The bound with the exact Δ_F is non-vacuous, and SNFS-like families come within about a factor 20 of it
at E = 40. This is a finite-box diagnostic, not a proof of asymptotic sharpness.

### E2: construction comparison (toy sizes; measured, not asymptotic; `src/rank3/e2.py`, 16 s)

**Setup.** n = log₂N ∈ {40, 50, 60}, 3 seeds each, B = 2^12. "Yield" counts the candidates with both norms B-smooth.

| n (median over seeds) | R2-rat | R2-nn | R3-rat | R3-nn | R3-lll | best3/best2 (median) |
|---|---|---|---|---|---|---|
| 40: yield / median bit length of \|F\| | 3818 / 49 | 695 / 64 | 3711 / 57 | 1219 / 66 | 80 / 83 | 0.952 |
| 50 | 2795 / 54 | 536 / 69 | 1451 / 67 | 691 / 77 | 11 / 99 | 0.519 |
| 60 | 1330 / 59 | 206 / 74 | 444 / 77 | 117 / 86 | 0 / 113 | 0.305 |

**Predictions.**
- **P6** (best rank-3 ≤ 2 × best rank-2): **held**, 9/9.
- **P7** (rank-3 medians ≥ best rank-2 median − 2 bits): **held**, 27/27. P7 also holds if "best rank-2" means the
  lower median.
- No lead signal fired.

**Disclosed deviations and qualifications.**
- "Best rank-2" was operationalised as the larger rank-2 yield.
- The boxes have 131 328 (rank 2) and 135 200 (rank 3) points, not equal populations. There are 79 791 and 112 753
  primitive candidates, which favours rank 3.
- "median_log2F" is the median bit length.
- The base m is ⌊N^{1/4}⌋, not the pre-registered ⌈·⌉.
- Primes were generated by seeded `nextprime`, not `randprime`.
- R3-lll is built from the integral-order fibre product {(y, y′) ∈ Z[ω] × Z[ω′] : φ(y) = φ′(y′)}, not from the
  pre-registered diagonal polynomial image.
- In all 9 instances R3-lll's selected generators were identical polynomials on both sides. The search found only
  common basis changes, not P ≠ Q designs.
- R3-lll's weight search uses rounded embeddings, and it uses algebraic-integer norms (Res/lc^{deg}). The other
  designs use Res(f, A).

## 10. Reviews

**Codex adversarial round (codex/rank3_2.out): MAJOR REVISION; all items applied.**
- What it confirmed:
  - R1(i) and R1(ii), including ramified and index-divisor primes (N = 2 index-8 test, v₂ = 168).
  - R2 and R4.
  - The Bean–Thunder hypotheses, from the abstract.
  - The scores, which it reproduced.
- Fixes applied:
  - #1: GP is added to the aggregate exclusion, and the height reading is made explicit via Mahler measure.
  - #2: "best possible" is deleted from R5, and the example is replaced by its GP version.
  - #3: the O3 conclusion is relabelled as model-level under H′, with the transfer not supplied.
  - #4: "rank 2 optimal" and "G_r becomes a theorem" are reworded.
  - #5: the yield-step bookkeeping is repaired (removal, tail truncation, lower endpoints).
  - #6: 1.6803 → 1.67989, and 2.4229 → 2.42283.
  - #7: Thunder 2005's hypotheses are restored.
  - #8: the GP wording and LP are reformulated.
  - #9: special-q is restricted and the error control stated.
  - #10: P2 is labelled vacuous, with the post-hoc non-vacuous check and all E2 deviations disclosed.
  - #11: the setting is fixed (no projection injectivity, tagged union).

**Opus proof critic: ACCEPT WITH FIXES; all applied.**
- What it confirmed:
  - R1 (including ramified primes), R2, R3 (the scoped constants, by independent optimisation) and R4.
  - The prereg commit order.
  - The scores.
  - 82 adversarial instances, with 0 violations.
- Fixes applied:
  - #1 (MAJOR): the setting no longer excludes rational sides or r > d. It verified v_p = 12 at P(3,4) and 144 at
    P(4,6); see `test_rational_side_and_r_greater_than_d`.
  - #2 (MAJOR): GP for §4.2.
  - #3: the Δ_F = 0 split example is replaced.
  - #4: the constant 1.6799.
  - #5: R4 in elementary Hadamard form.
  - #6: wording on G_r and on scope.
  - #7: GP wording, and the semi-discriminant route added to §7(a).
  - #8: the proof wording s → s+1.
  - #9: R2 error rD.
  - #10: special-q errors.
  - #11: the §7(c) range.
  - #12: the E2 deviations.
  - #13: post-hoc small-prime and index-divisor tests.

**Opus literature critic: all statements verified; citation fixes applied.**
- The C_r asymptotic is attributed to Bean–Thunder and zbMATH.
- The sharpness clause.
- Bean Thm 1 / Cor 1.
- Thunder 2003 is split from Thunder 2005.
- Thunder 2005's hypotheses.
- Lemma 14's form.
- Venues for Evertse, Schmidt (both papers), Evertse–Győry and Liu.
- The 2D-optimal folklore is attributed to JLSV 2006 and Barbulescu–Pierrot 2014 §5.2 (DLP analogues). GGV and Grémy
  are cited as counter-context.

**Not done.** No second codex round on the revised note. The revisions only narrow claims, add hypotheses, or add
checks.

## 11. Reproduction

```bash
.venv/bin/python -m pytest -q tests/test_rank3.py      # 18 tests, < 1 s
.venv/bin/python -m src.rank3.e1                        # E1 (pre-registered), ~1 s
.venv/bin/python -m src.rank3.e1_posthoc                # post-hoc non-vacuous count check, ~30 s
.venv/bin/python -m src.rank3.e2                        # E2 (pre-registered), ~16 s
.venv/bin/python -c "from src.rank3.ledger import *; print([float(scoped_constant(r)) for r in (2,3,4)])"
```
