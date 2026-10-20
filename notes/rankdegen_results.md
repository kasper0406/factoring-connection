# rankdegen: degenerate linked modules, large boxes, and the split side (rank3's open items)

Agent `rankdegen`, phase 18, 2026-10-20. Revised after codex rankdegen_2 (MAJOR REVISION) and two Opus critics
(proofs: ACCEPT WITH FIXES; literature: §10). §10 lists every applied fix.

**Files.**
- Code:
  - `src/rankdegen/fh.py`: certified minors, frames, exact LP, FH bounds, exact box counts.
  - `src/rankdegen/ledger.py`: Φ program, constants, split-side bound, H_split program.
  - `src/rankdegen/e{1,2,3}.py`: experiments. `e2.py --corrected` is the corrected D-C1 design (§8).
- Tests: `tests/test_rankdegen.py` (12 tests).
- Pre-registration: `notes/rankdegen_prereg.md`, commit 4c40b18 (before E1–E3 ran).
- Raw data: `experiments/raw/rankdegen_E{1,2,3}.jsonl`, `rankdegen_E2_corrected.jsonl`, `rankdegen_E1_P5.txt`.
- GPT-6.1-sol: blind round `codex/rankdegen_1.{md,out}`; adversarial round `codex/rankdegen_2.{md,out}`.

**Labels.**
- **proved.**
- **proved (model)**: an exact consequence of the stated model (heuristic H, B² linear algebra, single N) plus the
  named scope.
- **heuristic.**
- **KNOWN**: cited.
- **measured.**
- **open.**

No novelty is claimed.

**Notation.** As rank3 §1.
- A linked rank-r module is a family (u_j, u′_j) ∈ O × O′ with φ(u_j) = φ′(u′_j), injective as a pair map.
- Rows: ℓ_σ = (σ(u_j))_j (d of them) and ℓ′_τ (d′ of them); D = d + d′.
- W, W′ are the two norm forms and F = WW′.
- x = log N, with N squarefree.
- A **frame** is a set S of r+1 rows any r of which are independent; Φ_S := ∏_{i∈S}|det L_{S∖i}|.

---

## 0. Verdict

**Plain language.**
- **No degenerate family we can analyse beats rank 2, but degeneracy is not closed.** The degenerate types sort into
  three groups:
  - **Useless:** two sides that coincide in characteristic 0, and decoupled sides whose splitting is integral at N.
  - **Exact reductions:** a subfield module reduces exactly to a smaller non-degenerate one.
  - **The rest:** rank-deficient sides, and accidental dependencies with few rich lines. These obey exactly the same
    lower-bound constants as rank3's general-position ledger.
- **The tool is an elementary volume inequality, "frame–Hölder".** It splits the norm form into (r+1)-row pieces
  and applies Hölder's inequality. It is a weighted, side-aware form of the Hölder step in Bean's 1994 proof, whose
  uniform case is essentially Bean's method (literature critic).
  - It needs no general-position hypothesis.
  - With uniform weights it gives Bean–Thunder's exponent.
  - When both sides have enough rows it is strictly stronger than R1.
- **Still open in task 1: "frame-less" configurations**, in which no r+1 rows are in general position. Codex found
  three kinds:
  - Weil restrictions: modules over a common subfield k ≠ Q that is not imaginary quadratic. They have infinite
    volume because a unit group acts.
  - A finite-volume example with no special structure.
  - Decoupled modules whose splitting is not integral at N.

  Accidental dependencies with many points on a line are also open, unless an explicit LP is solved for each
  instance. We found no exploitable family among any of these, but we also have no ledger for them.
- **Large boxes: the naive count is false, and present tools cannot give the corrected one.**
  - Points in excess of the volume bound really occur at the L[1/3] scale (codex's explicit example).
  - Wherever we can analyse them, they lie on lower-rank sub-families. These are sublattices, or products of
    lower-rank candidates whose relations add no rank.
  - A precisely stated "lattice points with exceptions" hypothesis (LP_exc, §6.4) would close large boxes.
  - The standard counting theorems do not supply it. Their error terms ignore either N or the irrationality of the
    zero hyperplanes.
- **ptesplit O3 (large split side).**
  - The frame inequality charges the other side up to split-side size (r−1)x/d′, provided d′ ≥ r and d ≥ 2. R5
    stopped at (r−1)x/(D−1). An explicit special-N construction reaches this threshold in the size model.
  - Above the threshold, linking cannot charge the other side.
  - Only a split-side counting hypothesis excludes the regime: that split-by-construction candidates are as rare as
    smooth ones.
  - Under that hypothesis O3 is smoothenum's one-sided output-sensitive listing lever (no prior statement found), still at L[1/3]. The formal
    constant is (8(r−1)/(9r))^{1/3}, i.e. 0.7631 at r = 2. Inside the box scope it gains nothing.
  - Going below L[1/3] there would need a generator of split elements that beats smoothness statistics. That is
    ptesplit's open O1, not a geometric loophole.
- **Experiments (pre-registered).**
  - The frame bound held at all 259 194 exact checks (certified), and the linking exponent at 30/30.
  - The best degenerate design reached only 0.05–0.17 of the best rank-2 yield at n = 40–60 bits, falling with n.
  - One design had to be corrected and rerun (disclosed, §8).

**Where the L[1/3] limit bites, after this note.** The setting is relation methods in the linked-order model, for a
single N, inside the box scope.
- **Proved:** the limit sits at the frame–Hölder volume bound, for general-position modules and for every
  degenerate type with a frame-coverable configuration.
- **Not proved:**
  - for frame-less configurations there is neither a counterexample nor a theorem;
  - beyond the box scope the limit rests on the unproved lattice-point statement LP_exc.

**Precisely.**
1. **Theorem S (structure) [proved].**
   - Within a side, two rows are proportional exactly when their embeddings agree on the subfield
     K₀ = Q(u_j/u_{j₀}).
   - Across sides it is all-or-nothing: the two projective configurations either coincide or are disjoint.
   - Instability can come only from a whole side of deficient Q-rank.
2. **Lemma CC [proved, extraction form of Lemma ND].** Coincident sides (type A) give only a finite set of splits,
   computable without relations.
3. **Lemma SUB [proved].** Subfield modules (type B) reduce exactly to the subfield.
4. **Theorem FH [proved].** No GP hypothesis is needed.
   - J_y ≤ C_r ∏_S Φ_S^{−w_S/(r+1)} for any frame distribution w with prescribed marginals. This uses Bean–Thunder
     per frame; the self-contained route gives 2^r C_r.
   - With unramified linking, J_y ≤ C_r N^{−((r−1)/(r+1))t(w)}.
5. **Corollary FL [proved (model)]: a two-norm ledger.**
   - **Proposition PHI [proved]:** the constant is Φ ≥ (r−1)/r whenever PHI's one-sided hypotheses hold.
   - So the GP, C1 and low-loss C2 floors equal R3's: 1.5263 at r = 2; 1.6799 formal and 2.4228 scoped at r = 3.
   - Two nonlinear sides of equal degree pay more: Φ_bal = 2/3 at r = 2, floor 1.6799.
6. **Lemma PR [proved].** A projection/fibre bound via Pólya's lemma. It closes the one spurious optimum of
   frames-only bounds.
7. **Large boxes.** The negative statements are proved; the positive statement LP_exc is open (§6).
8. **Corollary SPL [proved (model)]** and the O3 verdict (§7).

---

## 1. Setting and what was open

rank3 proved V_F ≤ C_r N^{−r(r−1)/(D(D−1))} when Δ_F ≠ 0 (general position, GP), and the ledger in boxes of side
≤ B^{1/(r−1)}. Codex showed that GP cannot be dropped as stated (rank3 §2).

Open:
- (a) Δ_F = 0.
- (b) Boxes larger than B^{1/(r−1)}.
- (d) The split side above (r−1)x/(D−1).

## 2. Structure of degenerate modules

**Theorem S [proved].** Fix j₀ with u_{j₀} ≠ 0. Put v_j = u_j/u_{j₀} and K₀ = Q(v₁,…,v_r) ⊆ K; K₀ does not depend
on j₀. Define v′_j and K′₀ likewise. Let P_σ = [ℓ_σ] ∈ P^{r−1}(C).
- **(S1)** P_σ = P_{σ′} ⇔ σ|_{K₀} = σ′|_{K₀}.
  - Each point has multiplicity [K:K₀].
  - W(a) = N_K(u_{j₀})·N_{K₀}(Σ a_j v_j)^{[K:K₀]}.
  - *Proof.* ℓ_σ = σ(u_{j₀})(σ(v_j))_j with σ(u_{j₀}) ≠ 0, so proportionality means σ(v_j) = σ′(v_j) for every j,
    i.e. the two embeddings agree on Q(v). Embeddings of K₀ extend in [K:K₀] ways, and N_{K/Q} = N_{K₀}^{[K:K₀]}
    on K₀. ∎
- **(S2)** Suppose some P_σ = P′_τ.
  - Then the zero patterns of (u_j) and (u′_j) agree, so a common j₀ exists.
  - There is a field isomorphism ι: K₀ → K′₀ with ι(v_j) = v′_j, and {P_σ} = {P′_τ}.
  - Otherwise the two point sets are disjoint.
  - *Proof.* Take ι = τ^{−1}σ on Q(v). Any other σ″ is gσ on K₀ for some automorphism g of C, and then
    P_{σ″} = [gσ(v_j)] = [gτ(v′_j)] = P′_{gτ}. ∎
- **(S3)** dim_C span{ℓ_σ} = ρ := dim_Q span_Q(u₁,…,u_r).
- **(S4)** Let t* be the largest ratio |A|/rk(A) over nonempty sets of rows, counted with multiplicity.
  - This is standard: supermodularity, cf. Mumford's semistability criterion (MFK, GIT Ch. 3) and Kempf's optimal
    destabilising flag (Ann. of Math. 108 (1978)).
  - f(A) = |A| − t*·rk(A) is supermodular and f ≤ 0 = f(∅). So its maximisers are closed under ∪ and ∩, and there is
    a unique maximal maximiser A*.
  - A* is Galois-stable, hence a union of sides.
  - If t* > D/r then A* ≠ all rows, so **A* is one whole side, with d/ρ > D/r**: instability comes only from a
    rank-deficient side.
  - Strict semistability (t* = D/r) can come from flats that are not Galois-stable.

**Classification of Δ_F = 0.** The classes are defined so that they are exhaustive: the last one is a catch-all.
- **A (coincident):** {P_σ} = {P′_τ}. Linking holds in characteristic 0 on K₀ ≅ K′₀.
- **B (subfield):** K₀ ⊊ K or K′₀ ⊊ K′, so rows repeat. Reduce first (Lemma SUB); the classes below refer to the
  reduced module.
- **C1 (rank-deficient side):** a side with ρ < r and more than ρ distinct rows.
- **Decoupled:** the two row spans are complementary, ρ + ρ′ = r.
- **C2 (the rest):** the rows are distinct and some r of them are dependent.
  - This includes sides with ρ < r and exactly ρ rows, combined with cross dependencies. Codex's mixed example:
    K = Q(√2), K′ = Q(√2, √3), u = (1, s, 0), u′ = (1, s, t − 7).
  - It also includes accidental dependencies. Codex's example: f = x⁶ − 2 with modules (1, α, α³) and
    (1, α+N, α³). The 12 rows lie on two lines {z = ±√2·x} in P², six on each.
  - C2 splits into **frame-coverable** configurations, where FH applies, and **frame-less** ones.
  - Known frame-less C2 configurations:
    - **D (relative modules).** The whole module is a module over a common subfield k ≠ Q that is not imaginary
      quadratic. Codex's example: s = √2, z = √3, with modules (1, s, z, sz) and (1, s, z+N, s(z+N)) at r = 4. The
      units of Z[√2] act, so V_F = ∞, and there are no frames.
    - **A finite-volume frame-less example (codex rankdegen_2 #3).** The modules (1, s, t, t²) and (1, s+23, t, t²)
      with s = √2, t = ∛3, r = 4. It has 12 rows in three groups of rank 2 that share a direction, 288 nonzero
      4-minors and zero frames, and its volume is finite. It is not of type D.
  - Schmidt degeneracy means a *nonzero* subspace V^J stable under such a subfield J; type D is the case V^J = V.
    - When V^J ⊊ V it is a sub-family phenomenon, treated in §6.
    - References: Schmidt, Math. Ann. 191 (1971) 1–20, Satz 2; Ann. of Math. 96 (1972) 526–551. Both concern one norm
      form.
    - For products WW′, the relevant criteria are Thunder 2001's "finite type" and Evertse–Győry, Acta Arith. 50
      (1988) 357–379.

**Decoupled modules [proved, corrected after codex #1].** Let X = π₁(Z^r) and Y = π₂(Z^r) be the projections to the
two spans. Let M ⊆ X × Y be the image of Z^r, and put I = [X × Y : M].
- **If gcd(I, N) = 1**, every candidate's K-element lies in ker φ and its K′-element in ker φ′. Hence |W|, |W′| ≥ N,
  and the module is useless.
  - *Proof.* For every x ∈ X, (Ix, 0) ∈ M, so I·φ(θ̄(x)) = 0. Symmetrically for Y.
- **If 1 < gcd(I, N) < N**, then N splits.
- **Otherwise rad(N) | I**, so I ≥ N, and the module is a *fibre product* whose candidates satisfy congruences mod N.
  - Codex's example: u = (1, α, 0, 0), u′ = (1, 128, N, Nα), with W(e₁) = W′(e₁) = 1, has I = N².
  - If counts follow volume/I, its ledger pays the full x: ℓ ≥ min(d/ρ, d′/ρ′)(λ + x). That is an L[1/2]-type cost
    [heuristic]. It is frame-less, so there is no theorem.

## 3. Types A and B

**Lemma CC [proved, in the extraction form of resbarrier Lemma ND].**

*Hypotheses.*
- **The module is of type A**, with j₀ chosen so that c = φ(u_{j₀}) = φ′(u′_{j₀}) is a unit mod N.
  - If no j has φ(u_j) a unit, then either some gcd(φ(u_j), N) splits N, or all φ(u_j) ≡ 0 and the module is
    decoupled.
- **The maps extend to the square roots.** φ and φ′ extend to overorders containing the square roots used. As in ND,
  the relevant indices are coprime to N: a proper gcd splits N, and gcd = N is excluded.
- **The square roots are units.** The square roots √e below are units at the linking primes.

*Statement.*
- Then φ|_{O_{K₀}} = φ′∘ι|_{O_{K₀}}, by ND (i).
- Let C = ker(K₀^×/K₀^{×2} → K^×/K^{×2}) ∩ ι^{−1}ker(K′₀^×/K′₀^{×2} → K′^×/K′^{×2}).
- Let T = {[e] : u_{j₀}e ∈ K^{×2}, u′_{j₀}ι(e) ∈ K′^{×2}}. T is empty or a coset of C.
- For a relation set S with ∏θ(a) = δ² and ∏θ′(a) = δ′², the ratio φ(δ)/φ′(δ′) lies in a finite set computable
  without relations:
  - **|S| even:** {±φ(√e)/φ′(√ι e) : e ∈ C};
  - **|S| odd:** {±φ(√(u_{j₀}e))/φ′(√(u′_{j₀}ι e)) : e ∈ T}.

*Proof.*
- Write θ(a) = u_{j₀}γ_a and θ′(a) = u′_{j₀}ι(γ_a) with γ_a ∈ K₀.
- For |S| even, Γ = ∏γ_a is a square in K and ι(Γ) is a square in K′. So Γ = eη² with e ∈ C and η ∈ K₀.
- The u-parts and η-parts cancel under the agreeing maps.
- For |S| odd, apply the same argument to u_{j₀}Γ. ∎

Codex's identical-maps example (Φ₅, the same module on both sides) has relations, and all of them are of this
kind.

**Lemma SUB [proved].**
- *Form of the norms.* In type B, W = N_K(u_{j₀})·G^{[K:K₀]} with G = N_{K₀}(Σ a_j v_j), and likewise for W′.
- *Hypothesis.* Let c, c′ clear the denominators of the v_j and v′_j. Assume
  gcd(c c′ N_K(u_{j₀}) N_{K′}(u′_{j₀}), N) = 1. Otherwise N splits, or one of these numbers is divisible by N, which
  is a case we keep separate.
- *Smoothness.* Then a candidate's norms are smooth iff G and G′ are, up to the fixed scalar factors, which form a
  single fixed event.
- *Linking.* The reduced module (c v_j, c′ v′_j) is projectively linked: its rows at a linking prime are congruent
  up to a 𝔓-unit scalar. R1 and FH apply unchanged, after subtracting a scalar multiple of one row.
- *Conclusion.* Type B is therefore a module of total degree D₀ = [K₀:Q] + [K′₀:Q] with the same candidates. This is
  codex's bound (1) in `rankdegen_1.out`.

## 4. Theorem FH (frame–Hölder)

**Theorem FH [proved].** Let L₁..L_D be complex linear forms on R^r, y_i ≥ 0 with Σ y_i = r, and w a probability
distribution on frames with P_w(i ∈ S) = (r+1)y_i/r. Then

J_y := (1/r)∫_{S^{r−1}}∏_i |L_i(θ)|^{−y_i} dθ ≤ 2^r C_r ∏_S Φ_S^{−w_S/(r+1)},

and vol{x ∈ R^r : ∏|L_i(x)|^{y_i} ≤ T} = T·J_y.

**Improvement to C_r.** The constant 2^r C_r may be replaced by C_r:
- for real L_i, by step 2 below;
- in general, by applying Bean–Thunder's theorem to each frame.
  - Here d = r+1 > n = r, and |Δ_{F_S}|^{1/(r+1)!} = Φ_S^{1/(r+1)}.
  - The theorem is stated for complex coefficients: Zbl 0854.11019 (review by Bean); Bean, Compositio Math. 92 (1994)
    115–131, Thm 1, for r = 2.
- The proof below is the self-contained route, which loses the factor 2^r.

*Proof.*
1. **Hölder.**
   - The marginal condition gives ∏_i|L_i|^{−y_i} = ∏_S (∏_{i∈S}|L_i|^{−r/(r+1)})^{w_S}.
   - Generalised Hölder gives J_y ≤ ∏_S J_S^{w_S}, with J_S = (1/r)∫|F_S|^{−r/(r+1)} = vol{|F_S| ≤ 1}.
2. **Real frames, exact.**
   - Let S consist of real forms in GP. Write L₀ = Σ c_iL_i with every c_i ≠ 0, and substitute y_i = c_iL_i(x).
   - F_S becomes (∏c_i)^{−1}y₁⋯y_r(y₁+⋯+y_r). Hence
     vol{|F_S| ≤ 1} = V₀·|det L_{1..r}|^{−1}∏|c_i|^{−1/(r+1)} = V₀·Φ_S^{−1/(r+1)}.
   - The last step uses Φ_S = |det L_{1..r}|^{r+1}∏|c_i|, from Cramer's rule (`test_real_frame_identity`).
   - V₀ = vol{|x₁⋯x_r(x₁+⋯+x_r)| ≤ 1} = C_r is Bean–Thunder's extremal constant (Stewart, arXiv 2410.01097,
     eq. (19)). C₂ = 15.8997 and C₃ = 213.71.
3. **Complex frames.**
   - For phases φ_j put R_j = Re(e^{iφ_j}L_j). On R^r, |R_j| ≤ |L_j|, so vol{|F_S| ≤ 1} ≤ vol{|∏R_j| ≤ 1}.
   - With w_j = e^{2iφ_j}, det R_{S∖i} = (unimodular)·p_i(w). Here p_i is multilinear, with coefficient
     2^{−r}det L_{S∖i} at ∏w_j (`test_complex_to_real_top_coefficient`).
   - Iterated Jensen (∫log|a + bz| ≥ log|b|) gives ∫_T log|p_i| ≥ log(2^{−r}|det L_{S∖i}|).
   - Sum over i ∈ S. The average of the sum is finite, so some φ attains
     Σ_i log|det R_{S∖i}| ≥ Σ_i log(2^{−r}|det L_{S∖i}|).
   - For that φ the R_j are real and in GP, and step 2 gives J_S ≤ 2^r C_r Φ_S^{−1/(r+1)}.
4. Homogeneity of degree r gives the T-scaling. ∎

**Remarks.**
- **Uniform w reproduces Bean–Thunder's exponent.** In GP each r-subset lies in D − r frames, so uniform w over all
  frames gives Bean–Thunder's exponent.
  - For totally real F this is Bean–Thunder's inequality with their constant.
  - In the proof critic's r = 3 checks, FH(uniform)/2³ equals Bean–Thunder's bound numerically.
- **The factor 2^r.** Bean–Thunder holds for complex coefficients (Zbl 0854.11019; Bean 1994 Thm 1). Applied per
  frame, it removes the 2^r, so step 3 is an optional self-contained route. At r = 2 the lossless counterpart of
  step 3 is Bean's Thm 3, a maximum-principle argument credited to Bombieri.
- **Prior art.** Bean 1994 §4 applies generalised Hölder to the drop-one subforms; iterated, this is uniform Hölder
  over (r+1)-subsets at r = 2. Bean–Thunder very likely use the same induction (full text not seen). So uniform FH is
  Bean's method. What is not in Bean: non-uniform and side-aware weights, and the linking exponent t(w).
- **Context.** With Σy_i = r and real rows, J_y < ∞ iff y lies in the relative interior of Edmonds' base polytope of
  the row matroid. Frame points (r/(r+1))1_S lie in that interior, and FH interpolates between them. This parallels
  Barthe's rank-one Brascamp–Lieb theory (Invent. Math. 134 (1998)). Codex's frame-less finite-volume example shows
  that the hull of the frame points is strictly smaller than the interior.
- **Proof-critic numerics.**
  - r = 2: quadrature on 300 random mixed real/complex configurations gives J_y/FH ≤ 0.27.
  - Complex single frames give at most 0.70 of C₂Φ_S^{−1/3}.
  - r = 3: J = 10.04 against FH(uniform) = 356.8.

**Corollary FN (linking) [proved].**

*Hypotheses.*
- N is squarefree.
- At each p | N the linking primes are unramified: e(𝔭_K|p) = e(𝔭′_{K′}|p) = 1, in R1's notation.
- w is invariant under the side-preserving automorphisms of the row configuration. Example: uniform on the frames
  with exactly s K-rows.

*Statement.* Let t(w) be the minimum over cross pairs (σ, τ) of P_w({σ, τ} ⊂ S). Then
Σ_S w_S log Φ_S ≥ (r−1)t(w)log N, and hence J_y ≤ 2^r C_r N^{−((r−1)/(r+1))t(w)}.

*Proof.*
- For a Galois orbit O of frames, Ψ_O := ∏_{S∈O}∏_{i∈S}det(L_{S∖i})² is a Galois-invariant nonzero algebraic
  integer. So Ψ_O ∈ Z∖{0}, and v_p(Ψ_O) = v_𝔓(Ψ_O) for any 𝔓 | p.
- As in R1, at 𝔓 a K-row ℓ_{σ₀} and a K′-row ℓ′_{τ₀} are both congruent to the integer lift ṽ. So every minor
  containing both has v_𝔓 ≥ 1.
- A frame containing {σ₀, τ₀} has r − 1 such minors S∖i.
- Weight by w_O (w is constant on orbits) and sum. ∎

*Notes.*
- **Ramified linking primes.** The same argument gives at least 1/max(e, e′) per minor. R1's inclusion–exclusion
  recovers the full exponent only in GP. If some but not all p | N ramify, gcd(N, disc O·disc O′) already splits N.
- **Exact LP values are valid exponents.** Averaging an LP optimum over the symmetry group gives an optimal
  invariant w whenever the prescribed marginals are invariant. So the exact LP value max_w t(w)
  (`fh.maxmin_coverage`) is a valid exponent without knowing the Galois group.
- **When FN beats R1** (unramified case). The gain comes from side-balanced frames. It needs the hypergeometric law of
  the number of K-rows in a uniform frame to have at least three support points:
  min(r+1, d) − max(0, r+1−d′) ≥ 2.
  - A GP pair of cubics at r = 3 has t = 4/9, against 2/5 for uniform w (E1, all seeds).
  - There is no gain at r = 2 with d = d′ = 2, nor with a rational side. The latter is consistent with the SNFS
    saturation seen in rank3 E1.

## 5. Ledgers per type

**Corollary FL (two-norm frame ledger) [proved (model)].**
- *Weights.* Fix a Galois-invariant w. Galois acts transitively on each side, so w has constant marginals per side.
  The per-row weights y_K, y_{K′} then satisfy ∏|ℓ_σ|^{y_K} = |W|^{y_K}.
- *Count in a box of side 2E.*
  #{a : 0 < |W| ≤ Y₁, 0 < |W′| ≤ Y₂} ≤ J_y Y₁^{y_K} Y₂^{y_{K′}} + r(d+d′)(2E+1)^{r−1}.
  Each line meets the set in at most d + d′ − 1 intervals; this is Lemma R2's argument.
- *Yield step.* Use 2-dimensional dyadic bins, poly(x)² of them. Require the lattice error to be at most
  B/(4·#bins), which in L-units is still H5_r.
- *Constraints.* Every symmetric w then gives
  y_K ℓ₁ + y_{K′} ℓ₂ ≥ λ + ((r−1)/(r+1))t(w)x − log(2^r C_r), with ℓ_i = log Y_i.
  All these constraints hold simultaneously in the yielding bin.

**The Φ program.**
- *Setting.* The configuration has **complete classes**: every admissible (s, r+1−s)-subset is a frame. Use pure-s
  distributions, uniform on the frames with s K-rows.
- *The constraint, in raw units.*
  (r/(r+1))(s ℓ₁/d + (r+1−s) ℓ₂/d′) ≥ λ + ((r−1)/(r+1)) s(r+1−s) x/(dd′).
- *Rescaled.* Put q₁ = ℓ₁/(dλ), q₂ = ℓ₂/(d′λ) and R = λdd′/x. The constraint becomes
  (r/(r+1))(s q₁ + (r+1−s) q₂) − 1 ≥ ((r−1)/(r+1)) s(r+1−s)/R.
- *Φ.* Let Φ be the infimum of q₁q₂R over the constraint set. Then **ℓ₁ℓ₂ ≥ Φλx**, hence **ℓ₁ + ℓ₂ ≥ 2√(Φλx)**.
  The constraints are linear in w, so pure s suffice for complete classes.
- *Constants.* In rank3's L[1/3] units, the size budget optimises to 9β²(γ−β)² ≥ 4γΦ, which is H. So:
  - formal: **c³ = 64Φ/9**;
  - scoped (γ ≤ rβ/(r−1)): **c = 2(4r(r−1)Φ/9)^{1/3}**.

  rank3's R3 is Φ = (r−1)/r.
- **Scope of the reduction (proof critic #1).**
  - Φ > 0 needs the one-sided constraints in PHI. Each side must have ≥ r+1 rows in GP (supplying s = 0 or
    s = r+1 frames), or be covered by Lemma PR.
  - For a side of degree ≤ r with neither, the program's infimum is 0. This is spurious: the optimum needs that
    side's degree to grow.
  - For GP these cases are covered directly by R3: ℓ ≥ (r−1)x/(D−1) + (D/r)λ ≥ 2√((r−1)xλ/r) for any degrees.
  - For C1 with d′ = r at r = 3, the critic minimised ℓ₁ + ℓ₂ in raw units subject to PR, s = 1 and s = 2. It again
    found 2√(2λx/3), i.e. Φ = (r−1)/r.
  - Other bounded-degree C1 cases are not treated.

**Proposition PHI [proved].** If the admissible s form a range [s_lo, s_hi] ⊆ [1, r], and the constraints include
s_hi q₁ ≥ 1 and (r+1−s_lo) q₂ ≥ 1, then Φ ≥ (r−1)/r.

*Proof.*
- For each s, the bound q₁q₂R ≥ (r−1)/r follows from constraint s iff (sq₁ − 1)((r+1−s)q₂ − 1) ≥ −1/r. This is an
  algebraic identity (checked symbolically by the critic).
- Suppose the condition fails for every s. The endpoint constraints then give a consecutive pair s*, s*+1 with
  s*q₁ < 1 < (s*+1)q₁ and (r+1−s*)q₂ > 1 > (r−s*)q₂.
- Put α = 1 − s*q₁, α′ = (s*+1)q₁ − 1, β = (r+1−s*)q₂ − 1 and β′ = 1 − (r−s*)q₂. Failure gives αβ > 1/r and
  α′β′ > 1/r.
- Cauchy–Schwarz gives √(q₁q₂) ≥ √(αβ) + √(α′β′) > 2/√r.
- But q₁q₂ < 1/(s*(r−s*)) ≤ 1/(r−1). That forces 4(r−1) < r, which is false for r ≥ 2. ∎

*Checks and equality.*
- `test_phi_crossing_lemma_random` checks the proof; the critic's random search over r = 2..7 finds a ratio ≥ 1.017.
- The value (r−1)/r is approached as q₂ → ∞ in the large-degree program. That limit is not a genuine rational side,
  for which the program degenerates.
- Pre-registered P5 found (r−1)/r numerically for 9/9 types. The r = 4 (2, 3) case should use s = 2 only; with that
  correction it gives 3/4 (codex #5).

**Lemma PR (projection/fibre bound for a deficient side) [proved].**
- *Setting.* The K side has Q-rank ρ = r − 1, with kernel Λ₀ = Z z₀ and θ′(z₀) ≠ 0. (If θ′(z₀) = 0, both norms
  factor through Z^r/Zz₀ and the module has rank r − 1.)
- *The kernel element.* θ′(z₀) ∈ ker φ′, so |W′(z₀)| = [O′ : θ′(z₀)O′] ≥ N.
- *Fibres.* On each fibre g + Zz₀, W′(g + kz₀) is a real polynomial in k of degree d′ with leading coefficient
  W′(z₀). By Pólya's lemma (Pólya, Sitzungsber. Preuss. Akad. Wiss. Berlin 1928; Borwein–Erdélyi, GTM 161, Remez-type
  inequalities; sharp for Chebyshev polynomials) the fibre has at most d′ + 4(Y₂/(2N))^{1/d′} points with |W′| ≤ Y₂.
- *The count.* For Y₂ < 2N:
  count ≤ (d′+4)·#{g ∈ π(box) : |W̄(g)| ≤ Y₁} ≤ (d′+4)(V_{W̄}Y₁^{ρ/d} + E_π).
  - E_π is Davenport's error for the projected box π([−E,E]^r). It depends on z₀ and on the lift basis, and can be
    large when z₀ has large coordinates. **Its smallness is an extra scope condition**, not implied by H5_r.
  - V_{W̄} ≤ 2^ρC_ρ requires the K side's rows to be frame-coverable within their span, uniformly on its classes. It
    uses FN's orbit integrality with N replaced by 1.
- *In the program* this is the constraint ρ q₁ ≥ 1 (`test_kernel_norm_and_polya_fiber`).
- *Deficiency ≥ 2.* The fibre bound is then Davenport's, so it holds only inside the box scope.
- *Why PR is needed.* Without it, frames alone admit a spurious optimum: a deficient side of huge degree and zero
  size, at which the program collapses to c → 0.

**Accidental dependencies (C2), corrected [proved (model); proof critic #2].**
- *Marginals.* Galois symmetry makes the uniform distribution on class-s frames have marginals exactly s/d and
  (r+1−s)/d′.
- *Loss.* Let ε′ be the largest fraction, over cross pairs and the classes used, of class-s subsets containing the
  pair that fail to be frames. Assume every class used, including the one-sided ones, contains a frame.
- *Bound.* Then t_s ≥ (1−ε′)s(r+1−s)/(dd′), and rescaling R gives **Φ ≥ (1−ε′)(r−1)/r**.
- *Consequences.*
  - The C2 floor stays at least the rank-2 floor 1.5263 if ε′ ≤ 1/4 (formal, r = 3) or ε′ ≤ 3/4 (scoped, r = 3).
  - It equals the GP constant only as ε′ → 0.
  - For r = 3, ε′ = O(k/min(d, d′)), with k the largest number of rows on a line.
- *Rich lines.* Configurations with rich lines need the LP for each instance. For α⁶ = 2 it gives t = 1/9, against
  GP's 1/11, with uniform marginals.

**Ledger table** (formal = lattice error neglected; scoped = H5_r).

| type | Φ (proved lower bound) | formal c | scoped c (r = 3) | status |
|---|---|---|---|---|
| A coincident | — | — | — | useless (Lemma CC) |
| B subfield | as the reduced module (degree D₀) | ≥ R3 | ≥ R3 | proved (reduction) |
| GP, any degrees | (r−1)/r | (64(r−1)/(9r))^{1/3} | 2.4228 | proved (model) = R3 |
| GP, two sides of equal degree | 2/3 (r = 2), 8/9 (r = 3); (r²−1)/r² (r odd), (r−1)(r+2)/(r(r+1)) (r even) | 1.6799 (r = 2), 1.8490 (r = 3) | — | proved (model) |
| C1, one side of rank r−1, other side ≥ r+1 rows in GP, complete classes | (r−1)/r | as GP | as GP + projected-error condition | proved (model), FH + PR |
| C1, both sides of rank r−1, complete classes | (r−1)/r | as GP | as GP + projected-error conditions | proved (model) |
| C1, d′ = r (r = 3) | (r−1)/r | as GP | — | raw-unit check (proof critic) |
| C1, deficiency ≥ 2 | (r−1)/r | — | scoped only | proved (model), in scope |
| C2, loss ε′ | ≥ (1−ε′)(r−1)/r | ≥ 1.5263 if ε′ ≤ 1/4 (r = 3) | ≥ 1.5263 if ε′ ≤ 3/4 | proved (model) |
| C2 rich lines | per-instance LP (α⁶ = 2: t = 1/9) | | | computed |
| decoupled, gcd(I, N) = 1 | — | — | — | useless (proved) |
| decoupled, rad N \| I | — | — | — | L[1/2]-type under a volume heuristic; open |
| D relative; codex's finite-volume frame-less example; other frame-less C2 | — | — | — | **open** |
| rational side + rank-deficient K (frame-less) | exact reduction to the rank-(r−1) restriction when E_Γ‖v_Γ‖₁ + Y₂ < N (a = g + kz₀, R = v_Γ·g + kM, M ≡ 0 mod N, so only k = 0 survives) | | | proved in that box range; otherwise §6 |

**Balanced degrees.**
- With δ₁ = δ₂ the budget is √(γR)(q₁+q₂), so Φ_bal = inf (q₁+q₂)²R/4.
- Averaging the two middle constraints gives R ≥ c_r M/(ru/2 − 1), with u = q₁ + q₂ and M = max_s s(r+1−s). This is
  minimised at u = 4/r [proved].
- The bound is attained at q₁ = q₂ = 2/r, where every constraint holds, so it is exact (critic).
- So a single pair of two *equal-degree nonlinear* polynomials has floor 1.6799 in this model, above SNFS's 1.5263.

**Conclusion for task 1.**
- **Useless or exact reductions:**
  - types A and B;
  - decoupled modules with gcd(I, N) = 1.
- **Exactly the GP constants** ((r−1)/r, increasing in r, never below the rank-2 floor):
  - C1 under the stated conditions;
  - C2 in the limit ε′ → 0.
- **Open:**
  - frame-less configurations: type D, codex's finite-volume example, decoupled fibre products;
  - C2 with large losses, except by a per-instance LP;
  - bounded-degree C1 sides beyond the r = 3, d′ = r check.
- **No beater found.** Neither the theory nor the experiments (§8) produced a degenerate family that beats rank 2.
- **Type D, conjecturally.**
  - We expect it to behave as a module of rank r/[k:Q] and degree D/[k:Q], with the k-units contributing
    exp(o(log B)).
  - Its linking needs a ring map O_k → Z/N, for example √2 mod N. For general N that is as hard as factoring.

## 6. Boxes larger than B^{1/(r−1)} (task 2)

**6.1 The volume-only statement is false at the L[1/3] scale [proved; codex rankdegen_1].**
- *Construction.* Take two Eisenstein polynomials f and g = f + 4(X − m), with f(m) = g(m) = N, and the power bases
  (1, α, Tα²) and (1, α′, Tα′²). Both sides are nonlinear and the module is in GP.
- *Volume.* V_{F_T} = V_{F_1}/T. With T = N^{O(1)}, V_{F_T}X^{3/D} < 1.
- *Points.* Yet in a box of side B^{2/3} the plane c = 0 carries ≫ B^{4/3} primitive candidates with
  log|F| = O(x^{2/3}(log x)^{1/3}).
- *Where the excess sits.* On the rank-2 sublattice c = 0, whose binary restriction is in GP and obeys the rank-2
  ledger.
- *Even unscaled, the excess is of bulk order [proved (model), arithmetic].* At rank3's formal rank-3 optimum
  (log X = 4x/(D−1)), the c = 0 plane's main term N^{−2/(D(D−1))}X^{2/D} equals the rank-3 main term
  N^{−6/(D(D−1))}X^{3/D}. So a rank-2 sublattice can contribute an exception of the same order as the bulk.

**6.2 Exceptions need not lie on sublattices.**
- *Products.* Products A = L₁L₂ of two linear candidates (rank 3, power basis) satisfy F(A) = F(L₁)F(L₂) exactly
  [proved].
- *Shape.* Their coefficient map is bilinear in four parameters. The image is an arithmetic subset of coefficient
  space, {b² − 4ac = z²}, not a sublattice.
- *Size.* If counts follow volumes, their number up to X is ≍ A₂²X^{2/D} log X. That is below the rank-3 main term
  exactly when log X > 2x/(D−1), so at all yield-relevant thresholds [heuristic].
- *Rank.* They are harmless for relation rank: their relation vectors lie in the span of their constituents'
  vectors (codex rank3_1) [proved].

**6.3 Why these estimates do not supply the needed uniform bound [proved, narrowed after codex #9].**
- **Davenport and Barroero–Widmer, and Thunder's Lemma 14.**
  - Davenport (J. London Math. Soc. 26 (1951) 179–183; corrigendum 39 (1964) 580) and Barroero–Widmer (IMRN 2014,
    no. 18, 4932–4957, Thm 1.3) bound the discrepancy by projection volumes.
  - Thunder 2001 Lemma 14 is the uniform line-by-line bound dn(2B₀+1)^{n−1}, and has the same defect.
  - If some factor L_i is real, then for every X > 0 the set {|F| ≤ X} ∩ [−E,E]^r contains a full-width slab
    {|L_i| ≤ X(rEH)^{−(D−1)}}.
  - That slab projects to a coordinate hyperplane with volume ≍ E^{r−1}. So the error term is ≍ E^{r−1} > B beyond
    the scope, whatever N is.
  - Complex factors have real zero sets of codimension 2, where the analogous term is smaller; we did not pursue this.
  - These estimates do not see the irrationality of the hyperplanes, which is what keeps lattice points away from
    them.
- **Thunder 2005, Thm 3** (finite type, a′(F) < D/r): N_F(X) ≪_{r,D} (X/m(F))^{r/D} + X^{(r−1)/D}.
  - The second term is uniform, and ignores the N-saving of lower-rank restrictions.
  - At the formal rank-3 optimum it is X^{2/D} = L[1/3, 8β/3] > B (codex).
  - The genuine rank-2 contribution is N^{−2/(D(D−1))} times smaller.
- **Evertse, Publ. Math. Debrecen 56 (2000) 337–374.**
  - **Thm 1.** For a non-degenerate norm form of degree D in r variables,
    #{x ∈ Z^r : |F(x)| ≤ h} ≤ (16D)^{(r+11)³/3}h^{(r+Σ_{m=2}^{r−1}1/m)/D}(1+log h)^{r(r−1)/2}.
    - It counts all integer solutions.
    - It applies to W or W′ alone, not to WW′.
    - It is uniform in the form, so it has no N-saving, and its exponent exceeds r/D.
  - **Thm 2** bounds the number of solution subspaces. That bound decreases with the invariant height H*(F), but it
    still concerns a single norm form.
  - The "(r+7)³/3" quoted elsewhere is a misquotation, found in Thunder 2001's introduction.
- **Effective Diophantine input for algebraic hyperplanes.**
  - Liouville gives exponent d − 1, which is useless.
  - Schmidt's subspace theorem gives r − 1 + ε, but it is ineffective in the exceptional subspaces and not uniform
    in D.

**6.4 What would suffice (corrected after codex #8) [proved conditional].**

**LP_exc.** For every box, threshold X and examined set, the points with 0 < |F| ≤ X split into three parts:
- **(i) Bulk.** At most exp(o(log B))·J X^{r/D} points.
  - The loss is measured on the algorithm's own cost scale. An L[1/3, o(1)] allowance would be too coarse to exclude
    exponents below 1/3.
- **(ii) Linear exceptions.** Points on exp(o(log B)) proper sublattices.
  - Each induced module, in its induced box, must satisfy the hypotheses of an already-proved ledger: GP or
    frame-coverable, and in scope.
- **(iii) Multiplicative exceptions.** Points whose relation vectors lie in the span of those of lower-rank
  candidates.
  - These are charged by bounding the generation of the constituents, not by a lower-rank count.

**Consequence (proved, given LP_exc).** The yield step of §5 then runs in every box. Exceptions contribute their own
lower-rank yield or no new rank. Hence α ≥ 1/3, and no rank beats rank 2 [proved (model), given LP_exc].

**Status of LP_exc.** It is unproved.
- Restrictions need not inherit GP, frame feasibility or a safe box after a change of lattice basis. That is why the
  hypothesis in (ii) is stated explicitly.
- Norm-equation and unit-orbit families are further natural tests of (ii)–(iii). We found no useful family of either
  kind.

**Tiny boxes** (a_j ∈ {0, ±1}, large r).
- Davenport's error 3^{r−1}rD is the whole box, so nothing here applies.
- What remains is R4, which bounds generator heights below by N^{(r−1)/(D−1)} up to a factor r^{D/2}, together with a
  typical-size heuristic.

**Power-basis 3-dimensional sieving at its optimum** (e = 2β/3).
- The algorithm sieves the whole box at typical sizes, which min(e^C, main) covers.
- The open point is only whether a *selection* inside that box could harvest atypically small values beyond LP_exc.
- In E3, excess appeared only on Z·e₁, where F(a, 0, 0) = a⁶ and the two sides coincide (type A restricted to that
  line), and on c = 0. The toy box there is not large in the L[1/3] sense.

## 7. The split side, ptesplit O3 (task 3)

**Corollary SPL [proved (model), under the charged-side heuristic H′ as in rank3 §6].**
- *Setting.* The charged side W′ is smooth by chance, and the split side has log-size ℓ₁.
- *Available classes.* In GP, pure-s frames exist for max(1, r+1−d′) ≤ s ≤ min(r, d).
- *The bound.* For each available s:
  ℓ₂ ≥ (s/d)[(r−1)x/r − d′ℓ₁/(r+1−s)] + ((r+1)d′/(r(r+1−s)))(λ − log(2^rC_r)).
  When d′ ≥ r+1, the s = 0 class also gives ℓ₂ ≥ (d′/r)(λ − log(2^rC_r)).
- *Where the N-part vanishes.* At ℓ₁ = (r−1)(r+1−s)x/(r d′).
  - With d′ ≥ r, the s = 1 class is available and gives the largest threshold, **(r−1)x/d′**.
  - Otherwise the largest threshold is (r−1)min(r, d′)x/(r d′).
- *Comparison with R5.* R5's range ends at (r−1)x/(D−1). For d ≥ 2 the new range is strictly larger
  (`test_split_threshold`).

**The threshold is attained (size model).**
- *Setting.* Fix a small monic K′ = Q(β), with g of degree d′ ≥ r and a known root μ of g mod N (special N).
- *Construction.* Let ρ be an LLL-short vector of the lattice {ρ : deg ρ < d′, ρ(μ) ≡ 0 (mod N)}. That lattice has
  determinant N, so ‖ρ‖ ≲ N^{1/d′}.
- *The polynomial.* f = x^{d−d′}g + ρ is monic of degree d, with f(μ) ≡ 0 and Res(f, g) = ±Res(ρ, g), which is
  nonzero for generic ρ. Irreducibility is generic and can be checked per instance.
- *Sizes.* M(f) ≲ N^{1/d′}, so the power-basis rank-r module has ℓ₁ ≤ (r−1)x/d′ + d log E + O(d) and
  ℓ₂ ≤ d′ log E + O(d′).
- *Optimality.* N | Res(f, g) forces every such f to have height ≳ N^{1/d′}, so the construction is optimal.
- *Above the threshold.* Linking cannot charge W′ at all: codex's shifted-cyclotomic GP example keeps W′ fixed
  (rank3 §6).

**Exploitability.**
- *Products add no rank [proved].* When the charged lifts multiply, products of dictionary pieces add no relation
  rank. This is codex rankdegen_1, extending rank3's (s+tT)(u+vT) remark.
- *Repeated presentations [proved].* Different words for the same element give characteristic-0 identity rows, so
  presentations are not independent trials.
- *Non-multiplicative residues [open].* General exact residues whose charged sides do not multiply are not excluded.
  This is ptebarrier O3/O1.
- *Under H_split, O3 is one-sided output-sensitive listing.* H_split: among the box points of sizes (Y₁, Y₂), at
  most exp(o(log B))·ρ(u₁) are split by construction.
  - *The program.* The count entropy satisfies γ ≥ β + (η₁+η₂)/(3β). The cost is max(β + η₂/(3β), 2β). And
    η₁η₂ ≥ γΦ, from ℓ₁ℓ₂ ≥ Φλx.
  - *Formal infimum* (box free): **c³ = 8Φ/9 = 8(r−1)/(9r)**. This is a singular limit (γ → ∞), not an attained
    optimum.
    - At r = 2 it is 0.7631, smoothenum's infimum for an algebraic-side-only output-sensitive list.
    - At r = 3 it is 0.8399.
  - *Scoped floor:* 2(4r(r−1)Φ/9)^{1/3}, the same as the doubly-charged one.
  - `ledger.hsplit_floor` agrees numerically, approaching from above.
- **Exponent [proved (model)].**
  - ℓ₁ℓ₂ ≥ Φλx and λ ≥ log B + u₁log u₁ give ℓ₂ ≳ Φx/log B, whether ℓ₁ is small or large.
  - Then u₂ log u₂ ≳ x^{1−2α} ≫ x^α whenever log B ≲ x^α(log x)^{1−α} with α < 1/3.
  - Hence α ≥ 1/3 under H_split.

**Verdict on O3.**
- **Not excluded by norm volume.** No norm-volume argument excludes it. The frame bound shows exactly where linking
  stops charging, and that threshold is attained.
- **Excluded under H_split.** There it is a known constant-factor lever, not a change of exponent.
- **What going lower would need.** Getting below L[1/3] would require generating split elements of a linked module
  more densely than smoothness statistics allow. That is ptesplit's O1 (non-exact compressions), and no candidate is
  known.

## 8. Experiments (pre-registered, commit 4c40b18)

### E1: FH and linking, exact (`src/rankdegen/e1.py`, 35 s, 124 MB RSS)

**Instances.** 14 in all, box E = 30:
- C2 (α⁶ = 2) at N = 31, 31·41 and 31·41·47;
- B (x⁶ + x² + 1, N = 2759);
- A (Φ₅, N = 31·41);
- GP, C1 and FP, at three seeds each (N ≈ 2^20).

**Certification fix (disclosed).** The first run compared counts with the *upper* end of the certified main-term
ball, which does not certify the inequality (codex #12). The rerun compares with the lower end. Results are
identical; the committed file is the certified rerun.

| prediction | result |
|---|---|
| P1: 0 violations of the count bound (N-only and actual-minor versions) | **held**: 0/259 194; 2 923 informative (count > error) |
| P2: Σw log Φ ≥ (r−1)t log N for symmetric pure-s w | **held**, 30/30 |
| P3: type A LP value 0 | **held** |
| P4: C2 LP value 1/9 at every N | **held** (3/3) |
| P5: Φ program = (r−1)/r within 10⁻³ | **held**, 9/9 (r = 4 (2,3) recomputed with s = 2 only: 3/4) |

**Structure found.**
- C2: 40 zero minors of 220, and 225 frames (two rich lines).
- B: 60/220 zero minors, LP 1/9.
- GP pairs of cubics: LP 4/9, against 2/5 for the Bean–Thunder uniform distribution.
- C1 and FP keep t = 4/9, with 12 and 9 frames.

### E2: degenerate designs vs rank 2 (`src/rankdegen/e2.py`; measured, toy sizes)

**Disclosed deviation (codex #11).**
- The pre-registration requires D-C1's quadratic kernel element ν to be "not a multiple of a linear one".
- The implementation checked only that the quadratic coefficient was nonzero. All nine original D-C1 runs used
  ν = x(x − m), which the pre-registration excludes.
- `e2.py --corrected` reruns the registered design with an LLL-selected ν irreducible over Q
  (`rankdegen_E2_corrected.jsonl`). The other designs are unchanged. The original D-C1 runs are kept as exploratory.

**Norm convention.**
- As in rank3 E2, norms are resultants Res(f, A).
- The lattice polynomials are non-monic, so degree-2 candidates carry one more leading-coefficient factor than
  degree-1 ones. This slightly favours rank 2 and changes no conclusion.

| n | best rank-2 yield (R2-rat) | D-C1 original / corrected | D-FP | D-RK | best degenerate / best rank-2 |
|---|---|---|---|---|---|
| 40 (3 seeds) | 3818 / 3476 / 4059 | 42, 21, 20 / 1, 0, 0 | 654 / 551 / 488 | 467 / 480 / 458 | 0.171 / 0.159 / 0.120 |
| 50 | 3359 / 1573 / 2795 | 5, 0, 3 / 0, 0, 0 | 232 / 218 / 228 | 58 / 49 / 58 | 0.069 / 0.139 / 0.082 |
| 60 | 1330 / 1457 / 898 | 0, 0, 0 / 0, 0, 0 | 70 / 137 / 114 | 6 / 3 / 2 | 0.053 / 0.094 / 0.127 |

**Results.**
- **P6 held** (0/9 violations) and **P7 held** (0/9), in both the original and the corrected run.
- Median bit lengths of |F|:
  - rank 2: 49–59;
  - D-FP: 73–99;
  - D-C1: 83–113 original, 112–158 corrected;
  - D-RK: 66–91.
- No lead signal.

### E3: excess over the volume bound, by exceptional type (`src/rankdegen/e3.py`, 10 s)

- **Disclosed deviation.** The pre-registration says "the nine (B) rank-3 (3,3) instances". The rank3 file has three
  (B), three (A) and three (C) instances in that cell. E3 ran the three (B) and three (A), as the selection code
  specifies; the "nine" was a miscount.
- **P8 (excess in every instance): failed as stated**, holding in 3/6.
  - The (A) instances (GNFS-like, N ≈ 2^20) have T(X) > main(X) up to log₂X = 41, 32 and 38.
  - The (B) instances (SNFS-like, N ≈ 1100–1250) have no excess at any X.
- **P9 (G(X) ≤ 2 main(X) + 20): held**, 6/6.
  - Irreducible-quadratic points first appear at log₂X = 77, 72, 71 (A) and 5, 8, 5 (B).
  - The largest G/main is 0.135 (B); for A it is ≤ 4·10⁻⁴.
  - The draft's "≈ 80" and "10⁸" were misreadings, corrected after codex #12.
- **Where the excess lies.** Where there is excess, it lies entirely on Z·e₁ up to log₂X ≈ 38, and then on the plane
  c = 0.
- This is consistent with §6. At these sizes E3 is a consistency check, not evidence for LP_exc at scale.

## 9. GPT-6.1-sol blind round (codex/rankdegen_1.out)

**Agreements, and material we use:**
- the subfield reduction (its (1));
- the characteristic-0 coincidence analysis;
- the α⁶ = 2 and (√2, √3) examples;
- the large-box counterexample (§6.1);
- the result that products add no relation rank.

**Its nonzero-minor ledger is weaker than FH.** It proves v_p(S) ≥ r!b*/(D−r+1) and m(F) ≥ N^{1/(D−r+1)}, then
applies Thunder Thm 4 with an unspecified C_{r,D}. That gives exponent 1/(D−r+1) per unit instead of ≍ r²/D², with
non-explicit constants.

**Where we go further:**
- Theorem FH: a weighted, side-aware form of Bean's Hölder step; explicit; no GP needed;
- Proposition PHI and Lemma PR;
- the instability lemma (S4);
- Lemma CC in ND form;
- Corollary SPL with an attaining construction;
- the H_split costing;
- the experiments.

## 10. Reviews

**Codex adversarial round (codex/rankdegen_2.out): MAJOR REVISION; all items applied.**
- *Confirmed:*
  - FH, including the single phase vector;
  - FN;
  - S1–S4;
  - the Φ algebra and PHI;
  - the constants;
  - the scaled large-box counterexample;
  - the H_split constants;
  - the E1 re-scores, also with lower-end certification.
- *Fixes:*
  - #1: "decoupled ⇒ useless" was false. It is now conditional on gcd(I, N) = 1, and the fibre-product case is open.
  - #2: the classification is now definitional with a catch-all, includes the mixed example, and separates Schmidt
    degeneracy.
  - #3: the finite-volume frame-less example is added to the open cases.
  - #4: PR's projected-error condition and its side-frame hypothesis are stated.
  - #5: the complete-class hypothesis is stated, and P5's last case is corrected.
  - #6: the C2 claim is restated; it reaches the GP constants only as ε′ → 0.
  - #7: CC now has the torsor for odd |S| and the local hypotheses; SUB has its scalar conditions.
  - #8: LP_exc is reformulated (linear vs multiplicative exceptions; exp(o(log B)) losses).
  - #9: the quadric is in four variables, the slab argument needs a real factor, the Evertse citation is corrected,
    and the claim is narrowed.
  - #10: SPL's availability range, the d′ ≥ r condition, and an attaining construction.
  - #11: the E2 deviation is disclosed, and the design rerun.
  - #12: the E3 prose is corrected, and E1 certified against the lower end.

**Opus proof critic: ACCEPT WITH FIXES; all applied.**
- *Verified:*
  - every FH step, with numerics;
  - FN and S4;
  - the PHI identity and crossing step;
  - exactness of Φ_bal;
  - Pólya's constant;
  - SPL and H_split.
- *Fixes:*
  - #1 (MAJOR): the Φ reduction requires the one-sided constraints. The bounded-degree cases now go through R3 or the
    raw-unit check, and the "rational-side limit" wording is corrected.
  - #2 (MAJOR): the C2 overclaim is replaced by the proved Φ ≥ (1−ε′)(r−1)/r with explicit thresholds.
  - #3: the classification is exhaustive by definition.
  - #4: the yield-step bins.
  - #5: when FN beats R1, and the unramified hypothesis.
  - #6–#8: the hypotheses of CC, SUB and PR.
  - #9: ε′ = O(k/min(d, d′)).
  - #10: SPL needs d ≥ 2, and the scale-free argument for α ≥ 1/3 under H_split is added.
  - #11–#12: the per-row weights, and the index condition for decoupled modules.

**Opus literature critic: citation corrections applied.**
- **Bean–Thunder covers complex coefficients** (Zbl 0854.11019; Bean 1994 Thm 1). FH is now stated with C_r via
  Bean–Thunder per frame; the 2^r route is kept as an option.
- **Uniform FH is Bean's 1994 Hölder method.** FH is now presented as its weighted, side-aware form.
- **Evertse.** Thm 1 is now quoted exactly: (16D)^{(r+11)³/3}, all integer solutions, norm forms only. Thm 2 is
  added.
- **Schmidt degeneracy** means a nonzero stable subspace; type D is the case V^J = V. Schmidt 1971/1972 and
  Evertse–Győry 1988 are now cited.
- **Thunder's Lemma 14** is a line-by-line bound, not a projection bound.
- **Citations added or corrected:** the Barroero–Widmer venue; Davenport 1951 with its corrigendum; Pólya 1928 and
  Borwein–Erdélyi.
- **(S4)** is marked as standard (Mumford, Kempf).
- **"Known lever"** is replaced by "smoothenum's lever".
- **Context added:** Edmonds' base polytope and Barthe's Brascamp–Lieb parallel.
- **No prior art found** for the N-divisibility of frame minors or for the split threshold (r−1)x/d′.

## 11. Reproduction

```bash
PYTHONPATH=. .venv/bin/python -m pytest -q tests/test_rankdegen.py     # 12 tests, ~2 s
PYTHONPATH=. .venv/bin/python -m src.rankdegen.e1                       # E1, ~35 s (certified)
PYTHONPATH=. .venv/bin/python -m src.rankdegen.e2                       # E2 as originally run, ~16 s
PYTHONPATH=. .venv/bin/python -m src.rankdegen.e2 --corrected           # E2 with the registered D-C1 design
PYTHONPATH=. .venv/bin/python -m src.rankdegen.e3                       # E3, ~10 s
```
