# audit2L blind claims: Theorem 2L-final (rigorous Las Vegas factoring in L_N[α], 3/8 ≤ α < 2/5)

This file is a blind excerpt prepared for independent review. It contains the statements and proofs of a chain of
internal (unrefereed) research notes, with all review sections, reviewer verdicts, probability estimates and reviewer
attributions removed. Some sentences still describe how an earlier draft was repaired ("v1 used …"); these are
history, not verdicts. Parenthetical attributions were stripped mechanically, so a few sentences may read abruptly.

**Reviewer instruction.** You may read the underlying source notes (paths below) for statements this excerpt cites
but does not reproduce, but do **not** read the review sections of those notes (sections titled "Reviews", the
"Verdict"/"0. Verdict" sections, and any probability estimates) before you have formed and written down your own
verdict.

## The claim under review

Notation: λ := log log N, D := log N/λ, L[α, c] = L_N[α, c] = e^{(c+o(1))D^αλ}; "polylog" = e^{O(λ)}.

- **Theorem 2L-final(α).** For every fixed α ∈ (3/8, 2/5), a Las Vegas algorithm (FACTOR_2L in design 2L_DW(α))
  outputs the complete prime factorisation of every integer N ≥ 2 in expected time L_N[α, c(α)], where
  c(α) = max(2a, a + s + max(a, σ₁ + σ₂), (a + s)/(1 − 2^{−α})), s = 2κ(1 + C_γ), σ₁ = α/(ad₀),
  σ₂ = (α/a)(4/d₀ + d₀b₂), for any admissible rational design constants a, κ, d₀, b₂ > 0, and C_γ is the polylog
  exponent of the level-2 relation supply γ₂ ≥ (log N)^{−C_γ}.
- **Theorem 2L-final(3/8).** The same at α = 3/8 under the design constraint a > 2c₂, c₂ := 3d₀/(32κb₂), with
  σ₁(3/8) = 3/(8ad₀) + c₂.

Status claimed by the authors: proved relative to the project's internal unrefereed chain (paper2/main.tex,
notes/lv5unif_results.md = Theorem U at L[2/5], and the notes excerpted below) plus published results, of which the
essential one is Drappeau, Compositio Math. 151 (2015), Lemme 5 (stated there without proof as a consequence of
§3.3 of Harper's unrefereed preprint arXiv:1208.5992).

## Logical structure (as claimed)

1. Theorem U's one-level chain (paper2/main.tex, notes/lv5unif_results.md) supplies the split/field-selection
   machinery, LV6 (square roots), fruitfulness, recursion (Lemma F_α) and the second-moment tools.
2. **Hyp U′** (notes/uprime_results.md, excerpt E; with the DW replacements of c12alg §2.7, excerpt F): level-1
   moments μ₁ ≥ (log N)^{−C}, E[N₁²] ≤ (log N)^C(μ₁² + μ₁).
3. **Prop UB** (notes/multidescent_results.md, excerpt D): if the dead parent–child incidence mass is ≤ μ₁/2 then
   Pr[N₂ ≥ 1] ≥ μ₁²/(4E N₁²).
4. **Prop DC / DC_e** (notes/pc3_results.md, excerpt C): dead mass ≤ W(A + √(AB)), A, B the incidence-weighted
   relative variances of the child count and of the label's parent count.
5. **Lemma η, Prop PB, Lemma VD** (notes/c12alg_results.md, excerpt F): W/μ₁ ≤ 3u₂; dead mass ≤ μ₁/2 if
   A ≤ D^{−α}λ^C and B ≤ D^{7α−3}λ^C; A, B split into one-point, rational and algebraic blocks.
6. **Theorem RAT** (notes/c12rat_results.md, excerpt B): the rational blocks at second order.
7. **Theorem AB/AB₂** (notes/thmAB_results.md, excerpt A): the algebraic blocks at second order, and the final
   assembly (§13 of that note, Theorem 2L-final).
8. Prop CLB (excerpt G): one-point child lower bound.

Excerpts F and G (c12alg) predate excerpts A and B: where they call Theorem AB/AB₂ or Hyp RAT a "sketch" or a
hypothesis, excerpts A and B are the claimed proofs, and excerpt A §13 is the claimed final assembly.

Other source notes that may be consulted: notes/lv5unif_results.md, paper2/main.tex, notes/nested_results.md,
notes/escapeR2_results.md, notes/lem85_results.md.



---

# A. Theorem AB/AB2 and the assembly (source: notes/thmAB_results.md §§1-14)

## 1. Setting

### 1.1 The design

c12alg's design 2L_DW(α) (c12alg §1.1), for a fixed rational α ∈ [3/8, 2/5), with constants a, κ, d₀, b₂ > 0:
- d := least odd ≥ max(3, d₀D^{1−2α}); k := max(⌈κD^α⌉, d + 1); B = ⌊e^{aD^αλ}⌋; B₂ = ⌊e^{b₂D^{4α−1}λ}⌋;
  η₀ = 1/⌈D^{3/2−3α}⌉, so η₀ log N ≍ D^{3α−1/2}λ;
- the field law DW: g uniform on primitive points of [𝔐/2, 𝔐]²; w uniform on the integers of (W_m, 2W_m];
  c uniform on C_H; F := wF̂_N + G·R_c;
- level 1: Acc_ℓ(P) for P ∈ 𝒫 = {(a, b) : b prime in ((1 − 1/d)V₁, V₁], |a| ≤ V₁, gcd(a, b) = 1}, with the filter
  |G(P)| ∈ [𝔐V₁/32, 2𝔐V₁];
- level 2: child labels 𝔮 = (q, r), q ∈ (B, B₂] prime, q ∤ N; conservative boxes Box(𝔮) (Gauss-reduced basis
  e₁, e₂, coordinates |i|, |j| ≤ S₂, primitive, mod ±); Good(P′) :⇔ G(P′) ∈ 𝒮_B, the child filter
  |G(P′)| ≥ 𝔐|P′|_∞N^{−η₀/4}, and F(P′)/q ∈ 𝒮_B. The skewness is s(𝔮) := |e₂|/|e₁|.
- **Endpoint constraints, only at α = 3/8:** (SQ) a·b₂ > 3d₀/(32κ), and a > 2c₂ with c₂ := 3d₀/(32κb₂) (c12alg O4),
  which implies (SQ).

We assume N has no prime factor ≤ B: trial division to B costs B^{1+o(1)} = L[α, a]. (Not to B₂, which is
L[4α − 1] ≫ L[α].) The special-q primes ℓ ∈ J dividing N are dropped (lv5unif FM-U), and labels have q ∤ N by the
definition of Acc. A prime of N in (B, B₂] that divides a value is a factor exit (c12alg §2.8). Then
gcd(N, G(P)) = 1 for every B-smooth G(P) and every G(P) ∈ J·𝒮_B.

### 1.2 The objects (c12alg Lemma VD, cell by cell)

Fix a label 𝔮 and an environment cell e (a square Ω_e of side δ_e𝔐, δ_e = D^{−A₀}, A₀ ≥ 1; pc3 Prop DC_e).
Write E_e[·] := E[· | R_𝔮, g ∈ Ω_e] for expectations over g. By c12alg Lemma RE, the law of g under this conditioning
has density 1 + O(q/H + q/W_m) with respect to the uniform law on the primitive points of Ω_e for g with
q ∤ G_g(r, 1). The set of g with q | G_g(r, 1) keeps its mass ≈ 1/q (on it the root event is q | w). On that set every
point of Λ_𝔮 has q | G(P′), so Rat(P′; g) = 0 (q > B) and R(P; g) = 0 (q ∉ J, q > B). Such g contribute zero to every
sum below. If the surviving component has conditional mass t = 1 − O(1/q), then ϕ^alg_full = ϕ^alg_good/t and
ϕ^rat_full = ϕ^rat_good/t + (1 − t)/t. We work with the surviving component; the
corrections are a factor 1 + O(1/q) and an additive O(1/q) ≤ O(1/B), negligible against every budget.

*Children.* For P′ ∈ Box(𝔮):
- Rat(P′; g) := 1[G_g(P′) ∈ 𝒮_B, |G_g(P′)| ≥ 𝔐|P′|_∞N^{−η₀/4}];
- A_{P′}(g) := Pr_{w,c}[F(P′)/q ∈ 𝒮_B | g, R_𝔮], and A₁₂(g) := Pr[both | g, R_𝔮] for a pair;
- m(g) := Σ_{P′}Rat(P′; g)A_{P′}(g), and m_{𝔮,e} := E_e m(g) = E[N_𝔮 | R_𝔮, e];
- ϕ^alg_{𝔮,e} := E_e[Σ_{P₁≠P₂}Rat₁Rat₂(A₁₂ − A₁A₂)]/m²_{𝔮,e} and ϕ^rat_{𝔮,e} := Var_e(m(g))/m²_{𝔮,e}.

*Parents.* For P ∈ 𝒫 ∩ Λ_𝔮:
- R(P; g) := Σ_{ℓ∈J}1[Rat_ℓ(P; g)] = 1[G_g(P) ∈ J·𝒮_B, |G_g(P)| ∈ [𝔐V₁/32, 2𝔐V₁]]. At most one ℓ occurs, since two
  primes > B cannot both divide a number of the form ℓ·(B-smooth).
- a_P(g) := Pr_{w,c}[acc(P) | g, R_𝔮], where acc(P) means: F(P) is B₂-smooth, its (B, B₂]-part is squarefree, and no
  prime in (B, B₂] dividing F(P) divides b·f_d·disc(f)·N. Given R_𝔮, q | F(P) automatically (P ∈ Λ_𝔮, q ∤ b).
- z(g) := Σ_P R(P; g)a_P(g), z_{𝔮,e} := E_e z(g), and σ^alg_{𝔮,e} := E_e[Σ_{P₁≠P₂}R₁R₂(a₁₂ − a₁a₂)]/z²_{𝔮,e}.

This is c12alg's σ^alg: in Lemma VD's parent version the pairs are ((ℓ₁, P₁), (ℓ₂, P₂)). P₁ = P₂ forces ℓ₁ = ℓ₂, and
for P₁ ≠ P₂ we have Σ_{ℓ₁,ℓ₂}Rat_{ℓ₁}(P₁)Rat_{ℓ₂}(P₂) = R(P₁)R(P₂). The algebraic block therefore has no
ℓ-structure. The same-ℓ pairs that c12alg treated separately enter only through the rational weight E_e[R₁R₂]
(Lemma W-P below).

The weights are w_{𝔮,e} := Pr(R_𝔮, e)z_{𝔮,e}/#J and W := Σw_{𝔮,e} (pc3 Prop DC_e).

### 1.3 Parameters of the proof

Put θ := 2α − 1/2 and θ_p := max(1 − 2α, 2 − 5α). On [3/8, 2/5), θ = max(3α − 1, 1 − 2α, 2α − 1/2) is c12alg's loss
exponent and c12rat's θ. Fix a large constant C_T (chosen after all the constants C_Ξ below), and set:
- **child thresholds** (c12rat §7): T₀ := e^{C_TD^θλ²}, ϑ₀ := T₀³;
- **parent thresholds** (separate, obligation O3): ϑ₀^p := e^{C_TD^{θ_p}λ²};
- **cell parameters:** h := D^{−10}; h_* := h/(64(d + 1)ϑ₀) for children and h_*^p := h/(64(d + 1)ϑ₀^p) for parents;
  universal grids 𝒳 := {±(1 + h_*)^j : j ≥ 0} and 𝒳^p likewise;
- **cutoffs:** C_X := a fixed integer with B^{C_X − 1} ≥ D^{20}ρ̂_B(8dX_{F2})^{−2} and
  B₂^{C_X − 1} ≥ D^{20}ρ̂_{B₂}(8dX_F)^{−2}. This is possible since log(1/ρ̂) = O(log y) for both densities (c12alg §1.1
  sizes). Every centering constant c_i below is ≥ ρ̂_B(8dX_{F2}) (children) or ≥ ρ̂_{B₂}(8dX_F) (parents), so
  D^{10}(d + 1)²y^{−C_X} ≤ y^{−1}c₁c₂;
- **gcd-class thresholds:** P₀ := ⌈D^{7/10}⌉; S₀ := e^{D^{ζ}} with ζ := (θ + α)/2 (c12rat's S₀) for children, and
  S₀^p := e^{D^{ζ_p}} with ζ_p := (θ_p + 4α − 1)/2 for parents.

Define the losses Ξ := e^{C_ΞD^θλ²} (children) and Ξ_p := e^{C_ΞD^{θ_p}λ²} (parents). C_Ξ is large enough that Ξ
dominates:
- c12rat's Ξ (Lemma PUB);
- the box spread of Lemma SPR;
- the filter spread (N^{η₀/4})^γ = e^{O(D^{2α−1/2}λ²)};
- and Ξ_p dominates the parent losses of c12rat (PUB parents, Ξ_p = e^{O(D^{1−2α}λ)}) and the algebraic parent
  modulus loss (Δ′)^{γ₂} = e^{O(D^{2−5α}λ)}.

**Margins [proved at leading order; `src/thmAB/params.py`].** For every α ∈ [3/8, 2/5), with the first four
consequences used repeatedly:
- θ < α and θ_p < 3 − 7α. So Ξ^{O(1)}, T₀^{O(1)}, h_*^{−O(1)} = e^{O(D^θλ²)} are beaten by the child saving
  e^{−cD^α/λ²}, and Ξ_p^{O(1)}, (h_*^p)^{−O(1)} by the parent saving e^{−cD^{3−7α}/λ²}.
- h_* ≥ B^{−5/12} and h_*^p ≥ B₂^{−5/12} (Hildebrand's range).
- θ < ζ < α, θ_p < ζ_p < 4α − 1, ζ > α − 7/20 and ζ_p > 1 − α − 7/20. These are needed for the gcd-class tails,
  Lemma EXC-A.
- ω(Δ′) = o(P₀) for every child and parent modulus (α < 7/10 and 1 − α < 7/10).
- Every modulus Δ′ is x^{o(1)} at its value scale x (Fact D's level): child log Δ′ ≍ D^αλ against log x ≍ D^{2α}λ;
  parent D^{1−α}λ against D^{2−3α}λ.
- N^{η₀/2} beats every L[α]^{O(1)} and B₂^{O(1)}.
- The outputs beat the budgets: D^{−2α} ≪ D^{−α}; D^{2−8α} + D^{−2(3−7α)} ≪ D^{7α−3}; D^{2−8α} ≤ D^{−2α}. The AB₂
  budget holds on [3/8, 2/5).
- The parent acceptability subtraction is negligible for α > 3/8; at α = 3/8 it needs a > 2c₂.

---

## 2. Statements

**Theorem AB₂ (algebraic blocks at second order) [proved; this note].** In design 2L_DW(α), α ∈ [3/8, 2/5) (at
α = 3/8 with the endpoint constraints of §1.1), there are C, c, N₀ such that for N ≥ N₀:
- (i) Σ_{(𝔮,e): s(𝔮)≤T₀}(w_{𝔮,e}/W)|ϕ^alg_{𝔮,e}| ≤ CD^{−2α}λ^C(1 + Σ_{(𝔮,e): s(𝔮)≤T₀}(w_{𝔮,e}/W)ϕ^rat_{𝔮,e}) + T₀^{−1};
- (ii) Σ_{(𝔮,e)}(w_{𝔮,e}/W)|σ^alg_{𝔮,e}| ≤ C(D^{−2α} + D^{−2(3−7α)})λ^C(1 + Σ_{(𝔮,e)}(w_{𝔮,e}/W)σ^rat_{𝔮,e}) + e^{−cD^{3−7α}/λ²};
- (iii) Σ_{(𝔮,e): s(𝔮)>T₀}(w_{𝔮,e}/W)ϕ_{𝔮,e} ≤ T₀^{−1/2}, where ϕ_{𝔮,e} is the whole root- and cell-conditioned
  relative variance of N_𝔮.

**Theorem AB (c12alg §3.3) [proved; Corollary of AB₂].** For α ∈ (3/8, 2/5), with the threshold T := T₀:
- (i) Σ_{s≤T}(w/W)|ϕ^alg| ≤ D^{−α}λ^C(1 + Σ_{s≤T}(w/W)ϕ^rat);
- (ii) Σ(w/W)|σ^alg| ≤ D^{7α−3}λ^C(1 + Σ(w/W)σ^rat);
- (iii) Σ_{s>T}(w/W)ϕ ≤ D^{−α}.

*Proof of the corollary.* D^{−2α} ≤ D^{−α}. D^{−2(3−7α)} ≤ D^{7α−3}, and D^{−2α} ≤ D^{7α−3} since α > 1/3. The
additive terms T₀^{−1} and e^{−cD^{3−7α}/λ²} are ≤ D^{−α} and ≤ D^{7α−3} respectively, and the factor
(1 + Σϕ^rat) is ≥ 1. Finally T₀^{−1/2} ≤ D^{−α}. ∎

c12alg's statement used T = Ξ³(log N)^CD^α. The assembly (c12alg Theorem 2L-DW(α), step 1) needs (i) and (iii)
with one common threshold, and the same threshold in Hyp RAT(R1). c12rat's (RAT-1) holds for s(𝔮) ≤ T₀, so T₀ is
the natural common choice.

**Proof plan.** Fix (𝔮, e), g and a pair. §3 reduces A₁₂ − A₁A₂ to a model law. §4 gives its geometry. §5 gives the
cell factorisation. §6 gives the one-point moments (O1). §7 handles the exceptional sets. §8 gives the per-pair
covariance bound (Prop PD-A, O5). §9 handles the character errors (O2). §10 handles thin pairs and skew. §11 covers
parents (O3, O4). §12 assembles.

---

## 3. The model law of a pair [proved]

**Lemma ML.** Let (𝔮, e) be a child label and cell, g ∈ 𝒢_e, and P₁ ≠ P₂ ∈ Box(𝔮) with Rat(P₁; g) = Rat(P₂; g) = 1.
Put:
- Δ := |det(P₁, P₂)| (so q | Δ) and Δ′ := Δ/q;
- ι := λ_LAT^d mod Δ′, where λ_LAT is the unit of pc3 Lemma LAT(a) with P₁ ≡ λ_LATP₂ (mod Δ);
- L′ := {x ∈ Z² : x₁ ≡ ιx₂ (mod Δ′)};
- (j, k) a pair of columns chosen by c12alg Lemma MM, and Π₀ := q^{−1}diag(G₁, G₂)M_{jk}[−H, H + 1)², with
  G_i := G_g(P_i);
- Π_y := t_y + Π₀, where t_y := q^{−1}(W_m v + diag(G₁, G₂)y) and v = (F̂_N(P₁), F̂_N(P₂));
- π := Law(Σ_{i∉{j,k}}c_iw_i) for c uniform on C_H.

Define ν := Σ_yπ(y)Unif(L′ ∩ Π_y). Then:
- (a) the law of (F(P₁), F(P₂))/q given (g, R_𝔮) is within TV 7qN^{−η₀/2} of ν;
- (b) with p_i := ν(x_i ∈ 𝒮_B) and p₁₂ := ν(x₁, x₂ ∈ 𝒮_B), |(A₁₂ − A₁A₂) − (p₁₂ − p₁p₂)| ≤ 21qN^{−η₀/2}.

The same holds for parents P₁ ≠ P₂ ∈ 𝒫 ∩ Λ_𝔮 with R(P₁; g) = R(P₂; g) = 1, with Δ′ = Δ/q ≤ 2V₁²/q, **for the
smoothness events** s_P(g) := Pr[F(P)/q ∈ 𝒮_{B₂} | g, R_𝔮] and s₁₂ (in (b), A is replaced by s and 𝒮_B by 𝒮_{B₂}).
The acceptability events a_P are not value events; §11 transfers from s to a.

*Proof.*
- (a) is c12alg Cor TP-q. Its hypotheses hold:
  - P_i ∈ Λ_𝔮 are primitive and non-proportional;
  - q ∤ b_i, since a Z²-primitive point of Λ_𝔮 with q | b would have q | a;
  - q ∤ G_i, since G_i is B-smooth (children) or ℓ·(B-smooth) with ℓ > B₂ ≥ q (parents);
  - q ∤ N, and gcd(N, G₁G₂) = 1 (§1.1);
  - the filter |G_i| ≥ 𝔐|P_i|_∞N^{−η₀/4} (children: part of Rat; parents: |G| ≥ 𝔐V₁/32);
  - T = |G₁G₂| ≤ W_mN^{−η₀}.

  L′ is the divided lattice L_q/q of Cor TP-q. ι is a unit mod Δ′ because
  λ_LAT is a unit mod Δ (pc3 LAT(a)).
- (b) The marginals of TV-close laws are TV-close. So |A₁₂ − p₁₂| ≤ 7qN^{−η₀/2}, |A_i − p_i| ≤ 7qN^{−η₀/2}, and
  |A₁A₂ − p₁p₂| ≤ |A₁ − p₁| + |A₂ − p₂|. ∎

Note that the true one-point probability A_i(g) does not depend on the partner, while p_i does, through Π₀ and π.
Both are within 7qN^{−η₀/2} of each other, which is all that is used.

---

## 4. Geometry of the model [proved; elementary]

Fix (g, pair) as in Lemma ML and put s_i := (2H + 1)|G_i||P_i|_∞^{d−1}/q and X_i := (d + 1)s_i.

**Lemma PG (parallelogram geometry).** For N ≥ N₀:
- (a) **Widths.** The x_i-width W_i of Π₀ satisfies s_i/2 ≤ W_i ≤ 2s_i.
- (b) **Extent.** supp ν ⊂ [−X₁, X₁] × [−X₂, X₂].
- (c) **Chords.** Let ℓ_r and ℓ_c be the longest chords of Π₀ parallel to the x₁- and x₂-axes. Then
  ℓ_r ≥ s₁/(8ϑ) and ℓ_c ≥ s₂/(8ϑ), with ϑ = ϑ(P₁, P₂). Also |Π₀| ≥ ℓ_rW₂/2 and |Π₀| ≥ ℓ_cW₁/2.
- (d) **Fibre masses and anti-concentration.** Every x₁-fibre of every Unif(L′ ∩ Π_y) has mass ≤ 3/W₁. Hence
  ν(|x₁| ≤ tX₁) ≤ 6(d + 1)(2t + 1/X₁) for all t > 0, and likewise for x₂.
- (e) **Local law.** For every residue r mod Δ′, ν(x₁ ≡ r mod Δ′) ≤ (1 + N^{−Ω(μ)})/Δ′ + N^{−Ω(μ)}. Hence
  ν(gcd(x₁, Δ′) = δ) ≤ (1 + N^{−Ω(μ)})π₀(δ) + φ(Δ′/δ)N^{−Ω(μ)}, with π₀(δ) := φ(Δ′/δ)/Δ′ (Δ′ ≤ V₁² = 𝔐^{o(1)}). Also, for every interval I,
  ν(x₁ ∈ I, δ | x₁) ≤ (|I|/δ + 1)·3/W₁ (by (d)).
- (f) The same statements hold for the enlarged parallelograms Π_y⁺ of Lemma BL, with X_i replaced by (1 + 2h)X_i.

*Proof.*
- (a) The x_i-width is (2H + 1)|G_i|(|v_j(P_i)| + |v_k(P_i)|)/q. c12alg Lemma MM's column claim gives a column of
  the chosen pair with |v(P_i)| ≥ |P_i|_∞^{d−1}/2, and every |v_l(P_i)| ≤ |P_i|_∞^{d−1}.
- (b) The c-part of a point of supp ν is q^{−1}diag(G₁, G₂)Σ_lc_lw_l, with |Σ_lc_lv_l(P_i)| ≤ dH|P_i|_∞^{d−1}. The
  w-part is W_m|F̂_N(P_i)|/q ≤ W_md2^d𝔐|P_i|^d/q (c12alg Lemma GB-DW). With the filter this is
  ≤ N^{−η₀/2}s_i, because H ≥ 2^{d+8}d⁴W_mN^{η₀}.
- (c) |Π₀| = (2H + 1)²|G₁G₂||m_{jk}|/q², and |m_{jk}| ≥ Δ(|P₁|_∞|P₂|_∞)^{d−2}/4 by Lemma MM. Cavalieri gives
  |Π₀| = ∫(horizontal chord) dx₂ ≤ ℓ_rW₂, and concavity of the chord function gives |Π₀| ≥ ℓ_rW₂/2. So
  ℓ_r ≥ |Π₀|/W₂ ≥ (2H + 1)|G₁|Δ|P₁|_∞^{d−2}/(8q|P₂|_∞), using W₂ ≤ 2s₂. Since Δ/(|P₁|_∞|P₂|_∞) ≥ 1/ϑ, this is
  ≥ s₁/(8ϑ). The x₂-statement is symmetric.
- (d) An x₁-fibre of L′ ∩ Π_y is an arithmetic progression of difference Δ′ in an interval of length ≤ ℓ_c. So it has
  ≤ ℓ_c/Δ′ + 1 points. Write Z := #(L′ ∩ Π_y). Then Z ≥ |Π₀|/Δ′ − (2ℓ_c/Δ′ + W₁ + 1) ≥ ℓ_cW₁/(2Δ′)(1 − o(1)), so
  the fibre mass is ≤ (ℓ_c/Δ′ + 1)/Z ≤ 3/W₁. This uses ℓ_c/Δ′ → ∞ and W₁ → ∞: ℓ_c ≥ s₂/(8ϑ), with s₂ ≥ N^{μ} and
  ϑ, Δ′ ≤ L^{O(1)} for non-thin pairs, and ϑ ≤ |P₁||P₂| ≤ V₁² in general. Then sum over the ≤ 2tX₁ + 1 integers with
  |x₁| ≤ tX₁, and use X₁ ≤ 2(d + 1)W₁.
- (e) mass(x₁) = (ℓ_y(x₁)/Δ′ + θ_{x₁})/Z with |θ| ≤ 1. ℓ_y(·) is concave, hence unimodal, and for such a function
  and any arithmetic progression of difference s, |Σ_{x₁≡r (s)}f(x₁) − s^{−1}Σf(x₁)| ≤ 2 max f. This gives the
  finite form (`cells.fibre_mass_check`, for r = 0 and every s | Δ′): ν_y(x₁ ≡ r mod s) ≤ (Z₀/s + E₁)/(Z₀ − E₂), with
  Z₀ = |Π₀|/Δ′, E₁ = 2ℓ_c/Δ′ + W₁/s + 1, E₂ = 2ℓ_c/Δ′ + W₁ + 1. Here E_i/Z₀ ≤ N^{−Ω(μ)}, since
  Z₀ ≥ ℓ_cW₁/(2Δ′) and ℓ_c/Δ′, W₁ ≥ N^{Ω(μ)}. Take s = Δ′ and sum over the φ(Δ′/δ) residues of gcd δ. The last
  claim is (d) summed over the ≤ |I|/δ + 1 multiples of δ in I.
- (f) Π_y⁺ is convex, contains Π_y, and its widths and chords are at least those of Π_y, while its extent grows by
  2z_i ≤ 2hX_i. ∎

`tests/test_thmAB.py::test_parallelogram_geometry_PG` checks (a) and (c), including Lemma MM, exactly on 63 pairs
(d ∈ {3, 4, 5}, nearly parallel pairs included). `test_fibre_masses_PG_c` checks the finite forms of (d) and (e).

**Lemma BL (blocks, cells, boundary).** Fix (g, pair), non-thin: ϑ ≤ ϑ₀ (children) or ϑ ≤ ϑ₀^p (parents). Put
z₁ := hℓ_r, z₂ := hℓ_c, X_i⁻ := X_iB^{−C_X} (children) or X_iB₂^{−C_X} (parents).
- *Blocks.* Round X_i⁻ down and 2X_i up to points of the grid 𝒳 (resp. 𝒳^p), so every block endpoint is a grid
  point and hence a library endpoint. A block of coordinate i is a union of consecutive grid intervals,
  of relative width ≤ min(h^{1/2}, 1/8) and length ≤ 2z_i. The blocks partition the rounded range
  {X_i⁻ ≤ |x_i| ≤ 2X_i}, and this partition depends only on (g, pair), not on y. Every block contains a whole grid
  interval, so its relative width lies in [h_*/2, h^{1/2}].
- *Rectangles* are products R = I₁ × I₂ of blocks. *Cells* are pairs (R, δ) with δ | Δ′.
- *Core.* Π_y^∘ := {x ∈ Π_y : the x₂-chord through x has length ≥ hℓ_c and the x₁-chord has length ≥ hℓ_r}. It is
  convex: it is an intersection of Π_y with a vertical and a horizontal strip, since the chord functions are concave.
  A rectangle is *interior* for y if R ⊂ Π_y^∘. It is *boundary* for y if it meets Π_y but is not interior.

Then:
- (i) The grid resolves the blocks. The spacing of 𝒳 at |x| ≤ 2X_i is ≤ 2h_*X_i, and
  z_i ≥ hs_i/(8ϑ) = hX_i/(8(d + 1)ϑ) ≥ hX_i/(8(d + 1)ϑ₀) = 8h_*X_i (Lemma PG(c)), so the spacing is ≤ z_i/4.
  Every block contains at least one grid interval, so it has length ≥ h_*X_i⁻/2 ≫ N^{μ/2}Δ′S₀B² (resp. with S₀^p,
  B₂²), because X_i ≥ N^{h}.
- (ii) The boundary rectangles for y cover area ≤ Ch|Π₀|, so their ν_y-mass is ≤ Ch(1 + N^{−Ω(μ)}). All of them lie
  in the convex set Π_y⁺ := Π_y + [−2z₁, 2z₁] × [−2z₂, 2z₂], and |Π_y⁺| ≤ (1 + 40h)|Π₀|.
- (iii) For every rectangle R with sides ≥ Δ′, #(L′ ∩ R) = |I₁||I₂|/Δ′ + O(|I₁| + |I₂|). Each Δ′ × Δ′ square contains
  exactly Δ′ points of L′ ⊇ Δ′Z².

*Proof.* (i) is c12rat Lemma CELL(i), with Lemma PG(c) in place of the chord bound there. For (ii), a boundary
rectangle either meets ∂Π^∘ or lies inside Π∖Π^∘ ∪ (Π∖Π^∘ + box):
- the set Π∖Π^∘ has area ≤ hℓ_cW₁ + hℓ_rW₂ ≤ 4h|Π|, by PG(c), and it is the union of ≤ 4 caps at the ends of the
  chord functions; its (2z₁, 2z₂)-box neighbourhood has area ≤ 4h|Π| + 4z₁W₂ + 4z₂W₁ + 16z₁z₂;
- the rectangles meeting ∂Π^∘ lie within box-distance (2z₁, 2z₂) of it, an area ≤ 4z₁W₂ + 4z₂W₁ + 16z₁z₂
  ≤ (16h + 32h²)|Π|. This is c12rat Lemma CELL(ii)'s convexity argument, which uses only convexity and
  |Π| ≥ ℓ_rW₂/2, |Π| ≥ ℓ_cW₁/2;
- in total ≤ (40h + 64h²)|Π|;
- |Π⁺| ≤ |Π| + 2z₁·2W₂ + 2z₂·2W₁ + 16z₁z₂.

(iii) is elementary. ∎ (`test_boundary_area_BL` checks the bound (16h + 32h²)|Π| exactly on DW parallelograms for
h ∈ {1/4, 1/10, 1/25}.)

---

## 5. Exact factorisation in a cell [proved; finite algebra]

**Lemma CL (gcd classes).** Let ι be a unit mod Δ′, δ | Δ′ and Δ_δ := Δ′/δ.
- (i) For x ∈ L′, gcd(x₁, Δ′) = gcd(x₂, Δ′).
- (ii) x ↦ (x₁/δ, x₂/δ) is a bijection from {x ∈ L′ : gcd(x₁, Δ′) = δ} onto
  {(m₁, m₂) ∈ Z² : (m₁m₂, Δ_δ) = 1, m₁ ≡ ιm₂ (mod Δ_δ)}.

*Proof.* (i): x₁ ≡ ιx₂ with ι a unit. (ii): gcd(δm, Δ′) = δ·gcd(m, Δ_δ), and δm₁ ≡ ιδm₂ (mod Δ′) iff
m₁ ≡ ιm₂ (mod Δ_δ). ∎ (`test_class_bijection_CL`: exhaustive, six moduli including prime powers.)

So prime powers and "full" classes need no special treatment. If p^k ∥ Δ′ and p^{k+5} | x₁, then v_p(δ) = k and
m₁ is unrestricted at p, because p ∤ Δ_δ. This replaces c12rat's types, whose small-prime part was needed only for the
primitivity of g, which the algebraic side does not have.

**Lemma CF (cell factorisation).** Let R = I₁ × I₂ be a rectangle with sides ≥ N^{μ/2}Δ′, and δ | Δ′. Put
J_i := I_i/δ and:
- U_i := #{m ∈ J_i : (m, Δ_δ) = 1};
- a_i(m) := 1[m ∈ J_i, (m, Δ_δ) = 1, |δm| ∈ 𝒮_y], A_i := Σa_i, and d_i(I_i, δ) := A_i/U_i (the one-point cell
  density);
- n(R, δ) := #{x ∈ L′ ∩ R : gcd(x₁, Δ′) = δ}.

Then:
- (i) n(R, δ) = U₁U₂/φ(Δ_δ)·(1 + O(2^{ω(Δ′)}Δ′/min|J_i|));
- (ii) #{x ∈ L′ ∩ R, class δ, x₁ ∈ 𝒮_y} = n(R, δ)d₁(1 + O(2^{ω(Δ′)}Δ′/|J₂|)), and symmetrically for x₂;
- (iii) #{x ∈ L′ ∩ R, class δ, x₁, x₂ ∈ 𝒮_y} = n(R, δ)d₁d₂(1 + O(2^{ω(Δ′)}Δ′/min|J_i|)) + Err(R, δ), with
  |Err(R, δ)| ≤ (A₁A₂/φ(Δ_δ))·(Ẽ₁ + Ẽ₂)/2 and Ẽ_i := A_i^{−2}Σ_{χ≠χ₀ (mod Δ_δ)}|Σ_ma_i(m)χ(m)|².

*Proof.*
- (i), (ii): by Lemma CL, each m₁ ∈ J₁ coprime to Δ_δ has #{m₂ ∈ J₂ : m₂ ≡ ι^{−1}m₁ (mod Δ_δ)} = |J₂|/Δ_δ + O(1)
  partners, all automatically coprime to Δ_δ. Also U_i = |J_i|φ(Δ_δ)/Δ_δ + O(2^{ω(Δ′)}).
- (iii): the count is Σ_{m₁,m₂}a₁(m₁)a₂(m₂)1_C(m̄₁, m̄₂) with C := {(r₁, r₂) ∈ ((Z/Δ_δ)^×)² : r₁ ≡ ιr₂}. C is
  invariant under the diagonal action and |C| = φ(Δ_δ). c12rat Lemma CP gives |count − A₁A₂/φ| ≤ φ^{−1}(E₁E₂)^{1/2},
  where E_i := Σ_{χ≠χ₀}|â_i(χ)|² = A_i²Ẽ_i, and √(Ẽ₁Ẽ₂) ≤ (Ẽ₁ + Ẽ₂)/2. Finally A₁A₂/φ = n(R, δ)d₁d₂ times the factor
  of (i). ∎ (`test_cell_factorisation_CF`: the Parseval form, exactly, on 120 cells including non-reduced classes.)

The singles (ii) carry **no** character error. Every m₁ has the same number of partners up to O(1); this is the
algebraic counterpart of c12rat's Lemma FIB. Both A₁A₂/φ and the singles are therefore normalised by the same n(R, δ).

---

## 6. One-point cell densities and their moments (obligation O1) [proved]

Fix (g, pair) and a coordinate i (say i = 1). Let y = B (children) or B₂ (parents), and β := β(X₁, y), γ := 1 − β.
For δ | Δ′ write Δ_δ = Δ′/δ, and define the **local factor**
  τ(δ) := ∏_{p|Δ′} τ_p(v_p(δ)), where for p^k ∥ Δ′:
  - τ_p(j) := p^{jγ}(1 − p^{−β})/(1 − p^{−1}) for j < k and p ≤ y;
  - τ_p(k) := p^{kγ} for p ≤ y;
  - τ_p(0) := (1 − p^{−1})^{−1} and τ_p(j) := 0 for j ≥ 1, when p > y.

Put c₁ := ρ̂_y(X₁) (the **centering constant**; it depends on g and P₁, not on y or on the partner beyond Δ′).

**Lemma LF-A (local factor moments).** Let π₀(δ) := φ(Δ_δ)/Δ′, the law of gcd(r, Δ′) for r uniform mod Δ′. Then:
- E_{π₀}τ = 1;
- E_{π₀}(τ − 1)² ≤ Cγ²λ² + Cω(Δ′)/y;
- E_{π₀}τ^{16} ≤ C.

This holds uniformly for children (Δ′ ≤ 2S₂²) and parents (Δ′ ≤ 2V₁²/q).

*Proof.* π₀ is multiplicative, with law p^{−j}(1 − p^{−1}) for j < k and p^{−k} for j = k. Hence
E τ_p = Σ_{j<k}p^{−jβ}(1 − p^{−β}) + p^{−kβ} = 1 (telescoping) for p ≤ y, and (1 − 1/p)(1 − 1/p)^{−1} = 1 for p > y.

For the moments write E τ^r = ∏_pEτ_p^r ≤ exp(Σ_p(Eτ_p^r − 1)).
- For p ≤ e^{1/(16γ)} we have p^{16γ} ≤ e, and E τ_p^r − 1 ≤ C_r(p^γ − 1)²/p + C_r·(second-order terms)
  ≤ C_rγ²log²p/p. The linear terms cancel because E τ_p = 1.
- For e^{1/(16γ)} < p ≤ y: E τ_p^r − 1 ≤ 2p^{−1+rγ}, and the number of such p dividing Δ′ is
  ≤ 16γ log Δ′ = O(λ + D^{2−5α}λ). Their total is ≤ e^{O(λ)}D·e^{−1/(32γ)}, which is negligible since 1/γ ≥ cD^α/λ
  (children) or cD^{4α−1} (parents).
- For p > y: E τ_p^r − 1 ≤ C_r/p, total ≤ C_rω(Δ′)/y.
- Σ_{p|Δ′}log²p/p ≤ Σ_{p≤ω(Δ′)log ω(Δ′)}log²p/p ≤ C(log ω(Δ′))² ≤ Cλ².

The variance is Eτ² − 1. ∎

**Lemma OPV-A (one-point moments).** Let x₁ have the first marginal of ν (resp. of ν⁺ := Σ_yπ(y)Z_y⁺^{−1}·counting
measure on L′ ∩ Π_y⁺, the mixture of the uniform laws on the enlarged parallelograms). Call x₁ *exceptional* if
|x₁| < X₁⁻ or gcd(x₁, Δ′) > S₀ (resp. S₀^p). For non-exceptional x₁ let Ψ₁(x₁) := d₁(I(x₁), δ(x₁)) be the density of
its cell, I(x₁) its block and δ(x₁) its class. Then
  E[(Ψ₁/c₁ − 1)²·1_{non-exc}] ≤ 𝒱,  E[(Ψ₁/c₁)⁸·1_{non-exc}] ≤ C,
where
- 𝒱 := C[γ²(log Cd)² + γ²λ² + E_FO² + u^{−2}] ≤ CD^{−2α}λ^C (children);
- 𝒱 ≤ C(D^{2−8α} + D^{−2(3−7α)})λ^C ≤ C(D^{−2α} + D^{−2(3−7α)})λ^C (parents).

Moreover, independently of the above, p₁ ≤ 12(d + 1)c₁ (from PG(d) alone), and c₁/2 ≤ p₁ (step 6 below, which uses
only steps 1–5 and Lemma EXC-A's near-zero and large-class bounds for singles, which do not use p₁).

Here E_FO = C(log P₀)³(log(u + 1)/log y + 1/u) is c12rat's first-order error at the value scale: ≍ D^{−α}λ^C for
children and ≍ (D^{1−4α} + D^{7α−3})λ^C for parents.

*Proof.*
1. **Every block, non-exceptional classes.** Let ξ be the end of I nearer 0, so X₁⁻ ≤ |ξ| ≤ 2X₁, and let δ ≤ S₀
   be y-smooth. Every block has |I| ∈ [h_*|ξ|/2, h^{1/2}|ξ|] ⊂ [y^{−5/12}|ξ|, |ξ|] (Lemma BL; h_* ≥ y^{−5/12}).
   - Apply c12rat **Lemma FO** with x = |ξ|, z = |I|, s = δ, and 𝔪 := rad(Δ_δ)_{≤y}. Its hypotheses hold:
     - log x ≥ (1/2)log X₁;
     - s ≤ S₀;
     - 𝔪 | P₀#·𝔪₁ with 𝔪₁ the primes of Δ_δ in (P₀, y], and Σ_{p|𝔪₁}p^{−β} ≤ eω(Δ′)/P₀ + ω(Δ′)e^{−1/(2γ)} = o(1),
       since p^γ ≤ e for p ≤ e^{1/γ}.
   - FO gives #{m ∈ I/δ : P⁺(|m|) ≤ y, (m, 𝔪) = 1} = (|I|/δ)δ^{γ(ξ)}ρ̂(|ξ|)g_𝔪(β(ξ))(1 + O(E_FO)).
   - Primes p > y dividing Δ_δ need no condition for y-smooth m. Dividing by U = (|I|/δ)φ(Δ_δ)/Δ_δ·(1 + N^{−Ω(μ)})
     (Lemma BL(i)) gives Ψ₁ = ρ̂(|ξ|)·τ(δ)·(1 + O(E_FO)) with τ evaluated at β(ξ, y).
   - **β-drift.** Replacing β(ξ, y) by β(X₁, y) in δ^γ and in g_𝔪 costs a factor
     1 + O(|β(ξ) − β(X₁)|·(log S₀ + Σ_{p|Δ′, p≤y}p^{−β}log p)). By HT Thm 2, |β(ξ) − β(X₁)| ≪ log(X₁/|ξ|)/(log X₁ log y)
     ≤ C_X/log X₁ (as |ξ| ≥ X₁y^{−C_X}), and Σ_{p|Δ′}p^{−β}log p ≤ Cλ. So the factor is 1 + O(D^{ζ−2α}/λ) for
     children and 1 + O(D^{4α−2}/λ) for parents. These are absorbed
     into e (step 4); their squares, D^{−α−1/2} and D^{8α−4}, are ≤ D^{−2α} on [3/8, 2/5).
   - For δ not y-smooth, Ψ₁ = 0 = c₁τ(δ) exactly.
2. **Archimedean factor.** Put Y(x₁) := ρ̂(|ξ|)/ρ̂(X₁).
   - Fact HT3 iterated ≤ C_X + 1 times (|ξ| ≥ X₁y^{−C_X}), with the drift of β controlled by HT Thm 2, gives
     Y = (X₁/|ξ|)^γ(1 + O(1/u)). The block-to-point replacement |ξ| ↔ |x₁| costs (1 + h^{1/2})^γ = 1 + O(γh^{1/2}).
   - By Lemma PG(d), Pr[X₁/|x₁| > s] ≤ 6(d + 1)(2/s + 1/X₁). Integrating,
     E[log^r(X₁/|x₁|)(X₁/|x₁|)^{σ}·1_{|x₁|≥X₁⁻}] ≤ C_r(log Cd)^r for r ≤ 16 and σ ≤ 16γ, since (Cd)^{16γ} ≤ 2.
   - Hence E(Y − 1)⁴ ≤ Cγ⁴(log Cd)⁴ + Cu^{−4} and EY^{16} ≤ C.
3. **Local factor.** By Lemma PG(e), E[f(δ(x₁))] ≤ (1 + N^{−Ω(μ)})E_{π₀}f + Δ′N^{−Ω(μ)}max f for every f ≥ 0
   (max f ≤ e^{O(λ)}Δ′^{16γ} for the f used). Lemma LF-A then gives E(τ − 1)² ≤ Cγ²λ² + Cω/y and Eτ^{16} ≤ C under
   the x₁-marginal of ν (or ν⁺). No conditioning on the block is needed: the archimedean and local factors are
   combined by Cauchy–Schwarz and Hölder.
4. **Assembly.** On every non-exceptional x₁, Ψ₁/c₁ − 1 = (Y − 1)τ(1 + e) + (τ − 1)(1 + e) + e, with
   |e| ≤ CE_FO + Cγh^{1/2} + C(D^{ζ−2α} + D^{4α−2}). Cauchy–Schwarz with the fourth moments gives E[(Ψ₁/c₁ − 1)²1_{non-exc}]
   ≤ C(γ²(log Cd)² + γ²λ² + E_FO² + u^{−2}). Hölder gives E[(Ψ₁/c₁)⁸1_{non-exc}] ≤ C(EY^{16}Eτ^{16})^{1/2} ≤ C.

5. **Pointwise upper bound** (every block, every δ, used only on exceptional sets and in §9). Ψ₁ ≤ e^{O(λ)}δ^{γ}(X₁/|ξ|)^γc₁.
   Indeed:
   - c12rat **Lemma BTS** with modulus 1, applied to I/δ (length ≥ 2y), gives #smooth ≤ e^{O(λ)}(|I|/δ)·(…)^γρ̂;
   - Fact dBT1 gives ρ̂(x/δ) ≤ Cδ^γρ̂(x);
   - U ≥ (|I|/δ)/(Cλ);
   - h_*^{−γ} = e^{o(1)}.
6. **p₁ versus c₁**.
   - *Upper:* p₁ ≤ max fibre mass × #{x₁ ∈ 𝒮_y : |x₁| ≤ X₁} ≤ (3/W₁)·2Ψ(X₁, y) ≤ 12(d + 1)c₁, by PG(d).
   - *Lower*(1,0) + [0,2)(1,1)). Use the interior cells of §8 and Lemma CF(ii) instead:
     p₁ ≥ (1 − tiny)E_ν̂Ψ₁ and, since Ψ₁(∗) = 0,
     E_ν̂[Ψ₁/c₁] ≥ ν̂(interior) − E_ν̂[|Ψ₁/c₁ − 1|1_{int}] ≥ 1 − ν̂(∗) − 𝒱^{1/2} ≥ 1 − Ch − e^{−λ} − 𝒱^{1/2}.
     Here ν̂(∗) ≤ Ch + e^{−λ} uses only Lemma BL(ii) and EXC-A's mass bound ν(ℰ) ≤ e^{−λ} (PG(d) and the gcd tails),
     and 𝒱 is the moment bound of steps 1–4, an expectation of a function of the cell of x₁ under ν (no identity
     between masses and densities is needed). Hence p₁ ≥ c₁/2; this uses neither PU-A nor boundary smoothness.
   - For parents the same, with γ₂ log d ≍ D^{1−4α}λ → 0.

**Exponents.**
- *Children:* γ ≤ C₁/(aD^α)·λ (Fact HT(α2): γ log y = O(λ)), E_FO ≍ D^{−α}λ^C, 1/u ≍ D^{−α}.
- *Parents:* γ₂ = O(λ/log B₂) = O(D^{1−4α}), and E_FO ≍ (log(u₂ + 1)/log B₂ + 1/u₂)λ^C = (D^{1−4α} + D^{7α−3})λ^C.
  D^{2−8α} ≤ D^{−2α} since α ≥ 1/3. ∎

This is the **joint moment lemma O1**. It is uniform over:
- large gcd classes and prime powers (Lemma CL makes them exact; Lemma LF-A carries prime powers);
- rectangles near zero (the cutoff X₁⁻, and FO on every block down to X₁⁻);
- the enlarged boundary region (the ν⁺ version, PG(f));
- the mixture components (c₁ does not depend on y, and the statement is about the full mixture marginal).

The parent local factor and the large gcd classes δ > S₀^p, which c12alg listed under O1, are in Lemma LF-A and
Lemma EXC-A.

---

## 7. Exceptional sets [proved]

**Lemma EXC-A.** Let ℰ := {x : |x₁| < X₁⁻ or |x₂| < X₂⁻ or gcd(x₁, Δ′) > S₀}, with S₀^p for parents. Then
- ν(ℰ ∩ {x₁, x₂ ∈ 𝒮_y}) ≤ e^{−λ}c₁c₂;
- ν(ℰ ∩ {x_i ∈ 𝒮_y}) ≤ e^{−λ}c_i;
- ν(ℰ) ≤ e^{−λ}.

The same holds for ν⁺.

*Proof.*
- **Near zero.** ν(|x_i| < X_i⁻) ≤ 6(d + 1)(2y^{−C_X} + 1/X_i) (Lemma PG(d)). By the choice of C_X (§1.3) this is
  ≤ y^{−1}c₁c₂ ≤ e^{−λ}c₁c₂.
- **Large classes.** δ | x₁ forces δ | x₂ (Lemma CL). For fixed δ | Δ′:
  - ν(δ | x₁, x₁ ∈ 𝒮_y) ≤ (3/W₁)·2Ψ(X₁/δ, y) ≤ Cdδ^{−1+γ}c₁, by PG(d) and Fact dBT1.
  - Given x₁, the partners x₂ form an AP of difference Δ′ on a fibre, an interval of length ℓ(x₁). **First** split
    the fibre into its intersections with the dyadic ranges ±[2^kX₂⁻, 2^{k+1}X₂⁻) (O(log X₂) pieces; the part
    |x₂| < X₂⁻ is covered by the near-zero bound). **Then** cover each nonempty piece by an interval J of length
    max(|piece|, λ_min), λ_min := X₂y^{−C_X}/ϑ ≤ X₂⁻, placed **inside the piece's own dyadic range** [t, 2t]
    (t = 2^kX₂⁻; reflect negative pieces). On J, c12rat **Lemma BTS** applies:
    modulus Δ′, W = |J|/Δ′ ≥ λ_min/Δ′ ≥ 2y, log W ≥ log X₂ − C_X log y − log ϑ − log Δ′ = (1 − o(1))log X₂, and
    gcd(a, Δ′)W ≤ |J| ≤ t ≤ inf J. It bounds the smooth partners in the piece by
    e^{O(λ)}(|J|/Δ′)(2tΔ′/|J|)^{γ}ρ̂(2t) ≤ e^{O(λ)}(y^{C_X}ϑΔ′)^γ(|J|/Δ′)c₂, using
    ρ̂(2t) ≤ C(X₂/X₂⁻)^γρ̂(X₂) ≤ e^{O(λ)}c₂. Parent fibres use the parent size bounds.
    For children, (B^{C_X}Δ′)^γ = e^{O(λ)} (Δ′ ≤ 2S₂² = L[α]^{O(1)}, γ ≍ D^{−α}λ), and ϑ^γ ≤ ϑ₀^γ = e^{o(1)} for
    non-thin pairs.
  - So ν(δ | x, both smooth) ≤ e^{O(λ)}δ^{−1+γ}c₁c₂ for non-thin child pairs. For parents the factor is
    (B₂^{C_X}Δ′ϑ₀^p)^{γ₂} ≤ Ξ_p.
  - Sum over δ | Δ′ with δ > S₀. Split δ = δ_sδ_m into its P₀-smooth part and the rest; then δ_s > S₀^{1/2} or
    δ_m > S₀^{1/2}. Rankin with exponent σ = 1/log P₀ on δ_s gives Σ_{δ_s>S₀^{1/2}}δ_s^{−1+γ} ≤ S₀^{−σ/2}(log P₀)^C.
    Rankin with exponent 1/2 on δ_m gives Σ_{δ_m>S₀^{1/2}}δ_m^{−1+γ} ≤ S₀^{−1/4}exp(2ω(Δ′)P₀^{−1/2+γ}).
  - By the margins (ζ > θ, ζ > α − 7/20, and their parent versions), both are ≤ e^{−2λ}(Ξ_p)^{−1}.
- **Singles** are the same computation without the partner. ∎

c12alg's PU (alg) has the same structure. The algebraic modulus is Δ′, not qΔ′, so **for non-thin child pairs the
algebraic pair upper bound has only polylog loss**. For thin pairs the loss is ϑ^γ ≤ (V₂²)^γ = e^{O(D^{3α−1}λ²)} ≤ Ξ,
the "large-modulus" term of c12alg's Ξ. In particular:

**Lemma PU-A (algebraic pair upper bound).** For every pair and every admissible g:
- p₁₂ ≤ e^{O(λ)}(1 + ϑ)^γc₁c₂ (children): ≤ e^{O(λ)}c₁c₂ for ϑ ≤ ϑ₀, and ≤ Ξc₁c₂ always;
- p₁₂ ≤ Ξ_pc₁c₂ (parents);
- consequently A₁₂ ≤ p₁₂ + 7qN^{−η₀/2}.

*Proof.* p₁₂ = Σ_yπ(y)Z_y^{−1}Σ_{x₁∈𝒮_y}#{smooth partners on the fibre of x₁}. Bound the partners piece by piece as
in the large-classes bullet with δ = 1. The total covering length of a fibre is ≤ ℓ(x₁) + O(log X₂)λ_min. The
x₁-sum is then ≤ Σ_{x₁∈𝒮_y}ℓ(x₁)/(Δ′Z) + #{x₁ ∈ 𝒮_y : |x₁| ≤ X₁}·O(log X₂)λ_min/(Δ′Z).
- The first term is ≤ p₁(1 + o(1)) + #{x₁}/Z, and #{x₁}/Z ≤ 2X₁Δ′/|Π₀| = N^{−Ω(μ)}.
- For the second, #{x₁ ∈ 𝒮_y, |x₁| ≤ X₁} ≤ 2Ψ(X₁, y) = 2X₁c₁, and Δ′Z ≥ |Π₀|/2 ≥ s₁s₂/(64ϑ) (Lemma PG(c)). So the
  second term is ≤ O(log X₂)·2X₁c₁·X₂y^{−C_X}ϑ^{−1}·64ϑ/(s₁s₂) ≤ D^{2}(d + 1)²y^{−C_X}c₁ ≤ y^{−1}c₁c₂, by the choice of
  C_X. The thinness cancels against the 1/ϑ in λ_min.
- p₁ ≤ 12(d + 1)c₁ by OPV-A step 6 (from PG(d) alone).

The partner bound carries the factor e^{O(λ)}(y^{C_X}ϑΔ′)^γ. ∎

---

## 8. The per-pair covariance bound (Prop PD-A; obligation O5) [proved]

Fix (g, pair) non-thin, and the blocks and cells of Lemma BL. Define the probability ν̂ on cells ∪ {∗}:
- ν̂(R, δ) := Σ_yπ(y)1[R interior for y]·n(R, δ)/Z, for δ ≤ S₀ (S₀^p), R in the block partition;
- ν̂(∗) := 1 − Σν̂(R, δ);
- Ψ_i(R, δ) := d_i(I_i, δ) (Lemma CF), and Ψ_i(∗) := 0.

Here Z = Z_y := #(L′ ∩ Π_y), kept with its index y throughout; Z_y/Z_{y′} = 1 + O(q/H) (c12alg TP-W, step 5:
translates contain equally many fundamental domains) is used only in error estimates.

**Prop PD-A.** With 𝒱 as in Lemma OPV-A and 𝒱* := 𝒱 + Ch^{1/2} + e^{−λ}:
  |p₁₂ − p₁p₂| ≤ C𝒱*c₁c₂ + 𝔈,  where  𝔈 := Σ_yπ(y)Z^{−1}Σ_{R meets Π_y}Σ_{δ≤S₀}|Err(R, δ)|.

*Proof.*
1. **Decomposition.** For each y, split L′ ∩ Π_y into:
   - the points in interior rectangles with δ ≤ S₀;
   - the points in boundary rectangles with δ ≤ S₀;
   - the exceptional points ℰ (which include everything with |x_i| < X_i⁻).

   By Lemma CF(iii) on interior cells,
   p₁₂ = E_ν̂[Ψ₁Ψ₂] + Z^{−1}Σ_{int}Err + P_∂ + P_ℰ + tiny,
   where:
   - P_ℰ := ν(ℰ ∩ both smooth) ≤ e^{−λ}c₁c₂ (Lemma EXC-A);
   - 0 ≤ P_∂ ≤ Σ_yπ(y)Z^{−1}Σ_{∂ for y}[n(R, δ)Ψ₁Ψ₂(1 + tiny) + |Err(R, δ)|]. By positivity, a boundary
     rectangle's points in Π_y are at most the points of the full rectangle, to which Lemma CF applies, since its
     sides are ≥ h_*X_i⁻ ≫ N^{μ/2}Δ′.
   - "tiny" = N^{−Ω(μ)}c₁c₂ collects the factors of Lemma CF(i).

   By Lemma CF(ii), p_i = E_ν̂[Ψ_i] + B_i + O(e^{−λ}c_i) + tiny, with 0 ≤ B_i the boundary single parts.
2. **Identity.**
   p₁₂ − p₁p₂ = Cov_ν̂(Ψ₁, Ψ₂) + Z^{−1}Σ_{int}Err + P_∂ + P_ℰ − (E_ν̂Ψ₁)B₂ − B₁(E_ν̂Ψ₂) − B₁B₂ + O(e^{−λ}c₁c₂) + tiny.
   This is c12rat Prop PD's identity, with the atom ∗ carrying boundary and exceptional mass.
3. **Covariance.** |Cov_ν̂(Ψ₁, Ψ₂)| ≤ (Var_ν̂Ψ₁·Var_ν̂Ψ₂)^{1/2}, and Var_ν̂Ψ₁ ≤ E_ν̂(Ψ₁ − c₁)². The first marginal of
   ν̂ restricted to interior cells is dominated by the law of (I(x₁), δ(x₁)) under ν: an interior cell's ν̂-mass is
   the ν-mass of its points. So E_ν̂(Ψ₁ − c₁)² ≤ E_ν[(Ψ₁/c₁ − 1)²1_{non-exc}]c₁² + ν̂(∗)c₁² ≤ (𝒱 + Ch + e^{−λ})c₁²,
   by Lemma OPV-A, Lemma BL(ii) and Lemma EXC-A.
4. **Boundary terms** (Hölder against fourth moments on the joint space (y, x); this is the transposition O5).
   Let ∂ denote the boundary rectangles (for the relevant y), and write the joint measure
   μ_∂(y, R, δ) := π(y)Z_y^{−1}n(R, δ)1[R boundary for y]. Its total mass is b := Σ_yπ(y)Z_y^{−1}Σ_∂n ≤ Ch(1 + tiny)
   (Lemma BL(ii)). Then, by Hölder on the joint space:
   - B₁ ≤ Σμ_∂Ψ₁(1 + tiny) ≤ b^{3/4}(Σμ_∂Ψ₁⁴)^{1/4};
   - Σμ_∂Ψ₁Ψ₂ ≤ b^{1/2}(Σμ_∂Ψ₁⁴)^{1/4}(Σμ_∂Ψ₂⁴)^{1/4}.
   Every boundary rectangle for y lies in the convex Π_y⁺ (Lemma BL(ii)), so
   Σμ_∂Ψ₁⁴ ≤ Σ_yπ(y)Z_y^{−1}Σ_{R⊂Π_y⁺}n(R, δ)Ψ₁⁴ ≤ (1 + 40h)(1 + tiny)E_{ν⁺}[Ψ₁⁴1_{non-exc}] ≤ Cc₁⁴,
   by Lemma OPV-A for the mixture ν⁺ (PG(f)); the moment is the averaged one over (y, x), which is what OPV-A bounds.
   Moments, not maxima: no spread factor enters. So B_i ≤ Ch^{3/4}c_i and P_∂ ≤ Ch^{1/2}c₁c₂ + (boundary part of 𝔈).
5. **Total.** Collect, using c_i ≤ 2p_i (OPV-A step 6) and E_ν̂Ψ_i ≤ p_i ≤ 12(d + 1)c_i; in the products with B_i the
   factor d is absorbed by h^{3/4} = D^{−15/2}. ∎

**Corollary PD-A.** Under the hypotheses of Lemma ML, for non-thin pairs,
  |A₁₂ − A₁A₂| ≤ 16C𝒱*A₁A₂ + 𝔈 + 60qN^{−η₀/2}.
*Proof.* Lemma ML(b), Prop PD-A, and c_i ≤ 2p_i ≤ 2A_i + 14qN^{−η₀/2}. Also A_i ≥ c_i/4 ≥ N^{−η₀/4}, by OPV-A and
L[α]^{−O(1)} ≫ N^{−η₀/4}. ∎

The only first-order quantities are γ log d, γλ (local) and E_FO. They enter only through 𝒱, i.e. squared. The
boundary costs h^{1/2} = D^{−5} (Hölder, step 4), and the exceptional sets e^{−λ}.
No pair main term is computed and no transversality is used. The thinness enters only through h_* (the grid) in §9.

---

## 9. Character errors: the weighted energy inequality (obligation O2) [proved]

**Lemma EN-A (library energies).** For ξ ∈ 𝒳 and integer 1 ≤ δ ≤ S₀ with ξ/δ ≥ X_min^{1/2}, put
u_{ξ,δ}(m) := 1[1 ≤ m ≤ ξ/δ, P⁺(m) ≤ y] and U_{ξ,δ} := Ψ(ξ/δ, y). If Q ≤ (ξ/δ)^{η_D}, then
  Σ_{𝔨≤Q}φ(𝔨)^{−1}Σ_{χ≠χ₀ (mod 𝔨)}|û_{ξ,δ}(χ)|² ≤ ε_lib·U²_{ξ,δ},  ε_lib := C_D(H(u_min/2)^{−δ_D} + y^{−δ_D}).

*Proof.* |û(χ)| ≤ U, so |û(χ)|² ≤ U|Ψ(ξ/δ, y; χ)|. Every non-principal χ mod 𝔨 ≤ Q has conductor in
(1, Q] ⊂ (1, (ξ/δ)^{η_D}]. Fact D (Drappeau, Compositio 151 (2015), Lemme 5, with A = 1/2) bounds
Σ_𝔨φ(𝔨)^{−1}Σ_χ|Ψ(ξ/δ, y; χ)| by C_DΨ(ξ/δ, y)(H(u)^{−δ_D} + y^{−δ_D}). Its hypotheses (log x)^{c_D} ≤ y ≤ x and Q ≤ x
hold. ∎ (This is c12rat Lemma EN without the P₀-free restriction, which was needed there only for primitivity.)

In our ranges:
- children: Q = 2S₂² = L[α]^{O(1)}, ξ/δ ≥ X⁻/S₀ = x^{1−o(1)}, u ≍ D^α, so ε_lib ≤ e^{−c_εD^α/λ²};
- parents: Q = 2V₁²/q = x^{o(1)} (§1.3), u₂ ≍ D^{3−7α}, so ε_lib^p ≤ e^{−c_εD^{3−7α}/λ²}.

**From cells to the library.** For a cell (R, δ) with δ ≤ S₀ and a block I_i = [ξ_lo, ξ_hi] (or its negative):
- â_i(χ) = û_{ξ_hi,δ}(χ) − û_{ξ_lo,δ}(χ) (signs: χ(−1) factors have modulus 1). So
  Σ_{χ≠χ₀}|â_i|² ≤ 2Σ|û_{ξ_hi,δ}|² + 2Σ|û_{ξ_lo,δ}|².
- Lemma FO (as in OPV-A step 1) and Fact dBT1 give A_i ≥ (|I_i|/ξ_hi)U_{ξ_hi,δ}/(C(log P₀)^C) ≥ h_*U_{ξ_hi,δ}/λ^C.
  This uses g_𝔪(β) ≥ c(log P₀)^{−C}e^{−o(1)}, and U_{ξ_lo} ≤ U_{ξ_hi}.

Hence Ẽ_i ≤ λ^Ch_*^{−2}[𝔅_{ξ_hi,δ}(Δ_δ) + 𝔅_{ξ_lo,δ}(Δ_δ)], where 𝔅_{ξ,δ}(𝔨) := U_{ξ,δ}^{−2}Σ_{χ≠χ₀ (mod 𝔨)}|û_{ξ,δ}(χ)|².

**Lemma FAM-A (the family).** Let the terms be t = (𝔮, e, g, ordered non-thin pair, y, (R, δ) meeting Π_y, i ∈ {1, 2},
endpoint ε ∈ {lo, hi}), with library element j(t) := (ξ_ε(I_i), δ), modulus 𝔨(t) := Δ_δ = Δ′/δ, and coefficient
  c_t := (w_{𝔮,e}/W)·m_{𝔮,e}^{−2}·Pr_e(g)Rat₁Rat₂(g)·π(y)·Z^{−1}n(R, δ)Ψ₁Ψ₂·λ^Ch_*^{−2}.
Then:
- (a) Σ_{(𝔮,e): s≤T₀}(w/W)m^{−2}E_e[Σ_{non-thin}Rat₁Rat₂𝔈] ≤ Σ_tc_t𝔅_{j(t)}(𝔨(t));
- (b) 𝔇 := Σ_j max_𝔨 φ(𝔨)H_j(𝔨) ≤ e^{O(D^θλ²)}, where H_j(𝔨) := Σ_{t:j(t)=j,𝔨(t)=𝔨}c_t.

Consequently, by c12rat **Theorem WA** with ε = ε_lib, the left side of (a) is ≤ ε_lib𝔇 ≤ e^{−cD^α/λ²}.

*Proof.*
- (a) is Lemma CF(iii) with the previous display (A₁A₂/φ ≤ n(R, δ)Ψ₁Ψ₂(1 + tiny)).
- (b), step 1: **bounds per pair.** Fix (𝔮, e), j = (ξ, δ) and a pair with Δ′ = 𝔨δ. For fixed (g, y, i), the
  rectangles whose i-th block has ξ as an endpoint lie in Π_y⁺ and in the ≤ 2 blocks adjacent to ξ. With the pointwise
  bound Ψ_i ≤ e^{O(λ)}δ^γc_i (OPV-A step 6; (X_i/|ξ|)^γ ≤ B^{C_Xγ} = e^{O(λ)}) and the last claim of PG(e)
  (ν⁺(x_i ∈ I, δ | x_i) ≤ (|I|/δ + 1)·3/W_i with |I| ≤ 4h^{1/2}ξ and ξ ≤ 4X_i ≤ 8(d + 1)W_i):
  Σ_{such R}Z^{−1}n(R, δ)Ψ₁Ψ₂ ≤ e^{O(λ)}δ^{−1+2γ}c₁(g)c₂(g)·dh^{1/2}·1[X_{i,min}⁻/2 ≤ ξ ≤ 4X_{i,max}].
  Here X_{i,min} and X_{i,max} are the extremes of X_i(g) over the cell e. Since |G_g(P_i)| ∈
  [𝔐|P_i|N^{−η₀/4}, 2𝔐|P_i|], the range [X⁻_min/2, 4X_max] has logarithmic length O(η₀ log N + log B + log d)
  = O(D^{3α−1/2}λ).
- Step 2: **the g-integral.** E_e[Rat₁Rat₂c₁(g)c₂(g)] ≤ Ξ^Cϖ_{P₁}ϖ_{P₂}, with the per-point weight
  ϖ_P := ϱ_B(𝔐|P|)ρ̂_B((2d + 1)H𝔐|P|^d/q). This uses:
  - c_i(g) = ρ̂(X_i(g)) ≤ C(N^{η₀/4})^γρ̂(X_i,max) ≤ Ξρ̂_B((2d + 1)H𝔐|P_i|^d/q) (Fact dBT1; X_i ≤ 4(d + 1)H𝔐|P_i|^d/q,
    and a factor ≤ 2 in the argument of ρ̂ is absorbed);
  - Pr_e[Rat₁Rat₂] ≤ ΞPr_e[Rat₁]Pr_e[Rat₂] (c12rat **Lemma PUB** with weights ≡ 1; its proof is the line/BTS argument
    and holds verbatim for unit weights);
  - Pr_e[Rat_i] ≤ Ξϱ_B(𝔐|P_i|).

  There is no pointwise claim about Rat₁Rat₂ in g; only its integral is used.
- Step 3: **normalisation.** For s(𝔮) ≤ T₀, c12rat **Lemma SPR** (children) gives max_Pϖ_P ≤ Ξ min_Pϖ_P over Box(𝔮).
  This covers the spread of ϱ_B(𝔐|P|) by (4S₂T₀)^γ, and of the algebraic factor by ((4S₂T₀)^d)^γ =
  e^{O(D^{1−2α}λ² + D^{1−3α+θ}λ³)} ≤ Ξ. Prop CLB, per point (c12alg §4 steps 2–3, which bound each P′'s contribution
  by ϱ(𝔐|P′|)ϱ(X_{F2}(P′))(log N)^{−C}), gives m_{𝔮,e} ≥ Σ_Pϖ_P(log N)^{−C}. Hence ϖ_{P₁}ϖ_{P₂}/m² ≤ Ξ³/n², with n = #Box(𝔮) ≥ S₂².
- Step 4: **multiplicity.** For each Δ′ ≥ 1 the box has ≤ 10S₂² + 2S₂ ordered pairs with coordinate determinant Δ′
  (c12rat **Lemma MULT**), and 𝔨 ≤ Δ′/δ ≤ 2S₂²/δ. So φ(𝔨)·#{pairs : Δ′ = 𝔨δ}/n² ≤ 24/δ.
- Step 5: **sum.** Since max_𝔨 Σ_{(𝔮,e)} ≤ Σ_{(𝔮,e)} max_𝔨,
  𝔇 ≤ Σ_{(𝔮,e)}(w/W)Ξ^Ce^{O(λ)}dh^{1/2}λ^Ch_*^{−2}·Σ_δ24δ^{−2+2γ}·#(𝒳 ∩ ℛ_{𝔮,e}),
  where ℛ_{𝔮,e} is the **label-wide** scale range: the union over the label's pairs and over g ∈ Ω_e of
  [X⁻_min/2, 4X_max]. Its logarithmic length is O(d log(S₂T₀) + η₀ log N + C_X log B) = D^{O(1)}λ^{O(1)}
  (|P|^{d−1} varies over a box with s(𝔮) ≤ T₀ by (4S₂T₀)^{d}). Here Σ_δδ^{−2+2γ} ≤ 2 and
  #(𝒳 ∩ ℛ) ≤ D^{O(1)}/h_*. Also Σ_{(𝔮,e)}w/W = 1. So 𝔇 ≤ Ξ^Ch_*^{−3}D^{C} = e^{O(D^θλ²)}, since
  h_*^{−1} = 64(d + 1)ϑ₀D^{10} = e^{O(D^θλ²)}. ∎

Why this is simpler than c12rat's Lemma FAM: the modulus 𝔨 = Δ′/δ is a **coordinate** quantity (≤ 2S₂²) that does
not contain q. So for fixed 𝔨 every label can contribute, and the crude bound max_𝔨 Σ_𝔮 ≤ Σ_𝔮 max_𝔨 already
suffices. No comparison of label weights across residues r or primes q is needed (Lemma DOM is not used here), and
the g-dependence of the windows costs only the logarithmic length of the scale range (step 1).

**Parents.** The same, with y = B₂, the grid 𝒳^p and h_*^p, 𝔨 ≤ 2V₁²/(qδ), and the following normalisation package
.
Put R_typ := #{n ≤ 𝔐V₁ : n ∈ J·𝒮_B}/(𝔐V₁) (a purely rational density) and ρ_a(q) := ρ̂_{B₂}(X_F/q).

**Lemma W-P (parent weights).** For every label and cell, and parents P, P₁ ≠ P₂ ∈ 𝒫 ∩ Λ_𝔮, n_p := #(𝒫 ∩ Λ_𝔮):
- (P1) E_e[R(P₁; g)R(P₂; g)] ≤ Ξ_pD^CR_typ²;
- (P2) c_i(g) ≤ Cρ_a(q) for every admissible g;
- (P3) z_{𝔮,e} ≥ n_pR_typρ_a(q)D^{−C}.
Hence E_e[R₁R₂c₁c₂]/z²_{𝔮,e} ≤ Ξ_pD^C/n_p².

*Proof.*
- (P1) c12rat Lemma PUB's line argument with **unit (rational-only) weights** and target set J·𝒮_B ∩ [𝔐V₁/32, 2𝔐V₁].
  Along a line g* + jP₁ (G(P₁) = n₁ fixed), G(P₂) runs over an AP of difference K = |det(P₁, P₂)| ≤ 2V₁² with
  ≍ δ_e𝔐/|P₁| terms. For ℓ ∈ J with ℓ ∤ K, the n = ℓm in the AP have m in an AP of difference K, and BTS bounds the
  B-smooth m by e^{O(λ)}(K·L^{O(1)})^γ times the expected count; summing over ℓ gives Ξ_p·R_typ per unit length.
  At most one ℓ ∈ J divides K (K ≤ 2V₁² < Y² — since V₁² ≍ YS₁² and S₁² = L[α] < Y). For that ℓ, either ℓ ∤ G(P₂)
  on the whole line (contribution 0), or ℓ | G(P₂) throughout. In the latter case ℓ | n₁ as well (n₁ ≡ ιG(P₂) mod K).
  Since R_typ ≍_{D^C} ϱ_B(𝔐V₁/Y)/log Y, the lines with ℓ | n₁ carry a fraction ≤ D^C(log Y)/ℓ of the first-point
  mass R_typ, and on them the density of ℓ𝒮_B in the AP is ≤ Ξ_pϱ_B(𝔐V₁/Y) ≤ Ξ_pD^CR_typ log Y (BTS). So the same-ℓ
  total is ≤ Ξ_pD^CR_typ²(log Y)²/Y.
  Integrating n₁ against the line lengths (Hölder, as in PUB) gives (P1).
- (P2) c_i(g) = ρ̂_{B₂}(X_i(g)) with X_i(g) = (d + 1)(2H + 1)|G_g(P_i)||P_i|_∞^{d−1}/q. Here |G| ∈ [𝔐V₁/32, 2𝔐V₁] and
  |P_i|_∞ ∈ ((1 − 1/d)V₁, V₁] (b-window), so X_i(g) is within a factor 64e(d + 1) of X_F/q, which costs
  (64ed)^{γ₂} = 1 + o(1) (Fact HT3).
- (P3) At least half the points of 𝒫 ∩ Λ_𝔮 are *active*: |G_g(P)| ≥ 𝔐V₁/16 on all of Ω_e (c12rat Lemma SPR (parents),
  angle argument). For active P:
  - Pr_e[R(P)] ≥ R_typD^{−C}: G_g(P) over Ω_e ranges over an interval of relative length ≍ δ_e ≥ B^{−5/12} inside
    [𝔐V₁/16, 2𝔐V₁]. Hildebrand's short-interval theorem on each slice ℓ·(interval/ℓ) and Lemma FIB (Möbius for the
    primitivity of g, factor ≥ 1/(C log log N)) give a J·𝒮_B-density ≥ R_typD^{−C}.
  - a_P(g) ≥ ρ_a(q)D^{−C}: DW(a)⁺ with the root coset, Hildebrand at B₂ (relative length ≫ 1/d at scale X_F/q), and
    the root-conditioned acceptability subtraction of §11 (≤ D^C/B ≪ ρ_a(q); at α = 3/8 this uses (SQ)). Hildebrand
    is applied **after division by q**: at the scale X_F/q the DW(a)⁺ intervals have relative length ≥ 1/(64d(d + 1))
    ≫ B₂^{−5/12}.
  Sum over active P. Inactive points are only bounded above; no individual lower bound is claimed for them. ∎

Multiplicity: c12rat Lemma MULT (parent) gives ≤ 5 partners P₂ of each P₁ with a given **signed** det, so ≤ 10 for a
given |det|, and φ(𝔨)·#{pairs : Δ′ = 𝔨δ}/n_p² ≤ 20V₁²/(qδn_p) ≤ 80d log V₁/δ, as
n_p ≥ V₁²/(4dq log V₁). The scale range has log-length O(C_X log B₂ + log d).

So 𝔇_p ≤ Ξ_p^C(h_*^p)^{−3}D^C = e^{O(D^{θ_p}λ²)}, and the parent character total is ≤ ε^p_lib𝔇_p ≤ e^{−cD^{3−7α}/λ²}.
This is where the **separate parent thresholds** (O3) matter. With the child h_* the grid count would be
e^{O(D^{2α−1/2}λ²)}, which exceeds the parent saving for α ≥ 7/18 (`params.theta(7/18) ≥ 3 − 7α`, tested).

---

## 10. Thin pairs, skewed labels [proved]

**Thin child pairs** (ϑ > ϑ₀) in a label with s(𝔮) ≤ T₀ form a fraction ≤ 60s(𝔮)/ϑ₀ ≤ 60T₀^{−2} of the ordered
pairs (c12rat **Lemma TF**). For them,
|A₁₂ − A₁A₂| ≤ A₁₂ + A₁A₂ ≤ Ξc₁c₂ + C(d + 1)²c₁c₂ + tiny (Lemma PU-A for all ϑ; A_i ≤ p_i + tiny ≤ 12(d + 1)c_i
by OPV-A step 6, which needs no non-thin partner). So their contribution to
|ϕ^alg_{𝔮,e}| is ≤ m^{−2}Σ_{thin}E_e[Rat₁Rat₂·2Ξc₁c₂] ≤ 60T₀^{−2}·Ξ^C ≤ T₀^{−1}, by steps 2–3 of Lemma FAM-A.

**Thin parent pairs** (ϑ > ϑ₀^p). For every P, the fraction of thin partners in 𝒫 ∩ Λ_𝔮 is ≤ 3/ϑ₀^p + O(q/V₁)
(c12rat Lemma TF, parents). Their contribution to |σ^alg| is ≤ (3/ϑ₀^p + O(q/V₁))Ξ_p^C·polylog ≤ e^{−cD^{θ_p}λ²}.

**Skewed labels (iii).** For every child label and cell, without restriction on s(𝔮):
- E_e[N_𝔮(N_𝔮 − 1)] = Σ_{P₁≠P₂}E_e[Rat₁Rat₂A₁₂]
  ≤ Σ_{P₁≠P₂}ΞE_e[Rat₁Rat₂c₁c₂] + n²·tiny
  ≤ Ξ^C(Σ_Pϖ_P)²,
  by Lemma PU-A and step 2 of FAM-A, which is **per point** and uses no comparison across the box.
- m_{𝔮,e} ≥ Σ_Pϖ_P(log N)^{−C} (Prop CLB per point).
- So ϕ_{𝔮,e} ≤ 1/m + E_e[N(N − 1)]/m² ≤ Ξ^{C}(log N)^C.

c12rat **Lemma DOM**(i)–(ii) and **Lemma LS** (Prop SK(i)) give Σ_{(𝔮,e): s>T₀}w_{𝔮,e} ≤ (3/T₀ + 2/B)(log N)^Cμ₁, and
DOM(iii) gives W ≥ (1 − L^{−c})μ₁. Hence Σ_{s>T₀}(w/W)ϕ ≤ 2(3/T₀ + 2/B)Ξ^C(log N)^{2C} ≤ T₀^{−1/2}, once C_T is
large. ∎

This replaces c12alg §3.4(iii)'s "max/typical factor Ξ", which is false for very skewed labels: |P′|^{d−1} varies by
up to (S₂q)^{d−1}, a factor e^{O(D^αλ²)}. The per-point form avoids it.

---

## 11. Parents: the differences (obligations O3, O4) [proved]

The parent block (ii) is §§3–10 with:
- **Smoothness y = B₂**, u₂ ≍ D^{3−7α} → ∞. Hildebrand's theorem applies at B₂ (log B₂ ≍ D^{4α−1}λ ≫ λ^{5/3+ε}).
  The block widths satisfy h_*^p ≥ B₂^{−5/12} (θ_p < 4α − 1). Lemma FO holds verbatim at y = B₂ with S₀^p in place
  of S₀. Its proof uses y only through:
  - Hildebrand's range;
  - the k ≤ 3 HT3 steps at ratios ≤ y²;
  - BTS (any y);
  - the Rankin tails y^{−0.9/log P₀}(log P₀)^C and y^{−1/2}e^{ω(Δ′)P₀^{−0.49}}. These are smaller than every power of
    D at y = B₂, because ω(Δ′)P₀^{−0.49} ≤ D^{(1−α) − 0.343} and (1 − α) − 0.343 < 4α − 1.
- **Moduli** Δ′ = |det(P₁, P₂)|/q ≤ 2V₁²/q = L[1 − α]^{O(1)}, which is ≫ B₂^{O(1)} but x^{o(1)}. Nothing in §§5–9 iterates
  HT3 over divisors of Δ′: classes δ ≤ S₀^p use Lemma FO, larger ones use Lemma EXC-A. So c12alg's truncation at
  B₂^{1/2} is not needed.
- **Thresholds** ϑ₀^p and h_*^p, built from the parent losses Ξ_p only (O3). The child filter spread
  e^{O(D^{2α−1/2}λ²)} never occurs: the parent filter fixes |G(P)| to within a factor 64.
- **No skew:** parent sets are Euclidean squares.
- **Acceptability (O4).** acc(P) is not a function of the value x alone. Put s_P(g) := Pr[F(P)/q ∈ 𝒮_{B₂} | g, R_𝔮]
  and s₁₂ likewise. **Root-conditioned subtraction**.
  Conditionally on (g, R_𝔮), with c uniform (c12alg Lemma RE; uprime Lemma SQF's proof):
  - p ∈ (B, B₂], p ≠ q: Pr[p² | F(P) | g, R] ≤ p^{−2}(1 + o(1)). Since p ≠ q and q ∤ b^{d−1}G(P), the coefficient
    c_{d−1} varies F(P) through all classes mod p²q with tiny endpoint error, so the events mod p² and mod q are
    independent up to (1 + p²q/(2H + 1)).
  - p = q: Pr[q² | F(P) | g, R] ≤ 2/q (F(P) is equidistributed mod q² by the same coefficient; conditioning on q | F(P)
    leaves q classes mod q²).
  - q′ ∈ (B, B₂], q′ ≠ q, discriminant/leading coefficient: Pr[q′ | F(P) and q′ | f_d disc f | g, R] ≤ 2e(d − 1)/q′²
    (uprime SQF(b), applied mod q′, with CRT on c mod qq′; **both** q′ | F(P) and q′ | f_d disc f are required for q′
    to be an excluded child).
  - q′ = q: Pr[q | f_d disc f | g, R] ≤ 2e(d − 1)/q. This is SQF(b)'s inner computation: given F(P) ≡ 0 (mod q),
    F = L_P·H with H uniform on an affine hyperplane, and the event lies in {H(1, 0) = 0} ∪ {H(P) = 0} ∪ {Disc H = 0}.
  - primes of N in (B, B₂] dividing F(P): ≤ Σ_{q′|N}(1/q′ + 1/(2H + 1)) ≤ log N/B.

  Total: Pr[s_P ∧ ¬acc(P) | g, R] ≤ Σ_{p>B}(1 + 2ed)p^{−2} + (2 + 2ed)/q + log N/B ≤ D^C/B. Hence
  |(a₁₂ − a₁a₂) − (s₁₂ − s₁s₂)| ≤ 6D^C/B.
  Relative to s₁s₂ ≥ ρ̂_{B₂}(X_F)²D^{−C} (OPV-A step 6 lower bound, and Hildebrand at B₂), this is
  ≤ D^C/(Bρ̂_{B₂}(X_F)²) = exp(−aD^αλ + 2c₂D^{3−7α}λ(1 + o(1))).
  That is L^{−a+o(1)} for α > 3/8. At α = 3/8 it is ≤ D^{−3/4} iff **a > 2c₂**, the endpoint constraint. Then
  Prop PD-A is applied to (s₁₂, s₁, s₂). The relative bound transfers via a_i ≥ s_i − D^C/B ≥ s_i/2. Summed against
  ΣR₁R₂a₁a₂ ≤ z(g)², this involves no Ξ_p.

Everything else is word for word. The output is
|σ^alg_{𝔮,e}| ≤ C𝒱*_p(1 + σ^rat_{𝔮,e}) + (thin) + (character)/z² + tiny, with 𝒱*_p ≤ C(D^{−2α} + D^{−2(3−7α)})λ^C.

---

## 12. Proof of Theorem AB₂ [proved]

**(i).** Fix (𝔮, e) with s(𝔮) ≤ T₀. Split the ordered pairs into thin and non-thin.
- **Non-thin pairs** (Corollary PD-A):
  m²|ϕ^alg| ≤ 16C𝒱*E_e[Σ_{P₁≠P₂}Rat₁Rat₂A₁A₂] + E_e[ΣRat₁Rat₂𝔈] + n²·60qN^{−η₀/2}.
  - The first term is ≤ 16C𝒱*E_e[m(g)²] = 16C𝒱*m²(1 + ϕ^rat_{𝔮,e}). The full square dominates the off-diagonal
    sum, since all terms are ≥ 0.
  - The last term is ≤ N^{−η₀/3}m², since n²/m² ≤ (S₂²/m₂)²(log N)^C = L^{O(1)} (Prop CLB).
- **Thin pairs:** ≤ T₀^{−1}m² (§10).

Multiply by w_{𝔮,e}/W and sum over s(𝔮) ≤ T₀. The character terms total ≤ e^{−cD^α/λ²} (Lemma FAM-A and Theorem WA).
So
  Σ_{s≤T₀}(w/W)|ϕ^alg| ≤ 16C𝒱*(1 + Σ_{s≤T₀}(w/W)ϕ^rat) + T₀^{−1}·(1 + o(1)) + e^{−cD^α/λ²} + N^{−η₀/3},
with 𝒱* ≤ CD^{−2α}λ^C. Absorb the last two terms into T₀^{−1} (α > θ).

**(ii).** The same over all (𝔮, e), with §11.

**(iii).** §10. ∎

**Status of c12alg's obligations.**

| obligation | where | status |
|---|---|---|
| O1 joint moment lemma | Lemmas LF-A, OPV-A, EXC-A (incl. ν⁺, mixture, prime powers, large classes, near zero) | proved here |
| O2 weighted energy inequality | Lemma EN-A, Lemma FAM-A + c12rat Theorem WA; parents via Lemma W-P | proved here; parent normalisation rewritten in v2 |
| O3 separate parent thresholds | §1.3 (ϑ₀^p, h_*^p from Ξ_p), §9 parents, §11 | proved |
| O4 parent acceptability at 3/8 | §11 | discharged by the constraint |
| O5 boundary transposition | Lemma BL(ii) (convex Π^∘, Π⁺), Prop PD-A step 4 | proved |

**What Theorem AB needed that the sketch did not have.** c12alg's route (Lemma FB with Harper's BV on fibres, and
Lemma PP computing the pair main term to first order) is **not** used. Its grid/step bookkeeping, the short-fibre
corners and the parent truncation are all replaced by the cell device. Harper's Theorems 1–2 are no longer inputs;
Drappeau's Lemme 5 is the only BDH-type input, as in c12rat.

---

## 13. Theorem 2L-final(α) [proved, internal; conditional only on the unrefereed chain]

**Theorem 2L-final(α).** Let α ∈ (3/8, 2/5) be fixed, and let a, κ, d₀, b₂ > 0 be rationals. FACTOR_2L in design
2L_DW(α) (c12alg §1.1, with the search of multidescent Theorem 2L and POSTPROCESS-DW) is a Las Vegas algorithm. It
outputs the complete prime factorisation of every integer N ≥ 2 in expected time L_N[α, c(α)], where
  c(α) = max(2a, a + s + max(a, σ₁ + σ₂), (a + s)/(1 − 2^{−α})),  s = 2κ(1 + C_γ),
  σ₁ = α/(ad₀),  σ₂ = (α/a)(4/d₀ + d₀b₂),
and C_γ = C_γ(a, κ, d₀, b₂) is the polylog exponent of the level-2 relation supply γ₂ ≥ (log N)^{−C_γ}.

*Proof.* c12alg's proof of Theorem 2L-DW(α) (c12alg §5), whose two hypotheses are now supplied:
1. Lemma VD gives A ≤ Σ(w/W)/m + Σ_{s≤T₀}(w/W)(|ϕ^alg| + ϕ^rat) + Σ_{s>T₀}(w/W)ϕ, and likewise for B with
   1/z, σ^alg and σ^rat. Then:
   - the one-point terms are ≤ e^{−λ²}(log N)^C (Prop CLB) and ≤ e^{b₂D^{4α−1}λ}(log N)^C/L[1 − α] (the identity
     w/z = Pr(R)/#J, c12alg §5 step 1);
   - **Hyp RAT(α)** is c12rat's Corollary RAT(α), from Theorem RAT (RAT-1) for s ≤ T₀ and (RAT-2):
     Σ_{s≤T₀}(w/W)ϕ^rat ≤ CD^{−2α}λ^C + Σ(w/W)/m and Σ(w/W)σ^rat ≤ C(D^{−2α} + D^{−2(3−7α)})λ^C + Σ(w/W)/z;
   - **Theorem AB** (§2) gives the algebraic parts and (iii).

   **Weight bridge to RAT**. c12rat states RAT with the idealised weights a_{P′}(|G_g(P′)|) (c12alg
   §5), while Lemma VD's ϕ^rat uses the true A_{P′}(g). By DW(a)⁺, |A_{P′}(g) − a_{P′}(|G_g(P′)|)| ≤ 4qN^{−η₀/2} for
   every g, which is ≤ N^{−η₀/3}a_{P′}. Hence m(g) = m̃(g)(1 + O(N^{−η₀/3})) pointwise, and
   ϕ^rat_true ≤ (1 + N^{−η₀/4})ϕ^rat_ideal + N^{−η₀/4}.

   *Parents*. RAT's idealised parent weight is a_P^{id}(t) := the B₂-smoothness probability of
   the DW(a)⁺ model law at |G| = t, restricted to qZ and divided by q, i.e. the value-only weight; c12rat §1.2's
   a_P(t) := Pr[F(P) acceptable, q child | g, R_𝔮] is read in this value-only form, as c12alg §5 (R2) requires. The
   true weight a_P(g) satisfies |a_P(g) − a_P^{id}(|G_g(P)|)| ≤ |a_P − s_P| + |s_P − a_P^{id}| ≤ D^C/B + 4qN^{−η₀/2}
   (§11's root-conditioned subtraction, and the TV of DW(a)⁺ with the root coset). Relative to a_P^{id} ≥ ρ_a(q)D^{−C}
   (Lemma W-P, P3) this is ε_p ≤ D^C/(Bρ_a(q)) + N^{−η₀/3}. The elementary perturbation estimate gives
   σ^rat_true ≤ (1 + O(ε_p))σ^rat_ideal + O(ε_p). In the open window ε_p = L^{−a+o(1)}; at α = 3/8,
   ε_p ≤ L[3/8, −(a − c₂) + o(1)] ≤ D^{−1}, by a > 2c₂.

   **(A1) is not needed**. c12alg §5 made a² > 4α/d₀ conditional on c12rat using pc3 LAT(d). c12rat's
   Theorem RAT does not: Lemma MULT uses K = qΔ only, which holds without (A1), and Lemma TYPE treats q as a "large"
   prime of K (classes with q | n₁n₂ are in the exceptional set of Lemma EXC(a) and contribute 0), so v_q(K) ≥ 2 is
   covered. On the algebraic side q is divided out (Cor TP-q). So the illustrative parameters below, which violate
   (A1), are admissible.

   So A ≤ D^{−α}λ^C and B ≤ D^{7α−3}λ^C. (In fact A ≤ D^{−2α}λ^C and B ≤ (D^{−2α} + D^{−2(3−7α)})λ^C.)
2. c12alg Prop PB (with Lemma η, W/μ₁ ≤ 3u₂) gives a dead-incidence mass ≤ μ₁/2, so E N₂ ≥ μ₁/2 (pc3 Prop DC_e,
   multidescent Prop UB).
3. U′-DW(α) (c12alg §2.7: uprime's FM/LD*/D-U′, CS′, SQF, child scan, recursion and IR, with DW(a)⁺ in place of
   GM-U′) gives E[N₁²] ≤ (log N)^C(μ₁² + μ₁) and μ₁ ≥ (log N)^{−C}. Paley–Zygmund then gives γ₂ ≥ (log N)^{−C_γ}.
4. Theorem K-U and Prop select-U (lv5unif §8), with Lemma reject-DW (c12alg §2.4), give the split acceptance
   (γ₂/(16 log N))^{2k}/(32k²) = L[α, −2κ(1 + C_γ)].
5. Rows (nested Lemma RW), POSTPROCESS-DW (c12alg §2.8; lv5unif Prop LV6-U), Theorem DD with Lemma T_multi, the search
   cost L[α, σ₁ + σ₂ + o(1)] (uprime §6; escapeR2 §7.3), and Lemma F_α (lv5unif §10) finish the proof, as in
   multidescent Theorem 2L. ∎

**Theorem 2L-final(3/8).** The same at α = 3/8, under the constraint a > 2c₂, c₂ := 3d₀/(32κb₂), which implies (SQ).
The time is L_N[3/8, c(3/8)] with σ₁(3/8) = 3/(8ad₀) + c₂ (the parent B₂-smoothness cost is leading at the endpoint;
uprime §4.1).

*Proof.* As above, with:
- Theorem AB₂ in place of AB;
- c12rat's Corollary RAT (RAT₂ form);
- Prop PB's second-order budget: A/η ≤ D^{3−9α}λ^C → 0 and √(AB)/η ≤ D^{−α}λ^C → 0 on [3/8, 2/5).

The endpoint-specific inputs are:
- uprime's U′ at 3/8 (proved there in 2L′; its hypotheses log z/log B₂ ≥ D^{1/4} hold at 3/8 since u₂ ≍ D^{3/8}),
  with the §2.7 DW replacements. **Re-check at DW's η₀ = D^{−3/8}**,
  leading order only (`params`: "TV N^{−η₀/2} << every L[α]^C"): under DW, level 1 uses DW(a)⁺ instead of GM-U′,
  so no good-modulus/level-gap condition is needed. The η₀-dependent quantities are (1) the TV terms 2N^{−η₀/2} of
  DW(a)⁺/TP-W against ρ̂_{B₂}(X_F) = L[3/8, −c₂] and ρ̂_B(X_r) = L[3/8]^{−O(1)}: η₀ log N = D^{5/8}λ ≫ D^{3/8}λ;
  (2) Lemma W-DW (C_W = 2, no η₀); (3) the classes (B), (C) of L2-U′, whose ratios are N^{−Θ(h)}B^{O(1)} with
  h ≥ 2μ + 2ν + 2η₀; (4) Kern-U's bb′ ≤ V₁² ≤ (2H + 1)/d. All hold. FM-U′'s subtraction (log N)²/B against ρ̂_{B₂}
  needs (SQ), implied by a > 2c₂. This is a leading-order re-check, not a line-by-line re-proof of uprime §3 at the
  new η₀;
- (SQ) for FM-U′'s subtraction and for Lemma η's denominator (c12alg §5, 2L-DW₂);
- Lemma η (needs α > 4/11);
- the recursion on parent values L_{X_F}[2/5] = L^{o(1)} (needs α > 4/11);
- Lemma reject-DW;
- a > 2c₂ for §11. ∎

**c(α).** C_γ is not quantified here: it contains the unspecified implied constants of Hildebrand–Tenenbaum's
saddle-point estimates and of Drappeau's Lemme 5, as for Theorem U. Illustrations with C_γ = 0 (`src/c12alg/params.cost_constant`):
- c ≈ 5.39 at α = 7/18, (a, κ, d₀, b₂) = (3/5, 1/8, 1, 2);
- at α = 3/8 a small grid under a > 2c₂ gives c ≈ 4.71 at (4/5, 1/8, 1, 2).

The window formula is **degenerate**. In the open window the level-1 algebraic cost is L^{o(1)}, and
σ₁ + σ₂ = (α/a)(5/d₀ + d₀b₂) can be made small (d₀ large, b₂d₀² small). So for fixed C_γ the infimum of c(α) over
parameters is 0. We do **not** claim this: C_γ depends on the parameters and may diverge as they degenerate, and every
N₀ moves. The honest statement is: **for each admissible parameter choice, a finite c(α)**; numerically
c(α) ≤ 5.4 + O(C_γ) at the standard point.

### 13.1 Complete dependency list

**This note (proved here):** Lemmas ML, PG, BL, CL, CF, LF-A, OPV-A, EXC-A, PU-A, EN-A, FAM-A; Prop PD-A and its
corollary; Theorem AB₂; Theorem AB; the margins (`src/thmAB/params.py`).

**c12alg (internal, proved there):**
- Lemmas MM, BZ, TP, TP-W, Cor TP-q (the pair law);
- DW(a)⁺, GB-DW, CU, RE, CT, IR-DW, reject-DW, W-DW;
- Lemma VD, Lemma η, Prop PB, Prop CLB;
- §2.7 (U′-DW(α)), §2.8 (POSTPROCESS-DW);
- Lemma WD's root-conditioned SQF bound (used in §11);
- c12alg's Theorem 2L-DW(α)/2L-DW₂ proofs (the assembly).

c12alg's Lemmas FB, PP, LF, AR are **not** used (superseded). Lemma PU (alg) is replaced by Lemma PU-A.

**c12rat:**
- Theorem RAT and Corollary RAT(α);
- Lemmas used directly here: CP, FO, BTS, PUB (rat and parent-absolute), SPR, MULT, TF, LS, DOM, Theorem WA,
  Lemma EXC(b), CELL(i)–(ii) (the convexity arguments, re-proved here for parallelograms).

**pc3:** Lemma 0, Prop DC, Prop DC_e, Cor DC, Lemma LAT, Lemma DW-L.

**uprime:** Lemma 8.5H, GM/FM/DS/L2/LD*/D-U′, CS′, SQF, the child scan, recursion and IR at γ₂ (with c12alg's
DW replacements; Lemma 8.5H is not needed at level 1 under DW, but Theorem U's FACTOR_U still uses Lemma 8.5′ for the
recursive value factorisations).

**multidescent:** Prop UB, Theorem 2L's proof (search, rows, LV6, fruitfulness, recursion).

**nested:** Lemmas RW, T_multi. **escapeR2:** §7.3 (bounded-depth search recurrence).

**Theorem U's chain** (lv5unif §10 list): [pp] = paper2/main.tex §§3, 5–9 (the older Harper-based FM proofs in
pcest §2, hypD and thmA_proof §7 are superseded by paper2 + lem85 and not used), Lemma 8.5′(R′) (lem85),
Lemma S″, Lemmas OS, W, F, T_rec (recurse), Theorem DD (lemmaM), and lv5unif's Lemmas P, U1, CS, Kern-U, GM-U, W-U,
DS-U, L2-U, Theorems FM-U, LD*-U, D-U, K-U, Prop select-U, Prop LV6-U, Lemma F_α.

**Published inputs.** Theorem numbers as cited in the source notes (c12alg §10, c12rat §12.1,
uprime §10, lv5unif §14); new ones here are flagged.
- **Drappeau**, Compositio Math. 151 (2015) 828–862, doi 10.1112/S0010437X14007933, **Lemme 5** (Fact D; A = 1/2,
  effective for A < 1; checked verbatim in arXiv 1307.7554v3; also restated as Pascadi, Compos. Math. 161 (2025),
  Lemma 3.2). It is stated
  without proof as a consequence of §3.3 of **Harper's unrefereed preprint arXiv:1208.5992**. **This is the one
  input whose justification rests on an unrefereed source.** Also **Lemme 1(i)**.
- **de la Bretèche–Tenenbaum**, Ramanujan J. 9 (2005) 139–202, Thms 2.4, 2.1 and (4.1), via **Drappeau Lemme 1(i)**
  (Fact dBT1; upper bound only; primary not re-opened).
- **Hildebrand**, J. Number Theory 22 (1986) 289–307, **Theorems 1 and 3**: Ψ = xρ(u)(1 + O(log(u + 1)/log y)), and
  the short-interval form uniformly for xy^{−5/12} ≤ z ≤ x, y ≥ exp((log log x)^{5/3+ε}). These are restated in
  S. Jain, arXiv:2502.10530 (QJM 2026), Lemma 3.2, and in the HT survey JTNB 5 (1993) Thm 5.1. The primary was not
  opened; the statement is checked against both restatements.
- **Hildebrand–Tenenbaum**, Trans. AMS 296 (1986) 265–290, **Thms 1, 2(ii), 3** (checked via the survey's Cor 1.7 and
  Drappeau's quotations; the label "2(i)" for the α asymptotic is unconfirmed: survey (1.26)). The β-drift uses the
  defining equation of α plus Thm 2(ii).
- **Hildebrand–Tenenbaum survey**, JTNB 5 (1993) 411–484, **Cor 2.4** and **Thm 5.1** (Lemma η).
- **Alladi** (log-concavity of ρ; Lemma η), attribution per the Encyclopedia of Mathematics ("Dickman function");
  likely J. reine angew. Math. 335 (1982) 180–196, page not pinned.
- **Canfield–Erdős–Pomerance** 1983.
- **Tenenbaum**'s Ψ_q theorem (Phil. Trans. R. Soc. A 345 (1993)).
- **Iwaniec–Kowalski Thm 7.13**, **Montgomery–Vaughan** (large sieve, Brun–Titchmarsh), **Rosser–Schoenfeld**, the
  PNT with de la Vallée Poussin error, **Trudgian** (Ramanujan J. 39 (2016) Thm 1).
- **Lagarias–Odlyzko** 1977 Thm 1.3, **Stark** 1974 Lemma 8, **Landau** (LV6).
- **Schwartz–Zippel**; **AKS**; **Wiedemann/Kaltofen–Saunders**; **Cantor–Zassenhaus**; **Bernstein** (batch smooth parts).
- Elementary: Brunn–Minkowski in dimension 2, Cavalieri, Gauss reduction, inclusion–exclusion, Parseval on (Z/𝔨)^×.

**Not used:** Harper's Theorems 1–2 (BV/BDH) themselves; Saias; de la Bretèche–Drappeau; Lee–Venkatesan's Lemma 8.5.

---

## 14. AB₂ and L[3/8]: status, and what would still fail

- **AB₂ is proved here** (§12) on [3/8, 2/5), discharging O1, O2, O3, O5, with O4 discharged by the constraint.
  Therefore L[3/8] follows (Theorem 2L-final(3/8)), with the same status as the window theorem plus the endpoint
  checks listed there.
- **The exact estimates that carry the endpoint.** Every exponent margin used is strict at α = 3/8
  (`params.margins(3/8)`: all "ok"), except the parent acceptability subtraction. That one is exactly borderline
  (3 − 7α = α) and is closed by a > 2c₂. Without that constraint, the failing estimate would be the mixed bound
  Pr(bad₁, smooth₂ | g, R) ≪ (log N)^C B^{−1}s₂, which the crude subtraction does not give (c12alg O4).
- **Below 3/8.** Nothing here changes multidescent §3.1: the two-level design is infeasible below 3/8 (the child
  bound window 2 − 4α ≤ ζ₁ ≤ 4α − 1 closes). Three levels are borderline at 4/11 (pc3 §7.5).



---

# B. Theorem RAT (source: notes/c12rat_results.md §§1-8)

## 1. Setting

### 1.1 Designs

Two designs are covered; both use the DW field law of pc3 §6.2 as repaired by c12alg §2 (w uniform on all integers of
(W_m, 2W_m], c uniform on C_H, F = wF̂_N + G·R_c, g uniform on primitive points of [M/2, M]²).
- **(W)** c12alg's window design 2L_DW(α), α ∈ (3/8, 2/5) (c12alg §1.1): d ≍ D^{1−2α}, k ≍ D^α, B = e^{aD^αλ},
  B₂ = e^{b₂D^{4α−1}λ}, η₀ log N ≍ D^{3α−1/2}λ, special-q ℓ ∈ J = (Y, 2Y], Y = L[1 − α], level-1 box radius V₁ with
  V₁² ≍ YS₁², child boxes of radius S₂ (S₂² = L[α, σ₂]) in Gauss-reduced coordinates, child primes q ∈ (B, B₂].
- **(P)** pc3/uprime's 2L′ at α = 3/8 with DW: B₂ = L[1/2, b₂], η₀ = D^{−1/2}, otherwise as uprime §1. **Proviso**
: here N^{η₀/2} = L[1/2, 1/2] and q ≤ L[1/2, b₂], so the TV terms qN^{−η₀/2} of DW(a)⁺ are small only
  for b₂ < 1/2, and Lemma DOM(ii)'s sum over q and r needs **b₂ < 1/4**; moreover c12alg's DW lemmas were written for
  (W)'s larger η₀. So (P) is covered only conditionally on b₂ < 1/4 and on re-checking DW(a)⁺, Cor TP-q and Lemma RE at
  η₀ = D^{−1/2}. Design (W) at α = 3/8 (which c12alg's Theorem 2L-DW₂ uses) needs no proviso.

The same proof applies verbatim to (W) at α = 3/8 (whose other problems are on the algebraic side, c12alg §3.5).

Throughout, y := B. For a label 𝔮 = (q, r), Box(𝔮) is the conservative child box (uprime §4.2): primitive
P′ = ie₁ + je₂, max(|i|, |j|) ≤ S₂, gcd(i, j) = 1, P′ primitive in Z², modulo ±1, for a Gauss-reduced basis
|e₁| ≤ |e₂| of Λ_𝔮 = {a ≡ rb (mod q)}. The skewness is s(𝔮) := |e₂|/|e₁| ≥ 1. The parent set of 𝔮 is 𝒫 ∩ Λ_𝔮
(𝒫: b prime in ((1 − 1/d)V₁, V₁], |a| ≤ V₁, gcd(a, b) = 1).

**Environment cells.** As in pc3 Prop DC_e, e(A) is the cell of g in a fixed partition of [M/2, M]² into closed squares
Ω_e of side δ_eM, δ_e := D^{−A₀}. All statements are conditional on R_𝔮 ∩ {g ∈ Ω_e}. By c12alg's Lemma RE the law of g
given R_𝔮 and e is uniform on 𝒢_e := {g ∈ Ω_e primitive, q ∤ G_g(r, 1)} up to a density factor 1 + O(q/H + q/W_m)
(relative error N^{−Ω(η₀)}; we absorb it in "tiny" below), and every g with q | G_g(r, 1) has N_𝔮 = Z_𝔮 = 0.

### 1.2 What is to be proved (the interface)

**Lemma VD.** Conditional on R_𝔮 ∩ e,
Var(N_𝔮) = Var_g(T_𝔮(g)) + E_g Var(N_𝔮 | g), with T_𝔮(g) := Σ_{P′∈Box(𝔮)} f_{P′}(G_g(P′)) and
f_{P′}(n) := 1[|n| is y-smooth]·1[|n| ≥ X⁻_{P′}]·a_{P′}(|n|), X⁻_{P′} := M|P′|_∞N^{−η₀/4} (the DW child filter, c12alg
§1.1), a_{P′}(t) := Pr_{w,c}[F(P′)/q ∈ 𝒮_B | g, R_𝔮] for g with |G_g(P′)| = t (by DW(a)⁺ a function of t only, up to a
relative error N^{−η₀/3}, which we fold into the weight). Likewise for the parent count Z_𝔮 with
f_P(n) := 1[|n| ∈ J·𝒮_B]·1[MV₁/32 ≤ |n| ≤ 2MV₁]·a_P(|n|), a_P(t) := Pr[F(P) acceptable, q child | g, R_𝔮].

We call ϕ^rat_{𝔮,e} := Var_g(T_𝔮)/m_𝔮² (m_𝔮 := E_g T_𝔮) the **rational block of C₁**, and σ^rat the same for Z_𝔮 (C₂).
The algebraic blocks E_g Var(·| g) are c12alg's. Expanding the variance,
  ϕ^rat_𝔮 = Σ_{P′} Var_g f_{P′}/m² + Σ_{P′≠P″} Cov_g(f_{P′}(G′), f_{P″}(G″))/m²  ≤  1/m_𝔮 + Φ_𝔮,
  Φ_𝔮 := Σ_{P′≠P″}|Cov_g(f′(G′), f″(G″))|/m_𝔮²
(the diagonal is ≤ Σ E f² ≤ m since 0 ≤ f ≤ 1). The 1/m term is the one-point term, ≤ e^{−λ²}(log N)^C for every label
by c12alg's Prop CLB. **The task is Φ_𝔮 (and its parent analogue Φ^par_𝔮), weighted by w_{𝔮,e}.**

The precision needed: by c12alg's Prop PB, in (W) A ≤ D^{−α}λ^C and B ≤ D^{7α−3}λ^C suffice; at α = 3/8 one needs
A + √(AB) ≤ η/4 with η ≍ D^{−3/8} (pc3 Cor DC), so A, B ≤ D^{−3/8−ε} suffice. We prove D^{−2α}λ^C for both.

### 1.3 Exponents and the loss exponent θ (`params.margins`, 14 checks, exact rationals)

Leading order, in both designs: log y = aD^αλ; γ = 1 − β(x, y) ≤ C₁/(aD^α) for every relevant x (HT (α2): (1 − β) log y
≤ log(u log y) + C_α log log y = O(λ)); child and parent rational values have log x ≍ D^{2α}λ and u ≍ D^α; the BDH saving
of Fact D is ε := e^{−c_εD^α/λ²}. The moduli are K ≤ 4Q₂S₂² (children; log K ≲ D^{4α−1}λ in (W), D^{1/2}λ in (P)) and
K ≤ 2V₁² (parents; log K ≍ D^{1−α}λ). Define
  θ(α) := max( γ-exponent + child modulus exponent, γ-exponent + (1 − α), γ-exponent + η₀-exponent )
       = max(3α − 1, 1 − 2α, 2α − 1/2) in (W),   θ = 1/4 in (P);
it is the exponent of every "loss" Ξ := e^{O(D^θλ)} below ((|P′||P″|)^γ, (V₁²)^γ, (X/X⁻)^γ, and the spread
(|P′|^{d−1})^{γ_a} of the child algebraic weights over a box, of exponent (1 − 2α) + o(1); Lemma SPR). On [3/8, 2/5)
θ = 2α − 1/2. **On [3/8, 2/5), θ < α**,
so Ξ^{O(1)} = ε^{−o(1)}, and also Ξ^{O(1)} ≤ y^{1/12}. Further, with P₀ := ⌈D^{7/10}⌉:
- every modulus is a small conductor: log K ≤ (η_D/4) log x;
- ω(K) = o(P₀) for every child and parent modulus, so Σ_{p | K, p > P₀} p^{−β} → 0;
- P₀# = e^{(1+o(1))P₀} ≤ x^{o(1)}, far below every cell length.
`params.margins` checks these (and the window/3/8 comparisons of §1.2) for α = 3/8 (both designs) and on a grid of
(3/8, 2/5); at α = 2/5 only "B₂ ≪ V₁" becomes borderline.

### 1.4 The algebraic weights (regularity used)

We use only the weights a_P of DW(a)⁺ (children, with the root coset) and the parent weights a_P of §1.2 — not an
abstract class. Their regularity:

**Lemma ADM [proved].** Write γ_P(t) := 1 − β(t·L_P, 𝔶) for the saddle point at the target scale t·L_P (L_P := the
fixed scale factor of DW(a)⁺: L_P = (2H + 1)|P|_∞^{d−1}/q for children with 𝔶 = B, L_P = (2H + 1)|P|_∞^{d−1} with
𝔶 = B₂ for parents). Then:
- (local, two-sided) for t ≤ t′ ≤ y²t in the range of |G_g(P)|: a_P(t′)/a_P(t) = (t/t′)^{γ_P(t)}(1 + O(E_a)), with
  E_a ≤ Cλ³D^{−α} (children) and E_a ≤ Cλ³(D^{7α−3} + D^{1−4α}) (parents), i.e. the first-order errors 1/u + E_H of
  the respective target sets;
- (global, one-sided) for all t ≤ t′ in the range: a_P(t)/a_P(t′) ≤ C(t′/t)^{γ_a}, γ_a := max_t γ_P(t)
  ≤ C/(aD^α) (children), ≤ C/(b₂D^{4α−1}) (parents).
*Proof.* By DW(a)⁺/Cor TP-q (c12alg §§2.6, 3.2), given g with |G_g(P)| = t, F(P) (resp. F(P)/q) is within TV
ε_TV ≤ 7qN^{−η₀/2} of the law of t·L_P(Z + U), Z with a law independent of t and U uniform on [0, 1); the w-shift is
≤ N^{−η₀/4}·scale and is inside ε_TV. So a_P(t) = E_Z[Ψ_𝒯(tL_P(Z + [0, 1)))/(tL_P)] + O(ε_TV) for the target set 𝒯
(𝒮_B, resp. the acceptable B₂-smooth set). For each Z the interval has relative length ≥ 1/(2d) ≥ y^{−5/12}; by
Facts H1, H2 its 𝒯-density is ρ̂_𝒯(tL_P|Z|₊)(1 + O(E_H)). Local: Fact HT3 iterated ≤ 3 times gives
ρ̂_𝒯(t′x)/ρ̂_𝒯(tx) = (t/t′)^{1−β(tx)}(1 + O(1/u)) for t′/t ≤ y², and the drift of β between tx and tL_P is negligible
(|log|Z|| = O(log d) except on a ϑ-layer of mass ≤ dN^{−η₀/8}, absorbed into E_a). The curvature of the saddle-point
expansion over a ratio ≤ y² is O(log²(y²)/(log x log y)) = O(D^{−α}), inside E_a. Global: Fact dBT1.
For the acceptable set, the squarefree/discriminant subtraction is a relative O((log N)²/B) (uprime Lemma SQF);
ε_TV/a_P ≤ N^{−η₀/3}. ∎

Ratios t′/t > y² occur only on sets of small one-point mass (Lemma LC), where only the one-sided bound is used (OPV,
PUB). An arithmetic factor (e.g. a design with gcd(b′, g₁) = 1, making the weight depend on gcd(G′, b′)) could be
added if it is constant on the P₀-types of §4 and within 1 ± 1/P₀ otherwise; the DW weights need none.

---

## 2. Main statements

For a label 𝔮, a cell e and child points P′ ≠ P″ write S′ := E_{𝒢_e} f′(G′), Pair := E_{𝒢_e} f′(G′)f″(G″),
κ := |P′||P″|/K (K = |det(P′, P″)|; κ ≥ 1, large for nearly parallel pairs).

**Theorem RAT.** In design (W) for α ∈ [3/8, 2/5) (and in (P) under the proviso of §1.1), with the DW weights, there are constants C, c > 0 and
parameters T₀ = κ₀^{1/3} = e^{C_TD^θλ²} (§7; C_T arbitrary large, fixed) such that:
- (RAT-1) Σ_{𝔮,e: s(𝔮)≤T₀} w_{𝔮,e}Φ_{𝔮,e} ≤ C D^{−2α}λ^{C}·W + e^{−cD^α/λ²}μ₁;
- (RAT-2) Σ_{𝔮,e} w_{𝔮,e}Φ^par_{𝔮,e} ≤ C (D^{−2α} + D^{−2(3−7α)})λ^{C}·W + e^{−cD^α/λ²}μ₁;
- (RAT-0) Σ_{𝔮,e: s(𝔮)>T₀} w_{𝔮,e} ≤ T₀^{−1/2}μ₁.

Here W = Σ w_{𝔮,e} is the expected number of parent–child incidences, μ₁ = E N₁, and W ≥ (1 − o(1))μ₁ (Lemma DOM(iii)).

**Corollary RAT(α).** For every fixed α ∈ (3/8, 2/5), Hyp RAT(α) of c12alg §5 holds: (R1) Σ_{s(𝔮)≤T}(w/W)ϕ^rat
≤ D^{−2α}λ^C + Σ(w/W)/m and (R2) Σ(w/W)σ^rat ≤ (D^{−2α} + D^{14α−6})λ^C + Σ(w/W)/z, for every T ≤ T₀ (in particular
c12alg's T = Ξ_alg³(log N)^CD^α, once C_T is large); note D^{14α−6} ≪ D^{7α−3} for α < 3/7. The same bounds give Hyp
RAT₂ in the form of c12alg's Theorem 2L-DW₂ on [3/8, 2/5). [The one-point terms are c12alg's: Prop CLB for 1/m, and
the identity w_𝔮/z_𝔮 = Pr(R_𝔮)/#J for 1/z.]
*Proof.* Divide Theorem RAT by W ≥ μ₁/2. ∎

The proof occupies §§3–8. Its skeleton, for one non-thin pair (κ ≤ κ₀) and one cell, is **Prop PD** (§8.1):
  |Pair − S′S″| ≤ √(V′V″) + Σ_c |Err(c)| + O(h^{1/4}S′S″) + tiny,
where V′, V″ are one-point variances of cell densities (≤ CD^{−2α}λ⁶S′², pointwise in the pair: Lemma OPV) and Err(c)
are character errors (Theorem WA); primitivity at primes > P₀ adds O((γ²/P₀)λ² + P₀^{−2})S′S″ (Lemma PRIM). Thin pairs (κ > κ₀) use Lemma PUB, skewed labels Prop SK.

---

## 3. Geometry [all proved; elementary]

**Lemma LC (log-concave marginals).** Let Ω be a convex body in R² containing a square of side δM/2 and contained in a
square of side 2δM, with δM ≥ 16|P| (cells, their dilates, and the enlarged regions Π⁺ of Prop PD pulled back to g-space
qualify), P = (a, b) primitive, and
let 𝒢 ⊂ Ω ∩ Z² contain at least a fraction 1/2 of Ω ∩ Z² (true for the P₀-primitive points of Ω). Put
n(g) := det(P, g) = G_g(P), x_max := max_Ω |n| and H := length of n(Ω) ≥ δM|P|/2. Then for every interval I,
  Pr_{g∈𝒢}[n(g) ∈ I] ≤ 8(|I| + 1)/H,
and for 0 < t ≤ 1/2, Pr_{g∈𝒢}[|n(g)| ≤ t·x_max] ≤ 40t + 16/H. Consequently, for k ≤ 8 and 0 ≤ σ ≤ 1/2,
E_{g∈𝒢}[log^k(x_max/|n|)(x_max/|n|)^σ·1(|n| ≥ X⁻)] ≤ C_k + (16/H)(log(x_max/X⁻))^k(x_max/X⁻)^σ.
*Proof.* The g ∈ Z² with n(g) = n form a coset of ZP; its intersection with Ω has ≤ c(n)/|P| + 1 points, where c(n) is
the length of the chord {det(P, ·) = n} ∩ Ω. Lines with consecutive n are at distance 1/|P|, so |Ω| = |P|^{−1}∫c. The
function c is concave on n(Ω) (Brunn–Minkowski in dimension 2), so ∫c ≥ c_max·H/2. Hence
#{g ∈ Ω∩Z² : n(g) ∈ I} ≤ (|I| + 1)(c_max/|P| + 1) ≤ 2(|I| + 1)|Ω|(1 + |P|/c_max)/H; with c_max ≥ δM/2 ≥ 8|P|,
#(Ω∩Z²) ≥ |Ω| − O(δM) ≥ |Ω|/2 and #𝒢 ≥ #(Ω∩Z²)/2 this gives Pr ≤ 8(|I| + 1)/H; H ≥ δM|P|/2. For the second: if 0 ∈ n(Ω), H ≥ x_max; if not, the event is empty unless the
nearest endpoint is ≤ tx_max, in which case H ≥ (1 − t)x_max. The moment bound follows by integrating the tail. ∎

**Lemma LS (label shapes).** For a prime q and T ≥ 1, #{r mod q : s(q, r) > T} ≤ (π/2)(√(2q/(√3T)) + 1)²
≤ 1.82q/T + 3.4(q/T)^{1/2} + 1.6, hence ≤ 3q/T + 2 when q ≥ 16T.
*Proof.* s > T and |e₁||e₂| ≤ (2/√3)q (reduced basis) give |e₁|² < 2q/(√3T). The shortest vector ±e₁ is primitive in Z²
(a proper divisor would be a shorter vector of Λ_𝔮, as |e₁| < q), has b ≢ 0 (mod q), and determines r ≡ ab^{−1}.
Count lattice points in the disc. ∎ (E1 P1a checks the equivalent σ-form exactly for all r at q = 10007, 100003.)

**Lemma MULT (multiplicities).**
- (child) For a label and S := S₂ ≥ 8, the box has n ≥ S² points, every two distinct box points have coordinate
  determinant Δ = |ij″ − ji″| ∈ [1, 2S²] and K = qΔ, and for every Δ ≥ 1 the number of ordered pairs with this Δ is
  ≤ 10S² + 2S. Hence for every k ≥ 1: Pr_{pairs}[k | Δ] ≤ 50/k for k ≤ 2S² (and 0 otherwise).
- (parent) For P ∈ 𝒫 and K ≠ 0, #{P′ ∈ 𝒫 ∩ Λ_𝔮 : det(P, P′) = K} ≤ 5; and for every prime p ≠ q with p ≤ y and every
  P = (a₀, b₀) ∈ 𝒫, Pr_{P′∈𝒫∩Λ_𝔮}[p | det(P, P′)] ≤ (1/p)(1 + O(pq/V₁)).
*Proof.* (child) For a fixed primitive coordinate vector (i, j) with r := max(|i|, |j|), the (i″, j″) with
ij″ − ji″ = Δ form a coset of Z(i, j), so at most 2S/r + 1 lie in [−S, S]² (one representative per ± class). There are
≤ 4r coordinate classes with max-norm r, and Σ_{r≤S}4r(2S/r + 1) ≤ 10S² + 2S. Distinct ± classes of primitive vectors
are non-proportional, so Δ ≥ 1; |Δ| ≤ 2S² is trivial; K = qΔ is LAT(d). For k | Δ sum over the ≤ 2S²/k multiples and
divide by n(n − 1) ≥ S⁴/2. (n ≥ S² for S ≥ 8: the number of primitive points of [−S, S]² is (24/π²)S² + O(S log S),
and at most a fraction O(1/q) of them give Z²-imprimitive P′.) (parent) On the line {det(P, ·) = K} consecutive integer
points differ by P, |P| ≥ b > 2V₁/3, and the box has diameter 2√2V₁. For the second claim fix the prime b′ ∈ window.
Since p ≤ y < b₀, p ∤ b₀, and p | a₀b′ − a′b₀ fixes a′ mod p; with a′ ≡ rb′ (mod q) it fixes a′ mod pq, so at most
2V₁/(pq) + 1 values |a′| ≤ V₁ remain, against ≥ 2V₁/q − 2 values of a′ ≡ rb′ (mod q) with gcd(a′, b′) = 1. ∎

**Lemma TF (thin pairs).** κ(P′, P″) ≤ 2√3·s(𝔮)S²/Δ for child pairs, so the fraction of ordered pairs of Box(𝔮) with
κ > κ₀ is ≤ 2√3s(𝔮)S²(10S² + 2S)/(κ₀ n(n − 1)) ≤ 60 s(𝔮)/κ₀. For parents, the fraction of pairs P, P′ ∈ 𝒫∩Λ_𝔮 with
κ > κ₀ is ≤ 3/κ₀ + O(q/V₁).
*Proof.* For a reduced basis |ie₁ + je₂|² ≤ (3/2)(i²|e₁|² + j²|e₂|²) (the cross term is ≤ |ij||e₁|²), so
|P′|² ≤ 3S²|e₁||e₂|·s and |P′||P″| ≤ 3S²s|e₁||e₂| ≤ 2√3 qS²s; divide by K = qΔ and use MULT. Parents: |P||P′| ≤ 2V₁², so
κ > κ₀ forces P′ into the strip {|det(P, ·)| < 2V₁²/κ₀} of width ≤ 3V₁/κ₀ about RP; it holds ≤ 9V₁²/(qκ₀) + O(V₁)
points of Λ_𝔮 in the box, against ≍ 4V₁²/q. ∎ (`tests/test_c12rat.py::test_box_det_mult_thin_bounds`; E1 P1c.)

**Lemma CELL (tilings).** Fix a non-thin pair (κ ≤ κ₀) and a convex region Ω as in Lemma LC (a cell or a dilate), and
let Π := A(Ω), A(g) := (det(P′, g), det(P″, g)), a convex set of area K|Ω|. Let ℓ_r (ℓ_c) be the longest chord of Π
parallel to the n₁-axis (n₂-axis), H₁, H₂ the widths of Π in the two coordinates, x′_max := max_Π|n₁|, and put
z₁ := hℓ_r, z₂ := hℓ_c with h := D^{−A₁}. Let 𝒳 := {±(1 + h_*)^k : k ≥ 0} with h_* := hδ_e/(16κ₀) (the **universal
grid**, the same for all labels, cells, dilates and pairs). A **block** of coordinate i is a union of consecutive grid
intervals, with relative width ≤ min(h^{1/2}, 1/8) and length ≤ 2z_i; the blocks of a pair partition {X⁻_i ≤ |n_i|
≤ 2x_max} (the filter thresholds are made block endpoints, see Lemma EXC(d)). Then:
- (i) z₁ ≥ hδ_e x′_max/(8κ) ≥ 2h_*x′_max and z_i ≫ K·P₀#·B; every block has length ≥ min(z_i, |n|/16)·h_*-scale ≫
  K·P₀#·B (near the filter, |n| ≥ X⁻ ≫ K L^{C}); and blocks of length < z_i only occur at |n_i| ≤ 16z_i/min(h^{1/2},1/8).
- (ii) Every rectangle R = I₁ × I₂ contains (|I₁||I₂|/K)(1 + O(K/min|I_i|)) points of L := A(Z²) = {n₁ ≡ ιn₂ (mod K)}
  (L ⊇ KZ², and each K×K square contains exactly K points of L). The rectangles meeting Π but not contained in Π
  ("boundary rectangles") cover an area ≤ 20h|Π|, and all lie in the convex set Π⁺ := Π + [−2z₁, 2z₁] × [−2z₂, 2z₂],
  whose preimage A^{−1}(Π⁺) ⊂ Ω + B(0, 4h·diam Ω) is a convex region as in Lemma LC with |A^{−1}Π⁺| ≤ (1 + 20h)|Ω|.
- (iii) The n with |n₁| ≤ 128z₁h^{−1/2} have μ-mass ≤ Ch^{1/2} (μ the law of A(g) for g ∈ 𝒢(Ω)), by Lemma LC (H₁ ≥ ℓ_r).
*Proof.* (i) A chord of Π parallel to the n₁-axis is the image of a chord of Ω in direction P″, along which n₁ moves by
K per unit of the parameter t in g = g₀ + tP″; for Ω = Ω_e/e the central chord has length ≥ δ_eM/(2e), so
ℓ_r ≥ Kδ_eM/(2e|P″|), while x′_max ≤ 2√2M|P′|/e. So z₁/x′_max ≥ hδ_e/(8κ), **independently of the dilation e**. The
grid spacing at |n| ≤ 2x_max is ≤ 4h_*x_max ≤ z₁/2, and the relative spacing h_* is ≪ min(h^{1/2}, 1/8). And
z_i ≥ hδ_eM/(2y²κ₀) ≫ K·L^{C} since M = N^μ and e ≤ y². (ii) Convexity gives
|Π| ≥ ℓ_rH₂/2, |Π| ≥ ℓ_cH₁/2 and ℓ_rℓ_c ≤ ℓ_rH₂. Every rectangle has sides ≤ 2z_i, so the boundary rectangles lie within
box-distance (2z₁, 2z₂) of ∂Π: area ≤ 4z₁H₂ + 4z₂H₁ + 16z₁z₂ ≤ (16h + 32h²)|Π|. A^{−1} maps a displacement of z₁ in n₁
(at fixed n₂) to a displacement h·(chord) ≤ 2h·diam Ω along P″, and similarly for n₂. (iii) Lemma LC with
|I| ≤ 256z₁h^{−1/2} and H₁ ≥ ℓ_r = z₁/h. ∎
(Since z_i/x_max is dilation-invariant, one geometric grid serves all cells and all dilates Ω/e of Lemma PRIM; only
the range of scales grows, by log(y²), in the count of grid points.)

---

## 4. Types, coupled pairs and exact factorisation [proved; finite algebra]

Fix a non-thin pair P′, P″ with K = |det(P′, P″)| and unit ι (P′ ≡ ιP″ mod K, pc3 LAT), so (G_g(P′), G_g(P″)) = A(g)
runs over L = {n₁ ≡ ιn₂ (mod K)} as g runs over Z², bijectively. Write k_p := v_p(K). Let P₀ := ⌈D^{7/10}⌉ and split
the primes dividing K into **small** (p ≤ P₀), **mid** (P₀ < p ≤ y) and **large** (p > y; q is large).

**Types.** For n ∈ L with n₁n₂ ≠ 0 and gcd(n₁n₂, K_large) = 1, define:
- e_p(n_i) := v_p(n_i) for p ≤ P₀ (all small primes, dividing K or not);
- j_p(n) := min(v_p(n₁), k_p) = min(v_p(n₂), k_p) for mid p | K (equal by the lattice relation, ι being a unit);
- the type τ(n) := ((e_p(n₁), e_p(n₂))_{p≤P₀}, (j_p(n))_{p mid});
- s_i := ∏_{p≤P₀}p^{e_p(n_i)}·∏_{p mid}p^{j_p(n)}, m_i := n_i/s_i, δ := gcd(n₁, K) = gcd(n₂, K) (determined by τ);
- M_τ := (K/δ)·∏_{p≤P₀, p|δ} p.

n is **P₀-primitive** if the content of g = A^{−1}n has no prime factor ≤ P₀. The type τ is **admissible** if
(a) min(e_p(n₁), e_p(n₂)) = 0 for every p ≤ P₀ with p ∤ K, and (b) not both e_p(n_i) ≥ k_p + 1 for p ≤ P₀ with p | K.

**Lemma TYPE.** Let τ be admissible and put ι_τ := ι·(s₂/δ)(s₁/δ)^{−1} mod K/δ and, for p ≤ P₀ with p | δ and
e_p(n₁) = e_p(n₂) =: j ≤ k_p, ι_{τ,p} := ι·(s₂/p^j)(s₁/p^j)^{−1} mod p^{k_p−j+1}. Let
  C_τ := {(r₁, r₂) ∈ ((Z/M_τ)^×)² : r₁ ≡ ι_τ r₂ (mod K/δ), and r₁ ≢ ι_{τ,p}r₂ (mod p^{k_p−j+1}) for each such p}.
Then:
- (i) n ↦ (m₁, m₂) is a bijection from {n ∈ L : type τ, gcd(n₁n₂, K_large) = 1, P₀-primitive} onto
  {(m₁, m₂) ∈ (Z∖0)² : m_i coprime to P₀#, (m_i, M_τ) = 1, (m̄₁, m̄₂) ∈ C_τ};
- (ii) C_τ is invariant under the diagonal action (r₁, r₂) ↦ (ur₁, ur₂) of (Z/M_τ)^×; hence every r₁ has the same
  number |C_τ|/φ(M_τ) of partners.

*Proof.* (i) Given the type, n_i = s_im_i with m_i P₀-free (exact valuations at p ≤ P₀), coprime to the mid primes p
with j_p < k_p and to the large primes, and unrestricted at mid primes with j_p = k_p (which do not divide M_τ). The
lattice relation n₁ ≡ ιn₂ (mod p^{k_p}) at p | K: if v_p(δ) = k_p it is automatic; otherwise both valuations equal
j := v_p(δ) < k_p and it is m₁ ≡ ι(s₂/p^j)(s₁/p^j)^{−1}m₂ (mod p^{k_p−j}); together, r₁ ≡ ι_τr₂ (mod K/δ).
P₀-primitivity at p ≤ P₀: p | g iff n ∈ pL (A is a bijection Z² → L). If p ∤ K, n ∈ pL iff p | n₁ and p | n₂, excluded
by (a). If p | K: n ∈ pL iff n/p ∈ L iff (at p) p^{e₁−1}s′₁m₁ ≡ ιp^{e₂−1}s′₂m₂ (mod p^{k_p}). If e₁, e₂ ≥ k_p + 1 this
holds (excluded by (b)); if exactly one of them is ≥ k_p + 1 and the other equals k_p it fails; if e₁ = e₂ = j ≤ k_p
it is the congruence mod p^{k_p−j+1} excluded in C_τ. At primes ≠ p the condition n/p ∈ L is equivalent to n ∈ L. The
inverse map is (m₁, m₂) ↦ (s₁m₁, s₂m₂). (ii) All conditions are of the form r₁ ≡ c·r₂ or r₁ ≢ c·r₂ to some modulus
dividing M_τ, preserved by common scaling. ∎

(The sign matters: the relation is between *signed* values. `test_type_coupling_sets_TYPE` enumerates g exhaustively for
four pairs and checks that fibre sizes depend only on the type and that the realised sets are diagonal-invariant.)

**Lemma CP (coupled pairs).** Let M ≥ 1, C ⊂ ((Z/M)^×)² be diagonal-invariant, and a, b: Z → C finitely supported on
integers prime to M, with class sums a(r) := Σ_{m≡r}a(m), A := Σa, and â(χ) := Σ_m a(m)χ(m). Then
  |Σ_{m₁,m₂} a(m₁)b(m₂)1_C(m̄₁, m̄₂) − (|C|/φ(M)²)AB| ≤ (|C|/φ(M)²)·(Σ_{χ≠χ₀}|â(χ)|²)^{1/2}(Σ_{χ≠χ₀}|b̂(χ)|²)^{1/2}.
*Proof.* Write a(r) = A/φ + α(r), b(r) = B/φ + β(r), Σα = Σβ = 0. Since every r₁ has |C|/φ partners and every r₂ has
|C|/φ partners, Σ_C α(r₁) = Σ_C β(r₂) = 0. So the left side is |Σ_C α(r₁)β(r₂)| ≤ (Σ_Cα(r₁)²)^{1/2}(Σ_Cβ(r₂)²)^{1/2}
= (|C|/φ)‖α‖₂‖β‖₂, and by Parseval on (Z/M)^×, Σ_{χ≠χ₀}|â(χ)|² = φ‖α‖₂². ∎
(`test_coupled_pairs_lemma_CP`, `test_parseval_cyclic`.) This is the exact form of pc3's scaling-invariance
cancellation (Lemma SI): no cross terms, and the error is a product of two one-coordinate L² discrepancies.

**Corollary CPW (slowly varying weights).** If moreover a(m) = ā·u(m)·(1 + ϵ(m)), b(m) = b̄·v(m)(1 + ϵ′(m)) with
u, v ≥ 0, |ϵ|, |ϵ′| ≤ ε′ ≤ 1/4, then with Ẽ_u := Σ_{χ≠χ₀}|û(χ)|²/(Σu)²,
  |Σ a b 1_C − (|C|/φ²)AB| ≤ (|C|/φ²)AB·[3(Ẽ_u + Ẽ_v) + 4ε′²].
*Proof.* Class sums: a(r) = ā u(r) + e(r) with |e(r)| ≤ ε′āu(r) (positivity of u). So ‖α‖ ≤ ā(1 + ε′)‖u − Σu/φ‖
+ ε′āΣu/√φ, and similarly for β; A ≥ (1 − ε′)āΣu. Insert in Lemma CP and use xy ≤ (x² + y²)/2. ∎
*The point:* the weights' oscillation enters **squared**, and only multiplied by the main term.

**Lemma FIB (non-smooth counts).** For intervals J₁, J₂ of lengths ≥ M_τ P₀#·B and units r_i mod M_τ:
#{m ∈ J_i : m ≡ r_i (M_τ), m P₀-free} = (|J_i|/M_τ)∏_{p≤P₀, p∤M_τ}(1 − 1/p) + O(2^{π(P₀)}). Hence the number of
n of type τ (P₀-primitive) in a rectangle R with J_i = I_i/s_i is N(R, τ) = (|C_τ|/φ(M_τ)²)U₁U₂(1 + O(t)), where
U_i := #{m ∈ J_i : P₀-free, (m, M_τ) = 1} and t := 2^{π(P₀)}M_τ/|J_min| = N^{−Ω(μ)} ("tiny"); and for each fixed m₁ the
fibre has (|C_τ|/φ²)U₂(1 + O(t)) points, independent of m₁.
*Proof.* Inclusion–exclusion over e | P₀#/(P₀#, M_τ); Lemma TYPE(ii) for the fibres. 2^{π(P₀)} = e^{O(D^{7/10})}. ∎

---

## 5. Analytic inputs

### 5.1 Cited results (published; statements as used)

- **Fact D** (Drappeau, Compositio 151 (2015), Lemme 5; restated as Pascadi's Lemma 3.2 — paper2/uprime cite it as Compositio
  161 (2025), Jain cites it as a 2025 preprint; publication details UNVERIFIED by us; = paper2 Fact `fact:D`).
  There are η_D, c_D, δ_D > 0 such that for (log x)^{c_D} ≤ y ≤ x and Q ≤ x,
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ mod k, 1<cond χ≤x^{η_D}}|Ψ(x, y; χ)| ≪ Ψ(x, y){H(u)^{−δ_D} + y^{−δ_D}} (with A = 1/2 and the
  factor (log x)^{−A} discarded; the implicit constant is effective for A < 1, and we take A = 1/2), H(u) = e^{u/log²(u+1)}.
  **Drappeau states Lemme 5 without proof**, as a consequence of §3.3 of Harper's preprint arXiv:1208.5992 ([Har12b],
  pré-publication), which is unrefereed; so Fact D's justification rests on that preprint (paper2 discusses this). We use
  **only** this statement: no BDH theorem of Harper is used.
- **Fact dBT1** (de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), Thm 2.4(i), restated as Drappeau Lemme 1(i)). For 2 ≤ y ≤ x and
  1 ≤ d ≤ x, Ψ(x/d, y) ≤ C Ψ(x, y)d^{−β(x,y)}; equivalently ρ̂(x/d) ≤ C d^{γ(x)}ρ̂(x), ρ̂(x) := Ψ(x, y)/x.
- **Fact HT3** (Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Thm 3). Ψ(cx, y) = c^{β(x,y)}Ψ(x, y)(1 + O(1/u + log y/y))
  for x ≥ y ≥ 2, 1 ≤ c ≤ y. Iterating at most k times (with HT Thm 2's control of the drift of β, as in lv5unif Lemma
  U1): Ψ(x/d, y) = d^{−β(x,y)}Ψ(x, y)(1 + O_k(1/u)) for 1 ≤ d ≤ y^k, x/d ≥ x^{1/2}.
- **Fact HT(α2)** (a consequence of HT Thm 2, as paper2 (α2)): (1 − β(z, y)) log y ≤ log(u_z log y) + C_α log log y; in our ranges
  γ = 1 − β ≤ C₁/(aD^α), and y^{γ} = e^{O(λ)}.
- **Fact H1, H2**. For y ≥ exp((log log x)^{5/3+ε}):
  Ψ(x, y) = xρ(u)(1 + O_ε(E_H)) and, **uniformly** for xy^{−5/12} ≤ z ≤ x, Ψ(x + z, y) − Ψ(x, y) = zρ(u)(1 + O_ε(E_H)),
  the implied constant depending only on ε,
  E_H := log(u + 1)/log y. **Consequence used:** (Ψ(x + z) − Ψ(x))/z = ρ̂(x)(1 + O(E_H)) on that range ("the
  short-interval density equals the global density to first order"). In our ranges E_H ≍ λ/(aD^αλ) = O(D^{−α}).
- Iwaniec–Kowalski Thm 7.13 is **not** needed (all moduli are small conductors, §1.3).

### 5.2 Lemma EN (energies of the library sequences) [proved]

For real ξ and integer s with ξ/s ≥ X^{1/2} (X the relevant value scale), let 𝒜_{ξ,s} := {1 ≤ m ≤ ξ/s : m P₀-free,
P⁺(m) ≤ y} and A_{ξ,s}(χ) := Σ_{m∈𝒜_{ξ,s}}χ(m). For parents let 𝒜^J_{ξ,s} := {ℓm : ℓ ∈ J, m ∈ 𝒜_{ξ/ℓ,s}} with the
analogous A^J. Then for Q ≤ (ξ/s)^{η_D/2},
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀ (k)}|A_{ξ,s}(χ)|² ≤ ε_lib Ψ(ξ/s, y)²,  ε_lib := C(log P₀)³(H(u/3)^{−δ_D} + y^{−δ_D}),
and the same for A^J with Ψ(ξ/s, y) replaced by Σ_{ℓ∈J}Ψ(ξ/(sℓ), y). In our ranges ε_lib ≤ e^{−c_εD^α/λ²}.
*Proof.* Möbius: A(χ) = Σ_{e|P₀#}μ(e)χ(e)Ψ(ξ/(se), y; χ), so |A(χ)| ≤ Σ_e|Ψ(ξ/(se); χ)| and |A(χ)| ≤ A₀ :=
Σ_eΨ(ξ/(se)) ≤ CΨ(ξ/s)∏_{p≤P₀}(1 + p^{−β}) ≤ CΨ(ξ/s)(log P₀)^{3/2} (Fact dBT1; Σ_{p≤P₀}p^{−β} ≤ log log P₀ + C because
P₀^{γ} = O(1)). Hence |A|² ≤ A₀Σ_e|Ψ(ξ/(se); χ)|, and Fact D at x = ξ/(se) ≥ (ξ/s)^{1/2} (P₀# = e^{O(P₀)} =
(ξ/s)^{o(1)}) bounds Σ_kφ(k)^{−1}Σ_{χ≠χ₀}|Ψ(ξ/(se); χ)| by CΨ(ξ/(se))(H(u/3)^{−δ_D} + y^{−δ_D}): every non-principal
character mod k ≤ Q has conductor ≤ Q ≤ x^{η_D}, and (log x)^{c_D} ≤ y. Sum over e. For A^J add the ℓ-sum to the
triangle inequality. ∎
*Remark.* This is the L¹ → L² trick |Ψ(χ)|² ≤ Ψ·|Ψ(χ)| (pc3 §4.1). It only ever uses Fact D's 1/φ(k)-weighted form, which
is why the weighted application (§6) needs the moduli to be **spread** (multiplicity domination).

### 5.3 Lemma FO (first-order cell densities) [proved]

Let x lie in the relevant range (log x ≥ (1/2)log X), y^{−5/12}x ≤ z ≤ x, J := (x, x + z], 1 ≤ s ≤ S₀ := e^{D^{(θ+α)/2}}
(so S₀ ≤ y^{o(1)}, γ log S₀ → 0, and S₀^{1/log P₀} ≥ Ξ^{ω(1)}),
and 𝔪 squarefree with P⁺(𝔪) ≤ y, 𝔪 | P₀#·𝔪₁ and Σ_{p|𝔪₁}p^{−β} ≤ 1 (𝔪₁: mid primes of a modulus). Then
  #{m ∈ J/s : P⁺(m) ≤ y, (m, 𝔪) = 1} = (z/s)·s^{γ}ρ̂(x)·g_𝔪(β)·(1 + O(E_FO)),  E_FO := C(log P₀)³(E_H + 1/u),
with β = β(x, y), γ = 1 − β, g_𝔪(β) := ∏_{p|𝔪}(1 − p^{−β}).
*Proof.* Inclusion–exclusion: the count is Σ_{e|𝔪}μ(e)Ψ(J/(se)). For se ≤ y²: Facts H1, H2 at x/(se) (the relative
length z/x ≥ y^{−5/12} is unchanged) give Ψ(J/(se)) = (z/x)Ψ(x/(se))(1 + O(E_H)), and Fact HT3 (k ≤ 3 steps) gives
Ψ(x/(se)) = (se)^{−β}Ψ(x)(1 + O(1/u)). For se > y²: Lemma BTS below (modulus 1) and Fact dBT1 give
Ψ(J/(se)) ≤ (z/x)(se)^{−β}Ψ(x)·e^{O(λ)}(x/z)^{γ} ≤ (z/x)(se)^{−β}Ψ(x)·e^{O(λ)} (as (x/z)^γ ≤ y^{5γ/12}), and by
Rankin. Write e = e_s·e_m with
e_s | P₀# and e_m | 𝔪₁. Since s ≤ S₀ = y^{o(1)}, se > y² forces e_s > y^{0.9} or e_m > y. With σ := 1/log P₀ (p^σ ≤ e for
p ≤ P₀): Σ_{e_s>y^{0.9}}e_s^{−β} ≤ y^{−0.9σ}∏_{p≤P₀}(1 + e p^{−β}) ≤ e^{−0.9 log y/log P₀}(log P₀)^{C}, and
Σ_{e_m}e_m^{−β} ≤ exp(Σ_{p|𝔪₁}p^{−β}) ≤ e; with exponent 1/2 instead: Σ_{e_m>y}e_m^{−β} ≤ y^{−1/2}∏_{p|𝔪₁}(1 + p^{−1/2+γ})
≤ y^{−1/2}e^{ω(K)P₀^{−0.49}} ≤ y^{−1/2}e^{D^{0.3}}, and Σ_{e_s}e_s^{−β} ≤ (log P₀)^{C}. Both tails are smaller than every power of D (log y/log P₀ ≍ D^α and aD^αλ/2 ≫ D^{0.3}). So the count is
(z/x)s^{−β}Ψ(x)[g_𝔪(β) + O((E_H + 1/u)∏_{p|𝔪}(1 + p^{−β}))], and ∏(1 + p^{−β})/∏(1 − p^{−β}) ≤ exp(3Σ_{p|𝔪}p^{−β})
≤ C(log P₀)³. Finally (z/x)s^{−β}Ψ(x) = (z/s)s^{γ}ρ̂(x). ∎

### 5.4 Lemma BTS (smooth numbers in a progression in a short interval) [proved; elementary]

Let k ≥ 1, a ∈ Z, I = (X₀, X₀ + Z] with W := Z/k ≥ 2y, and suppose gcd(a, k)·W ≤ X₀. Then
  #{n ∈ I : n ≡ a (mod k), P⁺(n) ≤ y} ≤ 4W log(2y)·max_{W/y≤w≤2W}ρ̂(w),
and if log W ≥ (1/2)log(X₀ + Z), this is ≤ e^{O(λ)}·W·((X₀ + Z)/W)^{γ}ρ̂(X₀ + Z) (Facts dBT1, HT(α2)).
*Proof.* Put δ := gcd(a, k); n = δn′ with n′ ≡ a/δ (mod k/δ) a unit, n′ ∈ I/δ (if δ is not y-smooth there is nothing
to count). Each y-smooth n′ > W has a divisor m₁ ∈ (W/y, W] (take the largest product of its smallest prime factors
that is ≤ W). For fixed m₁ (smooth, coprime to k/δ), m₂ := n′/m₁ lies in an interval of length Z/(δm₁) in a fixed
class mod k/δ: ≤ W/m₁ + 1 values. Sum over m₁: W Σ_{W/y<m₁≤W}1_S(m₁)/m₁ + Ψ(W), and dyadic blocks give the bound.
(The exact inequality, before the dyadic step, is `bts.rhs`; E2 checks it on 2000 instances.) ∎
This is the large-modulus case of [pp] Lemma DS; c12alg's Lemma PU uses the same device.

### 5.5 Lemma PUB (universal pair upper bound) [proved]

For every label, region Ω (cell or dilate, Lemma LC), and every pair P′ ≠ P″ (thin or not):
  Pair ≤ Ξ·S′S″,  Ξ := e^{O(λ)}(|P′||P″|/δ_e)^{γ}·(a″_max/a″_min)·y^{γ} ≤ e^{C_ΞD^θλ},
with a″_max/a″_min the ratio of the weight over the range of |G″| on Ω (≤ C(X″/X″⁻)^{γ_a} by Lemma ADM, global part).
*Proof.* Σ_{g∈𝒢}f′(G′)f″(G″) ≤ Σ_{n₁}f′(n₁)·#{g ∈ Ω∩Z² : G′(g) = n₁, f″(G″(g)) > 0}. The g with G′(g) = n₁ form
g* + ZP′, along which G″ = G″(g*) ± jK: an arithmetic progression of difference K with W(n₁) := (chord of Ω along
the line)/|P′| + 1 terms, values ≤ X″ := max_Ω|G″|. Lines with W(n₁) < X″^{1/2} carry ≤ X″^{1/2}|P′|/(δ_eM) ≤ M^{−1/3}
of the g-mass, which is below S′S″·N^{−μ/4} (S′, S″ ≥ L^{−O(1)}). For the others Lemma BTS gives
≤ e^{O(λ)}a″_max W(n₁)^{1−γ}X″^{γ}ρ̂(X″). By Hölder, Σ_{n₁}f′W^{1−γ} ≤ (Σf′W)^{1−γ}(Σf′)^{γ}; here Σ_{n₁}f′(n₁)W(n₁) =
Σ_{g∈Ω∩Z²}f′(G′(g)) ≤ e^{O(λ)}#𝒢·S′ (the non-P₀-primitive g are a fraction ≤ 1/2 and carry density ≤ e^{O(λ)} times
the typical one, Fact dBT1), and Σ_{n₁}f′(n₁) ≤ e^{O(λ)}a′_max H′(X′/H′)^{γ}ρ̂(X′) (Lemma BTS with k = 1 on the range
of G′, of length H′), while #𝒢S′ ≥ e^{−O(λ)}a′_min|Ω|ρ̂(X′) and |Ω| ≥ H′·W_avg/2 (Lemma LC). So
Σf′W^{1−γ} ≤ e^{O(λ)}#𝒢S′·W_avg^{−γ}·(a′_max/a′_min)^{γ}·…, and with S″ ≥ e^{−O(λ)}a″_minρ̂(X″) (Lemma FO on the
interior blocks, which carry mass ≥ 1/2 by Lemma CELL) we get Pair ≤ e^{O(λ)}(a″_max/a″_min)(X″/W_avg)^{γ}y^{γ}S′S″.
Finally X″/W_avg ≤ 4|P′||P″|/δ_e. The exponent: γ log(|P′||P″|/δ_e) = O(D^{−α}·(log B₂ + log V₂ + log V₁²))λ = O(D^θλ),
and log(a″_max/a″_min) ≤ γ_a·(η₀/4) log N + O(1) = O(D^θλ) (§1.3). ∎
**Parent version (absolute).** For parent pairs, S″ can vanish (inactive points, Lemma SPR), so the relative form
is replaced by Pair ≤ Ξ_p·S′·S_typ, Ξ_p := e^{O(D^{1−2α}λ)}: the same line argument, with the density of f″ along a line
bounded by Ξ_p·S_typ (BTS on the line; the values lie in [MV₁/32, 2MV₁] and the weight at scale [X_F/(64d), X_F]).
c12alg's Lemma PU (alg and rat parts) is the same bound with the same exponent; we need only its rational part.

### 5.6 Lemma OPV (one-point variance of cell densities) [proved]

Fix a non-thin pair, a region Ω (cell or dilate), its blocks (Lemma CELL), types relative to K, and the first
coordinate. For a block I and a type t let ψ′(I, t) := Σ_{m∈U(I,t)}f′(s_t m)/#U(I, t), U(I, t) := {m ∈ I/s_t :
P₀-free, (m, M_t) = 1}. Let ν₁ be the law of (block(G′), t(G′)) for g ∈ 𝒢(Ω) (P₀-primitive points), restricted to
non-exceptional n (Lemma EXC). Let x₁ := max_Ω|G′|, β := β(x₁, y), γ := 1 − β, and
  c := a′(x₁)ρ̂(x₁)·∏_{p≤P₀}(1 − p^{−β})/(1 − p^{−1})·∏_{p|K, mid}(1 − p^{−β})/(1 − p^{−1})·∏_{p|K, large}(1 − p^{−1})^{−1}.
Then
  E_{ν₁}[(ψ′ − c)²] ≤ 𝒱₀c²,   E_{ν₁}[(ψ′/c)⁸] ≤ C,   and   c/2 ≤ S′ ≤ 2c,
  𝒱₀ := C[(γ + γ_a)²(log P₀)² + E_FO² + E_a² + γ_a²h + h^{1/4}] + C·γ²(log P₀)²ω(K)/P₀ ≤ C D^{−2α}λ⁶,
the last pointwise in K: ω(K)/P₀ ≤ D^{max(4α−1, 1/2, 1−α)−7/10+o(1)} = o(1) for child and parent moduli (§1.3), so the
mid-prime term is ≤ D^{−2α}λ²·o(1).
*Proof.*
- *Ordinary blocks* (relative width ≤ h^{1/2}; Lemma CELL(iii) puts all others — width ≤ 1/8, at |n₁| ≤ 128z₁h^{−1/2} —
  into a set of mass ≤ Ch^{1/2}) and types with s_t ≤ S₀: Lemma FO with 𝔪 = P₀#·rad(K/δ)_mid, Lemma FIB for #U, and
  Lemma ADM (local) across the block give ψ′(I, t) = c·X(ξ)·Y(t)·Z(t)·(1 + e), |e| ≤ C(E_FO + E_a + γ_ah^{1/2}), where
  X(ξ) := a′(ξ)ρ̂(ξ)/(a′(x₁)ρ̂(x₁)) (ξ the block's end nearer 0), Y(t) := s_t^{γ}, and
  Z(t) := ∏_{p mid, j_p = k_p}(1 − 1/p)/(1 − p^{−β}): for a mid prime with j_p = k_p the cofactor is no longer coprime
  to p (p ∤ M_t), so the factor (1 − p^{−β})/(1 − p^{−1}) present in c is cancelled. Thus Z ≥ 1 and
  Z − 1 ≤ 2Σ_{p mid, j_p=k_p}(p^{−β} − p^{−1}) ≤ 2Σ_{p mid, p^{k_p} | n₁}(p^γ − 1)/p. (Replacing β(ξ, y) by β(x₁, y) inside s^γ costs O(log s·log(x₁/ξ)/(log x·log y)), negligible.)
- *Archimedean factor.* For ξ ≥ x₁y^{−2}: X = (x₁/ξ)^{γ+γ_{P′}(ξ)}(1 + O(1/u + E_a)) (Fact HT3, Lemma ADM local); for
  all ξ: X ≤ C(x₁/ξ)^{γ+γ_a} (Fact dBT1, Lemma ADM global). Lemma LC gives E[(X − 1)⁴] ≤ C(γ + γ_a)⁴(1 + …) +
  C(1/u + E_a)⁴ and E[X¹⁶] ≤ C (the part ξ < x₁y^{−2} has mass ≤ 40y^{−2} and its moments are ≤
  Σ_{k≥2log y}40e^{−k(1−16(γ+γ_a))} ≤ y^{−1}).
- *Local factor.* Under ν₁, Pr[s | n₁] ≤ 2/s for every s ≤ S₀·y composed of small primes and mid prime powers: G′ is
  equidistributed mod s on Ω up to 1 + O(s|P′|/(δ_eM)), and conditioning on P₀-primitivity only lowers these
  probabilities (Pr[p^k | n₁] = p^{−k}·p/(p + 1) for p ≤ P₀). Hence, splitting s_t = s_small·s_mid:
  E[(s_small^γ − 1)⁴] ≤ γ⁴E[log⁴s_small·s_small^{4γ}] ≤ Cγ⁴(log P₀)⁴ (a finite Euler product, P₀^{4γ} = O(1));
  E[(s_mid^γ − 1)⁴] ≤ C Σ_{p|K mid}p^{−1}(p^γ − 1)⁴ + (lower order) ≤ Cγ⁴(log P₀)⁴ω(K)/P₀, because for P₀ < p ≤ e^{1/γ},
  (p^γ − 1)⁴/p ≤ γ⁴log⁴p/p ≤ γ⁴(log P₀)⁴/P₀ (decreasing in p), and for p > e^{1/γ}, (p^γ − 1)⁴/p ≤ p^{4γ−1} ≤ e^{4−1/γ};
  the second moments are bounded likewise by Cγ²(log P₀)²(1 + ω(K)/P₀). Similarly E[(Z − 1)²] ≤ CΣ_{p|K mid}(p^γ−1)²/p³
  ≤ Cγ²ω(K)/P₀³, and E[Y¹⁶], E[Z¹⁶] ≤ C. Types with s_t > S₀ are exceptional (Lemma EXC).
- *Assembly.* ψ′/c − 1 = (X − 1)YZ(1 + e) + (Y − 1)Z(1 + e) + (Z − 1)(1 + e) + e on ordinary blocks; Cauchy–Schwarz
  with the fourth moments gives E[(ψ′/c − 1)²·1_ord] ≤ 𝒱₀/2. On the remaining mass ≤ Ch^{1/2},
  E[(ψ′/c − 1)²·1_rest] ≤ (Ch^{1/2})^{1/2}E[(ψ′/c + 1)⁴]^{1/2} ≤ Ch^{1/4} (≤ D^{−2α} once A₁ ≥ 8α + 1). The eighth
  moment follows from the
  pointwise bound ψ′/c ≤ C·X·Y·Z (Fact dBT1/BTS for every block). Finally S′ = E_{ν₁}ψ′·(1 − O(h + e^{−λ})) + B′ with
  0 ≤ B′ ≤ C(h^{1/2} + e^{−λ})c (boundary and exceptional parts, Prop PD), so |S′ − c| ≤ c(𝒱₀^{1/2} + Ch^{1/4}) ≤ c/2. ∎

### 5.7 Lemma SPR (spread of one-point densities in a box) [proved]


- (children) For a label with s(𝔮) ≤ T₀, a cell e and every P′ ∈ Box(𝔮), let c_{P′} be OPV's c. Then
  max_{P′}c_{P′}/min_{P′}c_{P′} ≤ Ξ, and ψ′ ≤ Ξc_{P′} on every cell of every pair; hence S′S″/m_𝔮² ≤ 4Ξ³/n² and
  ψ′ψ″/m_𝔮² ≤ Ξ⁴/n².
- (parents) For every label and cell, call P ∈ 𝒫 ∩ Λ_𝔮 *active* if |G_g(P)| ≥ MV₁/16 on all of Ω_e. Then at least half
  of the points are active, every active P has S_P ≥ S_typ/polylog, and every P has ψ_P ≤ S_typ·polylog on every cell,
  where S_typ := ϱ̂_J(MV₁)·ρ₂(X_F)-scale. Hence m^par_𝔮 ≥ n·S_typ/polylog and ψψ/(m^par)² ≤ polylog/n²; **inactive P can
  have S_P = 0**, so no lower bound for individual S_P is claimed or used.
*Proof.* (children) x₁(P′) ∈ [cδ_eM|P′|, 2√2M|P′|] (the cell is a square of side δ_eM, so max|G′| over it is
≥ δ_eM|P′|/4), and |P′| varies over the box by a factor ≤ 4S₂s(𝔮) ≤ 4S₂T₀ (Lemma TF's proof). So ρ̂(x₁) varies by
≤ C(16S₂T₀/δ_e)^{γ} (Fact dBT1) = e^{o(D^θλ)}, and a′(x₁) — whose argument is x₁L_{P′} with L_{P′} ∝ |P′|^{d−1} — by
≤ C((16S₂T₀/δ_e)^{d})^{γ_a} = e^{O(dγ_a(log S₂ + log T₀))} = e^{O(D^{1−2α}λ + D^{1−3α+θ}λ²)} ≤ Ξ. The local factors in c differ by ≤ e^{O(ω(K)/P₀)} = 1 + o(1). The pointwise bound
ψ′ ≤ C·X·Y·Z·c with X ≤ C(x₁/X⁻)^{γ+γ_a} ≤ Ξ^{1/2}, Y ≤ S₀^γ = 1 + o(1), Z ≤ 2. m ≥ n·min S′ ≥ n min c/2.
(parents) For g ∈ Ω_e, |G_g(P)| = |P||g||sin(θ_P − θ_g)|, |P| ≥ 2V₁/3, |g| ≥ M/√2, so P is active when
|sin(θ_P − θ_e)| ≥ 1/4 + O(δ_e); the angles θ_P of 𝒫 ∩ Λ_𝔮 are equidistributed on [π/4, 3π/4] (for each prime b in the
window, a ≡ rb (mod q) runs over an AP of difference q ≪ V₁), and θ_e ∈ [0.46, 1.11], so at least half the points are
active. For active P, the values lie in [MV₁/16, 2MV₁] and F(P) at scale [X_F/(64d), X_F] (c12alg Lemma η), so
ρ̂_J and ρ₂ vary by (64d)^{γ+γ₂} = 1 + o(1); S_P ≥ (pass fraction 1)·S_typ/polylog by Lemma FO. Upper bounds hold for all P by BTS/dBT1. ∎

---

## 6. The weighted application theorem (pc3 §7 item 1)

pc3 §5 asked for "one theorem stating the exact measure, its multiplicity and domination bounds".
Given Lemma CP, the needed statement is linear in the error energies, and it is elementary.

**Theorem WA.** Let 𝔏 be a finite *library* of
finitely supported sequences u_j: Z → [0, 1], independent of the moduli, with normalisers U_j > 0 such that
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀ (k)}|û_j(χ)|² ≤ ε U_j²   for every j ∈ 𝔏.
Let (c_t, j(t), k(t))_{t∈𝒯} be a finite family of terms, c_t ≥ 0, j(t) ∈ 𝔏, k(t) ≤ Q, and H_j(k) := Σ_{t: j(t)=j,
k(t)=k}c_t. Then
  Σ_t c_t·Σ_{χ≠χ₀ (k(t))}|û_{j(t)}(χ)|²/U_{j(t)}²  ≤  ε·𝔇,   𝔇 := Σ_{j∈𝔏}max_{k≤Q}φ(k)H_j(k).
*Proof.* Σ_jΣ_k[φ(k)H_j(k)]·φ(k)^{−1}Σ_χ|û_j(χ)|²/U_j² ≤ Σ_j max_k[φH_j]·ε. ∎

The content is in bounding 𝔇 for the true family (Lemma FAM), which needs only an **upper** bound for the incidence
weights (Lemma DOM) and the lattice multiplicities (Lemma MULT); errors are measured against the dominating mass
Σŵ ≤ (log N)^Cμ₁, not against W. Small label means are harmless because every
error is *relative* to the cell's own main term (Cor CPW), so a label's normalisation m_𝔮² cancels up to the spread
ratio Ξ of one-point densities inside a box (Lemma PUB's proof: S_max/S_min ≤ Ξ for s(𝔮) ≤ T₀).

**Lemma DOM (weight domination) [proved].** Let α̂(q) := C(log N)^C q^{−1}ρ₂(X_F/q) + 7qN^{−η₀/2} (ρ₂ := ρ̂ at B₂) and
ϱ̂ := C(log N)^Cρ̂_J(2MV₁), ρ̂_J(x) := #(J·𝒮_B ∩ [1, x])/x. Then:
- (i) for every label (q, r) and cell e: w_{(q,r),e} ≤ ŵ_{q,e} := Pr_g[Ω_e]·α̂(q)ϱ̂·#𝒫(1 + O(q/V₁))/(q·#J), which does
  not depend on r;
- (ii) Σ_{q,r,e}ŵ_{q,e} ≤ (log N)^{C}μ₁;
- (iii) W ≥ (1 − L^{−c})μ₁.
*Proof.* (i) w_{𝔮,e} = Pr(R_𝔮, e)z_{𝔮,e}/#J = #J^{−1}Σ_{ℓ∈J}Σ_{P∈𝒫∩Λ_𝔮}Pr[g ∈ Ω_e, Acc_ℓ(P), q ∥ F(P)] (Lemma 0 of
pc3: the incidence implies R_𝔮). For fixed P and g with the rational filter, c12alg's DW(a)⁺ with the root coset and
an upper bound for B₂-smooth numbers in intervals of relative length ≥ 1/(2dq) (BTS with k = 1) give
Pr_{w,c}[q | F(P), F(P)/q B₂-smooth | g] ≤ α̂(q), pointwise in g (this is c12alg's Lemma η, first bullet). Summing over
ℓ, Σ_ℓPr_g[g ∈ Ω_e, Rat_ℓ(P)] = Pr_g[g ∈ Ω_e, G(P) ∈ J·𝒮_B, filter] ≤ Pr[Ω_e]·ϱ̂ (Lemma LC and BTS with k = 1 on the
range of G(P) over Ω_e). Finally #(𝒫 ∩ Λ_𝔮) ≤ (#𝒫/q)(1 + O(q/V₁)) uniformly in r (for each prime b in the window,
a ≡ rb mod q leaves 2V₁/q + O(1) values), and V₁ ≫ q·L^{C} on [3/8, 2/5) (§1.3). (ii) Σ_r has q terms, and
Σ_{B<q≤B₂}α̂(q) ≤ (log N)^{C}ρ₂(X_F)·Σ_q q^{−β₂} ≤ (log N)^{C}u₂ρ₂(X_F) (HT3 and B₂^{1−β₂} ≍ u₂ log u₂, as in c12alg
Lemma η), while μ₁ ≥ Σ_{ℓ,P}Pr[Rat_ℓ(P)]ρ₂(X_F)(log N)^{−C} (FM-U′-DW, c12alg §2.7) with Σ_{ℓ,P}Pr[Rat_ℓ(P)] ≥
#𝒫ϱ̂/(#J(log N)^{C}) (FIB and U1); u₂ ≤ D. (iii) An accepted parent with no child has F(P) B-smooth; by the same
one-point upper bound with B in place of B₂ their expected number is ≤ μ₁(log N)^Cρ_B(X_F)/ρ₂(X_F) ≤ μ₁L^{−c}, since
u_B := log X_F/log B ≍ D^{2−4α} ≫ u₂ ≍ D^{3−7α} for α > 1/3. ∎
(c12alg's Lemma WD proves the two-sided version w_{(q,r)} ≍_{polylog} w_{(q,s)}; we need only (i).)

**Lemma FAM (the accounting) [proved].** For the C₁ family of one window (q ∈ (Q₂, 2Q₂]) — terms t = (label with
s(𝔮) ≤ T₀, cell e, dilation e′ ∈ {1} ∪ {e′ ≤ y² : P⁻(e′) > P₀} (Lemma PRIM), ordered non-thin pair, cell c of the
pair in Ω_e/e′, coordinate i, endpoint), library j = (endpoint ξ ∈ 𝒳, s ≤ S₀) with u_j = 1_{𝒜_{ξ,s}},
U_j = Ψ(ξ/s, y), modulus k(t) = M_c, and weights
  c_t := ŵ_{q,e}·ϖ_{e′}·m_{𝔮,e}^{−2}·(N(c)/Z)ψ′(c)ψ″(c)·C(ξ/|I_i|)²(log P₀)⁴,  ϖ₁ := 1, ϖ_{e′} := e′^{−2+2γ},
— one has 𝔇 ≤ Ξ⁴·κ₀³h^{−1}δ_e^{−1}D^{4}·Σ_{q∈win,r,e}ŵ_{q,e}. The same holds for the C₂ family (parent pairs, library
𝒜^J_{ξ,s}) with Ξ replaced by polylog·Ξ_p.
*Proof.* *Size of a term.* N(c)/Z = ν(c); by Lemma SPR, ψ′ψ″/m_𝔮² ≤ Ξ⁴/n² (children), ≤ polylog/n² (parents, using the
active-point lower bound for m^par — individual parent densities may vanish, so only the aggregate is used); and
(ξ/|I_i|)² ≤ h_*^{−2} for every block (Lemma CELL: relative width ≥ h_*). So c_t ≤ ŵ_{q,e}ϖ_{e′}Ξ⁴n^{−2}ν(c)h_*^{−2}(log P₀)⁴.
*H_j(k) for fixed j = (ξ, s) and k.* The window primes q dividing k: ≤ 1 + log(4S₂²)/log B = O(1) for children
(k ≤ 4Q₂S₂², q > Q₂ ≥ B) and ≤ log(2V₁²)/log B ≍ D^{1−2α} for parents. The pair's K is
one of the ≤ τ(s)2^{ω(s)} numbers kδ/r_δ with δ | s, r_δ | rad(δ) (M_τ = (K/δ)r_δ and δ = gcd(s, K)); each gives
Δ = K/q, with ≤ 10S² + 2S ordered pairs per label (Lemma MULT). For one pair and one dilation, the cells with endpoint ξ
in coordinate i and s_i = s have ν-mass ≤ 2Pr[s | n_i] ≤ 4/s (Lemma OPV's local law). Summing over labels (all r, e)
and dilations (Σϖ_{e′} ≤ 2): H_j(k) ≤ O(1)·max_qŴ_q·Ξ⁴n^{−2}(10S² + 2S)·τ(s)2^{ω(s)}(8/s)·h_*^{−2}(log P₀)⁴, with
Ŵ_q := Σ_{r,e}ŵ_{q,e}.
*The max over k.* φ(k) ≤ k ≤ 4Q₂S₂² (M_τ ≤ K). So max_kφ(k)H_j(k) ≤ CQ₂·max_qŴ_q·Ξ⁴h_*^{−2}(log P₀)⁴τ(s)2^{ω(s)}/s.
*The sum over j.* The grid points in the relevant range number ≤ log(2x_max/X⁻_min)/h_* ≤ CD^{2α}λ/h_*, and Σ_{s≤S₀}τ(s)2^{ω(s)}/s ≤ (log S₀)⁴ = D^{2(θ+α)}. So 𝔇 ≤ CD^{2α+2θ+2α}λQ₂max_qŴ_qΞ⁴h_*^{−3}(log P₀)⁴.
*Comparison.* Σ_{q∈win}Ŵ_q ≥ (Q₂/(2 log Q₂))min_qŴ_q, and max_q Ŵ_q/min_q Ŵ_q ≤ C on a dyadic window (ŵ_q is the
explicit function α̂(q)ϱ̂#𝒫/(q#J)·Pr[Ω_e] of Lemma DOM, and α̂(q) = q^{−1}ρ₂(X_F/q)(log N)^C varies by ≤ 2^{1+γ₂} there).
With h_*^{−1} = 16κ₀/(hδ_e): 𝔇 ≤ C D⁴λ log Q₂·Ξ⁴κ₀³h^{−3}δ_e^{−3}(log P₀)⁴·Σ_{q∈win}Ŵ_q, which is Ξ^{O(1)}κ₀^{O(1)}·Σ ŵ.
For C₂ the computation uses the parent multiplicity ≤ 5#(𝒫 ∩ Λ_𝔮) per K, K ≤ 2V₁², #(𝒫 ∩ Λ_𝔮)² pairs and
#(𝒫∩Λ_𝔮) ≍ #𝒫/q: the ratio Q·5#(𝒫∩Λ)/(#(𝒫∩Λ)²·#q_win) ≤ C d log V₁ log Q₂ is polylogarithmic. ∎

**Corollary WA′.** For each window, Σ over the family of ŵ_{q,e}m^{−2}·Σ_c|Err(c)| (the character part of Cor CPW) is
≤ ε_lib·Ξ^{C}κ₀^{C}h^{−3}δ_e^{−3}D⁵·Σ_{win}ŵ; summed over the ≤ log₂B₂ windows and with Lemma DOM(ii), ≤ e^{−cD^α/λ²}μ₁.
*Proof.* Cor CPW bounds |Err(c)| by N(c)ψ′ψ″·3(Ẽ′ + Ẽ″) with Ẽ for the interval sequence 1_{𝒜}·1_{I/s} mod M_c. An
interval sum is a difference of two endpoint sums, so |û_J|² ≤ 2|û_{ξ₁,s}|² + 2|û_{ξ₂,s}|², and the principal sum is
≥ (|I|/ξ₂)Ψ(ξ₂/s)(log P₀)^{−2}/C (Lemma FO with 𝔪 ⊃ rad(M_c)_y; Σ_{p|M_c mid}p^{−β} → 0). This gives the weights c_t
above; apply Theorem WA with ε = ε_lib (Lemma EN, Q = 4Q₂S₂² ≤ x^{η_D/2}) and Lemma FAM. Ξ^{C}κ₀^{C}h^{−3}δ_e^{−3}D⁵ = e^{O(D^θλ²)}
and ε_lib = e^{−c_εD^α/λ²}, θ < α. ∎

*What was missing in pc3 §5, and why it is no longer missing.* (a) The product S₁(χ)S̄₂(χ) of two different sequences
is split by Lemma CP/CPW into two L² energies, each attached to a sequence of a **common library** (indexed by grid
endpoints and the small-prime part s, not by the pair); (b) the weights (labels, geometry, cells, algebraic densities,
normalisations) enter only through the multiplicity function H_j(k), for which crude pointwise bounds suffice because
the BDH saving e^{−cD^α/λ²} beats every loss e^{O(D^θλ)}; (c) "maximal over scales" is replaced by a universal grid of
e^{O(D^θλ)} points, each a separate library element; (d) LSQ
is unnecessary because Q ≤ x^{η_D} (all conductors are small), so Fact D alone covers characters with q | cond χ.

---

## 7. Skewed labels and nearly parallel pairs (pc3 §7 item 2)

Put T₀ := e^{C_TD^θλ²} and κ₀ := T₀³ (C_T large and fixed; since θ < α, κ₀^{O(1)} = ε_lib^{−o(1)}, κ₀ ≤ y^{1/100}, and
T₀ ≥ Ξ^{10} as well as T₀ ≥ c12alg's threshold Ξ_alg³(log N)^CD^α, whose exponent is O(D^θλ²)).

**Prop SK [proved].**
- (i) Σ_{q,r,e: s(q,r)>T₀}w_{(q,r),e} ≤ (3/T₀ + 2/B)(log N)^{C}μ₁ (Lemmas LS, DOM; q > B ≥ 16T₀).
- (ii) [conditional on c12alg's Theorem AB(ii), a sketch; not used in Theorem RAT] The dead-incidence mass (Prop DC_e)
  of the skewed labels is ≤ W_sk + (W_sk·Σ_{all}wσ²)^{1/2}, W_sk the mass in (i);
  with Σ_{all}wσ² ≤ (log N)^{C}μ₁·(B_rat + B_alg + one-point) (C₂ for all labels: RAT-2 and c12alg's Theorem AB(ii)),
  this is ≤ T₀^{−1/3}μ₁.
- (iii) In a label with s(𝔮) ≤ T₀, the thin pairs (κ > κ₀) are a fraction ≤ 60T₀/κ₀ = 60T₀^{−2} of all pairs, and
  contribute ≤ 120T₀^{−2}·4Ξ⁴ ≤ T₀^{−1} to Φ_𝔮 (Lemma PUB and Lemma SPR).
- (iv) Parent sets 𝒫 ∩ Λ_𝔮 lie in a Euclidean square; their thin pairs are a fraction ≤ 3/κ₀ + O(q/V₁) for **every**
  label (skewed child lattices included), contributing ≤ polylog·Ξ_p/κ₀ to Φ^par_𝔮.
*Proof.* (i) By DOM(i) the weight of a residue r is ≤ ŵ_{q,e}, independent of r, and by LS at most 3q/T₀ + 2 residues
are skewed; sum over e and q and use DOM(ii) (Σ_{r}ŵ = qŵ). (ii) Prop DC_e: the dead mass of a set 𝔖 of (label, cell)
is ≤ Σ_𝔖w(p + σ√p) ≤ Σ_𝔖w + (Σ_𝔖w)^{1/2}(Σ_𝔖wσ²)^{1/2} by p ≤ 1 and Cauchy–Schwarz. (iii) Lemma TF with s(𝔮) ≤ T₀;
each thin pair has |Cov| ≤ Pair + S′S″ ≤ (Ξ + 1)S′S″ ≤ (Ξ + 1)·4Ξ³m²/n² (Lemma SPR). (iv) Lemma TF (parents, per
point P: the fraction of thin partners of each P is ≤ 3/κ₀ + O(q/V₁)) and the absolute PUB: Σ_{thin}Pair ≤
Ξ_pc_typ·Σ_PS_P·(3n/κ₀) ≤ polylog·Ξ_p(m^par)²/κ₀. ∎

So skewness is handled without any shape-dependent analytic estimate: a skewed child lattice costs nothing beyond its
(pointwise dominated) weight, because the companion statement C₂ is insensitive to the child lattice's shape, and
transversality is never used — nearly parallel pairs are either thin (upper bound PUB, loss Ξ) or handled by the cell
decomposition, whose only κ-dependence is the factor h_*^{−3} ∝ κ₀³ in Lemma FAM. c12alg's Theorem AB(iii) gives the same
conclusion for the whole ϕ via its Lemmas PU and WD; the two arguments are independent.

---

## 8. Assembly

### 8.1 Prop PD (one non-thin pair) [proved]

Let Ω be a square of side ≥ δ_eM/y² (cells and their dilates Ω/e, e ≤ y², for Lemma PRIM), 𝒢 its P₀-primitive points,
P′ ≠ P″ child (or parent) points with κ ≤ κ₀, f′, f″ the DW-weighted indicators of §1.2, S′ := E_𝒢f′(G′), Pair := E_𝒢f′(G′)f″(G″). Then (for parents, where S′ may
vanish, read every S′ on the right as the OPV constant c′, and centre at π′c′; §8.3)
  |Pair − S′S″| ≤ (V′V″)^{1/2} + 𝔈 + C(h^{1/4} + e^{−λ} + κ₀^{−1/2})S′S″,
where V′ := Var_ν(ψ′) ≤ 4(𝒱₀ + 20h + e^{−λ})S′² (Lemma OPV; 𝒱₀ includes the mid-prime term of K) and likewise V″, and
𝔈 := Z^{−1}Σ_c N(c)ψ′(c)ψ″(c)·3(Ẽ′(c) + Ẽ″(c)) (sum over interior and boundary rectangles c, Cor CPW).
*Proof.* Let Z := #𝒢 and use the cells c = (R, τ) of Lemma CELL × Lemma TYPE (admissible types; the n with
gcd(n₁n₂, K_large) > 1, with ℓ | n for a prime ℓ ∈ J dividing K (parents), or with s_i > S₀, form an exceptional set
ℰ, handled by Lemma EXC below).
- *Interior cells.* By Lemma TYPE(i) the n ∈ c ∩ 𝒢 are the (s₁m₁, s₂m₂) with (m̄₁, m̄₂) ∈ C_τ, m_i ∈ U(I_i, t_i). Lemma
  CP/CPW with a(m) := f′(s₁m)1_{U₁}(m), b likewise, and Lemma FIB give Σ_{n∈c}f′f″ = N(c)ψ′(c)ψ″(c)(1 + O(tiny))
  + Err(c), |Err(c)| ≤ N(c)ψ′ψ″[3(Ẽ′ + Ẽ″) + 4ε′²], with ε′ ≤ C(γ_ah^{1/2} + E_a) the weight oscillation over the
  block (≤ γ_a/8 + E_a ≤ 1/4 on blocks of relative width ≤ 1/8, which carry mass ≤ Ch^{1/2}); and Σ_{n∈c}f′(n₁) = N(c)ψ′(c)(1 + O(tiny)) (the fibre
  over each m₁ has the same size, Lemma FIB). Here ψ′(c) depends on c only through (I₁, t₁) (it is the ψ′ of Lemma OPV).
- *Boundary rectangles.* 0 ≤ Σ_{n∈c∩𝒢∩Π}f′f″ ≤ Σ_{n∈c∩𝒢}f′f″, and the full rectangle obeys the interior formula; the
  same for the singles.
- *Covariance identity.* Let ν be the probability measure on {interior cells} ∪ {∗} with ν(c) = N(c)/Z and the atom ∗
  carrying the rest (boundary and ℰ), and put ψ′(∗) = ψ″(∗) := 0. Write S′ = E_νψ′ + B′ + O(tiny) with
  0 ≤ B′ ≤ b′ := Z^{−1}Σ_{∂}N(c)ψ′(c) (boundary part), and similarly S″. Then
  Pair − S′S″ = Cov_ν(ψ′, ψ″) + Z^{−1}Σ_{int}Err(c) + P_∂ + P_ℰ − (E_νψ′)B″ − B′(E_νψ″) − B′B″ + tiny,
  where B′ now includes the exceptional singles (≤ e^{−λ}S′) and P_ℰ ≤ e^{−λ}S′S″ + κ₀^{−1/2}S′S″ (Lemma EXC),
  with 0 ≤ P_∂ ≤ Z^{−1}Σ_∂(N(c)ψ′ψ″ + |Err(c)|). And |Cov_ν(ψ′, ψ″)| ≤ (Var_νψ′·Var_νψ″)^{1/2}, Var_νψ′
  ≤ E_ν[(ψ′ − c)²] ≤ E_{ν₁}[(ψ′ − c)²] + ν(∗)c² ≤ (𝒱₀ + 20h + e^{−λ})c² (ν's first marginal is dominated by the one-point
  law ν₁ of Lemma OPV; Var ≤ E(· − c)² for any constant), and c ≤ 2S′ (Lemma OPV).
- *Boundary terms.* By Cauchy–Schwarz, b′ ≤ (Z^{−1}Σ_∂N)^{1/2}(Z^{−1}Σ_∂Nψ′²)^{1/2} and Z^{−1}Σ_∂Nψ′ψ″ ≤
  (Z^{−1}Σ_∂Nψ′²)^{1/2}(Z^{−1}Σ_∂Nψ″²)^{1/2}; and Z^{−1}Σ_∂Nψ′² ≤ (Z^{−1}Σ_∂N)^{1/2}(Z^{−1}Σ_{R⊂Π̃}Nψ′⁴)^{1/2}, where Π̃ ⊃ Π
  is the union of all rectangles meeting Π; Π̃ is not convex, but Π̃ ⊂ Π⁺ (Lemma CELL(ii)), which is convex with
  A^{−1}(Π⁺) a region as in Lemma LC. So the last sum is ≤ a one-point fourth moment over A^{−1}(Π⁺), ≤ CS′⁴ by the
  bounded relative moments of Lemma OPV (applied to A^{−1}(Π⁺))
  (moments, not maxima: the spread ratio Ξ never enters). With Z^{−1}Σ_∂N ≤ 20h (Lemma CELL(ii)) all boundary terms are
  ≤ Ch^{1/4}S′S″.
- *Weights.* On ordinary blocks 4ε′² ≤ C(γ_a²h + E_a²) ≤ 𝒱₀; on the others (relative width ≤ 1/8, so ε′ ≤ γ_a/8 + E_a ≤ 1/4
  and Cor CPW applies) 4ε′² ≤ Cγ_a², and these carry mass ≤ Ch^{1/2}. This part of Err(c) sums to ≤ C𝒱₀S′S″.
- *Root event.* 𝒢_e excludes the sublattice {q | G_g(r, 1)}, on which f′ = f″ = 0; only Z changes, by 1 − O(1/q), so the
  covariance moves by O(S′S″/q). ∎

**Lemma EXC (exceptional types) [proved].** Let ℰ be the n ∈ L ∩ Π with (a) gcd(n₁n₂, K_large) > 1, or (b) (parents) ℓ | n₁n₂
for one of the ≤ 2 primes ℓ ∈ J dividing K, or (c) s₁ > S₀ or s₂ > S₀, or (d) n_i in the single grid interval containing
a filter threshold. Then Z^{−1}Σ_{n∈ℰ}f′(n₁)f″(n₂) ≤ (e^{−λ} + κ₀^{−1/2})c′c″ and Z^{−1}Σ_{n∈ℰ}f′(n₁) ≤ e^{−λ}c′, and
ν(ℰ) ≤ e^{−λ} (c′, c″ the constants of Lemma OPV); (d) contributes ≤ Ξh_*·S′S″ ≤ κ₀^{−1/2}S′S″ (children: mass ≤ 8h_*X⁻/x₁ ≪ N^{−η₀/8}; parents:
mass ≤ Ch_*, times the PUB loss Ξ_p).
*Proof.* (a) contributes 0 to every f-sum (a y-smooth n has no prime > y). For (b) and (c), the n with σ | n₁ for a fixed
integer σ (σ = ℓ, or σ a P₀-smooth s > S₀ dividing the P₀-part) come from the g on the lines g* + ZP′ with σ | G′(g);
running Lemma PUB's proof on these lines only gives Z^{−1}Σ_{σ|n₁}f′f″ ≤ Ξ·S′_σS″, where S′_σ := Z^{−1}Σ_{σ|G′}f′ ≤
e^{O(λ)}σ^{−1+γ}S′ (Lemma FO applied to the class σ | n₁; for σ = ℓ ∈ J the density of J·𝒮_B among multiples of ℓ is
≤ (log Y)·(its global density)). Summing, (b) gives ≤ Ξ(log Y)²Y^{−1}S′S″ and (c) gives ≤ ΞΣ_{s>S₀}s^{−1+γ}S′S″ ≤ Ξ·S₀^{−1/(2log P₀)}(log P₀)^C S′S″
(Rankin, exponent 1/log P₀); both are ≤ e^{−λ}S′S″ because Y and S₀^{1/log P₀} = e^{D^{(θ+α)/2}/log P₀} exceed every
power of Ξ = e^{O(D^θλ)}. (d): the threshold interval has relative width h_*, so its one-point mass is ≤ 8h_*·(position)/H
(Lemma LC), and PUB on the sub-sum gives the factor Ξ; κ₀ = T₀³ ≥ Ξ^{30}. The singles are the
same computation without the partner. ∎

**Remark (why first order suffices).** Every first-order quantity (E_H, 1/u, E_a, the archimedean factor's deviation
γ log(x₁/|n|), the local factor's deviation γ log s) appears only inside one-point L² norms, which are multiplied
together. No pair main term is ever computed or compared with the product of marginals. This is
"common normalized model", with the model chosen per pair: since S′ = E_νψ′ holds exactly for *each* pair's own cell
structure, different partners may use different cells. For interior cells the marginal identity is exact;
boundary and exceptional cells change it by the terms B′, P_∂, P_ℰ above, which are bounded separately.

### 8.2 Lemma PRIM (primitivity at primes > P₀) [proved]

Let Ω, P′, P″, f′, f″ be as in Prop PD and let 𝒢^prim ⊂ 𝒢 be the primitive points. Then
  |Cov_{𝒢^prim}(f′, f″) − Cov_𝒢(f′, f″)| ≤ C[(γ²log²P₀ + E_FO²)/(P₀ log P₀) + P₀^{−2}]S′S″ + Σ_{1<e≤y², P⁻(e)>P₀}e^{−2+2γ}|R_e|,
where R_e is the Prop PD remainder for Ω/e with weights f(e·) (bounded by Prop PD at Ω/e: (V′V″)^{1/2}_{Ω/e} + 𝔈_{Ω/e} + …;
its 𝔈-part joins the family of Lemma FAM with the extra weight e^{−2+2γ}, which is summable).
*Proof.* g ∈ 𝒢 is non-primitive iff some prime p > P₀ divides g, and then p | G′(g), G″(g), so f′f″ = 0 unless p ≤ y.
Bonferroni with terms e ∈ {1} ∪ {p} ∪ {pp′} (P₀ < p < p′ ≤ y): for F ∈ {f′f″, f′, f″, 1},
Σ_{𝒢^prim}F = Σ_𝒢F − Σ_pΣ_{𝒢(Ω/p)}F(p·) + O(Σ_{p<p′}Σ_{𝒢(Ω/pp′)}F(pp′·)), where 𝒢(Ω/e) are the P₀-primitive points of Ω/e
(e has no prime ≤ P₀); this is the Bonferroni inequality 1 − Σ_p1_{p|g} ≤ 1_prim ≤ 1 − Σ_p1_{p|g} + Σ_{p<p′}1_{pp′|g}.
The last sums are ≤ (pp′)^{−2+2γ}(1 + 𝒱₀)#𝒢S′S″ + (pp′)^{−2}|R_{pp′}| (Prop PD at Ω/pp′ for F = f′f″, Lemma FO for the
singles), total O(P₀^{−2}) + Σ_e e^{−2}|R_e|. Put w_p := −#𝒢(Ω/p)/#𝒢 (≈ −p^{−2}), w₁ := 1, and a_e := E_{𝒢(Ω/e)}f′(eG′), b_e likewise, π_e :=
E_{𝒢(Ω/e)}f′(eG′)f″(eG″). Then, normalising, Cov_prim = Σ_e w̃_eπ_e − (Σ_e w̃_ea_e)(Σ_e w̃_eb_e) + O(P₀^{−2}S′S″) + Σ_e e^{−2}|R_e|, with
w̃ := w/Σw (signed, Σw̃ = 1, Σ|w̃ − δ_{e1}| ≤ 2/P₀). For signed weights summing to one,
Σw̃_eπ_e − (Σw̃a)(Σw̃b) = Σ_e w̃_e(π_e − a_eb_e) + Σ_e w̃_e(a_e − ā)(b_e − b̄). The first sum is Cov_𝒢 + Σ_{p}w̃_pR_p-type
terms (Prop PD at Ω/p). In the second, a_p = p^{γ}a₁(1 + O(E_FO + γ log p)) by Lemma FO (the density of f′ on
multiples of p), so |a_p − a₁| ≤ C(γ log p + E_FO)p^{γ}a₁, and |a₁ − ā| ≤ Σ_p|w̃_p||a_p − a₁|; the sum is
≤ CΣ_{p>P₀} p^{−2}(γ log p + E_FO)²p^{2γ}a₁b₁ + C(Σ_pp^{−2}…)². Here p^{2γ} ≤ e for p ≤ e^{1/(2γ)} and the remaining
primes contribute ≤ e^{−1/(4γ)}; and Σ_{p>P₀}log²p/p² ≍ log P₀/P₀, Σ_{p>P₀}p^{−2} ≍ 1/(P₀ log P₀). So the sum is ≤ C(γ²log²P₀ + E_FO²)/(P₀ log P₀)·S′S″ ≤
CD^{−2α−7/10}λ⁶S′S″. For F = 1 the primes p > y also enter the normalisation; they contribute O(1/y). For a_p with
P₀ < p ≤ y, Lemma FO is applied with s = p via Fact HT3 (p ≤ y), not via S₀. ∎
(At toy scale, where γ ≈ 0.3, this term is visible: E3's R ≈ −0.08 under the primitive law is exactly the
untyped-prime primitivity covariance, §10.3.)

### 8.3 Proof of Theorem RAT

Choose h := D^{−A₁} with A₁ := 8α + 1, δ_e := D^{−A₀} with A₀ ≥ 1, P₀ := ⌈D^{7/10}⌉, T₀ := e^{C_TD^θλ²} and κ₀ := T₀³
(§7).

*(RAT-1).* Fix a window (Q₂, 2Q₂], a label with s(𝔮) ≤ T₀ and a cell e; the OPV constants c_{P′} of all box points are
within a factor Ξ of each other and c_{P′}/2 ≤ S′ ≤ 2c_{P′} (Lemmas SPR, OPV). Split the ordered pairs:
- **Thin pairs** (κ > κ₀): by Prop SK(iii) their contribution to Φ_𝔮 is ≤ T₀^{−1}.
- **Non-thin pairs.** Prop PD and Lemma PRIM give, for each pair,
  |Cov(f′, f″)| ≤ (V′V″)^{1/2} + 𝔈 + Σ_{e}e^{−2+2γ}(…)_{Ω/e} + C[h^{1/4} + (γ²log²P₀ + E_FO²)/P₀ + P₀^{−2} + κ₀^{−1/2}]S′S″.
  Here V′ ≤ 4𝒱₀(K)S′² and V″ likewise, with the **same** K (𝒱₀ contains the mid-prime term of K), so
  (V′V″)^{1/2} ≤ 4𝒱₀(K)S′S″ ≤ C D^{−2α}λ^{6}S′S″ (Lemma OPV, pointwise in K). Summing over pairs and dividing by m² = (Σ S′)² gives at most
  CD^{−2α}λ⁶ + C(h^{1/4} + λ⁶D^{−2α}/P₀ + P₀^{−2} + κ₀^{−1/2}), i.e. ≤ C D^{−2α}λ⁶ by the choice of A₁
  (h^{1/4} ≤ D^{−2α−1/4}) and P₀^{−2} = D^{−7/5} ≤ D^{−2α} on α ≤ 7/10.
- **Character errors.** The 𝔈-terms of all non-thin pairs of all labels with s(𝔮) ≤ T₀ and all cells (including the
  Ω/e terms of Lemma PRIM, with weights e^{−2+2γ}), multiplied by w_{𝔮,e}/m_{𝔮,e}² ≤ ŵ_{q,e}/m², form the family of
  Lemma FAM; Cor WA′ bounds their total by e^{−cD^α/λ²}μ₁.
Multiplying the per-label bounds by w_{𝔮,e} and summing: Σ_{s(𝔮)≤T₀}w_{𝔮,e}Φ_{𝔮,e} ≤ (CD^{−2α}λ⁶ + T₀^{−1})W +
e^{−cD^α/λ²}μ₁, which is (RAT-1). ∎

*(RAT-2).* As (RAT-1), with: the parent library 𝒜^J_{ξ,s}; the parent multiplicities (Lemma MULT, ≤ 5 per K) and
Pr_{P′}[p | det(P, P′)] ≤ (1/p)(1 + O(pq/V₁)); no skew (Prop SK(iv); parent thin pairs via the absolute PUB);
the J-structure (the ≤ 2 primes ℓ ∈ J dividing K are exceptional, Lemma EXC(b)); Lemma FO applied for each ℓ ∈ J;
and the **hard filter** |G(P)| ∈ [MV₁/32, 2MV₁], which needs care:
- by Lemma SPR, m^par_𝔮 ≥ n·c_typ/polylog (active points), while every OPV constant satisfies c_P ≤ c_typ·polylog;
- for a parent P let π_P ∈ [0, 1] be its pass fraction on Ω_e. Centering at π_Pc_P instead of c_P in Lemma OPV gives
  Var_{ν₁}(ψ_P) ≤ c_P²(π_P(1 − π_P) + 𝒱₀), and π_P(1 − π_P) ≠ 0 only for the *straddling* P (threshold inside the range
  of |G(P)| on Ω_e), which form a fraction ≤ Cδ_e of 𝒫 ∩ Λ_𝔮 for each cell (Lemma SPR's angle argument);
- Prop PD's argument goes through with all relative bounds read against c_P (absolute form), so
  Σ_{P≠P′}|Cov| ≤ Σ_{P,P′}c_Pc_{P′}(r_P + r_{P′})/2 + (𝔈 and smaller terms), r_P := π_P(1 − π_P) + 𝒱₀, and
  Σ_{P,P′}c_Pc_{P′}r_P ≤ (n c_typ polylog)²(Cδ_e + 𝒱₀) ≤ (m^par)²·polylog·(δ_e + 𝒱₀);
- the parent weights' first-order error is E_a ≍ D^{7α−3} + D^{1−4α} (Lemma ADM), so 𝒱₀^par ≤ C(D^{−2α} + D^{−2(3−7α)})λ⁶,
  and A₀ ≥ 1 > 2α makes the δ_e term smaller.
pc3's O(1) environment term M1 (the filter's g-dependent pass fraction) is exactly what the cells remove. Hence
(RAT-2): Σ_{𝔮,e}w_{𝔮,e}Φ^par_{𝔮,e} ≤ C(D^{−2α} + D^{−2(3−7α)})λ^{C}W + e^{−cD^α/λ²}μ₁. ∎

*(RAT-0).* Prop SK(i). ∎

---



---

# B′. RAT precision margins (source: notes/c12rat_results.md §9.2)

### 9.2 How much precision the rational blocks leave

Theorem RAT gives A_rat ≤ D^{−2α}λ^C and B_rat ≤ (D^{−2α} + D^{−2(3−7α)})λ^C. The window's budget is D^{−α}λ^C for A and D^{7α−3}λ^C for B
(c12alg Prop PB). So the margins are D^{α} for C₁ and D^{min(9α−3, 3−7α)} for C₂ — positive on (3/8, 2/5), and at
α = 3/8, against pc3's η ≍ D^{−3/8}, D^{3/8} for both. In the language of the task: the rational blocks do not
constrain α on [3/8, 2/5); c12alg's Theorem 2L-DW₂ uses B ≪ (D^{−2α} + D^{−2(3−7α)})λ^C, exactly what RAT-2 gives.



---

# C. Prop DC / DC_e, Lemma 0, Lemma LAT (source: notes/pc3_results.md §§1-3)

## 1. Setting

Notation as multidescent §1, §3.2, and uprime §1 (design 2L′: d ≍ D^{1/4}, B = L[3/8, a], B₂ = L[1/2, b₂], ℓ ∈ J = L[5/8],
uprime's height floor (H) h = μ + 4ν + η₀ and acceptance rule). Write η := μ₁/W, where W/μ₁ = E[m(P)]. Typically
W/μ₁ ≍ D^{3/8}; the deterministic cap gives ≤ m̄ ≍ D^{1/2}.
- A **label** is 𝔮 = (q, r), with Λ_𝔮 = {(a, b) : a ≡ rb (mod q)} and Gauss-reduced basis e₁, e₂.
- **Box(𝔮)** is the set of Z²-primitive points P′ = ie₁ + je₂, max(|i|, |j|) ≤ S₂, modulo ±1 (uprime's conservative
  geometry; S₂ from uprime's child scan).
- N_𝔮(A) := #{P′ ∈ Box(𝔮) : Good(P′)}.
- Each incidence (ℓ, P, q) has label 𝔮(P, q) = (q, ab^{−1} mod q). The **parent count** is
  Z_𝔮(A) := #{(ℓ, P) : Acc_ℓ(P), q is a child of P, 𝔮(P, q) = 𝔮}.
- The **root event** is R_𝔮 := {F_A(r, 1) ≡ 0 (mod q)}.
- The saddle point of Ψ(·, y) is written β (not α, which is reserved for the L-exponent).

**Lemma 0 [proved].** Z_𝔮 = N_𝔮 = 0 outside R_𝔮.
*Proof.* A Z²-primitive P′ ∈ Λ_𝔮 has q ∤ b′, so F(P′) ≡ b′^dF(r, 1) (mod q). Both Good(P′) and "q is a child of P" need
q ∥ F at such a point. (A point primitive only in coordinates, e.g. q·(1, 0), has q^d | F and fails q ∥ F.) ∎

*Size of Z_𝔮 [heuristic].* z_𝔮 := E[Z_𝔮 | R_𝔮] ≍ (Y/q)·m₁/(d log²N) = L[5/8]^{1−o(1)}. Every label has many potential
parents.

## 2. Decoupling [proved]

**Proposition DC.** Put m_𝔮 := E[N_𝔮 | R_𝔮] and ϕ_𝔮 := Var(N_𝔮 | R_𝔮)/m_𝔮² (ϕ := 1 if m_𝔮 = 0). Put
σ_𝔮 := Var(Z_𝔮 | R_𝔮)^{1/2}/z_𝔮, w_𝔮 := Pr(R_𝔮)z_𝔮/#J, and W := Σ w_𝔮. Then
  Σ_{(ℓ,P,q)} Pr_{ℓ,A}[E_{ℓ,P,q} ∧ N_{𝔮(P,q)} = 0] ≤ Σ_𝔮 w_𝔮(p_𝔮 + σ_𝔮√p_𝔮),   p_𝔮 := Pr(N_𝔮 = 0 | R_𝔮) ≤ min(1, ϕ_𝔮).
*Proof.*
- Every incidence of a field is counted under exactly one label. 𝒫 is common to all ℓ, and Box(𝔮) and Good(P′) do not
  depend on ℓ or P. Hence the left side is (1/#J)Σ_𝔮 E[Z_𝔮1{N_𝔮 = 0}], which by Lemma 0 equals
  Σ_𝔮 (Pr(R_𝔮)/#J)E[Z1{N = 0} | R].
- Under the law given R: E[Z1{N=0}] = zp + E[(Z − z)1{N=0}] ≤ zp + sd(Z)·√(p(1−p)).
- Finally, Chebyshev gives p ≤ ϕ. ∎

**Corollary DC.** Put A := Σ(w_𝔮/W)ϕ_𝔮 and
B := Σ(w_𝔮/W)σ_𝔮². Then the dead-incidence mass is ≤ W(A + √(AB)). In particular, A ≤ η/4 and B ≤ η/4 imply
E N₂ ≥ μ₁/2.

**Proposition DC_e (environment conditioning) [proved].** Let 𝒢 be any σ-algebra generated by a function e(A) of the
field (an "environment cell") such that N_𝔮 and Z_𝔮 are defined cell by cell. Then the same bound holds with
Σ_{𝔮,e}w_{𝔮,e}(p_{𝔮,e} + σ_{𝔮,e}√p_{𝔮,e}). Here every quantity is conditional on R_𝔮 ∩ {e(A) = e}, and
w_{𝔮,e} := Pr(R_𝔮, e)z_{𝔮,e}/#J. Cauchy–Schwarz then runs over (𝔮, e).
*Proof.* Apply the proof above under each conditional law and sum. ∎
- *Use.* There is an O(1) field-level term in σ_𝔮². The parent filter |G(P)| ≥ MV₁/32 has a pass fraction
  depending on g's direction (it varies by ≈ 1/32), so E[Z_𝔮 | g]/z_𝔮 is not concentrated, and σ² ≥ c₀ > 0.
  - With an O(1) σ², Corollary DC needs A ≲ η²/c₀ ≍ D^{−3/4}, which has no margin.
  - Taking e(A) := (cell of g of relative size D^{−A₀} in direction and scale) makes the pass fraction constant to
    O(D^{−A₀}) within each cell. C₁ and C₂ are then needed cell by cell.
  - The BDH/LSQ arguments are local in g, since a cell contains M²D^{−2A₀} values of g.
  - Any further thresholds should be randomised **independently per point**. A per-field random child-box radius would
    make ϕ = Θ(1).
- *Alternative.* Hölder gives E[Z1{N=0}] ≤ ‖Z‖₃p^{2/3}. With ϕ ≈ D^{−3/4} log D, a bounded third relative
  moment of Z_𝔮 suffices; losses up to D^{1/8−ε} are allowed. That tolerates O(1) environments, but needs three-point
  *upper* bounds for parents, which we do not have.

**Theorem 2L-DC [proved, conditional on C₁, C₂ and on what Theorem 2L uses besides PC₃].**
- **Hyp C₁:** A ≤ η/4, the incidence-weighted root- and cell-conditioned relative variance of the child count.
- **Hyp C₂:** B ≤ η/4, the same for the label's parent count.

Assume C₁ and C₂ in design 2L′ with environment cells, together with Hyp U′ (proved by uprime at α = 3/8, apart from the
CM bullet, which DC does not use). Then FACTOR_2L factors every N in expected time L_N[3/8, c].
*Proof.* Theorem 2L's proof (multidescent §3.3) uses PC₃ only through Prop UB's E N₂ ≥ μ₁/2. Corollary DC (or DC_e)
supplies it. The Paley–Zygmund step uses E[N₂²] ≤ E[N₁²], which U′ bounds. ∎

*What changed relative to PC₃.*
- PC₃ conditions on E_{ℓ,P,q}, a configuration containing the parent and the child points: three linear forms in g with
  coefficients up to x^{O(D^{−1/8})}. C₁ and C₂ are second moments of sums over a lattice family, so every configuration
  has two points.
- Lemma MP still applies: precision O(η) is necessary, and it is all that is asked.
- The CM bullet of U′ is replaced by the one-point lower bound m_𝔮 ≥ 8/η inside C₁. It concerns child points under the
  root-conditioned law, and is a GM-type statement at the child level: uprime's nominal m₂ ≥ exp((log log N)²) plus a
  one-point lower bound with polylog loss. The latter is a re-run of FM/GM at the child level with B-smooth moduli
  G(P′) (Lemma 8.5′ Step C), or pointwise under DW for filtered points [sketch].

## 3. Pair lattices and cancellation [proved]

**Lemma LAT.** Let P′, P″ be primitive and non-proportional, and D := |det(P′, P″)|.
- (a) There is a unique unit λ (mod D) with P′ ≡ λP″ (mod D).
- (b) {(G_g(P′), G_g(P″)) : g ∈ Z²} = {n₁ ≡ λn₂ (mod D)}.
- (c) Φ(P′) ≡ λ^dΦ(P″) (mod D) for every integral binary form Φ of degree d.
- (d) For P′, P″ ∈ Λ_𝔮, D = qΔ with Δ = |ij′ − ji′| ≤ 2S₂² < q. Here λ mod Δ depends only on the coordinates, and
  λ ≡ b′/b″ (mod q).

*Proof.* (a)–(c) as in v1: extend P″ to a basis; the index count gives (b); homogeneity gives (c). (d): the coordinate
vectors are primitive, so apply (a) to them. Here Δ ≤ 2S₂² < q because
S₂² = L[3/8, σ₂] < B ≤ q once a is large. ∎

**Lemma SI (scaling invariance) [proved].** If Ω̄ ⊂ (Z/K)² is a subgroup, the multiset {gcd(a₁, K) : (a₁, a₂) ∈ Ω̄}
depends only on gcd(a₂, K). Consequently, in the class decomposition
  Σ_{(a₁,a₂)∈Ω̄} S₁(a₁)S₂(a₂),  S = M(gcd class) + E,
the cross terms Σ M₁E₂ and Σ E₁M₂ vanish exactly. Here M is the average over *all* residues of the gcd class, and both
projections of Ω̄ must be surjective; they are for LAT(b) and DW-L.
*Proof.* Residues with equal gcd differ by a unit factor, and multiplying both coordinates by that unit is a bijection of
the fibres preserving first-coordinate gcd's. ∎

**Pair multiplicities [proved].** For a fixed primitive coordinate vector (i, j), the number of
(i′, j′) in the box with |ij′ − ji′| = Δ is ≤ 2S₂/|(i, j)|_∞ + 2. So #{pairs with determinant quotient Δ} ≪ S₂²
uniformly in Δ. There is no divisor function: a τ(Δ) weight would have been exp(O(aD^{3/8})), which swamps BDH's saving.



---

# C′. Design DW, Lemma DW(a), DW-L (source: notes/pc3_results.md §6.2)

**6.2 Design DW.** The G·R_c randomisation is Lee–Venkatesan's (JNT 187 (2018) 92–159; arXiv 1805.08873):
f = f̂ + Σc_i(x − my)y^ix^{d−1−i}. Their "stochastic deepening" (Lemma 5.6) is the closest prior device for avoiding
second-moment bounds. What DW adds is the integer multiplier w of N, an analogue of the Knuth–Schroeppel kN. We found
no rigorous-NFS precedent for it; novelty is UNVERIFIED. Set F_A := w·F̂_{g,N} + G·R_c, with w uniform among integers in (W_m, 2W_m] coprime to g₁. The
multiplier range W_m is renamed so as not to clash with the incidence count W. Take W_m ≥ (2MV₁)²N^{η₀} and
H ≥ W_mN^{η₀}, i.e. h ≥ 2μ + 2ν + 2η₀.
- F_A(g₀, g₁) = wN, and the resultant is ±wN.
- At level 1, μ ≍ D^{−1/4} is lower order than dν ≍ D^{−1/8}, so θ is unchanged at leading order.
- At level 2 the extra μ changes only the constant of the child values (σ₂), as uprime's (H) does.

**Lemma DW(a) [proved for filtered points].** Fix g and P with t = |G(P)| ≥ M|P|N^{−η₀/4}, gcd(N, t) = 1 and
gcd(b, g₁) = 1. Then F(P) ≡ w·F̂_N(P) (mod t), with F̂_N(P) a unit mod t, so F(P) mod t is uniform over w up to
O(tτ(g₁)/W_m). The c-part has [pp] Lemma smooth1's law at scale t·H|P|^{d−1}; the shift |wF̂_N(P)| is
≤ N^{−η₀/4} × that scale by the filter. So Pr_{w,c}[F(P) ∈ 𝒮] is the density of 𝒮 against an explicit smooth unimodal
weight, pointwise in g.
- Without the filter (e.g. G(P) = 1) the claim is false.
- Parents satisfy the filter. Child points need it added to Good(P′); its g-dependent pass fraction is absorbed by DC_e.

**Lemma DW-L (exact support lattice) [proved].** For P₁, P₂ with gcd(N, G(P₁)G(P₂)) = 1 and
primitive P_i, the set {(F_A(P₁), F_A(P₂))} over all (w, c) ∈ Z × Z^d with F_A primitive-compatible equals
{F₁ ≡ λ^dF₂ (mod |det(P₁, P₂)|)}. The order of v := (F̂_N(P₁), F̂_N(P₂)) modulo (G(P₁), G(P₂))·I is |G(P₁)G(P₂)|.
*Proof.* Two forms with prescribed values at P₁, P₂ differ by L_{P₁}L_{P₂}R. Their values at g differ by
G(P₁)G(P₂)Z. Coprimality with N allows imposing F(g) ∈ NZ, and Gauss's lemma then gives F = wF̂_N + G·R_c. ∎ (Checked on 1 624 small configurations.)

**Lemma DW(b) (two-point law) [sketch].** One would want: the law of (F(P₁), F(P₂)) is within TV ε of
the smooth-weighted uniform law on the DW-L lattice. Kern-U gives only fibre upper bounds. What is needed is a
two-point analogue of [pp] smooth1: a sum of d independent segments in R², anisotropic for nearly parallel pairs. Under
the root event R_𝔮 one also needs smoothing over a coset mod q (H ≫ qdV) and the q² | det case. Not written.

**6.3 Consequences [sketch].** Given DW(b), the algebraic pair for fixed g is a det-modulus lattice smooth-pair problem.
It is handled by the same tools as §5:
- Harper's Thm 2 (BDH, any moduli Q ≤ Ψ) at y = B for children and y = B₂ for parents, both with u ≍ D^{3/8}.
  Drappeau's Lemme 5 alone covers only |det| ≤ x^η;
- no mod-q ratio arises for F/q;
- the moduli are |det|, not G(P′), so the endogenous family of 6.1 disappears.

For Theorem U itself (relevant to `uprime`, not our file): DW(a) makes FM's algebraic factor pointwise for filtered
points, so FM would not need def:good or Lemma 8.5′. uprime's Lemma 8.5H is then needed only in the original design.
DW is **unchecked** against the internal chain. Three steps it breaks, all repairable:
- with integer w, gcd(w, t) > 1 is possible, so lem:coprime fails;
- SAMPLE_FIELD's "content 1 or w" fails, since any divisor of w can occur, so non-primitive forms must be rejected;
- Lemma IR's certificate primes need r ∤ w.



---

# D. Size model, Prop UB, Theorem 2L (source: notes/multidescent_results.md §§3.1-3.3)

## 3. Exponents and what must be proved (Q2)

### 3.1 The size model [model; `model.py`, tests]

j-level chain: level-i targets L[ζ_i], ζ₀ = 1 − α
(split), ζ_j = α (factor base); d = D^δ; base-m heights (exponent 1 − δ); smoothness of a value of exponent v over
L[ζ′] at cost L[α] needs v − ζ′ ≤ α. Constraints: 1 − δ ≤ 2α and δ + ζ_i − ζ_{i+1} ≤ α. Summing, jδ + 1 − 2α ≤ jα,
so **α_j = (j+1)/(3j+2)**, δ_j = j/(3j+2), ζ_i = 1 − α − i(3α − 1). Every level has D^{δ + 3α − 1} = D^{α} large algebraic
children per node. Two levels: α = 3/8, d = D^{1/4}, ζ₁ = 1/2, level-1 algebraic value D^{7/8}, level-2 values
D^{3/4}, m ≍ D^{3/8} = d^{3/2} children. The model bounds the exponent of designs of this shape; it proves no supply.

*Per-child losses.* A factor (log N)^{−Cm} from requiring all m children to
succeed is affordable as **search work per target** (retry until a parent with good children is found:
exp(O(C D^α ln t)) = L[α, O(C)]). It is **not** affordable as a loss of **root density**: the split raises the
descendable density to the power 2k ≍ D^{α}, so (log N)^{−Cm·2k} = exp(−Θ(D^{2α} ln t)) is outside L[α]. navigate's
bookkeeping implicitly assumes the first reading, with independent children (nested Prop BT's model). Prop UB
avoids the issue: under PC₃ the two-level count keeps a constant fraction of the level-1 first moment.

### 3.2 Prop UB: the union-bound reduction [proved]

*Design 2L* (two levels at α = 3/8; the precise parameters are those of lv5unif §1 with d ≍ D^{1/4}, k ≍ D^{3/8},
B = e^{aΛ}, an extra child bound B₂ = e^{b₂D^{1/2}ln t}, and a child-box radius S₂). For ℓ ∈ J and a field A:
- Acc_ℓ(P): P ∈ 𝒫 ∩ Λ_ℓ(g), G(P)/ℓ ∈ 𝒮_B, and F(P) = s·q₁⋯q_m with s B-smooth and q_i distinct primes in (B, B₂] not
  dividing b·D_A (children 𝔮_i = (q_i, ab^{−1} mod q_i)); m(P) := m.
- Box(𝔮): the primitive points, modulo ±1, of the reduced-basis box of Λ_𝔮 of radius S₂.
- Good(P′) for P′ ∈ Box(𝔮): q ∥ F(P′), G(P′) ∈ 𝒮_B, F(P′)/q ∈ 𝒮_B, plus the gcd checks.
- N₁(ℓ, A) := #{P : Acc_ℓ(P)}, N_𝔮(A) := #{P′ ∈ Box(𝔮) : Good(P′)},
  N₂(ℓ, A) := #{P : Acc_ℓ(P) and N_𝔮(A) ≥ 1 for every child 𝔮 of P}.
- For a triple (ℓ, P, q) let E = E_{ℓ,P,q} := {A : Acc_ℓ(P) holds and q is a child of P} (an event in the field law,
  ℓ fixed). Probabilities are over ℓ uniform in J and A ∼ 𝒜.

**Hypothesis PC₃.**
  Σ_{(ℓ,P,q)} Pr[ℓ, E_{ℓ,P,q}]·ϕ_{ℓ,P,q} ≤ μ₁/2,  where ϕ := Var(N_𝔮 | E)/E[N_𝔮 | E]² if E[N_𝔮 | E] > 0, and ϕ := 1 otherwise,
and μ₁ := E N₁. Here Var(N_𝔮|E)/E[N_𝔮|E]² = 1/E[N_𝔮|E] + (E[N_𝔮(N_𝔮−1)|E] − E[N_𝔮|E]²)/E[N_𝔮|E]², and
E[N_𝔮(N_𝔮 − 1)|E] = Σ_{P′≠P″} Pr[Good(P′) ∧ Good(P″) | E] concerns the events E ∧ Good(P′) ∧ Good(P″).

*A sufficient pointwise form, PC₃(η, ε):* outside a set 𝔈 of triples with Σ_𝔈 Pr[ℓ, E] ≤ εμ₁, (a) E[N_𝔮|E] ≥ 4/η and
(b) E[N_𝔮(N_𝔮 − 1)|E] ≤ (1 + η/2)E[N_𝔮|E]²; then ϕ ≤ 3η/4 there, and PC₃ follows when (3η/4)W + εμ₁ ≤ μ₁/2. Note
that the pointwise form forces an exceptional incidence fraction ≍ 1/m, so the weighted form is the natural one.

**Proposition UB.** Let μ₁ := E N₁ > 0 and W := E[Σ_{P : Acc} m(P)] (expected number of parent–child incidences;
Σ_{(ℓ,P,q)} Pr[ℓ, E_{ℓ,P,q}] = W). If PC₃ holds, then E N₂ ≥ μ₁/2 and
  Pr_{ℓ,A}[N₂ ≥ 1] ≥ μ₁²/(4·E[N₁²]).

*Proof.* Pointwise, N₂ ≥ N₁ − Σ_{P}Σ_{𝔮 child of P} 1[N_𝔮 = 0] (an accepted parent that is not counted in N₂ has at least
one dead child). Taking expectations, E N₂ ≥ μ₁ − Σ_{(ℓ,P,q)} Pr[ℓ, E ∧ N_𝔮 = 0]. For each triple, Chebyshev in the
conditional law gives Pr[N_𝔮 = 0 | E] ≤ ϕ_{ℓ,P,q} (trivially if the conditional mean is 0). Hence
Σ Pr[ℓ, E ∧ N_𝔮 = 0] ≤ Σ Pr[ℓ, E]ϕ ≤ μ₁/2. Finally 0 ≤ N₂ ≤ N₁, so by Paley–Zygmund
Pr[N₂ ≥ 1] ≥ (E N₂)²/E[N₂²] ≥ (μ₁/2)²/E[N₁²]. ∎

*Well-definedness*: for fixed (ℓ, P, q), E_{ℓ,P,q} is an event in the field law; membership
P ∈ Λ_ℓ(g) is part of it, and on it the child label (q, ab^{−1} mod q), hence Box(𝔮), is determined by (P, q).

**Remarks.**
- No statement about two children at once is used: their failure events are *added*, not intersected. The
  configurations have three **geometric points** (P, P′, P″) and up to six evaluations (G and F at each), plus the
  global acceptance conditions (discriminant, excluded primes). The rational part is three linear forms in the two
  coordinates of g; the algebraic part three evaluations in (w, c), of full rank in c for d ≥ 3 and pairwise
  non-proportional points (Kern-U's kernel argument with Q₃ = L_P L_{P′} L_{P″}).
- The price is precision: the weighted average of ϕ must be ≤ μ₁/(2W). Under the typical-factorisation count
  W/μ₁ ≍ D^{3/8} (so precision ≍ D^{−3/8}); the deterministic cap m(P) ≤ m̄ = O(D^{1/2}) only guarantees that
  precision D^{−1/2} suffices. By Lemma MP (almostall13) any moment-based proof that a child
  fails with probability < η needs moments accurate to relative error η, so this precision is necessary for the
  route. It is polylogarithmic, not L-scale: the union is over the children of one parent, not over all e^{Θ(Q)}
  targets of a field (AD's union bound, nested §6 / almostall13, which needs L[2/3] precision).
- *Precedent.* Removal arguments that sum bad incidences have precedent in rigorous descent: GKZ's proof of their
  Theorem 3 removes the ≤ q^{kd/2+4} values producing trap roots from ≥ q^{kd−4} solutions; KW avoid a bounded-degree
  trap subvariety. Those counts are deterministic and hold pointwise for every
  good target — the AD-type property unavailable here — so the analogy is formal. We found no prior use of a conditional 3-point smooth-number reduction of this kind.
- E[N₂²] ≤ E[N₁²] means the second moment needed is the *one-level* one (LD*-type, upper bound with polylog loss).

### 3.3 Theorem 2L [proved conditional on PC₃, U′ and the internal chain]

**Hyp U′ [OPEN at the time of that note; supplied by uprime].** Theorem U's FM-U, LD*-U, D-U hold for Acc_ℓ at the
3/8 parameters, with computable box radii S₁, S₂, giving μ₁ ≥ (log N)^{−C}, E[N₁²] ≤ (log N)^{C}(μ₁² + μ₁), search
sizes W₁, W₂ = L[3/8, O(1)] and E[N_𝔮|E] large enough for PC₃. What has to be re-proved:
- *Mixed smoothness.* The rational cofactor is B-smooth, the algebraic value B₂-smooth. Lemma 8.5′(R′) uses one y for
  the values and for the smooth parts of the restricted moduli. With y = B₂ it controls a much larger modulus family
  whose error e^{−cD^{3/8}/log²D} does not beat ϱ(z, B₂)/ϱ(z, B) = e^{Θ(D^{3/8} log D)}; with y = B the algebraic
  u ≍ D^{1/2} violates (R′3) (u log(4u log B)/log B ≍ D^{1/8} → ∞). A hybrid restricted BV —
  B-smooth modulus cofactors, B₂-smooth values, normalised by the density of the actual modulus family — is needed;
  GM-U, L2-U, DS-U must be redone with two smoothness parameters.
- *Box sizes.* Lemma CS evaluates both Rankin factors at B; at the root algebraic size D^{7/8} that gives an
  L[1/2]-scale box. A mixed (B, B₂) scan replacing CS, and an explicit child radius S₂, are needed.
- *Squarefree large part.* The lower bound must subtract values with p² | F(P), p ∈ (B, B₂]: either smooth counts in APs
  mod t·p² (not given by `def:good`) or a union bound Pr_c[p² | F(P)] ≲ 1/p², usable only when 1/B ≪ ρ_F(log N)^{−C}
.
- *Value factorisations* up to B₂ = L[1/2] must be recursive ([rec] Lemmas W/OS; the bounded-depth search recurrence
  of escapeR2 §7.3), and Lemma IR must be re-checked at γ₂.

**Theorem 2L.** Assume U′ and PC₃. Then FACTOR_2L —
FACTOR_U with the descendable set 𝒟_A := {ℓ : N₂(ℓ, A) ≥ 1}, DESCEND_ONCE replaced by the bounded-depth search below,
and POSTPROCESS on two-level rows — is Las Vegas and factors every N in expected time L_N[3/8, c] for a finite c.

*Search.* Scan the level-1 box of ℓ; for each P with Acc_ℓ(P), scan the box of each child; return the first P all of
whose children have a good point, with those points. Given U′'s search sizes, the cost is ≤ W₁(1 + m̄W₂) box points
(m̄ = O(D^{1/2})), each with recursive factorisations of values ≤ N^{o(1)} (escapeR2 §7.3's recurrence): L[3/8, O(1)].

*Proof sketch.* (i) Prop UB and U′ give γ₂ := Pr[N₂ ≥ 1] ≥ (log N)^{−C}. (ii) Theorem K-U and Prop select-U hold for
an arbitrary field-dependent target set (Lemma S″), so with δ_A := #𝒟_A/#J they apply verbatim with γ₂ in place of
γ; split acceptance (γ₂/(16 log N))^{2k}/(32k²) = L[3/8, −O(κ)]. (iii) Rows: nested Lemma RW (each algebraic child is
cancelled by its level-2 element z_{P′}; rational sides are B-smooth; |exponents| ≤ (C log N)^3). (iv) LV6 at d ≍ D^{1/4}:
lv5unif Prop LV6-U (factorial-bit characters, cost e^{O(d log d)} = e^{o(Λ_{3/8})}); character primes exceed B₂.
(v) Fruitfulness: Theorem DD with nested Lemma T_multi (the search is deterministic given canonical factorisations).
(vi) Recursion: halves ≤ √(2N) and values ≤ N^{θ+o(1)}, θ → 0; lv5unif Lemma F_α at α = 3/8 suffices for finiteness,
and escapeR2 NR-q step 8's three-term recursion gives a better constant. ∎

**Dependencies.** Internal: Theorem U's chain (lv5unif §10 list, including [pp] Thm FM/LD*, Lemma 8.5′, Kern-U,
L2-U, Theorem DD, Lemma S″, Lemmas OS/W/F/T_rec, Lemma IR), Prop LV6-U, nested Lemmas RW, TS, T_multi, escapeR2 §7.3
(bounded-depth search recurrence); this note's Prop UB. Hypotheses: PC₃ (open), U′ (open; hybrid BV, mixed CS, S₂,
squarefree part). Published: as Theorem U.



---

# E. Hyp U′: Lemma 8.5H, level-1 moments, CS′, SQF, recursion (source: notes/uprime_results.md §§1-6)

## 1. Setting

**Design 2L′.** Theorem U's design (lv5unif §1) at d := least odd ≥ max(3, d₀D^{1/4}), k := max(⌈κD^{3/8}⌉, d + 1),
ν = 1/(4k), μ := μ_opt = (1 − ν)/d (since k < (d² + 1)/4), η₀ := 1/⌈D^{1/2}⌉, B := ⌊e^{aΛ}⌋, B₂ := ⌊e^{b₂D^{1/2}λ}⌋,
Y = N^{1/(2k)+o(1)} (special-q primes ℓ = L[5/8]), with two changes:
- **(H) height floor**: h := μ + 4ν + η₀. Since μ ≍ D^{−1/4} and dν ≍ D^{−1/8}, θ_F := μ + h + dν = (d₀/(4κ))D^{−1/8}(1 + o(1)) is
  unchanged at leading order, and h ≥ h_*(μ) + η₀ (all three of lv5unif's floors hold: 2ν, μ − (d − 2)ν < 0,
  1 − dμ + ν = 2ν). At level 1 it is free; at level 2 it is a leading-order cost.
- **Acc_ℓ(P)** (multidescent §3.2) with one addition: P ∈ 𝒫 ∩ Λ_ℓ(g), |G(P)| ≥ MV₁/32, G(P)/ℓ ∈ 𝒮_B, and
  F(P) = s·q₁⋯q_m with s B-smooth and q_i distinct primes in (B, B₂] dividing none of b·f_d·disc(f)·**N**. (A child
  q | N is a factor exit through gcd(F(P), N); excluding it from Acc only simplifies the moments.)

Sizes (leading order; `params.sizes`): log Y ∼ D^{5/8}λ/(2κ), log X_r ∼ D^{3/4}λ/d₀ (rational cofactor),
log X_F ∼ (d₀/(4κ))D^{7/8}λ (algebraic value), log H ∼ log M ∼ D^{3/4}λ/d₀, u_r := log X_r/log B ∼ D^{3/8}/(ad₀),
u_{F2} := log X_F/log B₂ ∼ (d₀/(4κb₂))D^{3/8}, log u ∼ (3/8)λ. Child cap m̄ := log X_F/log B ∼ (d₀/(4κa))D^{1/2}.

**What fails if one smoothness parameter is used**, confirmed:
- y = B₂ in Lemma 8.5′ controls the family of moduli ℓs with s B₂-smooth, of density ρ₂(Q/Y)/log Y. Restricted to the
  actual family (s B-smooth) the relative error becomes ε·ρ₂(X_r)/ρ₁(X_r) = e^{−cD^{3/8}/λ²}·e^{Θ(D^{3/8}λ)} → ∞.
- y = B for the values: u = log X_F/log B ≍ D^{1/2}, and (R′3) u log(4u log y) ≤ K₀ log y fails (ratio ≍ D^{1/8}).
  Independently, B-smooth values at size X_F cost L[1/2] (that is why the children exist).

## 2. Lemma 8.5H: the hybrid restricted BV [proved]

Notation of paper2 Appendix A: E(z, y; q, a) := Ψ(z, y; q, a) − Ψ_q(z, y)/φ(q); Δ_{z,y}(q) := φ(q)^{−1}Σ_{χ≠χ₀}|Ψ(z, y; χ)|;
H(u) := exp(u/log²(u + 1)); α(z, y) the Hildebrand–Tenenbaum saddle point; δ_D, η_D, c_D, C_D the constants of
Drappeau's Lemme 5 (Fact D), C_dBT that of Fact dBT, C_α that of (α2). E1–E8, Lemma LS and Corollary LS (= lv5unif's
LS′) are the regime-free lemmas of paper2 Appendix A, used as proved there.

**Statement.** Let 2 ≤ y₁ ≤ y₂ ≤ x, and Y, Q > 0. Put u₂ := log x/log y₂ and, for z ≤ x, u₁(z) := log z/log y₁. Let
J_Y := {ℓ prime : Y < ℓ ≤ 2Y}, 𝒯_Q := {t = ℓs ∈ (Q/2, Q] : ℓ ∈ J_Y, P⁺(s) ≤ y₁}, δ_𝒯 := ρ₁(Q/Y)/log Y. Fix K₀ ≥ 1,
K := 2K₀ + 6, and U₀ ≥ 8 such that Hildebrand–Tenenbaum's Theorem 3 error O(1/u + log y/y) is ≤ 1/100 for
u ≥ U₀ − 1. Assume, for x ≥ x₀(K₀):
- (H1) log y₁ ≥ max{4(log log x)², (8K(1 + C_α)/δ_D)·log log x}, and (log x)^{c_D} ≤ y₂;
- (H2) u₂ ≥ max{(log log x)⁴, 2K/η_D};
- (H3) u₂ log(4u₂ log y₂) ≤ K₀ log y₂ and u₁(Q) log(4u₁(Q) log y₁) ≤ K₀ log y₂;
- (H4) K(log y₂/log y₁)(log log x + C_α log log y₁) + log(60C_dBT²C_D log log(x + 16)) ≤ δ_D u₂/(8 log²(u₂ + 1));
- (H5) Y > y₂, Y ≥ y₂^K, and Γ := min{log(Q/Y), log(√x/Q)} ≥ (2K₀ + U₀)log y₂ + 4 log log x.

Then, with C, c > 0 depending only on K₀ and the absolute constants above:
- (i) Qδ_𝒯/10 ≤ #𝒯_Q ≤ 3Qδ_𝒯;
- (ii) Σ_{t∈𝒯_Q} max_{(r,t)=1}|E(x, y₂; t, r)| ≤ C·Ψ(x, y₂)·δ_𝒯·ε₂, ε₂ := exp(−cu₂/log²(u₂ + 1)) + y₂^{−1/2}.

At y₁ = y₂ this is Lemma 8.5′(R′) up to the value of c (δ_D·min(3/8, 3/(4K₀)) here, δ_D/(4 + 2K₀) there); (H1)'s
second clause and (H4) are then implied by (R′1)–(R′3) for x ≥ x₀. Besides (H4), the mixed density
floor (H3b) is the other genuinely two-parameter hypothesis.

**Proof.** We follow lv5unif §5 (Steps A–E), marking every change.

*Step A (saddle points), now at both parameters.* For 2y_i ≤ z ≤ x, (α2) gives 1 − α(z, y_i) ≤
(log(u log y_i) + C_α log log y_i)/log y_i with u log y_i = log z ≤ log x, so by (H1) 1 − α(z, y_i) ≤ (1 + C_α)/(4 log log x)
≤ 1/100 for i = 1, 2 and x ≥ x₀ (y₂ ≥ y₁). In particular α(z, y_i) ≥ 3/4. Moreover, by the second part of (H1),
  K(1 − α(z, y₁)) ≤ K(1 + C_α)log log x/log y₁ ≤ δ_D/8.   (A′)

*Step B (density floors), now separately.* E7 at y₂ and z = x: ρ₂(x) ≥ y₂^{−1}exp(−u₂ log(4u₂ log y₂)) ≥ y₂^{−1−K₀}
by (H3). E7 at y₁ and 2y₁ ≤ z ≤ Q: ρ₁(z) ≥ y₁^{−1}exp(−u₁(z)log(4u₁(z)log y₁)) ≥ y₂^{−1−K₀}, because
u ↦ u log(4u log y₁) is increasing, u₁(z) ≤ u₁(Q), and y₁ ≤ y₂. Hence for 2y₁ ≤ Q′ ≤ Q:
  Ψ(x, y₂)ρ₁(Q′) ≥ x·y₂^{−2−2K₀}.   (B′)

*Step C (Lemma SM-H: y₁-smooth moduli, y₂-smooth values, relative saving).* Let 2y₁ ≤ Q′ ≤ min{Q, √x·y₂^{−(2K₀+U₀)}(log x)^{−4}}
and 𝒮₁ := {s ∈ (Q′/2, Q′] : P⁺(s) ≤ y₁}. Claim: Σ_{s∈𝒮₁}Δ_{x,y₂}(s) ≤ C·Ψ(x, y₂)ρ₁(Q′)·ε₂.
By E1 this bounds Σ_s max_a|E(x, y₂; s, a)|. Split the characters by conductor r: (a) 1 < r ≤ R₁ := y₂^K; (b) r > R₁.
- (a) **[the only changed step].** E3 is applied with the value parameter y₂: the Möbius variable d runs over divisors
  of s_{y₂} = ∏_{p|s, p≤y₂}p = rad(s), because s is y₁-smooth and y₁ ≤ y₂; so d is squarefree, y₁-smooth, d ≤ Q′,
  (d, r) = 1. Exchanging sums, the modulus weights are Σ_{s∈𝒮₁, rd|s}1/φ(s), and E8 **with the modulus parameter y₁**
  gives ≤ 10C_dBT log log(Q′ + 16)·ρ₁(Q′)·(rd)^{1−β₁}/(φ(r)φ(d)), β₁ := α(Q′, y₁). (If rd is not y₁-smooth the
  sum is empty and the bound is trivial.) Bound r^{1−β₁} ≤ R₁^{1−β₁}. The d-sum is
  Σ_d d^{1−β₁}φ(d)^{−1}|Ψ(x/d, y₂; χ*)| summed over primitive χ*; by Fact D at (x/d, y₂, Q := R₁, A = 1/2) — its
  hypotheses hold since x/d ≥ √x, (log x)^{c_D} ≤ y₂ by (H1), R₁ = y₂^K ≤ x^{η_D/2} ≤ (x/d)^{η_D} by (H2) — and Fact dBT
  at y₂ (Ψ(x/d, y₂) ≤ C_dBTΨ(x, y₂)d^{−β₂}, β₂ := α(x, y₂)), and E6 (Σ_d d^{1−β₁−β₂}/φ(d) < 6, as β₁, β₂ ≥ 3/4):
    (a) ≤ 60C_dBT²C_D log log(x + 16)·R₁^{1−β₁}·{H(u₂/2)^{−δ_D} + y₂^{−δ_D}}·Ψ(x, y₂)ρ₁(Q′).
  *The transfer factor.* By (α2) at (Q′, y₁), log R₁^{1−β₁} = K log y₂·(1 − β₁) ≤ K(log y₂/log y₁)(log log x +
  C_α log log y₁). In Lemma 8.5′ (y₁ = y₂) this was (log x)^{O(1)}; here it is exp(O((log y₂/log y₁)log log x)).
  By (H4) and δ_D(u₂/2)/log²(u₂/2 + 1) ≥ δ_Du₂/(2 log²(u₂ + 1)), the H-term contributes ≤ exp(−3δ_Du₂/(8 log²(u₂+1))).
  By (A′), R₁^{1−β₁} ≤ y₂^{δ_D/8}, so the y₂-term contributes ≤ y₂^{−3δ_D/4} (the log log prefactor is ≤ y₂^{δ_D/8} for
  x ≥ x₀ by (H1)), and by (H3) (u₂ log u₂ ≤ K₀ log y₂) y₂^{−3δ_D/4} ≤ exp(−3δ_Du₂/(4K₀ log²(u₂ + 1))). So
    (a) ≤ C·Ψ(x, y₂)ρ₁(Q′)·exp(−cu₂/log²(u₂ + 1)), c := δ_D·min(3/8, 3/(4K₀)).
- (b) If Q′ ≤ R₁ the case is empty. Otherwise, unchanged in form: Corollary LS with value parameter y₂, z₀ = x, Q₀ = Q′, R₀ = R₁, 𝒳(s) = {χ mod s : cond χ > R₁}
  for s ∈ 𝒮₁ (any subset of moduli is allowed). Its hypotheses hold since x/Q′ ≥ √x ≥ max(4, y₂^{3/2}). It gives
  (b) ≤ 1.5·10³π(y₂)(log x)²(y₂^{3/2}x/R₁ + y₂^{1/2}x^{1/2}Q′). Relative to (B′) this is
  ≤ 1.5·10³(log x)²(y₂^{9/2+2K₀−K} + y₂^{7/2+2K₀}Q′/√x) ≤ 1.5·10³(log x)²(y₂^{−3/2} + y₂^{7/2−U₀}(log x)^{−4}) ≤ y₂^{−1/2}.

*Step D (count of 𝒯_Q), at y₁.* As lv5unif Step D with y = y₁: ℓ > y₂ ≥ y₁ ≥ P⁺(s), so t determines (ℓ, s);
z := Q/(2ℓ) ≥ e^Γ/4 ≥ y₂^{U₀−1} ≥ y₁^{U₀−1}, α(z, y₁) ≥ 0.99 (Step A), and HT Theorem 3 at y₁ gives (i).

*Step E (the error sum).* E2 (lift averaging) needs ℓ > y₂ (every y₂-smooth n is coprime to ℓ) — this is (H5).
- Part 1 (χ_ℓ principal): Σ_ℓ(ℓ − 1)^{−1}Σ_{s∈(Q/2ℓ,Q/ℓ], P⁺(s)≤y₁}max_a|E(x, y₂; s, a)|. Step C applies with Q′ = Q/ℓ
  (2y₁ ≤ Q/(2Y) by Γ; Q/ℓ ≤ √x e^{−Γ}). Since Ψ(·, y₁) is nondecreasing, ρ₁(Q/ℓ) ≤ 2ρ₁(Q/Y); Σ_{ℓ∈J_Y}1/(ℓ − 1) ≤ 2.6/log Y.
  So Part 1 ≤ 5.2C·Ψ(x, y₂)δ_𝒯ε₂.
- Part 2 (ℓ | cond χ, so cond χ > Y): Corollary LS with z₀ = x, Q₀ = Q, R₀ = Y gives
  ≤ 1.5·10³π(y₂)(log x)²(y₂^{3/2}x/Y + y₂^{1/2}x^{1/2}Q). By (B′) with Q′ = Q/Y, Ψ(x, y₂)δ_𝒯 ≥ x y₂^{−2−2K₀}/log Y, and
  by (H5) Part 2 ≤ Ψ(x, y₂)δ_𝒯·y₂^{−1/2}.
- The principal character of t contributes the main term Ψ_t/φ(t) exactly. ∎

**Corollary bad₂.** Call t *good₂* if for every x in the grid 𝒳 (§3) and every (r, t) = 1,
|Ψ(x, B₂; t, r) − Ψ_t(x, B₂)/φ(t)| ≤ Ψ_t(x, B₂)/(2φ(t)). In the setting of Lemma 8.5H with x ranging over 𝒳 and
y₁ = B, y₂ = B₂, #{t ∈ 𝒯_Q bad₂} ≤ 30C|𝒳|ε₂·#𝒯_Q.
*Proof.* As paper2 Corollary `cor:bad`: Tenenbaum's Ψ_q theorem (Fact Ten) at y = B₂ — its range ω(t) ≤ log N ≤
B₂^{c/log(1+u₂)} holds and its error O(log(1 + u₂)log(1 + ω)/log B₂) + O(ω/B₂) → 0 — gives
Ψ_t(x, B₂)/(2φ(t)) ≥ Ψ(x, B₂)/(3t) ≥ Ψ(x, B₂)/(3Q); Markov and (i). ∎

**Where the design sits** (`params.margins`, all hold for every (a, d₀, κ, b₂)). x ∈ 𝒳 has log x = θ_F log N(1 + o(1)),
u₂ ∼ (d₀/(4κb₂))D^{3/8}; Q ∈ [MV₁/32, 2MV₁], log(Q/Y) ∼ log X_r ∼ D^{3/4}λ/d₀. (H1): log B = aD^{3/8}λ ≫ λ². (H2):
D^{3/8} ≫ λ⁴. (H3): u₂ log u₂ ≍ D^{3/8}λ ≪ log B₂ ≍ D^{1/2}λ, and u_r log u_r ≍ D^{3/8}λ ≪ D^{1/2}λ. **(H4):** LHS
≍ (b₂/a)D^{1/8}λ, RHS ≍ D^{3/8}/λ² — margin D^{1/4}/λ³. (H5): log Y ≍ D^{5/8}λ ≫ log B₂; Γ ≥ min(log X_r,
(1/2)θ_F log N − log Q) ≍ D^{3/4}λ ≫ log B₂ (the level gap is automatic at 3/8: it does not use η₀). The saving
ε₂ = exp(−Θ(D^{3/8}/λ²)) beats every polylogarithm, as FM needs. *General form:* with y₁ = L_x[β₁], y₂ = L_x[β₂]
relative to the value size x, (H4) requires β₂ − β₁ < 1 − β₂; at 3/8, β₁ = 3/7, β₂ = 4/7.

## 3. Level 1 with two smoothness parameters [proved]

The grid is lv5unif's: ϖ := B²S₁², 𝒳 := {X_Fϖ^{−2}(3/2)^i} ∩ [X_Fϖ^{−2}, 16dX_F], X_F := (2d + 1)MHV₁^d, V₁ = ⌊√(2Y)S₁⌋,
X_r := 2MV₁/Y. ρ_r := ρ₁(X_r), ρ_{F2} := ρ₂(X_F).

**Lemma U1₂ (local behaviour at y₂) [proved].** For z₁ ≤ N with log z₁/log B₂ ≥ D^{1/4}/2 and z₁ ≤ z₂ ≤ z₁B₂^C:
ρ₂(z₂) ≤ 2ρ₂(z₁) and ρ₂(z₁) ≤ (log N)^{C_U(C+1)}ρ₂(z₂). *Proof.* HT Theorem 3 in ≤ C + 1 steps with c ≤ y = B₂;
(1 − α(z, B₂))log B₂ ≤ 2 log log N. ∎ Every shift used below is by a factor ≤ ϖ² = B⁴S₁⁴ ≤ B₂ or ≤ B₂² — one or two
steps — and in fact (1 − α₂)log ϖ² = O(D^{−1/8}λ) → 0, so these shifts cost 1 + o(1) at y₂.

**Lemma GM-U′ [proved].** Hypotheses of lv5unif GM-U, with t good₂. Then Pr_c[F_A(P) ∈ 𝒮_{B₂}] ≥ ρ_{F2}/(log N)^{C_g}.
*Proof.* paper2 `lem:goodmod` / lv5unif GM-U line by line with y = B₂ for the values: the small-layer bound, the
envelope (|F̂(P)| ≤ MHV^dN^{−η₀/4}; Lemma `lem:valuesizes` needs h ≥ 1 − dμ + ν + η₀, true under (H)), the choice of
x₁ ∈ 𝒳 ∩ (ty′/4, ty′/2], goodness₂ and Fact Ten at y = B₂ give Ψ(x₁; t, r′) ≥ Ψ(x₁, B₂)/(3t), and U1₂ gives
ρ₂(x₁) ≥ ρ_{F2}/(log N)^{2C_U}. The two error terms 33eϖN^{−η₀/4} and 3d²V/(2H + 1) are ≤ ρ_{F2}B^{−1}, since
ρ_{F2} = L[3/8]^{−O(1)} while N^{−η₀/4} = L[1/2]^{−Θ(1)} and V/H ≤ N^{−μ}. ∎

**Lemma SQF (local exclusions) [proved].** Let (ℓ, P, g, w) satisfy the hypotheses of GM-U′ (so t = |G(P)| = ℓn,
n B-smooth, gcd(t, bg₁Nw) = 1; b prime > (1 − 1/d)V₁ > B₂), and let c be uniform on C_H.
- (a) For a prime p ∈ (B, B₂]: Pr_c[p² | F_A(P)] ≤ 1/p² + 1/(2H + 1).
- (b) For a prime q ∈ (B, B₂], q ∤ N: Pr_c[q | F_A(P) and q | f_d·disc(f)] ≤ 2e(d − 1)/q².
- (c) For a prime q ∈ (B, B₂], q | N: Pr_c[q | F_A(P)] ≤ 1/q + 1/(2H + 1).
Hence Pr_c[F_A(P) ∈ 𝒮_{B₂} ∖ 𝒜cc] ≤ (log N)²/B for N ≥ N₀, where 𝒜cc is the algebraic part of Acc.

*Proof.* (a) p ∤ G(P) (its prime factors are ≤ B or equal to ℓ > B₂). F_A(P) = F̂(P) + G(P)⟨c, v(P)⟩, so p² | F_A(P)
fixes ⟨c, v(P)⟩ modulo p². The coefficient of c_{d−1} in ⟨c, v(P)⟩ is b^{d−1}, a unit mod p. For fixed c₀, …, c_{d−2},
c_{d−1} ∈ [−H, H] lies in one class mod p²: ≤ (2H + 1)/p² + 1 values. (c) The same with q in place of p².
(b) Since q ≤ B₂ ≤ (2H + 1)/d, every class of C_H mod q has between ⌊(2H+1)/q⌋^d and ⌈(2H+1)/q⌉^d elements, so the law
of c mod q is ≤ (1 + q/(2H + 1))^d ≤ e times uniform. Under the uniform law, F_A mod q is uniform on
𝒜_q := {F ∈ F_q[x, y]_d : F(g₀, g₁) = Nw}: forms vanishing at the point (g₀ : g₁) are exactly G·F_q[x, y]_{d−1}
(G ≢ 0 since gcd(g₀, g₁) = 1), and Nw ≢ 0 (q ∤ N, q ≠ w). P ≢ λ(g₀, g₁) mod q since q ∤ G(P); so F ↦ F(P) is a
nonconstant affine function on 𝒜_q and Pr[F(P) = 0] = 1/q. Given F(P) = 0, F = L_P·H with L_P := bx − ay and H uniform
on the affine hyperplane ℋ := {H ∈ F_q[x, y]_{d−1} : H(g₀, g₁) = κ₁}, κ₁ := Nw/L_P(g₀, g₁) ≠ 0 (L_P(g₀, g₁) = −G(P)).
Now f_d = F(1, 0) = b·H(1, 0), and if f_d ≢ 0 then q | disc(f) iff the binary form F has a repeated linear factor over
F̄_q, i.e. iff L_P | H or H has a repeated factor. So the event is contained in {H(1, 0) = 0} ∪ {H(P) = 0} ∪
{Disc(H) = 0}. The first two are nonconstant affine conditions on ℋ (or empty, when (1 : 0) or P equals (g₀ : g₁)
projectively), probability ≤ 1/q each. Disc(H) is a polynomial of degree 2(d − 2) in the coefficients, not identically
zero on ℋ when q > d (take H = κ∏_{i<d}(x − r_iy) with distinct r_i ∈ F_q ∖ {g₀/g₁} and κ fixing H(g₀, g₁) = κ₁), so
Schwartz–Zippel gives ≤ 2(d − 2)/q. Total ≤ e·(1/q)·(2(d − 1)/q). ∎ (Sum: Σ_{p>B}(1 + 2e(d−1))/p² + (log N/log B)(1/B) +
2π(B₂)/(2H + 1) ≤ (log N)²/B.)

**Theorem FM-U′ [proved].** μ₁ = E_{ℓ,A}N₁(ℓ, A) ≥ S₁²ρ_rρ_{F2}/(d(log N)^{C_FM}), provided (SQ) holds.
*Proof.* lv5unif FM-U / paper2 `thm:FM` with: (1) for each counted (ℓ, P, n, g) with t good₂ and every w,
Pr_c[F_A(P) ∈ 𝒜cc] ≥ ρ_{F2}/(log N)^{C_g} − (log N)²/B ≥ ρ_{F2}/(2(log N)^{C_g}) (GM-U′, Lemma SQF), where the last step
is **(SQ)**: log(1/ρ_{F2}) = u_{F2}log u_{F2}(1 + o(1)) = (3d₀/(32κb₂))Λ(1 + o(1)) < aΛ = log B. (2) Good rational mass
unchanged (B-smooth n, Lemma FIB, Lemma U1 at B). (3) Bad mass: the moduli t = ℓ|n| lie in 𝒯_Q for the 7 dyadic Q
covering [MV₁/32, 2MV₁], with y₁ = B; Corollary bad₂ gives π_J^{−1}Σ_ℓ𝔟_ℓ ≤ C|𝒳|ε₂(MV₁/Y)ρ_r(log N)^C, which is
(log N)^{−C′} times the good mass for every C′, since |𝒳| = O(Λ) and ε₂ = e^{−Θ(D^{3/8}/λ²)}. ∎

**If (SQ) fails** (a·b₂ < 3d₀/(32κ) strictly), the subtraction (log N)²/B exceeds ρ_{F2} by e^{Θ(Λ)} and the union
bound gives nothing; that is a failure of the elementary estimate, not evidence that acceptable values are scarce
. The alternative, B₂-smooth counts in APs mod t·p² for p ∈ (B, B₂], would need moduli t·p² — still in a
hybrid family (p² is y₂-smooth but not y₁-smooth), so it would need a three-parameter version; not pursued, as (SQ)
costs only a constraint on constants.

**Lemma DS-U′ [proved].** For fixed C₀, X ≤ N with log X/log B₂ ≥ D^{1/4}, 1 ≤ q ≤ B^{C₀}, an interval I ⊂ [XB^{−C₀}, X] of
length A ≥ XB^{−C₀}: #{n ∈ I : n ≡ r (q), P⁺(n) ≤ B₂} ≤ (A/q)ρ₂(X)(log N)^{C}.
*Proof.* paper2 `lem:DS` with divisors d₁ ∈ (D′/B₂, D′] (y = B₂); the dyadic sum has ≤ log₂(2B₂) ≤ log N blocks;
X/(D′/B₂) ≤ B₂B^{2C₀} ≤ B₂², so U1₂ with C = 2 applies. ∎

**Theorem L2-U′ [proved].** paper2 `thm:L2` / lv5unif L2-U with 𝒮* := 𝒮_{B₂} ∩ [−X_F, X_F] ∖ {0} and the rational
family 𝒩_ℓ (B-smooth n) unchanged:
  E_ℓΣ_{n∈𝒩_ℓ}𝔔(ℓ²|n|q)/φ(ℓ|n|) ≤ τ(q)(log N)^{C_L}Ψ(X_r, B)Ψ(X_F, B₂)²/(Yq(MV₁)²), uniformly in 1 ≤ q ≤ 4S₁².
*Proof.* R₀ := B^{C₀} ≥ 64S₁²B/(ρ_rρ_{F2}) with C₀ ≤ 2σ₁/a + 3 fixed. Class (A): N_c ≤ (X_F/m′)ρ_{F2}(log N)^C by DS-U′
(and R₀^{−1} ≤ ρ_{F2}); multiplicities mult(m′) ≤ (Ψ(X_r, B)/m₁)(log N)^{O(C₀)} by Lemma U1 at B (unchanged).
Class (B): large sieve with Σ|a_F|² = #𝒮* ≤ 2Ψ(X_F, B₂); ratios q/(R₀ρ_rρ_{F2}) ≤ 1 and X_rq²/(X_Fρ_rρ_{F2}) =
N^{−h−(d+1)ν}B^{O(1)} = o(1). Class (C): ratios log N/(π_Jρ_{F2}) and X_rq²Y⁴/(π_JX_Fρ_rρ_{F2}) = N^{(5−d)ν−h}B^{O(1)},
both o(1). ∎

**Theorem LD*-U′ [proved].** E_{ℓ,A}[N₁(N₁ − 1)] ≤ d·η₀^{−1}(log N)^{C_LD}S₁⁴ρ_r²ρ_{F2}².
*Proof.* An upper bound may drop the squarefree, discriminant and q ∤ N conditions: Acc ⊂ {rational filter, F(P) ∈ 𝒮_{B₂}}.
Then lv5unif LD*-U Steps 0–7 verbatim (Kern-U with bb′ ≤ V₁² ≤ (2H + 1)/d; W-U; Step 5 with DS-U at B; the pair bound
q ≤ 4S₁²), with L2-U′ in Step 7. The power count V⁴/Y² ≤ 4S₁⁴ and the constant (2d + 1)²/d ≤ 9d are unchanged. ∎

**Theorem D-U′ [proved].** γ₁ := Pr[N₁ ≥ 1] ≥ (log N)^{−C_γ}, and E[N₁²] ≤ (log N)^C(μ₁² + μ₁).
*Proof.* Paley–Zygmund with Lemma CS′ (§4.1): m₁ ≥ (log N)^{−2C_R}/4 and m₁ ≤ μ₁d(log N)^{C_FM}. ∎

These are exactly the level-1 moments Prop UB consumes (μ₁ and E[N₁²]; E[N₂²] ≤ E[N₁²]).

## 4. Box sizes

### 4.1 Lemma CS′ (mixed scan) [proved]

For y ∈ {B, B₂} and z with log z/log y ≥ D^{1/4}, let α̂_y(z) solve Σ_{p≤y}log p/(p^α − 1) = log z and
R_y(z) := z^{α̂−1}∏_{p≤y}(1 − p^{−α̂})^{−1}; R*_y(z) the exact value. As in Lemma CS, ρ_y(z) ≤ R*_y(z) ≤ (log N)^{C_R}ρ_y(z)
(Rankin above; Hildebrand–Tenenbaum Theorem 1 at the saddle point below, φ₂ ≍ log z log y).

*Computation without primes up to B₂.* Enumerating p ≤ B₂ costs L[1/2]. Instead: sum exactly over p ≤ B, and replace
the sums over p ∈ (B, B₂] of log p/(p^α − 1) and −log(1 − p^{−α}) by the integrals ∫_B^{B₂}dt/(t^α − 1) and
∫_B^{B₂}−log(1 − t^{−α})dt/log t. With a de la Vallée Poussin-type explicit error |θ(t) − t| ≤ C₁t exp(−c₂√log t)
(Trudgian, Ramanujan J. 39 (2016) 225–234, Thm 1: |θ(t) − t| ≤ t·√(8/(17π))X^{1/2}e^{−X}, X = √(log t/6.455), t ≥ 149;
hence |θ(t) − t| ≤ 0.5·t·exp(−0.35√log t), our weakening), integration by parts bounds both errors by
B₂^{1−α}(log B₂)²exp(−c₂√log B) = exp(O(λ) − c₂D^{3/16}√(aλ)) → 0, uniformly for α in the window
[1 − Cλ/log B₂, 1] that contains the saddle point. The approximate
saddle sum is decreasing in α, so bisection locates its root; the Rankin logarithm has zero derivative at the saddle
point, so an O((log N)^{−3}) error in α̂ costs O(1) in log R. Hence R_{B₂} is obtained to relative accuracy 2 at cost
B^{1+o(1)}. Certified quadrature to
(log N)^{O(1)} digits is routine and not written out. (φ₂ ≍ log z·log y is HT Theorem 2(ii).)

**Lemma CS′.** The scan S = 2^j, j = 0, 1, …, stopping at the first S with S²R_B(X_r(S))R_{B₂}(X_F(S)) ≥ 1, stops at some
S₁ with S₁² ≤ e^{(σ₁+o(1))Λ}, σ₁ := 3/(8ad₀) + 3d₀/(32κb₂), and m₁ := S₁²ρ_rρ_{F2} ∈ [(log N)^{−2C_R}/4, 64].
*Proof.* Lemma CS verbatim with the two Rankin factors at different y. Existence: for S² ≤ e^{2σ₁Λ},
d log S = O(D^{5/8}λ) = o(log X_F), so log X_r and log X_F keep their leading terms, and CEP (valid: u ≤ log z/(2 log log z))
gives log(R_BR_{B₂}) = −(u_r log u_r + u_{F2}log u_{F2})(1 + o(1)) = −σ₁Λ(1 + o(1)). Overshoot: X(S/2) ≤ X(S) ≤ 2^{d+1}X(S/2)
≤ BX(S/2), and Lemma U1(i) at B and U1₂(i) at B₂ give m₁ ≤ 64. ∎

*Single-y comparison.* Lemma CS with both factors at B gives log(1/(ρ_B(X_r)ρ_B(X_F))) ∼ (d₀/(8κa))D^{1/2}λ (u ∼ D^{1/2},
log u ∼ λ/2): an L[1/2]-scale box (`params.sigma1_single_y`), as multidescent said.

### 4.2 Child boxes and S₂ [proved]

For a child 𝔮 = (q, r) of P (r = a/b mod q), Λ_𝔮 := {P′ : a′ ≡ rb′ (mod q)}. Two geometries:
- **Conservative.** Box(𝔮) := primitive points, mod ±1, of {ie₁ + je₂ : |i|, |j| ≤ S₂} for a
  Gauss-reduced basis e₁, e₂ of Λ_𝔮. There are ≥ cS₂² box points. Envelopes: X_{r2} := 2MV₂, X_{F2} := (2d + 1)MHV₂^d/Q₂ ≥ |F(P′)|/q.
- *Square (sharper; needs CM's label-law bound).* Λ_𝔮 ∩ [−V₂, V₂]² with V₂ := ⌊√(2Q₂)S₂⌋. Regular lattices
  (λ₁ ≥ √q/(2S₂)) have ≍ S₂² points; irregular ones have ≤ 2 primitive points (Lemma `lem:reg`), so such children are
  dead. The number of residues r mod q with an irregular lattice is ≤ #{primitive v : |v|² < q/(4S₂²)} ≤ q/S₂², i.e. a
  fraction ≤ 1/S₂² of labels. Since the incidence weight is bounded pointwise in (ℓ, P, q) (CM1) and P is
  equidistributed mod q (V₁ ≫ q), this is also a negligible fraction of incidences; stated as part of CM.

**Scan S₂(Q₂)** (conservative geometry): the first 2^j with S²R_B(X_{r2}(S))R_B(X_{F2}(S)) ≥ 4(log N)^{2C_R}·
exp((log log N)²). Then m₂(Q₂) := S₂²ρ₁(X_{r2})ρ₁(X_{F2}) ≥ exp((log log N)²) and S₂² = e^{(σ₂(Q₂)+o(1))Λ} with
  σ₂(Q₂) = (3/(8a))(3/d₀ + d₀b(Q₂)) ≤ σ₂ := (3/(8a))(3/d₀ + d₀b₂),  b(Q₂) := log(2Q₂)/(D^{1/2}λ),
using log X_{r2} ∼ D^{3/4}λ/d₀ and, under (H), log X_{F2} ∼ (2/d₀ + d₀b(Q₂))D^{3/4}λ (square geometry: d₀b/2). One scan
per dyadic window that actually occurs (≤ log₂B₂ of them), each B^{1+o(1)}. The threshold exp((log log N)²) is
super-polylogarithmic but L^{o(1)}: it is chosen so that m₂ beats the unknown polylogarithmic loss of CM whatever its
exponent (the algorithm cannot know C_CM).

**Search sizes.** After paper2's regularity precheck, W₁ := #(Λ_ℓ(g) ∩ [−V₁, V₁]²) ≤ 400S₁² = L[3/8, σ₁]; W₂ := max (2S₂ + 1)² = L[3/8, σ₂].

**What this does and does not give.** m₂ is the *nominal* (product) main term of a child box. That the conditional
mean E[N_𝔮 | E] is ≥ m₂/(log N)^C is Hyp CM (§7) — not a consequence of anything in this section.

## 5. Squarefree large part: constants

Item 3 of multidescent's U′ asked for the constants of the union-bound option. They are in Lemma SQF and FM-U′:
- per prime p ∈ (B, B₂]: Pr_c[p² | F(P)] ≤ p^{−2} + (2H + 1)^{−1} (needs B₂² ≤ (2H + 1)/d; true, log H ≍ D^{3/4}λ);
- discriminant/leading coefficient: ≤ 2e(d − 1)q^{−2};
- children dividing N: ≤ (log N/log B)(1/B + 1/(2H + 1));
- total ≤ (log N)²/B; usable iff **a > 3d₀/(32κb₂)**, i.e. log B beats log(1/ρ_{F2}) at leading order. Since σ₁ contains
  3d₀/(32κb₂) as well, (SQ) says: the factor-base exponent a must exceed the algebraic half of σ₁.
- (SQ) is satisfiable with room (e.g. a = 3/5, d₀ = 1, κ = 1/8, b₂ = 2: 3d₀/(32κb₂) = 3/8 < 3/5). Increasing b₂ relaxes
  (SQ) and lowers σ₁, but raises σ₂ (larger children); the cost optimum (`params.optimise_2L`) balances these.

## 6. Recursion, search and Lemma IR [proved]

- **Value factorisations.** Batch smoothness (Bernstein) to B decides the rational side and the B-part of F(P); the
  cofactor c₀ := |F(P)|/(B-smooth part), 1 ≤ c₀ ≤ X_F, is factored by FACTOR_U (Theorem U) with fresh coins, and the
  canonical factorisation is checked against Acc. No sieving or trial division to B₂. Expected cost per call
  ≤ L_{X_F}[2/5, c_U + o(1)] = exp(O((θ_F log N)^{2/5}λ^{3/5})) = exp(O(D^{7/20}λ)) = L[3/8]^{o(1)}
  (7/20 < 3/8).
- **Lemma OS** (canonical outputs ⇒ the outer transcript's law is independent of the oracle) and **Lemma W** (Wald for
  adaptive fresh-coin sub-calls) apply verbatim; hence (T) and Theorem DD's fruitfulness (multidescent Theorem 2L (v)).
- **Search (escapeR2 §7.3 recurrence at depth 2).** For each of the ≤ W₁ level-1 candidates: batch work plus ≤ 1
  FACTOR_U call; for each accepted point and each of its ≤ m̄ children: a child scan (B^{1+o(1)}, cached per window)
  and ≤ W₂ candidates with batch work. Deterministic count ≤ W₁(1 + m̄W₂) candidate examinations; expected cost
  T_SQ ≤ W₁(1 + m̄W₂)·L^{o(1)} = L[3/8, σ₁ + σ₂ + o(1)], uniformly in the field and the target, including failed searches.
- **Lemma IR at γ₂.** With (H), log H/(13d log log H) ∼ D^{1/2}/(13d₀²·(3/4)) → the rejection probability
  exp(−Θ(D^{1/2})) + 3/(2H + 1) + 1/W is o(γ₂) for γ₂ ≥ (log N)^{−C}. This is used *globally* (as in
  select-U: the set 𝒢* = {A ∉ 𝓡 : δ_A ≥ γ₂/4} has measure ≥ γ₂/4), not conditionally on a rare incidence E. Lemma IR's
  primes r ≤ (1/2)log H are < B, and Q_IR = H^{1/2+o(1)} ≤ (2H + 1)/d.
- **Cost.** Lemma F_α at α = 3/8 with (10.1) of lv5unif, T_SQ as above: c₂L = max(2a, a + s + max(a, σ₁ + σ₂),
  (a + s)/(1 − 2^{−3/8})), s = 2κ(1 + C_{γ₂}); 1/(1 − 2^{−3/8}) ∈ [4.36881, 4.36883] (exact bracket via x⁸ vs 8).
  Illustration (rational grid, (SQ) imposed; C_{γ₂} not computable): c ≤ 4.23, 4.61, 5.35 for C_{γ₂} = 0, 1, 5.



---

# F. Window design, Lemma η, Prop PB, Lemma VD, DW lemmas, Lemma TP (source: notes/c12alg_results.md §§1-3.2)

## 1. The window design 2L_DW(α) and the precision budget

### 1.1 Parameters

Fix rationals α ∈ (3/8, 2/5) and a, κ, d₀, b₂ > 0.
- **Degrees and special-q.** d := least odd ≥ max(3, d₀D^{1−2α}); k := max(⌈κD^α⌉, d + 1); ν := 1/(4k); μ := (1 − ν)/d
  (μ_opt, since k < (d² + 1)/4); Y, J and ℓ ≈ N^{1/(2k)} = L[1 − α] as in lv5unif.
- **Smoothness bounds.** B := ⌊e^{aD^αλ}⌋; children are the primes in (B, B₂], B₂ := ⌊e^{b₂D^{4α−1}λ}⌋.
- **Height margin.** η₀ := 1/⌈D^{3/2−3α}⌉, i.e. η₀ log N ≍ D^{3α−1/2}λ. This lies strictly between max(α, 4α − 1) and
  2α, so N^{η₀} beats every L[α]^{O(1)} and every B₂^{O(1)}, and η₀ = o(μ).
- **Field law DW.** g uniform on primitive points of [M/2, M]², M := ⌊N^μ⌋; w uniform on the integers of (W_m, 2W_m],
  W_m := ⌈(2MV₁)²N^{η₀}⌉; c uniform on C_H, H := ⌊2^{d+8}d⁴W_mN^{η₀}⌋ (so h = 2μ + 2ν + 2η₀ + o(η₀)); F := wF̂_{g,N} + G·R_c.
- **Level 1.** Acc_ℓ(P) as uprime §1 (B-smooth G(P)/ℓ, filter |G(P)| ≥ MV₁/32, F(P) = s·q₁⋯q_m with s B-smooth and
  distinct q_i ∈ (B, B₂] not dividing b f_d disc(f) N). Box S₁ from Lemma CS′ with DW's envelope.
- **Level 2.** Conservative child boxes Box(𝔮) (uprime §4.2) with S₂ from the child scan run with DW's envelope
  X_{F2} = (2d + 1)MHV₂^d/Q₂. Then:
  - **Good(P′)** :⇔ G(P′) ∈ 𝒮_B, the child filter |G(P′)| ≥ M|P′|_∞N^{−η₀/4}, and F(P′)/q ∈ 𝒮_B (then q ∥ F(P′),
    since q > B).

**Sizes (leading order; `params.sizes`).**
- Logarithms of the main quantities:
  - log Y ∼ D^{1−α}λ/(2κ);
  - log M ∼ log X_r ∼ D^{2α}λ/d₀;
  - log X_F ∼ (d₀/(4κ))D^{2−3α}λ (the d·log V₁ term dominates, as 2 − 3α > 2α);
  - u₂ := log X_F/log B₂ ∼ (d₀/(4κb₂))D^{3−7α}.
- Child values: log X_{r2} ∼ D^{2α}λ/d₀ and log X_{F2} ∼ (3/d₀ + d₀b(Q₂))D^{2α}λ, so u_child ≍ D^α.
- Box exponents: σ₁ = α/(ad₀) (the algebraic level-1 cost u₂ log u₂ ≍ D^{3−7α}λ is L^{o(1)}) and
  σ₂ = (α/a)(4/d₀ + d₀b₂).
- Cost constant c = max(2a, a + s + max(a, σ₁ + σ₂), (a + s)/(1 − 2^{−α})), with s = 2κ(1 + C_γ)
  (`params.cost_constant`; for example 5.39 at α = 7/18, a = 3/5, κ = 1/8, d₀ = 1, b₂ = 2 and C_γ = 0, an illustrative
  value; it violates the constraint (A1) of §5, see there).

**Margins [proved at leading order; `params.margins`, 27 checks].** Every comparison used in §§2–5 holds for every
α ∈ (3/8, 2/5). At α = 3/8, two fail (first-order budget, √(AB)) and (SQ) is borderline. At α = 2/5, one fails (the
parent rational-pair loss e^{O(D^{1−2α}λ²)} against the B₂-saving e^{−cD^{3−7α}/λ²}) and three are borderline
(d log V₁ vs (μ + h) log N, the μ_opt branch, the level floor). The window is open at both ends. Among
the checks:
- (SQ) is automatic: u₂ log u₂ = o(log B);
- the thin-pair losses Ξ = e^{O(D^{3α−1}λ)} (child) and Ξ₂ = e^{O(D^{2−5α}λ)}, e^{O(D^{1−2α}λ)} (parent) are below the
  BV savings e^{−cD^α/λ²} and e^{−cD^{3−7α}/λ²}, and inside Hildebrand's range;
- the recursion on parent values costs L_{X_F}[2/5] = L[α]^{o(1)} (needs α > 4/11);
- IR rejection is e^{−Θ(D^{4α−1})} and content rejection 2^{2−d} = e^{−Θ(D^{1−2α})}.

### 1.2 Lemma η (typical number of children) [proved; rewritten in v2]

v1 combined an FM-type lower bound and a DS-type upper bound. Both lose (log N)^C = D^Cλ^C, which is fatal against the
D^{8α−3} margins. v2 bounds the ratio pointwise, with **no polylog loss**.

**Lemma η.** In design 2L_DW(α), for N ≥ N₀: W/μ₁ ≤ 3u₂, so η = μ₁/W ≥ (4κb₂/(3d₀))D^{−(3−7α)}(1 + o(1)).

*Proof.*
- **Reduction to a pointwise ratio.** W = Σ_{ℓ,P}E[m(P)1_{Acc_ℓ(P)}] and μ₁ = Σ_{ℓ,P}E[1_{Acc_ℓ(P)}]. The rational
  part of Acc is a function of g. A ratio of sums is at most the largest ratio of terms, so
  W/μ₁ ≤ sup_{(ℓ,P,g) filtered} E_{w,c}[m(P)1_{𝒮_{B₂}}(F(P)) | g]/E_{w,c}[1_{Acc}(P) | g].
- **Denominator.** E[1_Acc | g] ≥ E[1_{𝒮_{B₂}} | g](1 − L^{−a+o(1)}). This uses Lemma SQF's subtraction (log N)²/B
  against ϱ₂ = L^{−o(1)}.
- **The model law.** By Lemma DW(a)⁺, given g, F(P) is within TV 2N^{−η₀/2} of a mixture of uniform laws on
  I = Z ∩ (s + [0, ℓ_P)), with ℓ_P ≥ X_F/(64d) (b-window and the parent filter). These TV terms are ≤ N^{−η₀/3}
  relative to ϱ₂ = L^{−o(1)}. So it suffices to show, for every such interval I ⊂ [−X_F, X_F],
  R(I) := Σ_{B<q≤B₂}Ψ(I/q, B₂)/Ψ(I, B₂) ≤ 3u₂
  (n ∈ I with q | n is B₂-smooth iff n/q is).
- **Splitting I.** Split I into the part with |n| < X_F/K₁, K₁ := D⁴, and at most 2·64d·K₁ pieces J of relative length
  ≥ 1/(64dK₁) at scales x ∈ [X_F/K₁, X_F].
- **The pieces J.** For each J and q ≤ B₂, Hildebrand's short-interval theorem at y = B₂ gives
  Ψ(J/q, B₂) = (|J|/q)ρ(u_x − v_q)(1 + O(1/u₂)) and Ψ(J, B₂) = |J|ρ(u_x)(1 + O(1/u₂)), where v_q = log q/log B₂ ≤ 1.
  Hildebrand's error is log(u + 1)/log y. It is ≪ 1/u₂ because u₂ log u₂ ≪ log B₂, which holds for α > 4/11. The
  hypotheses hold: B₂ ≥ exp((log log x)^{5/3+ε}), and relative length 1/(64dK₁) ≥ B₂^{−5/12}.
- **Dickman.** −ρ′/ρ = ρ(u − 1)/(uρ(u)) is increasing (log-concavity of ρ, Alladi) and equals ξ(u)(1 + O(1/u))
  (Hildebrand–Tenenbaum survey, JTNB 1993, Cor 2.4: ρ(u − v) = ρ(u)exp{vξ(u) + O((1 + v²)/u)}, |v| ≤ u/2). So
  ρ(u − v)/ρ(u) ≤ e^{vξ(u)(1+O(1/u))} for v ≤ 1. Hence
  Σ_q(1/q)ρ(u_x − v_q)/ρ(u_x) ≤ (1 + o(1))∫_{log B}^{log B₂}e^{(ξ/log B₂)t}dt/t + O(λ) ≤ (1 + o(1))e^{ξ}/ξ + O(λ)
  ≤ (1 + o(1))u_x, using e^ξ = 1 + uξ. And u_x ≤ u₂.
- **The part |n| < X_F/K₁.** Its Ψ-mass relative to Ψ(I) is ≤ CdK₁^{−β₂} ≤ CdK₁^{−3/4}. On it m ≤ m̄ ≍ D^{2−4α}, so
  it contributes ≤ CdD^{−3}D^{2−4α} = o(1).
- **Total.** R(I) ≤ (1 + o(1))u₂ + o(1) ≤ 3u₂. ∎

The deterministic cap m̄ ≍ D^{2−4α} would not suffice: D^{−α} ≪ D^{−(2−4α)} fails for α < 2/5. The input
"Hildebrand at B₂" is added to the dependency list.

### 1.3 Prop PB (precision budget) [proved]

**Proposition PB.** Assume A := Σ_𝔮(w_𝔮/W)ϕ_𝔮 ≤ D^{−α}λ^{C_A} and B := Σ_𝔮(w_𝔮/W)σ_𝔮² ≤ D^{7α−3}λ^{C_B}, cell by cell
(Prop DC_e). Then, for N ≥ N₀(α, C_A, C_B), the dead-incidence mass is ≤ μ₁/2.

*Proof.* By Prop DC/DC_e and Cauchy–Schwarz (pc3 Cor DC), the dead mass is ≤ W(A + √(AB)). With Lemma η
(η ≥ cD^{−(3−7α)}, with no λ- or log-loss):
- A/η ≤ D^{−α+(3−7α)}λ^{C} → 0, because α > 3/8;
- √(AB)/η ≤ D^{(3−8α)/2}λ^{C} → 0, because (α + 3 − 7α)/2 − (3 − 7α) = (8α − 3)/2 > 0.

So W(A + √(AB)) ≤ Wη/2 = μ₁/2. ∎

**Lemma VD (variance decomposition) [proved].** Condition on R_𝔮 (and on a cell). Write:
- Rat(P′; g) for the rational part of Good(P′);
- A_{P′}(g) := Pr[F(P′)/q ∈ 𝒮_B | g, R_𝔮] and A₁₂(g) := Pr[both | g, R_𝔮];
- m(g) := Σ_{P′}Rat(P′; g)A_{P′}(g).

Then:
- ϕ_𝔮 ≤ 1/m_𝔮 + |ϕ^alg_𝔮| + ϕ^rat_𝔮;
- ϕ^alg_𝔮 := E[Σ_{P₁≠P₂}Rat₁Rat₂(A₁₂ − A₁A₂) | R]/m_𝔮²;
- ϕ^rat_𝔮 := Var(m(g) | R)/m_𝔮².

The same holds for the parent count, with σ^alg, σ^rat and 1/z_𝔮.

*Proof.* This is the law of total variance with respect to g. Given g, Var(N | g, R) = Σ Rat A(1 − A) + Σ_{≠}Rat₁Rat₂(A₁₂ − A₁A₂),
and the first sum is ≤ m(g). ∎

So C₁ splits into three parts:
- a one-point term (§4);
- an **algebraic block** ϕ^alg (§3, ours);
- a **rational block** ϕ^rat: the variance over g of a sum of rational indicators with *smooth* weights A_{P′}(g)
  (c12rat). By Lemma DW(a)⁺, A_{P′}(g) depends on g only through |G_g(P′)|, up to a relative error N^{−η₀/3}.

---

## 2. Design DW and the internal chain (task 2) [proved]

### 2.1 The field law and Lemma GB-DW

**Lemma GB-DW.** Let gcd(g₀, g₁) = 1, gcd(g₁, N) = 1 and g₀, g₁ ∈ [M/2, M]. Let F̂_N be paper2 lem:GB applied to T_d = N
(multiplier outside). Then:
- F̂_N(g₀, g₁) = N;
- a_d g₀^d ≡ N (mod g₁), so a_d is a unit mod g₁;
- ‖F̂_N‖_∞ ≤ d2^dM, because |a₀| ≤ 2^dN/M^d + d2^{d−1}M and N/M^d ≤ 2^dN^ν ≤ M at μ = μ_opt.

The DW form F := wF̂_N + G·R_c satisfies F(g₀, g₁) = wN, f(m) ≡ 0 (mod N) for m = g₀/g₁ mod N, and Res(f, g₁x − g₀) = ±wN.
Its leading coefficient is f_d = wa_d + g₁c₀.
*Proof.* paper2's proof with Nw replaced by N; f(m) ≡ g₁^{−d}F(g₀, g₁) = g₁^{−d}wN ≡ 0. ∎

*Degree drop.* Without gcd(w, g₁) = 1, f_d = 0 is possible. For fixed (g, w) it happens for at most one c₀, so with
probability ≤ 1/(2H + 1). SAMPLE_FIELD-DW **rejects** it. Also 0 < |f_d| < N, so the gcd check of f_d is a factor exit or
a pass.

*Why no coprimality with g₁.*
- An example shows that w coprime to g₁ restricts F(P) mod p for p | gcd(t, g₁). The law is then not
  uniform on the DW-L lattice.
- With unrestricted w, Lemma DW-L describes the support, and the law is (Lemma TP-W) uniform on it.

### 2.2 Lemma CU (unit residue) and the forced residue [proved]

**Lemma CU.** Let P = (a, b) be primitive and t := |G_g(P)| ≥ 2 with gcd(t, N) = 1. Then gcd(F̂_N(P), t) = 1, and
F(P) ≡ w·F̂_N(P) (mod t). Hence gcd(F(P), t) = gcd(w, t), and F(P) mod t is uniform over all residues as w runs over a
complete system mod t.
*Proof.* Let p | t.
- *If p ∤ g₁:* g₁P ≡ b(g₀, g₁) (mod p), so g₁^dF̂_N(P) ≡ b^dN. Here p ∤ b, since otherwise p | g₁a, so p | a,
  contradicting primitivity.
- *If p | g₁:* then p ∤ g₀, p | g₀b, so p | b and p ∤ a. Then F̂_N(P) ≡ a_da^d ≢ 0 (mod p) by Lemma GB-DW.

The congruence F ≡ wF̂_N (mod G(P)) holds because F − wF̂_N = G·R_c. ∎

This replaces [pp] lem:coprime, which used w prime > t. No hypothesis on b or g₁ is needed;
conclusion. E2 checked Lemma CU and the forced residue on 5 911 atoms: 0 violations.

### 2.3 Lemma RE (root event is uniform in g) [proved]

**Lemma RE.** For a prime q ∤ N with q ≤ B₂ and any r mod q: Pr_{w,c}[q | F(r, 1) | g] = 1/q + O(1/H + 1/W_m), uniformly
in g.
*Proof.*
- If q ∤ G(r, 1): F(r, 1) = wF̂_N(r, 1) + G(r, 1)R_c(r, 1), and R_c(r, 1) = Σ c_ir^{d−1−i} has coefficient 1 at c_{d−1}.
  Fix the other coordinates; c_{d−1} uniform on 2H + 1 consecutive integers is uniform mod q up to q/(2H + 1).
- If q | G(r, 1): then q ∤ g₁ (else q | g₀), and F(r, 1) ≡ wg₁^{−d}N (mod q), so the event is q | w. That has probability
  1/q + O(1/W_m).

Exhaustive check: `test_root_event_uniform_in_g`. ∎

If q | G(r, 1), every primitive P′ ∈ Λ_𝔮 has q | G(P′) and is not Good, and the label has no parent (q ∤ G(P)). Such
(g, 𝔮) contribute a dead environment of conditional probability ≈ 1/q ≤ 1/B, which is negligible against η.

### 2.4 Lemma CT (content) and Lemma IR-DW [proved]

**Lemma CT.** For a prime p ∤ N, p divides the content of F iff p | w and p | c_i for every i. In particular the content
e (coprime to N) divides w, and
Pr_{w,c}[content > 1, coprime to N] ≤ Σ_{p≤2W_m}(1/p + 1/W_m)(1/p + 1/(2H + 1))^d ≤ 2^{2−d}. The sum is over p | w, so
p ≤ 2W_m. For p > 2H + 1, "p | c" means c = 0, of total probability ≤ (2H + 1)^{−d}.
*Proof.*
- If F ≡ 0 (mod p), then wN = F(g₀, g₁) ≡ 0, so p | w (as p ∤ N).
- If p | w: F ≡ G·R_c, and G ≢ 0 (mod p), so F ≡ 0 iff R_c ≡ 0, i.e. iff p | c for every coordinate. ∎

SAMPLE_FIELD-DW: gcd-check the content e against N (a proper divisor is a factor exit, since e ≤ |f_d| < N). If e > 1,
reject. This replaces "content ∈ {1, w}; reject if w". E2 observed 5 392 content rejections against 5 398.6 predicted (z = −0.09). *This probability is not N^{−Ω(η₀)}*. Rejection is therefore used
**globally**, as Lemma reject is used in select-U, and all moment statements (DC, C₁, C₂, TV) are made under the raw law.

**Lemma IR-DW.** Let 𝒫_IR be the primes in (r₀(d), ½ log H], and 𝒫_w := {r ∈ 𝒫_IR : r | w}.
1. For r ∈ 𝒫_IR ∖ 𝒫_w, the law of f mod r over c is the one in [pp] lem:IR, so the irreducible fraction is ≥ 1/(4d):
   - if r ∤ g₁, f mod r ranges over the polynomials with f(m_r) = wNg₁^{−d} ≠ 0;
   - if r | g₁, over those with leading coefficient wa_d ≢ 0.
2. For r ∈ 𝒫_w, r cannot certify: if r ∤ g₁, f(m_r) ≡ 0; if r | g₁, deg(f mod r) < d.
3. Pr_{w,c}[no certificate] ≤ 2^{1−|𝒫_IR|} + e(1 − 1/(4d))^{|𝒫_IR|/2} ≤ exp(−log H/(26d log log H)).

*Proof of 3.* Write n := |𝒫_IR|, Q_IR := ∏𝒫_IR = H^{1/2+o(1)} ≤ W_m (since H ≤ W_m²).
- For each S ⊂ 𝒫_IR of size s := ⌈n/2⌉: Pr_w[∏_S r | w] ≤ 2/∏_S r ≤ 2·16^{−s}.
- So Pr_w[|𝒫_w| ≥ n/2] ≤ 2^n·2·16^{−n/2} = 2^{1−n}.
- On the complement, the CRT argument of [pp] lem:IR (c mod Q_IR is ≤ e × uniform) gives e(1 − 1/(4d))^{n/2}. ∎

A worst-case w divisible by every certificate prime exists, so averaging over w is necessary. The sketch is
the same. E2: 0 violations of claim 2.

**Lemma reject-DW.** Under the raw DW law, the total rejection probability is
Pr[ℛ] ≤ exp(−log H/(26d log log H)) + 3/(2H + 1) + 2^{2−d}.
It collects: no certificate; N | f′(m) (unchanged: coefficient g₁ of c_{d−1}); degree drop; content > 1. In the window
this is e^{−Θ(D^{4α−1})} + e^{−Θ(D^{1−2α})} ≤ γ₂/4.

### 2.5 Lemma W-DW (the second moment's multiplier step) [proved]

**Lemma W-DW.** For 1 ≤ t ≤ W_m and every Φ₀: Z/t → [0, ∞), E_w Φ₀(wu mod t) ≤ (2/t)Σ_{c mod t}Φ₀(c), for every u with
gcd(u, t) = 1.
*Proof.* w ↦ wu permutes Z/t. And w is uniform on W_m consecutive integers, so each class mod t has probability
≤ (⌊W_m/t⌋ + 1)/W_m ≤ 2/t. ∎

*Use in LD*-U′.* The steps of [pp] Thm LD / lv5unif LD*-U are:
- Step 3: 𝔮(P) ≤ 𝔔₀(t, |Δ|; r₀) with r₀ = F(P) mod t = wu (Lemma CU).
- Step 4: Lemma W-DW gives E_w𝔔₀ ≤ (2/t)Σ_{r mod t}𝔔₀(t, |Δ|; r) = (2/t)𝔔(t|Δ|). The last equality holds because
  F ≡ F′ (mod t|Δ|) implies equal residues mod t. Since 2/t ≤ 2/φ(t), this is the old (C_W/φ(t))𝔔(t|Δ|) with C_W = 2
  instead of 24/η₀.
- Steps 5–7 depend on w only through this bound, so they are verbatim.

The old proof summed over reduced classes only, a subset; the all-residue sum is exactly 𝔔(t|Δ|). So LD*-U′ is valid
under DW, with a better constant.

### 2.6 Lemma DW(a)⁺ (one-point law, with root coset) [proved]

**Lemma DW(a)⁺.** Fix g and a primitive P with t := |G(P)| ≥ M|P|_∞N^{−η₀/4}, gcd(t, N) = 1, and put
ℓ_P := t(2H + 1)|P|_∞^{d−1}. Then:
- the law of F(P) over (w, c) is within TV
  ε₁ := [3d²|P|_∞ + 4 + 4d²2^dW_mN^{η₀/4}]/(2H + 1) + 2t/W_m ≤ 2N^{−η₀/2}
  of a mixture of uniform measures on Z ∩ (s + [0, ℓ_P)), s ∈ R. That is, translates of **one** interval of length
  ℓ_P, inside [−X_F(P), X_F(P)], X_F(P) := (2d + 1)MH|P|_∞^d;
- for a prime q ∤ tN with q ≤ B₂, conditional on q | F(P), the law is within TV 2qε₁ + 4q/(2H + 1) of the same mixture
  restricted to qZ.

*Proof.*
- [pp] lem:smooth1: Law(R_c(P)) is within TV 3d²|P|/(2H + 1) of ϑ = Law(ΣY_i), with Y_i independent uniform on integer
  intervals of lengths (2H + 1)l_i. Take i₀ with l_{i₀} = |P|_∞^{d−1}. Then ϑ = Law(Σ_{i≠i₀}Y_i) ∗ Unif(interval of length
  (2H + 1)|P|^{d−1}).
- Multiply by G(P) and add wF̂_N(P). For w in a fixed class mod t, the values lie in one class mod t (Lemma CU). As w runs
  over a complete system mod t, the classes run over Z/t, so the union is uniform on Z ∩ (translate of [0, ℓ_P)).
- The archimedean spread of wF̂_N(P) over w ∈ (W_m, 2W_m] is ≤ 2W_m(d + 1)d2^dM|P|^d (Lemma GB-DW), which is
  ≤ 4d²2^dW_mN^{η₀/4}·ℓ_P/(2H + 1). This costs that much TV, and the w-classes cost 2t/W_m.
- If ab = 0 (P = (±1, 0) or (0, ±1)), smooth1 does not apply, but R_c(P) = ±c₀ or ±c_{d−1} is itself uniform on an
  interval of length 2H + 1. The same conclusion holds.
- This is pc3's DW(a) with explicit constants. Note the layer-cake intervals of smooth1 can be short; the
  one-long-factor decomposition above avoids "small layers". The q-coset statement is the one-point case of
  Cor TP-q below. ∎

### 2.7 Level 1 under DW (U′-DW(α)) [proved, internal-chain status]

uprime's §§3–6 hold in design 2L_DW(α), with these replacements (and only these).
- **GM-U′ → DW(a)⁺.** Pr_{w,c}[F(P) ∈ 𝒮_{B₂} | g] ≥ ϱ₂(X_F)(log N)^{−C} for **every** filtered (ℓ, P, g). No good
  modulus is needed. Lemma 8.5H and Cor bad₂ are therefore **not used**, and (H1)–(H5) need not be checked in the
  window.
- **FM-U′.** The rational good mass (FIB) times the pointwise algebraic factor, minus the (SQ) subtraction
  (log N)²/B. That subtraction is now ≤ ϱ₂ L^{−a+o(1)}, since ϱ₂ = L^{−o(1)}. There is no bad-mass step.
- **SQF under DW.** Cases (a) and (c) are unchanged. In (b), if q | w (probability 1/q), F ≡ G·R_c mod q, and the same
  Schwartz–Zippel argument applied to G·L_P·H′ gives ≤ 2e(d − 1)/q².
- **LD*-U′.** Lemma W-DW in Step 4. Kern-U needs bb′ ≤ V₁² ≤ (2H + 1)/d, which holds. DS-U′ and L2-U′ hold, with u₂ ≍
  D^{3−7α} → ∞. The class (B)/(C) ratios are N^{−Θ(h)}B^{O(1)}.
- **CS′, the child scans, recursion and IR.** Unchanged in form, with DW's envelopes and Lemma IR-DW. Lemma OS, W and
  T_rec are verbatim.
- **Range restatement.** uprime's U1₂, DS-U′ and CS′ are stated for log z/log B₂ ≥ D^{1/4}. In the window
  u₂ ≍ D^{3−7α}, which is < D^{1/4} once α > 11/28. Their proofs use only u₂ → ∞, i.e. (log log N)⁴ ≤ u₂, together with
  HT Theorems 1–3, CEP and Rankin. They are restated with the hypothesis log z/log B₂ ≥ u₂/2. This is a routine
  restatement, not re-proved line by line here.

"These replacements and only these" in v1 was an overclaim: the range restatement is a further, routine change.

### 2.8 POSTPROCESS and LV6-U under DW [proved]

Fix an accepted DW field: gcd(g₁f_dwN-checks), f primitive and irreducible, f_d^{d−2}f′(m) a unit mod N. Then:
- ψ: Z[β] → Z/N, β ↦ f_dm, is well defined (Lemma GB-DW).
- Norm(z_P) = f_d^{d−1}F(P) and ψ(z_P) = f_dG(P)/g₁.
- The row identity (pp eq:row) is unchanged, since it involves g₁, n_P and ψ(z_P) only.
- Lemma C needs r ∤ D = f_d disc(f) only.
- Lemma H uses ‖f‖_∞ ≤ 2W_m d2^dM + 2MH ≤ 3MH < N.
- Lemmas E and X are degree-only (lv5unif §9: factorial-bit characters).
- SQ uses Euler's lemma and exact verification.
- The congruence 𝒳² ≡ 𝒴² is unchanged.
- Two-level rows: nested Lemma RW, unchanged.

Composite w is **never factored** and appears in no column. The factor exits (gcd checks of g₁, w, f_d, the content,
G(P), F(P) and the child primes) and transcript invariance (T) are unchanged. E2 ran the thmA pipeline
with DW fields:
- 31/32 complete factorisations; the single failure is a field-selection cap, factored correctly post hoc;
- 0 violations of the row identity and the congruence;
- 120/120 tested kernel vectors were squares;
- DD fruitfulness 26/56 (ω = 2, z = −0.53) and 23/32 (ω = 3, z = −0.41).

### 2.9 If DW were rejected

DW-L needs a multiplier that is uniform modulo G(P₁)G(P₂) ≈ (MV₁)² and independent of c.
- *A prime w* would need primes in (W, 2W] to be equidistributed modulo individual moduli ≈ (MV₁)². That is not known,
  since BV holds only on average, and a pair needs the joint residue.
- *A coefficient-only randomisation* at fixed resultant cannot move the class of F(P) mod G(P) (paper2 Prop gauss).

We found no alternative that keeps the exact two-point lattice. DW's integer multiplier is the minimal change we know.

---

## 3. Two-point smoothing and the algebraic blocks (task 1)

### 3.1 Two elementary lemmas [proved]

For points P_i = (a_i, b_i), write v(P) = (a^{d−1−i}b^i)_{0≤i<d}, w_i := (v_i(P₁), v_i(P₂)), m_{jk} := det(w_j, w_k) and
Δ := |a₁b₂ − a₂b₁|.

**Lemma MM (a large minor).** For d ≥ 3 and non-proportional P₁, P₂ (zero coordinates allowed),
max_{j<k}|m_{jk}| ≥ Δ(|P₁|_∞|P₂|_∞)^{d−2}/4.
The maximum is attained, up to the factor 4, at one of (0, 1), (d − 2, d − 1), (0, d − 1), (0, d − 2), (1, d − 1). The
chosen pair contains, for each i, a column with |v(P_i)| ≥ |P_i|_∞^{d−1}/2.
*Proof.* Put X := a₁b₂ and Y := a₂b₁, so Δ = |X − Y|. Then m_{jk} = (a₁a₂)^{d−1−k}(b₁b₂)^j(X^{k−j} − Y^{k−j}), and
A := |P₁|_∞|P₂|_∞ = max(|a₁a₂|, |b₁b₂|, |X|, |Y|).
- If A = |a₁a₂|: |m₀₁| = A^{d−2}Δ.
- If A = |b₁b₂|: |m_{d−2,d−1}| = A^{d−2}Δ.
- If A = |X| (the case |Y| is symmetric): |a₁| ≥ |b₁| and |b₂| ≥ |a₂|, and |a₁a₂|, |b₁b₂| ≥ |Y| (from
  |X||Y| = |a₁a₂||b₁b₂| ≤ A·min(…)).
  - If |Y| < |X|/2: |m_{0,d−1}| ≥ |X|^{d−1} − |Y|^{d−1} ≥ |X|^{d−1}/2 ≥ A^{d−2}Δ/4.
  - If |Y| ≥ |X|/2 and X, Y have the same sign: m_{0,d−1} = (X − Y)ΣX^{d−2−i}Y^i, with all terms of one sign, so
    |m_{0,d−1}| ≥ Δ|X|^{d−2}.
  - If |Y| ≥ |X|/2 and the signs differ: then Δ = |X| + |Y| ≤ 2|X|.
    - For d − 2 odd: |m_{0,d−2}| = |a₁a₂|(|X|^{d−2} + |Y|^{d−2}) ≥ (|X|/2)|X|^{d−2} ≥ A^{d−2}Δ/4.
    - For d − 1 odd: |m_{0,d−1}| = |X|^{d−1} + |Y|^{d−1} ≥ A^{d−2}Δ/2.
- The column claim follows from the dominance relations in each case. ∎

(`test_max_minor_lemma_MM`: 3 000 random and 5 adversarial pairs, d ∈ {3, 4, 5, 7, 9}.)

**Lemma BZ (Bezoutian bound).** For binary forms A, K of degree n and any P₁, P₂:
- Δ_s := a₁b₂ − a₂b₁ divides A(P₁)K(P₂) − A(P₂)K(P₁);
- |A(P₁)K(P₂) − A(P₂)K(P₁)| ≤ 2n³‖A‖_∞‖K‖_∞Δ(|P₁|_∞|P₂|_∞)^{n−1}.

*Proof.* The quotient is the Bezoutian Σ_{i,j<n}B_{ij}m_i(P₁)m_j(P₂). Each Bezout-matrix entry is a sum of ≤ n terms
α_sκ_t − α_tκ_s, and each monomial is ≤ (|P₁||P₂|)^{n−1}. ∎ (`test_bezoutian_bound`.)

### 3.2 Lemma TP and its corollaries [proved]

**Lemma TP (coefficient part).** Let d ≥ 3, let P₁, P₂ be primitive and non-proportional, V := max|P_i|_∞, and let
(j, k) be any pair with m_{jk} ≠ 0. Put Φ(c) := (R_c(P₁), R_c(P₂)), Γ := Φ(Z^d), M_{jk} := (w_j | w_k),
Π := M_{jk}[−H, H + 1)², and π := Law(Σ_{i∉{j,k}}c_iw_i). For c uniform on C_H:
- ‖Law Φ(c) − π ∗ Unif(Γ ∩ Π)‖_TV ≤ d(1 + 4V²)/(2H + 1);
- #(Γ ∩ (s + Π)) = (2H + 1)²|m_{jk}|/Δ for every s ∈ R².

*Proof.*
1. ker Φ has the basis e_m = Q·x^{d−3−m}y^m (Kern-U), Q = (b₁x − a₁y)(b₂x − a₂y). Each coordinate of Σ_m s_me_m with
   s ∈ [0, 1)^{d−2} has absolute value ≤ |b₁b₂| + |a₁b₂ + a₂b₁| + |a₁a₂| ≤ 4V².
2. Let R := Γ ∩ M_{jk}[0, 1)², a set of representatives of Γ/Γ_{jk}, of size |m_{jk}|/Δ = [Γ:Γ_{jk}]. For z = M_{jk}s ∈ R:
   - p := (s at coordinates j, k; 0 elsewhere) satisfies Φ(p) = z and ‖p‖_∞ < 1;
   - for any integer preimage x₀, p − x₀ ∈ ker Φ ⊗ R;
   - rounding in the basis e_m gives x(z) ∈ Z^d with Φ(x(z)) = z and ‖x(z)‖_∞ ≤ 1 + 4V².
3. Let ξ := x(z) with z uniform on R, independent of c, and c′ := c + ξ. Then Φ(c′) = Σ_{i∉{j,k}}c_iw_i + [M_{jk}(c_j, c_k) + z].
   The bracket is uniform on Γ ∩ Π: unique decomposition, since R represents Γ/M_{jk}Z².
4. TV(Law c, Law c′) ≤ max‖ξ‖₁/(2H + 1), because translating the cube by y costs ≤ ‖y‖₁/(2H + 1).
5. Counting: s + Π is a disjoint union of (2H + 1)² translates of the half-open fundamental domain of Γ_{jk}. ∎

d(1 + 4V²) is pessimistic. (`dw.short_preimage`, `test_short_preimages_and_representatives`.)

**Lemma TP-W (with the multiplier).** In design DW, fix g. Take P₁, P₂ as above with G_i := G_g(P_i),
|G_i| ≥ M|P_i|_∞N^{−η₀/4}, gcd(N, G₁G₂) = 1. Put v := (F̂_N(P₁), F̂_N(P₂)), 𝔇 := diag(G₁, G₂), X := (F(P₁), F(P₂)) =
wv + 𝔇Φ(c), T := |G₁G₂|, and choose (j, k) by Lemma MM. Let L := Zv + 𝔇Γ, which equals {X₁ ≡ λ^dX₂ (mod Δ)} by
pc3's DW-L, with [L : 𝔇Γ] = T. Then
  ‖Law(X) − Σ_y π(y)·Unif(L ∩ (W_mv + 𝔇y + 𝔇Π))‖_TV ≤ [d(1 + 4V²) + 4 + 16d⁴2^dW_mN^{η₀/2}]/(2H + 1) + 2T/W_m ≤ 3N^{−η₀/2}.

*Proof.*
1. For fixed w, apply Lemma TP and the injective map 𝔇. Law(X | w) is within TV d(1 + 4V²)/(2H + 1) of
   κ_w := Σ_yπ(y)Unif(Λ_w ∩ (wv + 𝔇y + 𝔇Π)), where Λ_w := wv + 𝔇Γ depends only on w mod T.
2. *Archimedean shift.* Moving the window from wv + 𝔇y + 𝔇Π to W_mv + 𝔇y + 𝔇Π moves it by σ := (w − W_m)(𝔇M_{jk})^{−1}v in
   the parallelogram's own coordinates.
   - Every translate holds the same number of points of Λ_w (Lemma TP's count, applied to 𝔇Γ_{jk} ⊂ 𝔇Γ).
   - The symmetric difference consists of ≤ (2H + 1)(|σ|₁ + 4) unit cells of 𝔇Γ_{jk}.
   - So the TV cost is ≤ (|σ|₁ + 4)/(2H + 1).
3. *The size of σ.* By Cramer, |σ|_∞ ≤ W_m·max_{i∈{j,k}}|det(v, 𝔇w_i)|/|G₁G₂m_{jk}|. Here
   det(v, 𝔇w_i) = F̂(P₁)K_i(P₂) − F̂(P₂)K_i(P₁) with K_i := G·x^{d−1−i}y^i, a form of degree d with ‖K_i‖ ≤ M.
   - Lemma BZ: |det(v, 𝔇w_i)| ≤ 2d³·d2^dM·M·Δ(|P₁||P₂|)^{d−1}.
   - Lemma MM: |m_{jk}| ≥ Δ(|P₁||P₂|)^{d−2}/4.
   - The filter: |G₁G₂| ≥ M²|P₁||P₂|N^{−η₀/2}.
   - Hence |σ|_∞ ≤ 8d⁴2^dW_mN^{η₀/2}, **independently of the angle between P₁ and P₂**. Obtained
     by a different route.
4. *Mixing the classes.* Σ_wPr(w)κ′_w = Σ_{r mod T}p_rκ′_r, with Σ_r|p_r − 1/T| ≤ 2T/W_m. Moreover
   (1/T)Σ_rκ′_r = Σ_yπ(y)Unif(L ∩ (…)), because {rv + 𝔇Γ} runs once over L/𝔇Γ and each coset has the same number of
   points in each translate.
5. *Design bound.* With H ≥ 2^{d+8}d⁴W_mN^{η₀} and T ≤ (2MV₁)² ≤ W_mN^{−η₀}, the total is ≤ 3N^{−η₀/2}. ∎

**Corollary TP-q (root coset).** Suppose moreover P₁, P₂ ∈ Λ_𝔮 with q ∤ b₁b₂G₁G₂N, so that q | Δ. Then:
- R_𝔮 ⇔ q | X₁ ⇔ q | X₂;
- conditional on R_𝔮, X/q is within TV 7qN^{−η₀/2} of Σ_yπ(y)Unif(L′ ∩ q^{−1}(W_mv + 𝔇y + 𝔇Π)), where
  L′ := {x₁ ≡ λ^dx₂ (mod Δ/q)}.

*Proof.*
- F(P_i) ≡ b_i^dF(r, 1) (mod q).
- Under the model, the probability of q | X₁ is (1/q)(1 + O(q/H)) in every translate. Each translate contains
  ⌊(2H + 1)/q⌋² fundamental domains of q𝔇Γ_{jk} ⊂ L ∩ (qZ)² = L_q, and [L : L_q] = q.
- Conditioning multiplies the TV by ≤ 2/Pr(event).
- q² | Δ is possible unless constraint (A1) of §5 holds. The divided lattice L_q/q = {x₁ ≡ λ^dx₂
  (mod Δ/q)} is valid regardless, and q | Δ′ is then covered by retaining its gcd-class factor. v2's "q² ∤ Δ" was
  false without (A1). ∎

So the **DW(b) TV statement is proved**. pc3 §7 item 3 asked for this, and it holds uniformly,
nearly parallel pairs included.



---

# G. Prop CLB and the conditional assembly 2L-DW (source: notes/c12alg_results.md §§4-5)

## 4. The one-point child lower bound (task 3) [proved]

**Proposition CLB.** In design 2L_DW(α) (and at 3/8 with the analogous parameters), for every label 𝔮 = (q, r) and every
cell e of relative size D^{−A₀}:
m_{𝔮,e} := E[N_𝔮 | R_𝔮, e] ≥ m₂(q)(log N)^{−C} ≥ e^{λ²}(log N)^{−C} ≥ 8/η (for N ≥ N₀).
Here m₂(q) = S₂²ϱ(X_{r2})ϱ(X_{F2}) ≥ e^{λ²} is the nominal main term of uprime's child scan, **rerun with DW's envelope**
X_{F2} = (2d + 1)MHV₂^d/Q₂, H = H_DW.

*Proof.*
1. **Box points.** Box(𝔮) has ≥ (2S₂ + 1)²/4 primitive coordinate vectors (i, j) mod ±. P′ = ie₁ + je₂ fails to be
   Z²-primitive only if q divides its content. That happens for ≤ (2S₂ + 1)²/q + 1 of them. In all cases |P′| ≤ V₂.
   This holds for skewed labels too, since the conservative box is defined in coordinates.
2. **Rational factor.** Fix a primitive P′ ∈ Λ_𝔮. As in [pp] lem:fibres, g ↦ G_g(P′) has fibres g* + jP′ of length
   ≥ M/(8|P′|) ≥ N^{η₀/2} inside the cell (|P′| ≤ V₂ ≪ MN^{−η₀}D^{−A₀}).
   - Möbius over e | G(P′) gives the primitive points with loss φ(n)/n ≥ 1/(C log log N).
   - **No angular discard.** v1 discarded points within angle D^{−A₀} of the cell's
     direction. That set can be almost the whole box for skewed labels.
     - For each fixed P′, G_g(P′) over g in the cell ranges over an interval J(P′, e) of length ≍ M|P′|D^{−A₀}, possibly
       around 0. The filter removes only the part with |G| < M|P′|N^{−η₀/4}, a fraction ≤ N^{−η₀/4}D^{A₀}.
     - Restrict to an interior interval of the remaining values, of relative length ≥ 1/4, at scale between
       M|P′|N^{−η₀/4} and M|P′|. On it the fibres have length ≥ M D^{−A₀}/(8|P′|) ≥ N^{η₀/2}, and the B-smooth values
       have density ≥ ϱ(M|P′|)(log N)^{−C}, since ϱ is nonincreasing.
     - Smaller values only help the lower bound. On the algebraic side the scale |G|H|P′|^{d−1}/q is also smaller, so
       ϱ is again larger.
   - By Lemma U1 at B, Pr_g[Rat(P′; g) | e] ≥ ϱ(M|P′|)(log N)^{−C} ≥ ϱ(X_{r2})(log N)^{−C}, since ϱ is nonincreasing and
     M|P′| ≤ X_{r2}.
3. **Algebraic factor.** For g with Rat(P′; g), Lemma DW(a)⁺ with the root coset applies (q ∤ G(P′) since G(P′) is
   B-smooth; q ∤ N; q ∤ b′): F(P′)/q is within TV 2qN^{−η₀/2} of a mixture of uniform measures on integer intervals of
   length ≥ (2H + 1)|G(P′)||P′|^{d−1}/(2dq) in [−X_{F2}, X_{F2}].
   - Every such interval I of length ℓ_I: Ψ(I) ≥ ℓ_Iϱ(X_{F2})(log N)^{−C}, by Lemma U1(iii) (interval) and (ii) (shift
     by ≤ B^{C}). Points near 0 are a ≤ 1/ϖ fraction (small layers).
   - So A_{P′}(g) ≥ ϱ(X_{F2})(log N)^{−C} − 2qN^{−η₀/2} ≥ ϱ(X_{F2})(log N)^{−C−1}.
   - Lemma RE: conditioning on R_𝔮 changes the g-law by 1 + O(q/H).
4. **Sum.** m_{𝔮,e} ≥ Σ_{P′}Pr[Rat | e]·inf_gA ≥ cS₂²ϱ(X_{r2})ϱ(X_{F2})(log N)^{−2C−1}. The scan threshold
   4(log N)^{2C_R}e^{λ²} absorbs every fixed polylog (uprime §4.2). ∎

*Comparison with uprime's "routine".* (Step 1's "this holds for skewed labels too" is now true without the discard.) The claim is confirmed, with two necessary corrections:
- (a) The child filter must be part of Good(P′), otherwise DW(a) fails at points with tiny |G(P′)|. Its
  pass fraction is 1 − O(N^{−η₀/4}), so it introduces no O(1) environment term (unlike the parent filter, pc3 M1).
- (b) S₂ must be rescanned with DW's larger height. Keeping uprime's S₂ would lose an L[α]-scale factor.

Fields with q | G(r, 1) give N_𝔮 = 0. They are included in m_𝔮's average with conditional weight ≈ 1/q, which changes
m_𝔮 by a factor 1 − O(1/q). Hence (A's one-point term) Σ(w/W)/m_𝔮 ≤ e^{−λ²}(log N)^{C} ≪ η.

---

## 5. The achieved exponent

**Hyp RAT(α) (the rational blocks, first order).** In design 2L_DW(α), cell by cell:
- (R1) Σ_{𝔮: s(𝔮)≤T}(w_𝔮/W)ϕ^rat_𝔮 ≤ D^{−α}λ^{C}, where ϕ^rat_𝔮 = Var_g(Σ_{P′∈Box(𝔮)}Rat(P′; g)a_{P′}(|G_g(P′)|))/m_𝔮².
  Here a_{P′}(n) is the **idealised** weight: the probability that a variable with law (the DW(a)⁺ mixture for
  t = n, restricted to qZ) divided by q is B-smooth. It is explicit and depends only on (P′, n, q). It differs from
  A_{P′}(g) by ≤ 4qN^{−η₀/2} in absolute terms (TV), uniformly in g with |G_g(P′)| = n.
  - It is ≍_{polylog} ϱ(nH|P′|^{d−1}/q).
  - v1 claimed an identification with a single density value to relative error N^{−η₀/3}. That is false: the shape of
    the convolution changes its smoothness average at first order (1 − β)log d.
  - RAT must be stated with the idealised weight.
- (R2) Σ_𝔮(w_𝔮/W)σ^rat_𝔮 ≤ D^{7α−3}λ^{C}, the same for the label's parent count, with weights a_P(|G_g(P)|) at B₂.

These are pc3 §7 items 1 and 4 at **first order** in the window, for rational values with smooth weights. Skewed labels
(item 2) are excluded by Theorem AB(iii). By the same template as Theorem AB (with LSQ for the mod-q part), they are
plausible, but they are c12rat's. We did not attempt them.

**Theorem 2L-DW(α) [proved, conditional].** Let α ∈ (3/8, 2/5) be fixed. Assume:
- (i) the internal chain as in Theorem U and Theorem 2L;
- (ii) Hyp RAT(α);
- (iii) Theorem AB (§3.3, sketch).

Then FACTOR_2L in design 2L_DW(α) is Las Vegas and factors every N in expected time L_N[α, c(α)], with c(α) =
max(2a, a + s + max(a, σ₁ + σ₂), (a + s)/(1 − 2^{−α})) for any admissible (a, κ, d₀, b₂).

*Proof.*
1. Lemma VD gives A ≤ Σ(w/W)/m + A^alg + A^rat and B ≤ Σ(w/W)/z + B^alg + B^rat.
   - The one-point terms are ≤ e^{−λ²}(log N)^C (Prop CLB) and, for parents, use the exact identity
     w_𝔮/z_𝔮 = Pr(R_𝔮)/#J. Then Σ_𝔮(w_𝔮/W)/z_𝔮 = Σ_𝔮Pr(R_𝔮)/(#J·W).
     - Σ_r Pr(R_{(q,r)}) ≤ d + 1 + q·Pr(F ≡ 0 mod q) ≤ d + 2, averaged under the raw law. F ≡ 0 mod q is possible but
       has probability ≤ (1/q)(1/q + 1/H)^d, by Lemmas RE and CT.
     - W ≥ (1 − o(1))μ₁: an accepted parent with no child has F(P) B-smooth, a negligible fraction (cf. c12rat Lemma
       DOM(iii)).
     - So the sum is ≤ (d + 2)π(B₂)/((1 − o(1))#J·μ₁) ≤ e^{b₂D^{4α−1}λ}(log N)^C/L[1 − α] → 0.

     No lower bound on z_𝔮 is needed.
   - Theorem AB(i), (iii) and RAT(R1) give A ≤ D^{−α}λ^{C}. RAT(R1) is restricted to s ≤ T, matching AB(i).
   - Theorem AB(ii) and RAT(R2) give B ≤ D^{7α−3}λ^{C}.
2. Prop PB: E N₂ ≥ μ₁/2.
3. U′-DW(α) (§2.7) gives E[N₁²] ≤ (log N)^C(μ₁² + μ₁) and μ₁ ≥ (log N)^{−C}, so γ₂ ≥ (log N)^{−C} by Paley–Zygmund
   (Prop UB).
4. Theorem K-U and Prop select-U with Lemma reject-DW give the split cost (γ₂/(16 log N))^{2k}/(32k²) = L[α, −O(κ)].
5. Rows (nested Lemma RW), POSTPROCESS-DW (§2.8), Theorem DD with Lemma T_multi, and Lemma F_α at α finish the proof, as
   in multidescent Theorem 2L. ∎

**Theorem 2L-DW₂ [proved, conditional; v2.1].** For every fixed α ∈ [3/8, 2/5), FACTOR_2L in design 2L_DW(α) (at 3/8
with the constraint a > 2·3d₀/(32κb₂) of (O4)) factors every N in expected time L_N[α, c(α)]. At α = 3/8, c includes
σ₁(3/8) = 3/(8ad₀) + 3d₀/(32κb₂), since the parent B₂-smoothness cost is leading at the endpoint. No numerical c(3/8) is claimed. It assumes:
- (i) the internal chain;
- (ii) **Hyp RAT₂(α)**: (R1) with D^{−2α}λ^C, and (R2) with (D^{−2α} + D^{−2(3−7α)})λ^C on the right-hand side. This is the
  form c12rat's Theorem RAT/Corollary claims;
- (iii) **the statement of Theorem AB₂(i)–(iii)** (§3.6). Its proof is a sketch with the open obligations (O1)–(O3)
  and (O5); O4 is discharged by the constraint a > 2c₂.

*Admissible parameters at α = 3/8:* (SQ) a·b₂ > 3d₀/(32κ); a > 2c₂ with c₂ = 3d₀/(32κb₂), which implies (SQ); and (A1)
a² > 4α/d₀ if c12rat's rational block uses pc3 LAT(d). At 3/8, §2.7's FM-U′ line ("ϱ₂ = L^{−o(1)}") and Lemma η's
denominator step need (SQ) explicitly, because ϱ₂ = L^{−Θ(1)} at the endpoint.

*Proof.* As Theorem 2L-DW(α), with A ≪ D^{−2α}λ^C and B ≪ (D^{−2α} + D^{−2(3−7α)})λ^C. The budget (§3.6, "Consequence") holds on [3/8, 2/5), and the
remaining margins hold at 3/8 (`params.margins2`). The other inputs at 3/8 are uprime's U′ at 3/8 with the §2.7 DW
replacements, Lemma η (any α), Prop CLB and Lemma reject-DW. ∎

In particular, conditionally on (ii) and (iii), **L[3/8]**.

**Complete dependency list of Theorem 2L-DW(α).**
- *Hypothesis:* Hyp RAT(α) (open; c12rat), with the idealised weights.
- *Design constraint for c12rat (optional for us).* (A1) a > σ₂(q) for q near B, i.e. a² > 4α/d₀·(1 + o(1)), which
  makes Δ′ < q in every child box. pc3 LAT(d) and the rational mod-q ratio need it. The algebraic block does not:
  pairs with q | Δ′ are covered by Cor TP-q, since L_q/q = {y₁ ≡ λ^dy₂ (mod Δ/q)} holds even when q² | Δ. They are
  included in the BV averaging. v1's illustrative parameters (a = 3/5) violate (A1). With (A1), e.g. α = 7/18, a = 2,
  κ = 1/8, d₀ = 1, b₂ = 2, C_γ = 0 gives c ≈ 9.5. A rational grid search under (A1) gives c ≈ 3.91 at
  (a, κ, d₀, b₂) = (4/5, 1/16, 5/2, 1/2) (`params.cost_constant`; illustrative, C_γ not computable).
- *Sketch:* Theorem AB (this note; the bookkeeping in Lemmas FB and PP is not written to the last constant). For
  2L-DW₂, also Theorem AB₂ (§3.6), whose boundary bookkeeping transposes c12rat's Prop PD.
- *Proved here:* Lemma η, Prop PB, Lemma VD, Lemmas GB-DW, CU, RE, CT, IR-DW, reject-DW, W-DW, DW(a)⁺, MM, BZ, TP,
  TP-W, Cor TP-q, LF, AR, PU, TH, WD, Prop CLB, §2.7 (U′-DW), §2.8 (POSTPROCESS-DW), and the margins (`params.py`).
- *Internal (unrefereed):*
  - pc3: Prop DC, DC_e, Cor DC, Lemma DW-L, Lemma LAT;
  - uprime: §§3–6 (FM/LD*/D-U′, CS′, SQF, child scan, recursion and IR at γ₂), with the §2.7 replacements;
  - multidescent: Prop UB, Theorem 2L's proof;
  - nested: Lemmas RW and T_multi;
  - escapeR2: §7.3;
  - Theorem U's chain (lv5unif §10 list: [pp] §§3, 5–9; Kern-U, L2-U, Lemma S″, Lemmas OS, W, F, T_rec, Theorem DD,
    LV6-U).
- *Published:*
  - Harper arXiv 1208.5992 Theorems 1–2 (BV/BDH for smooth numbers; **unrefereed preprint**);
  - Hildebrand, J. Number Theory 22 (1986) 289–307. It is used at
    y = B (Lemma PP, child) and at y = B₂ (Lemma η, parent PP);
  - Dickman-function facts for Lemma η: log-concavity of ρ (Alladi) and the Hildebrand–Tenenbaum survey (JTNB 5 (1993)
    411–484), Cor 2.4 and Thm 5.1;
  - Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Theorems 1, 2, 3;
  - de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), Thm 2.4(i);
  - Canfield–Erdős–Pomerance;
  - Tenenbaum's Ψ_q theorem;
  - Iwaniec–Kowalski Thm 7.13;
  - Montgomery–Vaughan;
  - Lagarias–Odlyzko, Stark, Landau (LV6);
  - AKS, Wiedemann/Kaltofen–Saunders, Cantor–Zassenhaus, Bernstein;
  - Trudgian 2016 (CS′).

  Not needed: Lemma 8.5′/8.5H, which are replaced by DW(a). Drappeau's Lemme 5, if used as §3.6's alternative to
  Harper's Thm 2, carries the [H] dependence: its proof uses Harper 1208.5992 §3.3. Theorem U itself still uses them, for
  the recursive value factorisations (FACTOR_U).

**Achieved α:**
- none unconditionally;
- conditionally on RAT(α) and Theorem AB: every α > 3/8;
- conditionally on RAT₂ and Theorem AB₂: α = 3/8.

---