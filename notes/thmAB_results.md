# thmAB: Theorem AB written out in full (algebraic two-point blocks under DW) — results

Agent `thmAB`, phase 30, 2026-10-28. **v2**: v1 + the repairs from codex adversarial round `thmAB_2` (MAJOR REVISION,
no FATAL) and the Opus proof critic (ACCEPT WITH MINOR REVISIONS); §15 records the reviews.

Files:
- this note;
- code: `src/thmAB/params.py` (exact-rational exponent margins of every comparison below, 26 checks + the
  acceptability flag), `src/thmAB/cells.py` (exact checks of the finite statements: Lemmas CL, CF, PG, BL);
- tests: `tests/test_thmAB.py` (6 tests, exact arithmetic);
- codex: `codex/thmAB_1.{md,out}` (blind round, launched before anything here was written), `codex/thmAB_2.*`
  (adversarial round on v1).
- No experiments were run (none was needed: every finite claim is an exact identity or inequality, checked by the
  tests). No pre-registration was therefore filed.

Labels as in c12alg/c12rat: **proved** (from the project's internal, unrefereed chain plus cited published results;
anything resting on Drappeau's Lemme 5 rests, through it, on §3.3 of Harper's unrefereed preprint arXiv:1208.5992);
**sketch**; **heuristic**; **open**; **UNVERIFIED** (a literature statement not checked by us).

Notation: λ := log log N, D := log N/λ, L[α, c] = e^{(c+o(1))D^αλ}, "polylog" = e^{O(λ)}. ρ̂_y(x) := Ψ(x, y)/x; β = β(x, y)
is the Hildebrand–Tenenbaum saddle point and γ := 1 − β. To avoid the clashes between c12alg and c12rat:
- κ is the design constant (k = κD^α);
- ι is the unit of a pair lattice;
- ϑ(P₁, P₂) := |P₁||P₂|/|det(P₁, P₂)| is the thinness of a pair (Euclidean norms; c12rat's κ);
- 𝔐 = N^μ is c12alg's M (the size of g).

---

## 0. Verdict

**Plain language.** Theorem AB is now written out in full, and the proof gives more than was asked. "In full" means
every step cites a lemma proved in the project or a published theorem; it has had one adversarial round and one proof
critic, whose repairs are applied (v2), and it is not refereed. On the
algebraic side, the covariance of two child points (or two parents) is bounded by the **product** of two one-point
standard deviations, each of first order. So the algebraic blocks are of **second order**, on the whole range
3/8 ≤ α < 2/5. That is c12alg's Theorem AB₂, and its four open obligations O1, O2, O3, O5 are discharged here.

Combined with the proved pieces (U′ under DW, Prop DC, Theorem RAT, the DW lemmas, the child lower bound), this gives:
- **Theorem 2L-final(α):** a Las Vegas algorithm that factors every N in expected time L_N[α, c(α)], for every fixed
  α ∈ (3/8, 2/5);
- **Theorem 2L-final(3/8):** the same at α = 3/8, under the design constraint a > 2c₂ (c12alg's O4).

Status: **proved, internal and unrefereed**, like everything in this project. It rests on the project's whole
internal chain, on c12rat's Theorem RAT, and on published inputs, of which the essential one is Drappeau's Lemme 5,
which itself rests on Harper's unrefereed preprint. The v2 repairs have been checked by one further codex round
(§15).

How the proof works:
- **The pair law is explicit.** Given the rational data g and the root event, (F(P₁), F(P₂))/q is within total
  variation 7qN^{−η₀/2} of a mixture of uniform measures on L′ ∩ Π_y. Here L′ = {x₁ ≡ ιx₂ mod Δ′} and the Π_y are
  translates of one parallelogram (c12alg Cor TP-q, proved). There is no primitivity constraint, and the modulus is
  the small coordinate determinant Δ′, not qΔ′.
- **Cells.** Cut the plane into rectangles (one block per coordinate) and gcd classes δ = gcd(x₁, Δ′) = gcd(x₂, Δ′).
  In a cell, the lattice condition becomes m₁ ≡ ιm₂ mod Δ′/δ on units. c12rat's coupled-pair Lemma CP then factorises
  the smooth-pair count exactly into the product of the two one-point cell densities, plus a character error.
- **Cauchy–Schwarz.** The pair covariance becomes the covariance of the two cell densities under a cell measure. That
  covariance is at most the product of the two one-point standard deviations, each of relative size
  D^{−α}λ^C (children) or (D^{−α} + D^{7α−3})λ^C (parents). First-order errors therefore enter only squared.
- **Character errors.** They are summed over the whole incidence-weighted family with c12rat's Theorem WA and
  Drappeau's Lemme 5. For the algebraic block this bookkeeping is simpler than c12rat's: the modulus Δ′/δ does not
  involve the child prime q, so no weight domination across residues is needed.
- **Thin pairs, skew labels, parents.** Thin pairs are bounded by a pair upper bound with only polylog algebraic loss.
  Skewed labels are negligible by weight (c12rat Lemma DOM). Parents use their own thresholds, built from parent
  losses only.

What remains, honestly:
- the review of this note (§15);
- the v2 repairs of c12rat's Theorem RAT, which have not been re-reviewed;
- the internal chain, which no outside referee has read;
- C_γ (hence c(α)) is not computable, as for Theorem U.

**Probabilities (ours, after the v2 repairs; reviewers' in §15).**
- Theorem AB as proved here is correct: 75%. [codex round 2: 80% after repairs; proof critic: 80%]
- Theorem AB₂ (second order on [3/8, 2/5)): 68%. [codex 70%; critic 72%]
- Theorem 2L-final(α), α ∈ (3/8, 2/5), correct relative to the internal chain and Theorem RAT: 65%. [codex 70%;
  critic 65%]
- Theorem 2L-final(3/8), likewise: 52%. [codex 60%; critic 52%]

The residual risk is mostly in the glue (c12rat's unreviewed v2 repairs, c12alg §2.7/§5 bookkeeping, U′ under DW at
η₀ = D^{−3/8}), not in the new algebraic estimates.

---

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
sum below; they enter only through the normalisation, identically in numerators and denominators, so all ratios are
unaffected and absolute terms change by a factor 1 + O(1/q) (codex thmAB_2 #7).

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
  D^{10}(d + 1)²y^{−C_X} ≤ y^{−1}c₁c₂ (codex thmAB_2 #8);
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

  L′ is the divided lattice L_q/q of Cor TP-q (valid also when q | Δ′, c12alg v2.2). ι is a unit mod Δ′ because
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
  point and hence a library endpoint (critic item 2). A block of coordinate i is a union of consecutive grid intervals,
  of relative width ≤ min(h^{1/2}, 1/8) and length ≤ 2z_i. The blocks partition the rounded range
  {X_i⁻ ≤ |x_i| ≤ 2X_i}, and this partition depends only on (g, pair), not on y. Every block contains a whole grid
  interval, so its relative width lies in [h_*/2, h^{1/2}]. (v1 also split blocks into "ordinary" and "non-ordinary";
  that split was stated backwards — away from 0 the length cap binds — and is removed in v2: codex thmAB_2 #1,
  critic 4b.)
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
     ≤ C_X/log X₁ (as |ξ| ≥ X₁y^{−C_X}), and Σ_{p|Δ′}p^{−β}log p ≤ Cλ. So the factor is 1 + O(D^{−2α}) (critic 4a).
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
   |e| ≤ CE_FO + Cγh^{1/2} + CD^{−2α}. Cauchy–Schwarz with the fourth moments gives E[(Ψ₁/c₁ − 1)²1_{non-exc}]
   ≤ C(γ²(log Cd)² + γ²λ² + E_FO² + u^{−2}). Hölder gives E[(Ψ₁/c₁)⁸1_{non-exc}] ≤ C(EY^{16}Eτ^{16})^{1/2} ≤ C.
   (v1 used an e^{O(λ)} pointwise bound on a set it believed to have mass Cdh^{1/2}; that step is gone.)
5. **Pointwise upper bound** (every block, every δ, used only on exceptional sets and in §9). Ψ₁ ≤ e^{O(λ)}δ^{γ}(X₁/|ξ|)^γc₁.
   Indeed:
   - c12rat **Lemma BTS** with modulus 1, applied to I/δ (length ≥ 2y), gives #smooth ≤ e^{O(λ)}(|I|/δ)·(…)^γρ̂;
   - Fact dBT1 gives ρ̂(x/δ) ≤ Cδ^γρ̂(x);
   - U ≥ (|I|/δ)/(Cλ);
   - h_*^{−γ} = e^{o(1)}.
6. **p₁ versus c₁** (critic item 5; codex thmAB_2 #6 — this breaks the v1 loop OPV/EXC/PU).
   - *Upper:* p₁ ≤ max fibre mass × #{x₁ ∈ 𝒮_y : |x₁| ≤ X₁} ≤ (3/W₁)·2Ψ(X₁, y) ≤ 12(d + 1)c₁, by PG(d).
   - *Lower:* p₁ ≥ ν(x₁ ∈ 𝒮_y, x₁ non-exceptional) = c₁E[Yτ(1 + e)1_{non-exc}]. Now E[Yτ] ≥ (1 − O(1/u))Eτ,
     E[τ1_{non-exc}] ≥ 1 − o(1) − E[τ1_{exc}], and E[τ1_{exc}] ≤ (ν(exc))^{1/2}(Eτ²)^{1/2} = o(1) by EXC-A's mass bound
     (which uses only PG(d) and the gcd tails). Hence p₁ ≥ c₁/2.
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
  ≤ y^{−1}c₁c₂ ≤ e^{−λ}c₁c₂ (y = B for children, B₂ for parents; codex thmAB_2 #8: v1's (c₁c₂)⁸ was unjustified and
  is not needed).
- **Large classes.** δ | x₁ forces δ | x₂ (Lemma CL). For fixed δ | Δ′:
  - ν(δ | x₁, x₁ ∈ 𝒮_y) ≤ (3/W₁)·2Ψ(X₁/δ, y) ≤ Cdδ^{−1+γ}c₁, by PG(d) and Fact dBT1.
  - Given x₁, the partners x₂ form an AP of difference Δ′ on a fibre, an interval of length ℓ(x₁). **First** split
    the fibre into its intersections with the dyadic ranges ±[2^kX₂⁻, 2^{k+1}X₂⁻) (O(log X₂) pieces; the part
    |x₂| < X₂⁻ is covered by the near-zero bound). **Then** cover each nonempty piece by an interval J of length
    max(|piece|, λ_min), λ_min := X₂y^{−C_X}/ϑ ≤ X₂⁻, inside ±[2^{k−1}X₂⁻, 2^{k+2}X₂⁻) (codex thmAB_2 #4, critic 5).
    On J, c12rat **Lemma BTS** applies: modulus Δ′, W = |J|/Δ′ ≥ λ_min/Δ′ ≥ 2y (λ_min ≥ N^{h}/(V₂²S₂²y^{C_X})), and
    gcd(a, Δ′)W ≤ |J| ≤ 2^{k}X₂⁻ ≤ (left end of J)·2. It bounds the smooth partners in the piece by
    e^{O(λ)}(|J|/Δ′)(2^{k+2}X₂⁻Δ′/|J|)^{γ}ρ̂(2^{k+2}X₂⁻) ≤ e^{O(λ)}(y^{C_X}ϑΔ′)^γ(|J|/Δ′)c₂, using
    ρ̂(2^{k+2}X₂⁻) ≤ C(X₂/X₂⁻)^γρ̂(X₂) ≤ e^{O(λ)}c₂.
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
   (Lemma BL(ii)). Then, by Hölder on the joint space (codex thmAB_2 #5: not per component):
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

  There is no pointwise claim about Rat₁Rat₂ in g; only its integral is used (codex blind §5 makes the same point).
- Step 3: **normalisation.** For s(𝔮) ≤ T₀, c12rat **Lemma SPR** (children) gives max_Pϖ_P ≤ Ξ min_Pϖ_P over Box(𝔮).
  This covers the spread of ϱ_B(𝔐|P|) by (4S₂T₀)^γ, and of the algebraic factor by ((4S₂T₀)^d)^γ =
  e^{O(D^{1−2α}λ² + D^{1−3α+θ}λ³)} ≤ Ξ. Prop CLB, per point (c12alg §4 steps 2–3, which bound each P′'s contribution
  by ϱ(𝔐|P′|)ϱ(X_{F2}(P′))(log N)^{−C}), gives m_{𝔮,e} ≥ Σ_Pϖ_P(log N)^{−C}. Hence ϖ_{P₁}ϖ_{P₂}/m² ≤ Ξ³/n², with n = #Box(𝔮) ≥ S₂².
- Step 4: **multiplicity.** For each Δ′ ≥ 1 the box has ≤ 10S₂² + 2S₂ ordered pairs with coordinate determinant Δ′
  (c12rat **Lemma MULT**), and 𝔨 ≤ Δ′/δ ≤ 2S₂²/δ. So φ(𝔨)·#{pairs : Δ′ = 𝔨δ}/n² ≤ 24/δ.
- Step 5: **sum.** Since max_𝔨 Σ_{(𝔮,e)} ≤ Σ_{(𝔮,e)} max_𝔨 (legitimate because steps 1–4 give the domination label by
  label; codex thmAB_2),
  𝔇 ≤ Σ_{(𝔮,e)}(w/W)Ξ^Ce^{O(λ)}dh^{1/2}λ^Ch_*^{−2}·Σ_δ24δ^{−2+2γ}·#(𝒳 ∩ ℛ_{𝔮,e}),
  where ℛ_{𝔮,e} is the **label-wide** scale range: the union over the label's pairs and over g ∈ Ω_e of
  [X⁻_min/2, 4X_max] (critic item 7). Its logarithmic length is O(d log(S₂T₀) + η₀ log N + C_X log B) = D^{O(1)}λ^{O(1)}
  (|P|^{d−1} varies over a box with s(𝔮) ≤ T₀ by (4S₂T₀)^{d}). Here Σ_δδ^{−2+2γ} ≤ 2 and
  #(𝒳 ∩ ℛ) ≤ D^{O(1)}/h_*. Also Σ_{(𝔮,e)}w/W = 1. So 𝔇 ≤ Ξ^Ch_*^{−3}D^{C} = e^{O(D^θλ²)}, since
  h_*^{−1} = 64(d + 1)ϑ₀D^{10} = e^{O(D^θλ²)}. ∎

Why this is simpler than c12rat's Lemma FAM: the modulus 𝔨 = Δ′/δ is a **coordinate** quantity (≤ 2S₂²) that does
not contain q. So for fixed 𝔨 every label can contribute, and the crude bound max_𝔨 Σ_𝔮 ≤ Σ_𝔮 max_𝔨 already
suffices. No comparison of label weights across residues r or primes q is needed (Lemma DOM is not used here), and
the g-dependence of the windows costs only the logarithmic length of the scale range (step 1).

**Parents.** The same, with y = B₂, the grid 𝒳^p and h_*^p, 𝔨 ≤ 2V₁²/(qδ), and the following normalisation package
(v2; codex thmAB_2 #2: v1 used c12rat's S_typ, which already contains an algebraic factor, and would lose ρ_a(q)²).
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
  - a_P(g) ≥ ρ_a(q)D^{−C}: DW(a)⁺ with the root coset, Hildebrand at B₂ on intervals of relative length ≥ 1/(2dq), and
    the root-conditioned acceptability subtraction of §11 (≤ D^C/B ≪ ρ_a(q); at α = 3/8 this uses (SQ)).
  Sum over active P. Inactive points are only bounded above; no individual lower bound is claimed for them. ∎

Multiplicity: c12rat Lemma MULT (parent) gives ≤ 5 partners P₂ of each P₁ with a given **signed** det, so ≤ 10 for a
given |det| (codex thmAB_2), and φ(𝔨)·#{pairs : Δ′ = 𝔨δ}/n_p² ≤ 20V₁²/(qδn_p) ≤ 80d log V₁/δ, as
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
up to (S₂q)^{d−1}, a factor e^{O(D^αλ²)}. The per-point form avoids it. Codex's blind round found the same repair.

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
    D at y = B₂, because ω(Δ′)P₀^{−0.49} ≤ D^{(1−α) − 0.343} and (1 − α) − 0.343 < 4α − 1 (`params`: "parent FO tail";
    v1's params used P₀^{0.49} = D^{0.49} instead of D^{0.343}, corrected; critic 4a).
- **Moduli** Δ′ = |det(P₁, P₂)|/q ≤ 2V₁²/q = L[1 − α]^{O(1)}, which is ≫ B₂^{O(1)} but x^{o(1)}. Nothing in §§5–9 iterates
  HT3 over divisors of Δ′: classes δ ≤ S₀^p use Lemma FO, larger ones use Lemma EXC-A. So c12alg's truncation at
  B₂^{1/2} is not needed.
- **Thresholds** ϑ₀^p and h_*^p, built from the parent losses Ξ_p only (O3). The child filter spread
  e^{O(D^{2α−1/2}λ²)} never occurs: the parent filter fixes |G(P)| to within a factor 64.
- **No skew:** parent sets are Euclidean squares.
- **Acceptability (O4).** acc(P) is not a function of the value x alone. Put s_P(g) := Pr[F(P)/q ∈ 𝒮_{B₂} | g, R_𝔮]
  and s₁₂ likewise. **Root-conditioned subtraction** (v2; codex thmAB_2 #3 found v1's displayed bound false for p = q).
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
  ΣR₁R₂a₁a₂ ≤ z(g)², this involves no Ξ_p (critic item 8).

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

| obligation | where | status (v2) |
|---|---|---|
| O1 joint moment lemma | Lemmas LF-A, OPV-A (v2: FO on every block), EXC-A (incl. ν⁺, mixture, prime powers, large classes, near zero) | proved here; critic: "genuinely proved" after fix 4b (applied); codex round 2 kept it open pending the same fix |
| O2 weighted energy inequality | Lemma EN-A, Lemma FAM-A + c12rat Theorem WA; parents via Lemma W-P (v2) | proved here; parent normalisation rewritten in v2 (codex #2) |
| O3 separate parent thresholds | §1.3 (ϑ₀^p, h_*^p from Ξ_p), §9 parents, §11 | proved; both reviewers agree |
| O4 parent acceptability at 3/8 | §11 (a > 2c₂; v2 root-conditioned subtraction) | discharged by the constraint |
| O5 boundary transposition | Lemma BL(ii) (convex Π^∘, Π⁺), Prop PD-A step 4 (v2: Hölder on the joint (y, x) space) | proved; both reviewers accept |

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

   **Weight bridge to RAT** (critic item 1). c12rat states RAT with the idealised weights a_{P′}(|G_g(P′)|) (c12alg
   §5), while Lemma VD's ϕ^rat uses the true A_{P′}(g). By DW(a)⁺, |A_{P′}(g) − a_{P′}(|G_g(P′)|)| ≤ 4qN^{−η₀/2} for
   every g, which is ≤ N^{−η₀/3}a_{P′}. Hence m(g) = m̃(g)(1 + O(N^{−η₀/3})) pointwise, and
   ϕ^rat_true ≤ (1 + N^{−η₀/4})ϕ^rat_ideal + N^{−η₀/4}; likewise for σ^rat.

   **(A1) is not needed** (critic item 9). c12alg §5 made a² > 4α/d₀ conditional on c12rat using pc3 LAT(d). c12rat's
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
  with the §2.7 DW replacements. **Re-check at DW's η₀ = D^{−3/8}** (uprime used η₀ = D^{−1/2}; critic item 9),
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
saddle-point estimates and of Drappeau's Lemme 5, as for Theorem U (codex thmAB_2 #9: "not quantified", not
"not computable"). Illustrations with C_γ = 0 (`src/c12alg/params.cost_constant`):
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

**c12rat (internal; v2 repairs not re-reviewed):**
- Theorem RAT and Corollary RAT(α);
- Lemmas used directly here: CP, FO, BTS, PUB (rat and parent-absolute), SPR, MULT, TF, LS, DOM, Theorem WA,
  Lemma EXC(b), CELL(i)–(ii) (the convexity arguments, re-proved here for parallelograms).

**pc3:** Lemma 0, Prop DC, Prop DC_e, Cor DC, Lemma LAT, Lemma DW-L.

**uprime:** Lemma 8.5H, GM/FM/DS/L2/LD*/D-U′, CS′, SQF, the child scan, recursion and IR at γ₂ (with c12alg's
DW replacements; Lemma 8.5H is not needed at level 1 under DW, but Theorem U's FACTOR_U still uses Lemma 8.5′ for the
recursive value factorisations).

**multidescent:** Prop UB, Theorem 2L's proof (search, rows, LV6, fruitfulness, recursion).

**nested:** Lemmas RW, T_multi. **escapeR2:** §7.3 (bounded-depth search recurrence).

**Theorem U's chain** (lv5unif §10 list): [pp] §§3, 5–9 (hypD, pcest, ldbound, lv6, thmA), Lemma 8.5′(R′) (lem85),
Lemma S″, Lemmas OS, W, F, T_rec (recurse), Theorem DD (lemmaM), and lv5unif's Lemmas P, U1, CS, Kern-U, GM-U, W-U,
DS-U, L2-U, Theorems FM-U, LD*-U, D-U, K-U, Prop select-U, Prop LV6-U, Lemma F_α.

**Published inputs.** Theorem numbers are as checked by the project's literature critics (c12alg §10, c12rat §12.1,
uprime §10, lv5unif §14); new ones here are flagged.
- **Drappeau**, Compositio Math. 151 (2015) 828–862, **Lemme 5** (Fact D; A = 1/2, effective for A < 1). It is stated
  without proof as a consequence of §3.3 of **Harper's unrefereed preprint arXiv:1208.5992**. **This is the one
  input whose justification rests on an unrefereed source.** Also **Lemme 1(i)**.
- **de la Bretèche–Tenenbaum**, Ramanujan J. 9 (2005) 139–202, **Thm 2.4(i)** (Fact dBT1; upper bound only).
- **Hildebrand**, J. Number Theory 22 (1986) 289–307, **Theorems 1 and 3**: Ψ = xρ(u)(1 + O(log(u + 1)/log y)), and
  the short-interval form uniformly for xy^{−5/12} ≤ z ≤ x, y ≥ exp((log log x)^{5/3+ε}). These are restated in
  S. Jain, arXiv:2502.10530 (QJM 2026), Lemma 3.2, and in the HT survey JTNB 5 (1993) Thm 5.1. The primary was not
  opened by the project's critics; the statement is verified against both restatements.
- **Hildebrand–Tenenbaum**, Trans. AMS 296 (1986) 265–290, **Thms 1, 2(i), 2(ii), 3**.
- **Hildebrand–Tenenbaum survey**, JTNB 5 (1993) 411–484, **Cor 2.4** and **Thm 5.1** (Lemma η).
- **Alladi** (log-concavity of ρ; Lemma η).
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

## 15. Reviews

**Codex adversarial `thmAB_2`** (xhigh, on v1): **MAJOR REVISION, no FATAL.** "I found no FATAL obstruction to the
proposed method." It confirmed ML (for smoothness events), PG, the boundary geometry, CL, CF, FO and LF-A (incl. prime
powers and primes > y), EN-A against Drappeau's Lemme 5, child FAM-A ("the decisive summation is
Σ_δδ^{−2+2γ} = O(1)"), thin pairs, skew, the separate parent thresholds, and the assembly exponents at 3/8.

| # | sev. | finding | action (v2) |
|---|---|---|---|
| 1 | MAJOR | ordinary/non-ordinary mass claim reversed; OPV-A steps 5–6 fail | split removed; FO on every block; 8th moment by Hölder (§6) |
| 2 | MAJOR | parent FAM normalisation: S_typ already contains an algebraic factor (loses ρ_a(q)²) | Lemma W-P with rational R_typ and ρ_a(q), (P1)–(P3) proved (§9) |
| 3 | MAJOR | root-conditioned SQF false for p = q; discriminant event needs q′ \| F(P) too | subtraction rewritten case by case, total D^C/B; ML stated for smoothness events (§§3, 11) |
| 4 | MAJOR | fibre enlarged before dyadic splitting: pieces may be too short for BTS | split first, then cover pieces; thinness cancels (§7) |
| 5 | MINOR | boundary fourth moments per component vs mixture | Hölder on the joint (y, x) space; Z_y kept (§8) |
| 6 | MINOR | p_i comparison and the OPV/EXC/PU loop | p_i ≤ 12(d + 1)c_i from PG(d) alone; lower bound via non-exceptional mass (§6 step 6) |
| 7 | MINOR | R_𝔮 does not remove q \| G(r, 1) | wording; zero contribution, normalisation 1 + O(1/q) (§1.2) |
| 8 | MINOR | cutoff C_X vs the (c₁c₂)⁸ claim | C_X redefined for the needed y^{−1}c₁c₂ (§1.3, §7) |
| 9 | MINOR | "written out completely", "discharged", "not computable" overstate | wording (§§0, 12, 13) |

Probabilities (after the repairs): AB 80%, AB₂ 70%, 2L-final(α) 70%, 2L-final(3/8) 60%.

**Opus proof critic** (read-only, on v1; also asked the triage points T1–T4): **ACCEPT WITH MINOR REVISIONS**, no
FATAL, no MAJOR in the argument itself. "Transposing the device to the algebraic side is sound … simpler than in
c12rat." Required fixes, all applied:
1. the ordinary/non-ordinary step (its 4b, = codex #1): the only place where a polylog multiplied something merely
   D-polynomially small;
2. snap the cutoffs to grid points (§4 Lemma BL);
3. label-wide grid range in FAM-A step 5 (§9);
4. PG(e) needs + φ(Δ′/δ)N^{−Ω(μ)};
5. BL(ii) area accounting ("Ch");
6. enlarge dyadic pieces, not fibres (= codex #4);
7. the idealised-to-true weight bridge to RAT (§13);
8. settle (A1) (not needed; §13) and the P₀^{0.49} exponent in `params` (fixed).
Its triage answers: **(T1)** covered — non-reduced classes exactly via Lemma CL, endpoints divided by δ are the library
u_{ξ,δ}, the actual incidence weights are used with upper bounds and g-integrals only; **(T2)** yes — every modulus is
x^{o(1)} (children D^{α}λ vs D^{2α}λ; parents D^{1−α}λ vs D^{2−3α}λ), so Lemme 5 alone suffices; **(T3)** no unpaid
loss — modulus φ(𝔨) paid by MULT (net 24/δ), no τ(Δ′) since (δ, 𝔨) determines Δ′, Σδ^{−2+2γ} ≤ 2, grid count
e^{O(D^θλ²)} (understated in v1, fixed); **(T4)** O1 is genuinely proved by LF-A/OPV-A/EXC-A after fix 4b. It
disagreed with codex on the parent normalisation (it accepted c12rat's S_typ); v2 follows codex's stricter form.
Probabilities: AB 80%, AB₂ 72%, 2L-final(α) 65%, 2L-final(3/8) 52%; residual risk in the glue (RAT interface, (A1),
c12alg §2.7/§5) and, at 3/8, U′-DW at η₀ = D^{−3/8} and the exact ties θ = 1 − 2α, 3 − 7α = α closed only by
constants.

**Codex round 3 `thmAB_3`** (on v2): *(pending)*

**Opus literature critic**: *(pending)*

**Codex blind `thmAB_1`** (xhigh, launched before anything here was written). It chose the same route (b) for the
same two reasons: no primitivity, and modulus Δ′ instead of qΔ′. It reached the same structure:
- exact cell coupling;
- the joint moment lemma (M), with prime powers kept exactly;
- a weighted energy inequality (W), which "never assumes independence of Π_z and Rat₁Rat₂";
- separate thresholds;
- per-point skew bounds;
- acceptability via SQF with a > 2c₂ at 3/8.

Its probabilities: AB 80%, AB₂ 65%, end-to-end L[α] (α > 3/8) given RAT 70%. It flagged (W) for parents as the
estimate to attack hardest. Differences: it truncates Möbius at y^{1/4} instead of our S₀ (both work); it states its
thresholds in Λ_c, Λ_p (our Ξ, Ξ_p).

**Codex triage `r29_triage`** (relayed by the orchestrator). It ranked this lead first among new-theorem leads (30%
for a complete positive result in two months of specialist effort). It advised "close first-order AB before AB₂" and
warned that "covariance is second order" does not by itself discharge O1. The cell route discharges O1 by Lemma OPV-A,
and AB is obtained as a corollary of AB₂, not the other way round.

---

## 16. Reproduction

```
.venv/bin/python -m pytest -q tests/test_thmAB.py
.venv/bin/python -c "from fractions import Fraction as F; from src.thmAB import params as P; print(P.margins(F(3,8))); print(P.all_ok(F(7,18)))"
```

## 17. Candidate card

- **Mechanism.** Two-level descent, α ∈ [3/8, 2/5), design DW. The algebraic two-point blocks are bounded by the cell
  covariance device: exact coupling on the DW lattice L′ = {x₁ ≡ ιx₂ mod Δ′}, Cauchy–Schwarz against one-point cell
  variances, and character errors through Drappeau's Lemme 5 with a coordinate-modulus weighted energy inequality.
- **Result.**
  - Proved (internal): Theorem AB₂, hence AB.
  - Hence Theorem 2L-final(α): L_N[α, c(α)] for every fixed α ∈ (3/8, 2/5), and L_N[3/8, c] under a > 2c₂.
  - Status: conditional on the internal chain and c12rat's Theorem RAT, and pending this note's review.
- **Failure scope.** No claim below 3/8. c(α) is not numerically certified (C_γ).
- **Next step.** Independent review of §§6–9, and of c12rat v2 (Theorem RAT) jointly with this note, since both now
  carry the headline.
