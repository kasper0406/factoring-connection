# pc3: Hyp PC₃ (multidescent §3.2–3.4) — results

Agent `pc3`, phase 26, 2026-10-24. **v2.** It applies the codex adversarial round (MAJOR REVISION), the Opus proof
critic and uprime's results; the Opus literature critic is in §11.

Files:
- this note;
- `notes/pc3_prereg.md` (commit `1025bda`; code frozen there);
- code: `src/pc3/toy.py`, `experiments.py`, `analyze.py`, `dc.py`;
- tests: `tests/test_pc3.py` (7 tests);
- data: `experiments/raw/pc3_e{1,2}.jsonl` (+ `_summary.json`);
- codex: `codex/pc3_1.{md,out}` (blind round, launched before any analysis here) and `codex/pc3_2.{md,out}`
  (adversarial round on v1).

Labels:
- **proved**: from the project's internal, unrefereed chain plus the cited published results.
- **sketch**: the steps are named and the inputs identified, but the argument is not written. Reviewers found gaps in
  these, listed where they occur.
- **heuristic**, **open**, **experimentally supported** (toy), **UNVERIFIED** (literature statement not checked by us).

---

## 0. Verdict

**Plain language.** We did not prove Hyp PC₃, and so we have no rigorous algorithm below L[2/5]. But PC₃ turned out to
be a stronger statement than the two-level descent needs.

What multidescent's union bound (Prop UB) actually needs is that, summed over all parent–child incidences, the child's
search box is rarely empty. By an elementary inequality (Prop DC), this follows from two concentration statements that
never mention the parent and the child together:
- **(C₁)** the number of good points in a child's box varies little from field to field, *without conditioning on any
  parent*;
- **(C₂)** the number of potential parents sharing that child's label varies little.

Each is a statement about **pairs** of points: pairs of child points, and pairs of parents. The hard object multidescent
identified, a ternary equation in smooth numbers with coefficients up to x^{o(1)} and conditional weights (Hyp RB), is
no longer needed. Neither is uprime's Hyp CM, the lower bound on the parent-conditioned child mean.

Why this helps:
- The precision both statements need is only "first order": relative variance ≤ c·D^{−3/8}.
- Their natural size is "second order", ≈ D^{−3/4}.
- Their rational parts are smooth-number pair counts on explicit lattices. That is the territory of
  Barban–Davenport–Halberstam (BDH) estimates for smooth numbers (Harper) and of the large sieve over the child primes.

**We do not prove C₁ or C₂.**
- Our proof sketches failed review in specific places: the weighted application of BDH and the large sieve, skewed child
  lattices, the archimedean term for nearly parallel pairs, and two-point smoothing (§7).
- We propose a design change, the integer multiplier DW. It turns the algebraic two-point block into a lattice problem
  of the same kind, but it is unchecked against the internal chain.

**Proved.**
- **Theorem 2L-DC.** With uprime's proof of U′, the two-level algorithm runs in L_N[3/8, c] *conditional on C₁ and C₂*
  (in place of PC₃) and on the internal chain.
- The tools:
  - Lemma LAT (pair lattices);
  - Lemma LSQ (large sieve over child primes);
  - Lemma DW-L (the exact two-point support lattice under DW, found by codex);
  - the scaling-invariance cancellation;
  - Lemma DW(a) (the one-point algebraic law under DW is pointwise in g, for filtered points).

Precisely:
1. **Prop DC and Theorem 2L-DC [proved].** Let ϕ_𝔮 and σ_𝔮² be the root-conditioned relative variances of the child
   count and of the label's parent count. Let A := E_w ϕ and B := E_w σ² be their averages under the incidence weights
   w_𝔮. Then
   - Σ_{(ℓ,P,q)} Pr[ℓ, E ∧ N_𝔮 = 0] ≤ W(A + √(AB)).
   - Hence Prop UB's conclusion holds when A ≤ η/4 and B ≤ η/4, with η := μ₁/W.
   - An environment-conditioned version (Prop DC_e) removes constant-size field-level fluctuations, such as the
     g-dependent parent filter (proof critic M1).
   - Codex found the Cauchy–Schwarz inequality independently in its blind round.
2. **Heuristic sizes [heuristic].** A, B ≈ D^{−3/4} up to constants and logarithms; the requirement is ≍ D^{−3/8}.
   The margin is D^{3/8}, or D^{1/4} if the deterministic child cap D^{1/2} replaces the typical count D^{3/8}.
   - CM is not needed: the 1/m term of C₁ involves the *unconditional* (root-conditioned) child mean, a one-point
     quantity. uprime's child scan gives nominal m₂ ≥ exp((log log N)²).
3. **Rational two-point blocks [tools proved; application sketch].**
   - Lemma LAT: each pair lies on the lattice n₁ ≡ λn₂ (mod D).
   - The cross terms of the gcd-class decomposition vanish exactly, by scaling invariance.
   - The mod-q part is controlled on average over children by LSQ, because q ≥ B = L[3/8, a] ≫ 1/ρ = L[3/8, O(1/a)].
   - The mod-Δ part is controlled by BDH for smooth numbers (Harper 1208.5992 Thm 2) on average over the pair family.
   - The principal part needs Saias's smooth approximation of Ψ (proof critic M3) and [pp] Lemma DS for thin pairs.
   - Open in the assembly: the weighted application theorem (codex #1) and skewed labels (codex #4).
4. **Algebraic two-point blocks [design DW; DW-L, DW(a) proved; DW(b) sketch].**
   - With F = w·F̂_{g,N} + G·R_c and w a uniform integer, the pair (F(P₁), F(P₂)) ranges exactly over the lattice
     {F₁ ≡ λ^dF₂ (mod |det(P₁, P₂)|)}. The multiplier order is |G(P₁)G(P₂)| (codex).
   - So the algebraic block is again a det-modulus lattice problem. No BV over the endogenous moduli G(P′) is needed,
     and no τ-weighted product family (uprime's trap 6).
   - Missing: a two-point smoothing (TV) lemma, and DW's compatibility with Lemmas GB/coprime/IR and with POSTPROCESS
     (codex #7: composite w breaks three stated steps, all repairably).
5. **Deeper descents [model].** At depth 2 of a three-level chain (α = 4/11), DC is borderline (margin 0), like the
   direct route. Below 3/8 nothing changes.
6. **Experiments (pre-registered; toy; rational block only): 9/9 predictions held.**
   - E1: per-pair covariances in a child box are as large as ±0.3 at toy scale. A local × archimedean proxy tracks them
     (slope 1.18), and they cancel over the box to 5·10⁻⁴.
   - E2: the parent-conditioned relative variance (PC₃'s object) equals the unconditional one (C₁'s object) to 0.0014,
     at every parent coefficient scale K ∈ {20, 80, 320}. The label's parent count is nearly Poisson.

**Probability** that C₁ and C₂ are provable by known methods in a focused effort, in design 2L′ with DW:
- codex: 35% (that the outlined route can be repaired) and 80% (that DW survives the internal chain);
- ours: 25–35%.

This is higher than multidescent's 5–15% for PC₃, because no ternary problem remains. The limiting items are §7.1–7.4.

---

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
- Finally, Chebyshev gives p ≤ ϕ. ∎ (Codex blind §3; proof critic ACCEPT.)

**Corollary DC (Cauchy–Schwarz over labels; proof critic M2, codex #8).** Put A := Σ(w_𝔮/W)ϕ_𝔮 and
B := Σ(w_𝔮/W)σ_𝔮². Then the dead-incidence mass is ≤ W(A + √(AB)). In particular, A ≤ η/4 and B ≤ η/4 imply
E N₂ ≥ μ₁/2.

**Proposition DC_e (environment conditioning) [proved].** Let 𝒢 be any σ-algebra generated by a function e(A) of the
field (an "environment cell") such that N_𝔮 and Z_𝔮 are defined cell by cell. Then the same bound holds with
Σ_{𝔮,e}w_{𝔮,e}(p_{𝔮,e} + σ_{𝔮,e}√p_{𝔮,e}). Here every quantity is conditional on R_𝔮 ∩ {e(A) = e}, and
w_{𝔮,e} := Pr(R_𝔮, e)z_{𝔮,e}/#J. Cauchy–Schwarz then runs over (𝔮, e).
*Proof.* Apply the proof above under each conditional law and sum. ∎
- *Use.* The proof critic found an O(1) field-level term in σ_𝔮². The parent filter |G(P)| ≥ MV₁/32 has a pass fraction
  depending on g's direction (it varies by ≈ 1/32), so E[Z_𝔮 | g]/z_𝔮 is not concentrated, and σ² ≥ c₀ > 0.
  - With an O(1) σ², Corollary DC needs A ≲ η²/c₀ ≍ D^{−3/4}, which has no margin.
  - Taking e(A) := (cell of g of relative size D^{−A₀} in direction and scale) makes the pass fraction constant to
    O(D^{−A₀}) within each cell. C₁ and C₂ are then needed cell by cell.
  - The BDH/LSQ arguments are local in g, since a cell contains M²D^{−2A₀} values of g.
  - Any further thresholds should be randomised **independently per point**. A per-field random child-box radius would
    make ϕ = Θ(1) (proof critic M1).
- *Alternative (proof critic).* Hölder gives E[Z1{N=0}] ≤ ‖Z‖₃p^{2/3}. With ϕ ≈ D^{−3/4} log D, a bounded third relative
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
vectors are primitive, so apply (a) to them (proof critic, removing v1's caveat). Here Δ ≤ 2S₂² < q because
S₂² = L[3/8, σ₂] < B ≤ q once a is large. ∎ (`tests/test_pc3.py::test_pair_lattice_is_unit_graph`; proof critic, 2000
random pairs.)

**Lemma SI (scaling invariance) [proved].** If Ω̄ ⊂ (Z/K)² is a subgroup, the multiset {gcd(a₁, K) : (a₁, a₂) ∈ Ω̄}
depends only on gcd(a₂, K). Consequently, in the class decomposition
  Σ_{(a₁,a₂)∈Ω̄} S₁(a₁)S₂(a₂),  S = M(gcd class) + E,
the cross terms Σ M₁E₂ and Σ E₁M₂ vanish exactly. Here M is the average over *all* residues of the gcd class, and both
projections of Ω̄ must be surjective; they are for LAT(b) and DW-L.
*Proof.* Residues with equal gcd differ by a unit factor, and multiplying both coordinates by that unit is a bijection of
the fibres preserving first-coordinate gcd's. ∎ (`test_scaling_invariance_cross_terms`; proof critic, 3000 random
subgroups.)

**Pair multiplicities (proof critic m5) [proved].** For a fixed primitive coordinate vector (i, j), the number of
(i′, j′) in the box with |ij′ − ji′| = Δ is ≤ 2S₂/|(i, j)|_∞ + 2. So #{pairs with determinant quotient Δ} ≪ S₂²
uniformly in Δ. There is no divisor function: a τ(Δ) weight would have been exp(O(aD^{3/8})), which swamps BDH's saving.

## 4. Analytic tools

**4.1 BDH for smooth numbers [cited; codex adversarial confirmed the statement].** Harper, arXiv 1208.5992, Thm 2: for
(log x)^K ≤ y ≤ x and Q ≤ Ψ(x, y),
  Σ_{Δ≤Q} Σ_{(a,Δ)=1} |Ψ(x, y; Δ, a) − Ψ_Δ(x, y)/φ(Δ)|² ≪ Ψ(x, y)²ε_y + QΨ(x, y),  ε_y := e^{−cu/log²(u+1)} + y^{−c}.
In the range Q = x^{o(1)} this can also be derived from Drappeau's Lemme 5 (conductors ≤ R := (log x)^C/(ε_yρ), using
log(1/ρ) = o(log x)) together with the multiplicative large sieve. That derivation is not independent of Harper:
Lemme 5 rests on Harper's calculations (codex #10).
- Non-reduced classes reduce to reduced classes mod Δ/δ at x/δ, and contribute nothing unless δ is y-smooth.
- Intervals and twists require normalisation and variation bounds (codex #10).
- At the design's scales u ≍ D^{3/8}, so ε_y ≤ D^{−A} for every A.

**4.2 Lemma LSQ [proved].** Let 𝒬 be a set of primes in (Q, 2Q] and |a_n| ≤ 1 for n ≤ X. Then
  Σ_{q∈𝒬}Σ_{ψ≠ψ₀ (q)} max_t |Σ_{n≤t}a_nψ(n)|² ≪ (X + Q²)(log 2X)²Σ|a_n|²
(Iwaniec–Kowalski Thm 7.13 plus a dyadic maximal argument; proof critic checked it).
**Corollary (ratio classes).** For A₁, A₂ ⊂ [1, X] coprime to 𝒬 and c_q a unit,
  Σ_q w_q|#{n₁ ≡ c_qn₂ (q)} − |A₁||A₂|/(q−1)| ≪ (max w/Q)(X + Q²)(log X)²(|A₁||A₂|)^{1/2}.
*Size.* Relative to the pair main term the error is ≪ (max w/avg w)·κ·(1/(Q√(ρ₁ρ₂)) + Q/(X√(ρ₁ρ₂)))(log X)³. The
factor κ appears because the main term is thinner than A₁×A₂ (proof critic m1). For Q ≥ B = L[3/8, a] with
1/ρ, κ_max, Δ_max = L[3/8, O(1/a)], the error is super-polylogarithmically small once a² exceeds an absolute constant.
Classes mod Δ should be handled with characters mod qΔ (hybrid), which avoids a Δ-loss (proof critic m1, codex #6).

**4.3 The principal part [sketch; v1's PCV corrected].** For a pair lattice of index K with smooth weight ω, write the
principal part as Σ over joint gcd classes (δ₁, δ₂) of the products of class averages.
- **(i) Local terms.** These need the correctly normalised class density
  τ(δ) := Pr(n ∈ 𝒮 | gcd(n, K) = δ)/Pr(n ∈ 𝒮) = KΨ_{K/δ}(x/δ, y)/(φ(K/δ)Ψ(x, y))
  (codex #3; v1's normalisation was wrong).
  - τ(δ) − 1 = O((1−β)log δ + 1/u) for δ ≤ y^{1/log u}.
  - Larger p have τ up to u log u but weight ≤ y^{−1/log u} (proof critic m7).
  - Primitivity of g adds −Σ_{p∤K}p^{−2}(…). E1 sees this term.
  - Total, on average over the family: O(((1−β) + 1/u)²).
- **(ii) Archimedean term.** v1 claimed a κ¹ loss; the method gives κ² (proof critic M3). Repair:
  - Use Saias's smooth approximation Ψ(x, y) = Λ(x, y)(1 + O(exp(−(log y)^{3/5−ε}))), valid for
    y ≥ exp((log log x)^{5/3+ε}), which holds here [UNVERIFIED range]. The smooth part of the archimedean covariance then
    has no κ-loss, and the residual costs ε_S·κ², harmless for κ ≤ exp((log y)^{1/2}).
  - Thin pairs (κ up to B^{C₀}) are bounded above with polylog loss by [pp] Lemma DS (smooth numbers in short
    progressions, internal). Their weight is ≍ 1/κ.
- **(iii) Skewed labels (codex #4; uprime §7.3 item 3).** For labels whose reduced basis is very anisotropic (e.g.
  Λ_{(q,0)}), almost every pair has κ ≍ q. A label-shape estimate is needed.
  - Such labels form a fraction ≲ 1/T² of the r mod q (uprime §4.2: irregular lattices ≤ q/S₂² residues).
  - Their incidence weight is near-uniform in r, because P is equidistributed mod q (V₁ ≫ q); this needs the weight
    domination below.
  - Not written.

## 5. Rational blocks of C₁ and C₂ [sketch]

Chain: Lemma LAT → class decomposition mod qΔ → characters mod q split off (LSQ, averaged over q) → mod-Δ errors by BDH
summed over the pair family (multiplicities ≪ S₂² uniformly) → principal part §4.3.
- For C₂ the parent values are J·S_B, with character sums Σ_ℓχ(ℓ)Ψ(x/ℓ, B; χ) and u_s ≍ D^{3/8}.
- For parent pairs Δ₁₂ ≤ 2V₁²/q = L[5/8] ≪ Ψ.

**Missing (codex #1, FATAL to v1's claimed outline): a weighted application theorem.** BDH and LSQ bound sums of
discrepancies over moduli and classes. The quantities in C₁ and C₂ weight those discrepancies by:
- the label weights w_𝔮;
- the region geometry, which depends on (q, r) through e₁, e₂;
- the environment cell;
- the algebraic pushforward density;
- normalisation by each label's own first moment.

Separate control of Σ_t|e_t| and of the average of the weights r_t does not control Σ_t|e_t|r_t. A single theorem
stating the exact measure, its multiplicity and domination bounds, and the maximal (over scales) versions of BDH it
needs is the main missing piece. Harper's BDH has no supremum over x inside the Δ-sum, so scales must be discretised
(proof critic m9).

## 6. Algebraic blocks and design DW

**6.1 Obstruction in Theorem U's design.** With w prime and W ≈ t·N^{η₀}, a two-point algebraic statement needs BV or BDH
over the moduli G(P′), G(P″). These moduli are the rational values themselves, weighted by the partner's smoothness:
- this family is endogenous and exponentially sparse;
- dropping the partner's indicator costs 1/ρ ≫ ε^{−1};
- uprime's CM plan meets the same family as {s·k} with τ-weights (uprime §7.3, items 6–7).

We know no density-relative BV, and no Brun–Titchmarsh bound for smooth numbers to these moduli.

**6.2 Design DW.** Set F_A := w·F̂_{g,N} + G·R_c, with w uniform among integers in (W_m, 2W_m] coprime to g₁. The
multiplier range W_m is renamed so as not to clash with the incidence count W. Take W_m ≥ (2MV₁)²N^{η₀} and
H ≥ W_mN^{η₀}, i.e. h ≥ 2μ + 2ν + 2η₀ (proof critic m4: H ≥ W_mN^{η₀} suffices for the archimedean shift).
- F_A(g₀, g₁) = wN, and the resultant is ±wN.
- At level 1, μ ≍ D^{−1/4} is lower order than dν ≍ D^{−1/8}, so θ is unchanged at leading order.
- At level 2 the extra μ changes only the constant of the child values (σ₂), as uprime's (H) does.

**Lemma DW(a) [proved for filtered points].** Fix g and P with t = |G(P)| ≥ M|P|N^{−η₀/4}, gcd(N, t) = 1 and
gcd(b, g₁) = 1. Then F(P) ≡ w·F̂_N(P) (mod t), with F̂_N(P) a unit mod t, so F(P) mod t is uniform over w up to
O(tτ(g₁)/W_m). The c-part has [pp] Lemma smooth1's law at scale t·H|P|^{d−1}; the shift |wF̂_N(P)| is
≤ N^{−η₀/4} × that scale by the filter. So Pr_{w,c}[F(P) ∈ 𝒮] is the density of 𝒮 against an explicit smooth unimodal
weight, pointwise in g.
- Codex #2: without the filter (e.g. G(P) = 1) the claim is false; v1 stated it universally.
- Parents satisfy the filter. Child points need it added to Good(P′); its g-dependent pass fraction is absorbed by DC_e.

**Lemma DW-L (exact support lattice) [proved; codex adversarial #5].** For P₁, P₂ with gcd(N, G(P₁)G(P₂)) = 1 and
primitive P_i, the set {(F_A(P₁), F_A(P₂))} over all (w, c) ∈ Z × Z^d with F_A primitive-compatible equals
{F₁ ≡ λ^dF₂ (mod |det(P₁, P₂)|)}. The order of v := (F̂_N(P₁), F̂_N(P₂)) modulo (G(P₁), G(P₂))·I is |G(P₁)G(P₂)|.
*Proof (codex).* Two forms with prescribed values at P₁, P₂ differ by L_{P₁}L_{P₂}R. Their values at g differ by
G(P₁)G(P₂)Z. Coprimality with N allows imposing F(g) ∈ NZ, and Gauss's lemma then gives F = wF̂_N + G·R_c. ∎ Codex
checked 1 624 small configurations; v1's "order divides K" was wrong.

**Lemma DW(b) (two-point law) [sketch; proof critic M4].** One would want: the law of (F(P₁), F(P₂)) is within TV ε of
the smooth-weighted uniform law on the DW-L lattice. Kern-U gives only fibre upper bounds. What is needed is a
two-point analogue of [pp] smooth1: a sum of d independent segments in R², anisotropic for nearly parallel pairs. Under
the root event R_𝔮 one also needs smoothing over a coset mod q (H ≫ qdV) and the q² | det case. Not written.

**6.3 Consequences [sketch].** Given DW(b), the algebraic pair for fixed g is a det-modulus lattice smooth-pair problem.
It is handled by the same tools as §5:
- Lemme 5 at y = B for children and y = B₂ for parents, both with u ≍ D^{3/8};
- no mod-q ratio arises for F/q;
- the moduli are |det|, not G(P′), so the endogenous family of 6.1 disappears.

For Theorem U itself (relevant to `uprime`, not our file): DW(a) makes FM's algebraic factor pointwise for filtered
points, so FM would not need def:good or Lemma 8.5′. uprime's Lemma 8.5H is then needed only in the original design.
DW is **unchecked** against the internal chain. Codex #7 lists three steps it breaks, all repairable:
- with integer w, gcd(w, t) > 1 is possible, so lem:coprime fails;
- SAMPLE_FIELD's "content 1 or w" fails, since any divisor of w can occur, so non-primitive forms must be rejected;
- Lemma IR's certificate primes need r ∤ w.

## 7. What remains (the honest gap), ranked

1. **The weighted application theorem** (§5; codex #1). This is the main missing mathematics.
2. **Skewed labels** in the archimedean and transversality steps (§4.3(iii); codex #4; uprime §7.3 item 3).
3. **Two-point smoothing** (DW(b); proof critic M4; uprime §7.3 item 4).
4. **The archimedean second-order bound** via Saias plus Lemma DS (§4.3(ii)), with explicit treatment of Fourier tails:
   the smooth-pair upper bound must not be circular, and the weights must be differentiable enough (codex #4).
5. **DW's compatibility** with Lemmas GB/coprime/IR, LV6-U and POSTPROCESS (codex #7).
6. **The one-point child lower bound** m_𝔮 ≥ 8/η (FM/GM at the child level) [sketch, routine per uprime].

If DW is rejected, the algebraic two-point block in the original design remains open, with the obstruction of 6.1.

**Where the outlines are not a known open problem.** None of items 1–6 is of the type of Hyp RB, a ternary smooth
equation with growing coefficients. All are binary, lattice-based estimates at polylogarithmic precision in a regime
(moduli x^{o(1)}, u ≍ D^{3/8}) where BDH and the large sieve have super-polylogarithmic slack.

**7.5 Deeper descents [model; codex #9 and proof critic m3 corrected v1].** At depth i, DC needs ϕ_i ≤ η_i and
σ_i²ϕ_i ≲ η_i², with η_i ≍ D^{−iα}. Against the heuristic ϕ ≍ σ² ≍ D^{−2α}:
- the margin is positive at i = 1;
- it is borderline (0) at i = 2, the same as multidescent's direct route;
- it is negative for i ≥ 3.

So 4/11 is borderline and below that it is open.

**7.6 3/8 < α < 2/5 [model; codex blind §4].** A two-level chain with m ≍ D^{3−7α} children per parent needs precision
D^{−(3−7α)} > D^{−α}. The rest is unchanged.

## 8. Answers to the task, item by item

1. **Averaging.** Yes. Averaged over pairs inside a child box and over the child primes, BDH and the large sieve
   control the non-principal parts of the two-point blocks with super-polylogarithmic slack. Their weighted application
   (§7.1) is not written. The averaging does not yield the ternary RB; Prop DC makes RB unnecessary.
2. **Circle method.** For the two-point problems the principal arc is the gcd-class/archimedean principal part. It
   does **not** reproduce the product of marginals exactly. The residual is second-order (codex blind §2; §4.3). We did
   not carry out the dilated ternary circle method for RB. In it, effective denominators ≤ R correspond to
   θ-denominators up to qHR (codex blind §2).
3. **Weaker variants.** None needed at 3/8 if C₁ and C₂ hold. Above 3/8 the precision relaxes (7.6). DC does not help
   below 3/8 (7.5).
4. **Algebraic block.** This is the clean reduction under DW (DW-L proved, DW(b) sketch). The original-design
   obstruction is in 6.1.
5. **Obstruction.** These are the items of §7, ranked. The single most important one is the weighted application
   theorem for BDH and the large sieve under the incidence measure.

## 9. Experiments (pre-registered, `notes/pc3_prereg.md`, commit `1025bda`)

Toy setup: rational block only, exact smoothness; y = 150, M = 2000, S = 20, labels q ∈ (500, 1000], 20 labels per
experiment. Child values are ≈ 10⁶ (u ≈ 2.9) and parent values reach 8·10⁷. The smoke-run disclosure is in the prereg:
the principal proxy was introduced after seeing the primitivity effect.

### 9.1 E1: two-point structure of the child box (20 labels × 511 points, 3·10⁵ g per label)

| statistic | value | 95% CI (label bootstrap) | band | verdict |
|---|---|---|---|---|
| P1a ⟨C_emp − C_mod⟩ | −0.0002 | [−0.0008, 0.0005] | ±0.01, CI ⊂ ±0.02 | held |
| P1b slope cov/var | 1.18 | [1.16, 1.21] | [0.8, 1.3] | held |
| P1c excess dispersion | 0.0070 (var C_mod = 0.0150) | — | ≤ 0.5·var | held (barely) |
| P1d box mean ⟨C_emp⟩ | 0.0005 | [−0.0007, 0.0018] | ±0.01 | held |

Reading:
- Per-pair covariances are large at toy scale: sd 0.15, with binned means from −0.28 to +0.23. The proxy tracks them;
  it uses local factors, including the primitivity of g, times the size density.
- The proxy under-disperses by ≈ 18%, probably from non-multiplicativity at p ≥ 17 and the local × size interaction.
- The pair covariances cancel over the box, to 5·10⁻⁴, because the box is complete mod small p.
- Toy illustration of the C₁ premise only.

### 9.2 E2: parent conditioning vs the unconditional box (20 labels × K ∈ {20, 80, 320})

| statistic (pooled; label bootstrap CI) | K = 20 | K = 80 | K = 320 | band | verdict |
|---|---|---|---|---|---|
| m; z (parents per g) | 40.5; 80.5 | 40.5; 62.6 | 40.5; 39.5 | — | — |
| ϕ_ex = ϕ − 1/m | −0.0043 | −0.0043 | −0.0043 | — | — |
| P2a Φ_ex − ϕ_ex | −0.0013 [−0.0015, −0.0012] | −0.0014 [−0.0015, −0.0012] | −0.0013 [−0.0015, −0.0012] | abs ≤ 0.01 | held |
| P2c σ² − 1/z | −0.0018 [−0.0021, −0.0014] | 0.0024 [0.0020, 0.0029] | 0.0029 [0.0023, 0.0036] | [−0.02, 0.05] | held |
| P2d shift ⟨E[N∣E_P]/m − 1⟩ | 0.0043 | 0.0045 | 0.0046 | ≤ 0.03 | held |
| P2e dc − ϕ, dc = E[Z(N−m)²]/(EZ·m²) | −0.0012 | −0.0012 | −0.0012 | ≤ 0.01 | held |
| P2b Φ_ex(320) − Φ_ex(20) | −3·10⁻⁶ [−3·10⁻⁵, 3·10⁻⁵] | | | ≤ 0.01 | held |

Reading (fixed in advance):
- At toy scale, PC₃'s parent-conditioned variance and C₁'s unconditional variance agree to 0.0014.
- A 16-fold change in the parent's coefficients moves the former by 3·10⁻⁶.
- The parent count is nearly Poisson. The toy has no g-dependent parent filter, so it does not show the O(1) σ² term
  that the proof critic identified (M1).
- Nothing here is asymptotic, and the algebraic block is absent.

## 10. Literature

- **BDH for smooth numbers.**
  - Harper, arXiv 1208.5992, Theorems 1–2 (preprint; Thm 2 confirmed by codex).
  - Harper, J. London Math. Soc. (2025), doi:10.1112/jlms.70293: BDH variance asymptotics for general sequences and
    smooth numbers, at large moduli (codex adversarial; we have not read it; it does not directly address the
    small-modulus weighted problem).
  - Drappeau, Compositio 151 (2015), Lemme 5.
- **Large sieve.** Iwaniec–Kowalski Thm 7.13.
- **Local behaviour of Ψ.**
  - Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Thm 3.
  - de la Bretèche–Tenenbaum, Ramanujan J. (2005).
  - Saias, J. Number Theory 32 (1989): smooth approximation Λ(x, y) (UNVERIFIED range).
  - Short intervals: Hildebrand, J. Number Theory (1986), Thm 3, range x/y^{5/12} (codex; v1 cited the wrong paper).
  - Internal: [pp] Lemma DS.
- **Ternary smooth equations.** As multidescent §7: Harper, Compositio 152 (2016); Lagarias–Soundararajan; Drappeau,
  MPCPS 2013.
- **Correction.** Matthiesen–Wang Thm 1.3 assumes Ψ(X, y) > X/log₈X with log₈ the eighth iterated logarithm (codex
  blind).
- **Decoupling via the parent count of a label.** No prior use found by us. This is a search result, not a novelty
  claim.

## 11. Reviews

**Codex blind `pc3_1`** (xhigh, before any analysis here). Verdict OPEN.
- Found the Cauchy–Schwarz decoupling E[Z1_{D₀}]/z ≤ p + √(vp(1−p)), and warned that z ≫ 1 does not give v ≪ 1.
- Showed that the principal arc does not reproduce the conditional marginals exactly.
- Mapped the moduli sizes, and Harper's Thms 1–2 and Compositio minor arcs.
- Found the α ∈ (3/8, 2/5) slack.
- Probabilities: rational block 10%, algebraic block 20% (for PC₃ as stated).

**Codex adversarial `pc3_2`** (on v1): **MAJOR REVISION.**

| # | severity | finding | action |
|---|---|---|---|
| 1 | FATAL (to v1's outline) | weighted application theorem missing | accepted: §5, §7.1; RB2/AB2 demoted to sketch |
| 2 | FATAL | DW(a) false without a size filter (G(P) = 1) | accepted: DW(a) restricted to filtered points; filter added; DC_e |
| 3 | MAJOR | τ normalisation wrong | accepted: §4.3(i) |
| 4 | MAJOR | κ-tail not uniform over skewed labels; Hildebrand citation wrong; Fourier tails circular | accepted: §4.3(ii–iii), §7.2/7.4 |
| 5 | MAJOR | DW quotient argument wrong; exact lattice supplied | accepted: Lemma DW-L (codex's proof) |
| 6 | MAJOR | LSQ correct; weighted/hybrid extensions missing; w-domination not explicit | accepted: §4.2, §7.1 |
| 7 | MAJOR | composite w breaks lem:coprime, SAMPLE_FIELD content, IR certificates (repairable) | accepted: §6.3, §7.5 |
| 8 | MAJOR | product of averages; use A + √(AB) | accepted: Corollary DC |
| 9 | MAJOR | j = 3 is borderline, not failing | accepted: §7.5 |
| 10 | MINOR | BDH-S quote correct; Lemme 5 derivation not independent | accepted: §4.1 |
| 11 | MINOR | LAT and cancellation correct | — |
| 12 | MINOR | verdict overstated; add Harper JLMS 2025 | accepted: §0, §10 |

Codex's probabilities: 35% that Prop UB's conclusion in 2L′ can be justified along the repaired route; 80% that DW
survives the internal chain.

**Opus proof critic (read-only, on v1).**
- Accepted: Lemma 0, Prop DC and the corollary statement, LAT(a–c), LAT(d) (caveat unnecessary), scaling invariance,
  LSQ, DW(a) (for the right points).
- Major:
  - M1: an O(1) environment term in σ² from the g-dependent parent filter, and per-field threshold randomisation would
    give ϕ = Θ(1). Accepted: Prop DC_e and per-point randomisation.
  - M2: the product-average issue. Accepted, as codex #8.
  - M3: the PCV archimedean bound gives κ², not κ. Accepted: Saias plus [pp] Lemma DS.
  - M4: DW(b)/(c) labelled proved without a two-point smoothing lemma. Accepted: demoted.
- Minor: m1 (κ in LSQ size; hybrid characters mod qΔ), m3 (j = 3), m4 (sizes, notation), m5 (no τ(Δ) in
  multiplicities, which is a real trap; fixed in §3), m7, m8, m9 (supremum over scales; rational pair upper bound via
  Lemma DS), m10 (overclaims), m11 (α/β notation), m12. All applied.

**Opus literature critic:** [pending at the time of this draft; see §11 addendum].

**Coordination with `uprime`** (coordinator's message).
- U′ is proved at 3/8 except the bullet "E[N_𝔮 | E] large enough for PC₃", which uprime moved into PC₃ as Hyp CM.
- Prop DC does not use the parent-conditioned mean at all. CM is therefore not needed by Theorem 2L-DC: C₁'s 1/m term is
  the root-conditioned one-point mean.
- uprime's §7.3 traps reappear as follows:
  - multiplicity per |P′|-block: avoided, since we average over pairs in coordinate space, multiplicity ≪ S₂²;
  - fixed-q large sieve: we average over q (LSQ);
  - skewed labels: §4.3(iii), open;
  - two-point smoothing scales: DW(b), open;
  - non-reduced residues mod Δ: gcd-class decomposition, §3;
  - τ-weighted product family: avoided by DW, whose moduli are |det|.
- Used from uprime: the design 2L′ with (H), the child box geometry and S₂ from its child scan (σ₂, nominal
  m₂ ≥ exp((log log N)²)), Lemma CS′ for S₁, and Lemma 8.5H for the level-1 moments behind U′.
- DW's height floor h ≥ 2μ + 2ν + 2η₀ is larger than (H) h = μ + 4ν + η₀ only by lower-order terms at level 1. At level
  2 it changes σ₂'s constant.

## 12. Reproduction

```
.venv/bin/python -m pytest -q tests/test_pc3.py
.venv/bin/python -m src.pc3.experiments E1 experiments/raw/pc3_e1.jsonl 20
.venv/bin/python -m src.pc3.experiments E2 experiments/raw/pc3_e2.jsonl 20
.venv/bin/python -m src.pc3.analyze E1 experiments/raw/pc3_e1.jsonl
.venv/bin/python -m src.pc3.analyze E2 experiments/raw/pc3_e2.jsonl
```

## 13. Candidate card

- **Mechanism.** Two-level descent at 3/8 (multidescent Prop UB). The parent–child conditioning is removed through the
  label's parent count (Prop DC/DC_e). The field multiplier is made a linear integer (DW).
- **Result.**
  - Proved: Prop DC, DC_e, Corollary DC, Theorem 2L-DC (L[3/8] conditional on two two-point concentration hypotheses
    C₁, C₂, replacing PC₃ and CM); Lemmas 0, LAT, SI, LSQ, DW-L, DW(a) (for filtered points); the pair-multiplicity
    bound.
  - Sketch: C₁ and C₂ via BDH, LSQ and principal parts.
- **Failure scope.** Not proved: C₁, C₂ (§7). DW is unchecked against the chain. Depth-2 margins of three-level chains
  are 0.
- **Next step.** The weighted application theorem (§7.1), then skewed labels and two-point smoothing.
