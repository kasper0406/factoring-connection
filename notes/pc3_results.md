# pc3: Hyp PC₃ (multidescent §3.2–3.4) — results

Agent `pc3`, phase 26, 2026-10-24. Draft v1 (before reviews).

Files: this note; `notes/pc3_prereg.md` (commit `1025bda`, code frozen there); code `src/pc3/toy.py`, `experiments.py`,
`analyze.py`, `dc.py`; tests `tests/test_pc3.py`; data `experiments/raw/pc3_e{1,2}.jsonl` (+ `_summary.json`); codex
`codex/pc3_1.{md,out}` (blind, launched before any analysis here), `codex/pc3_2.*` (adversarial, §11).

Labels: **proved**; **proved conditional on X**; **outline** (argument written, not checked line by line;
UNVERIFIED); **heuristic**; **open**; **experimentally supported** (toy). Everything about the factoring algorithm
rests on the project's internal, unrefereed chain (Theorem U and its dependencies).

---

## 0. Verdict

**Plain language.** We did not prove Hyp PC₃ as multidescent stated it. We did find that it is the wrong target.
multidescent's union bound (Prop UB) does not need the hard three-point statement. That statement was a ternary
equation in smooth numbers with huge coefficients and conditional weights (Hyp RB). Prop UB needs a bound on how often
an accepted parent has a dead child box.

We show by an elementary argument (Prop DC) that this bound follows from two *two-point* statements:
- the child box of a label has small relative variance *without any conditioning on the parent*;
- the number of parents of a label (it is large, about Y/q) has small relative variance.

Codex found the same Cauchy–Schwarz decoupling independently in its blind round. Two-point statements about smooth
values of two linear forms are lattice-point problems for pairs of smooth numbers. They are of
Barban–Davenport–Halberstam type and need no ternary circle method.
- The rational parts reduce to two known tools: BDH for smooth numbers (Harper), averaged over the pair family, and the
  multiplicative large sieve over the child primes q. The large sieve works because the children (q ≥ B = L[3/8, a])
  are far larger than the inverse smooth density (L[3/8, O(1/a)]).
- The algebraic parts need one design change (DW): the field multiplier w becomes a uniform *integer* and enters
  linearly, F = w·F̂_{g,N} + G·R_c. Then, for fixed g, the algebraic values at two points are spread uniformly over the
  homogeneity lattice F(P′) ≡ λ^d F(P″) mod det(P′, P″). The algebraic block becomes the same kind of lattice
  smooth-pair problem. The "endogenous sparse modulus family" disappears.

The required precisions turn out to be first-order (≈ D^{−3/8}) for each statement. The natural size of both
covariances is second-order (≈ D^{−3/4}), which leaves a margin of D^{3/8} up to logarithms.

**What is proved and what is not.**
- Proved: Prop DC; the lattice structure (Lemma LAT); the q-large-sieve lemma (LSQ); the integer-multiplier
  distribution lemma (Lemma DW).
- BDH for smooth numbers is cited from Harper (arXiv 1208.5992, Thm 2; also derivable from Drappeau's Lemme 5 and the
  large sieve, §4.1).
- The assembly into the two-point covariance bounds (Props RB2, AB2) is a detailed **outline**, not checked line by
  line. Its archimedean and local main terms need the standard local estimates (Hildebrand–Tenenbaum Thm 3) and an
  upper bound for smooth numbers in short intervals (Hildebrand) for nearly parallel pairs.
- **No rigorous α < 2/5 is claimed.** Theorem 2L still needs Hyp U′ (sister agent `uprime`). If the outlines hold, the
  conditional L[3/8] of Theorem 2L rests on U′ alone, in the modified design. The same design change also removes
  Lemma 8.5′ from the algebraic half of Theorem U's first moment, which is relevant to U′ (§6.4).

Precisely:
1. **Prop DC (decoupling) [proved; codex blind independently].** Σ_{(ℓ,P,q)} Pr[ℓ, E ∧ N_𝔮 = 0] ≤
   Σ_𝔮 w_𝔮(ϕ_𝔮 + σ_𝔮√ϕ_𝔮). Here ϕ_𝔮 is the root-conditioned relative variance of the child count, σ_𝔮² the
   root-conditioned relative variance of the label's parent count, and Σ w_𝔮 = W. So Prop UB's conclusion E N₂ ≥ μ₁/2
   holds if, on w-weighted average, ϕ ≤ η/4 and σ²ϕ ≤ η²/16, with η = μ₁/W ≍ D^{−3/8} (§2). This replaces PC₃. PC₃
   (conditional variances) is sufficient for the same conclusion but is not what we prove.
2. **The two-point problems are lattice problems [proved].** For primitive non-proportional P′, P″ the pair
   (G_g(P′), G_g(P″)), g ∈ Z², fills exactly the lattice n₁ ≡ λn₂ (mod D), D = |det(P′, P″)|, λ a unit. For any binary
   form Φ of degree d, Φ(P′) ≡ λ^dΦ(P″) (mod D). For child pairs D = qΔ with Δ ≤ 2S₂², and for parent pairs of one
   label D = qΔ₁₂ with Δ₁₂ ≤ 2V₁²/q (§3).
3. **Rational block [BDH + LSQ proved/cited; assembly outline].** The relative pair covariance splits into three parts:
   - a principal (local × archimedean) part, which is second-order, O(((1−α) + 1/u)²·log D) on average over the pair
     family;
   - a mod-Δ non-principal part, bounded by BDH for smooth numbers on average over Δ with error e^{−cu/log²u} + y^{−c};
   - a mod-q non-principal part, bounded by Lemma LSQ on average over the child primes q, with error
     ≪ (log x)^{O(1)}/(qρ).

   This applies to child pairs (DC1) and to parent pairs (DC2) (§4–5). It does NOT give the ternary Hyp RB, which is
   no longer needed.
4. **Algebraic block [Lemma DW proved; AB2 outline].** Make the design change DW: F = w·F̂_{g,N} + G·R_c with w a
   uniform integer and W ≥ (2MV)²·2V²·N^{η₀}, which costs a constant factor in θ. Then for fixed g the algebraic pair
   law is the smooth-weighted uniform law on an explicit lattice of index D·(local factor), up to TV error N^{−η₀/4}.
   AB2 then follows as for the rational block, by BDH over the pair family. No BV over endogenous moduli and no
   individual-modulus input is used (§6).
5. **Precision bookkeeping [model].** With ϕ, σ² = O(D^{−3/4}log D) (what 3–4 give), the condition in 1 holds with
   margin D^{3/8}/log D. The binding constraints are the constants of the design, notably a large factor-base constant
   a, and the hypotheses listed in §7. For j = 3 levels (α = 4/11) the DC route needs σ ≤ D^{−2α}·(…) at depth 2.
   This fails against the heuristic σ ≍ D^{−α}, so DC is a two-level tool (§7.3).
6. **Experiments (pre-registered; toy; rational block only).** 9/9 predictions held.
   - E1: pair covariances in a child box range over ±0.3 at toy scale, but the principal proxy tracks them (slope 1.18,
     mean residual −0.0002 [−0.0008, 0.0005]) and they cancel over the box (mean C_emp = 0.0005).
   - E2: the parent-weighted conditional excess variance equals the unconditional one to within 0.0014 (systematically
     −0.0013, CI ±0.0002). This holds at every parent coefficient scale K ∈ {20, 80, 320}; the difference between
     K = 320 and K = 20 is −3·10⁻⁶. The parent count is nearly Poisson (σ²_ex ≤ 0.003).

---

## 1. Setting

Notation as multidescent §1 and §3.2 (design 2L at α = 3/8: d ≍ D^{1/4}, B = L[3/8, a], ℓ ∈ J = L[5/8],
children q ∈ (B, B₂], B₂ = L[1/2], η := μ₁/W, W/μ₁ = E[m(P)] ≍ D^{3/8}).
- A **label** is 𝔮 = (q, r). Its lattice is Λ_𝔮 = {(a, b) : a ≡ rb (mod q)}, with Gauss-reduced basis e₁, e₂ (so
  |e₁||e₂| ≤ (2/√3)q).
- Box(𝔮) is the set of primitive P′ = ie₁ + je₂ with max(|i|, |j|) ≤ S₂, modulo ±1.
- N_𝔮(A) = #{P′ ∈ Box(𝔮) : Good(P′)}.
- An incidence (ℓ, P, q) has label 𝔮(P, q) = (q, ab^{−1} mod q). The **parent count** of a label is
  Z_𝔮(A) := #{(ℓ, P) : Acc_ℓ(P), q a child of P, 𝔮(P, q) = 𝔮}. Each incidence is counted under exactly one label, and
  Σ_𝔮 E Z_𝔮/#J = W.
- The **root event** is R_𝔮 := {F_A(r, 1) ≡ 0 (mod q)}.

**Lemma 0 (supports) [proved].** Z_𝔮 = N_𝔮 = 0 outside R_𝔮.
*Proof.* A primitive P′ = (a′, b′) ∈ Λ_𝔮 has q ∤ b′, since otherwise q | a′ ≡ rb′ as well. So P′ ≡ b′(r, 1) (mod q)
and F(P′) ≡ b′^d F(r, 1) (mod q). Both Good(P′) and "q is a child of P" require q | F at a primitive point of Λ_𝔮. ∎

*Size of Z_𝔮 [heuristic, typical-count model].* There are #(Λ_𝔮 ∩ 𝒫) ≍ V₁²/(dq log V₁) candidates. Each is accepted
with probability ≍ ρ_r ρ′_F/log Y given R_𝔮, so z_𝔮 := E[Z_𝔮 | R_𝔮] ≍ (Y/q)·m₁/(d log² N) with m₁ = S₁²ρ_rρ_F
polylogarithmic. Since Y = L[5/8] and q ≤ L[1/2], z_𝔮 = L[5/8]^{1−o(1)} ≫ 1: **every label has many potential
parents**. This is what makes decoupling possible.

## 2. Prop DC: decoupling the parent from the child box [proved]

**Proposition DC.** Put m_𝔮 := E[N_𝔮 | R_𝔮] > 0, ϕ_𝔮 := Var(N_𝔮 | R_𝔮)/m_𝔮², σ_𝔮 := Var(Z_𝔮 | R_𝔮)^{1/2}/z_𝔮 and
w_𝔮 := Pr(R_𝔮) z_𝔮/#J. Then
  Σ_{(ℓ,P,q)} Pr_{ℓ,A}[E_{ℓ,P,q} ∧ N_{𝔮(P,q)} = 0] ≤ Σ_𝔮 w_𝔮 (ϕ_𝔮 + σ_𝔮 √ϕ_𝔮),   Σ_𝔮 w_𝔮 = W.
*Proof.* Every incidence of a field A is counted under one label, and N_𝔮 does not depend on ℓ, so the left side is
(1/#J)Σ_𝔮 E[Z_𝔮 1{N_𝔮 = 0}]. By Lemma 0 this is Σ_𝔮 (Pr(R_𝔮)/#J)·E[Z 1{N = 0} | R] (label index dropped). Under the
law given R, put p := Pr(N = 0). Then E[Z 1{N=0}] = zp + E[(Z − z)1{N=0}] ≤ zp + Var(Z)^{1/2}(p(1 − p))^{1/2}
≤ z(p + σ√p), and Chebyshev gives p ≤ ϕ. ∎

**Corollary DC.** If Σ_𝔮 w_𝔮 ϕ_𝔮 ≤ μ₁/4 and Σ_𝔮 w_𝔮 σ_𝔮²ϕ_𝔮 ≤ μ₁²/(16W), then E N₂ ≥ μ₁/2, which is Prop UB's
conclusion. The Paley–Zygmund step of Prop UB is unchanged.
*Proof.* Cauchy–Schwarz: Σ wσ√ϕ ≤ (Σ w)^{1/2}(Σ wσ²ϕ)^{1/2} ≤ μ₁/4. Then follow Prop UB's proof. ∎

*Remarks.*
- With η = μ₁/W, the conditions read: ϕ ≤ η/4 on average, and σ²ϕ ≤ η²/16 on average. Under the heuristics
  ϕ ≈ 1/m + C(1−α)², σ² ≈ 1/z + C′(1−α)², with (1 − α) ≍ D^{−3/8}, the product σ²ϕ ≍ D^{−3/2} sits D^{3/4} below the
  required η² ≍ D^{−3/4}, and ϕ sits D^{3/8} below η. The precision is **first-order** in each factor.
- Codex's blind round (`codex/pc3_1.out` §3) found exactly this inequality, E[Z1_{D₀}]/z ≤ p + √(vp(1−p)). It noted
  that z ≫ 1 alone does not give v ≪ 1, because a common random environment can make Z non-concentrated. That is the
  content of DC2 below.
- **What DC changes.** PC₃ asked for Var(N_𝔮 | E_{ℓ,P,q}). That is a statement about the parent point P together with
  two child points: three linear forms in g with coefficients up to x^{o(1)} (multidescent's Hyp RB), plus three
  algebraic evaluations. DC needs:
  - **DC1**: ϕ_𝔮, the child box with *no* parent conditioning;
  - **DC2**: σ_𝔮, the parent count of a label.

  Both are second moments of sums over a lattice family, so every configuration has **two** points. Lemma MP still
  applies: the precision is O(η), not better.
- Prop DC was tested on random finite joint laws (`tests/test_pc3.py::test_prop_dc_inequality_random`, together with
  the K-split variant of `src/pc3/dc.py`).

## 3. The two-point lattices [proved]

**Lemma LAT.** Let P′, P″ ∈ Z² be primitive and non-proportional, with D := |det(P′, P″)|.
- (a) There is a unique λ ∈ (Z/D)^× with P′ ≡ λP″ (mod D).
- (b) {(G_g(P′), G_g(P″)) : g ∈ Z²} = {(n₁, n₂) ∈ Z² : n₁ ≡ λn₂ (mod D)}.
- (c) Φ(P′) ≡ λ^dΦ(P″) (mod D) for every binary form Φ ∈ Z[x, y] of degree d.
- (d) If P′, P″ ∈ Λ_𝔮, then q | D. In coordinates P′ = ie₁ + je₂, P″ = i′e₁ + j′e₂ one has D = q|ij′ − ji′|.
  The class of λ mod D/q depends only on (i, j, i′, j′), and λ ≡ b′/b″ (mod q).

*Proof.*
- (a) Extend P″ to a basis (P″, Q) of Z². Write P′ = λP″ + μQ; then μ = ±det(P″, P′) = ±D, so P′ ≡ λP″ (mod D). If a
  prime p divides gcd(λ, D), then P′ ≡ 0 (mod p), contradicting primitivity.
- (b) The image is a lattice of index |det| = D. It lies in the right-hand side because G is linear and (a) holds.
  The right-hand side also has index D.
- (c) follows from (a) by homogeneity.
- (d) det(e₁, e₂) = ±q. Modulo D/q the basis is invertible (gcd(q, D/q) = 1 is the generic case; if q² | D, treat the
  q-part separately). Modulo q both points are multiples of (r, 1). ∎

`tests/test_pc3.py::test_pair_lattice_is_unit_graph` checks (b) on random pairs.

*Consequences.*
- **Rational pairs.** For g uniform on a box, (n₁, n₂) is uniform on the lattice of (b), intersected with a
  parallelogram. Primitive g are handled by Möbius over g = e·g′.
- **Local structure.** At p ∤ D the two values are independent mod p, except for primitivity: for primitive g they
  are never both ≡ 0 (mod p). This is a negative covariance; E1 found it, §9.1. At p | D, p ∤ q, the gcd's with p
  agree. The mod-q part is the ratio condition n₁/n₂ ≡ b′/b″ (q). Since q > B, it involves only units.
- **Size of the moduli.** Child pairs: Δ ≤ 2S₂² = L[3/8, O(1)] = y^{O(1)}. Parent pairs of a label:
  Δ₁₂ = |det(P₁, P₂)|/q ≤ 2V₁²/q = L[5/8, O(1)]. Both are x^{o(1)}, while Ψ(x, y) = x^{1−o(1)}, so there is no
  level-of-distribution problem (codex blind §1).
- **Lattice submodules are scaling-invariant.** If Ω ⊂ Z² is any sublattice containing K·Z², then its image
  Ω̄ ⊂ (Z/K)² is closed under (a₁, a₂) ↦ (ua₁, ua₂) for units u. Hence, for each a₂, the multiset
  {gcd(a₁, K) : (a₁, a₂) ∈ Ω̄} depends only on gcd(a₂, K). This is what makes the cross terms of §4.3 vanish exactly.

## 4. Analytic inputs

Throughout, S = 𝒮_y is the set of y-smooth numbers, Ψ = Ψ(x, y), u = log x/log y, α = α(x, y) the saddle point, and
ε_y(u) := e^{−cu/log²(u+1)} + y^{−c}. At the design's scales, u ≍ D^{3/8}, so ε_y(u) ≤ D^{−A} for every A.

### 4.1 BDH for smooth numbers

**Lemma BDH-S [cited: Harper, arXiv 1208.5992, Theorem 2, as quoted and checked by codex blind §1; exact constants
UNVERIFIED by us].** For (log x)^K ≤ y ≤ x and 1 ≤ Q ≤ Ψ(x, y),
  Σ_{Δ≤Q} Σ_{(a,Δ)=1} |Ψ(x, y; Δ, a) − Ψ_Δ(x, y)/φ(Δ)|² ≪ Ψ(x, y)² ε_y(u) + QΨ(x, y).
*Independent derivation, in the range we use [outline].* By orthogonality the left side is
Σ_Δ φ(Δ)^{−1}Σ_{χ≠χ₀ (Δ)}|Ψ(x, y; χ)|². Split by the conductor r of χ.
- (i) r ≤ R. Use |Ψ(χ)|² ≤ Ψ·|Ψ(χ)| and reduce induced characters to primitive ones (lv5unif §5, E3). The divisor sum
  Σ_{rd|Δ≤Q} 1/φ(Δ) ≪ log Q·log log x/(φ(r)φ(d)), and Drappeau's Lemme 5 (Compositio 2015, as used in lv5unif (D),
  conductors ≤ x^{η_D}) gives ≪ Ψ² ε_y(u)(log x)^{O(1)}.
- (ii) r > R. The multiplicative large sieve gives ≪ (log x)^{O(1)}(x/R + Q)Ψ.

With R := (log x)^{C}/(ε_y(u)ρ) ≤ x^{η_D/2}, which holds because log(1/ρ) ≍ u log u = o(log x), the total is
≪ Ψ²ε_y(u)(log x)^{O(1)} + QΨ(log x)^{O(1)}.

*Extensions used below [outline].*
- (a) **Intervals.** Replace Ψ(x, y; Δ, a) by Ψ(I; Δ, a) for I = (x₁, x₂], with x₁ ≥ x/y, by differencing.
- (b) **Non-reduced classes.** For a with gcd(a, Δ) = δ, n = δn′ reduces to reduced classes mod Δ/δ at x/δ.
- (c) **Smooth weights and twists.** Partial summation for Σ_{n∈S, n≡a} ω(n)e(ξn) with |ξ| ≤ D^A/x costs a factor
  D^{2A}.

### 4.2 The large sieve over the child primes

**Lemma LSQ [proved, from the multiplicative large sieve].** Let 𝒬 be a set of primes in (Q, 2Q] and (a_n)_{n≤X}
complex with |a_n| ≤ 1. Then
  Σ_{q∈𝒬} Σ_{ψ≠ψ₀ (q)} max_{t≤X} |Σ_{n≤t} a_nψ(n)|² ≪ (X + Q²)(log 2X)² Σ_n |a_n|².
*Proof.* Without the maximum this is Iwaniec–Kowalski Theorem 7.13, since characters mod a prime q are primitive. The
maximum over t costs (log 2X)², by the Rademacher–Menshov dyadic argument applied to the bilinear form. ∎

**Corollary LSQ (ratio classes).** Let A₁, A₂ ⊂ [1, X] be sets of integers coprime to every q ∈ 𝒬, weights
w_q ≥ 0, and c_q ∈ (Z/q)^×. Then
  Σ_{q∈𝒬} w_q |#{(n₁, n₂) ∈ A₁×A₂ : n₁ ≡ c_qn₂ (q)} − |A₁||A₂|/(q − 1)| ≤ (max_q w_q/Q)·C(X + Q²)(log 2X)²(|A₁||A₂|)^{1/2}.
*Proof.* The count is (q − 1)^{−1}Σ_ψ ψ̄(c_q)A₁(ψ)Ā₂(ψ), where A(ψ) := Σ_{n∈A}ψ(n). Remove ψ₀, apply Cauchy–Schwarz over
(q, ψ), then the lemma. ∎

*Size in the design.* Let |A_i| = ρ_iX. Relative to the main term Σ_q w_q|A₁||A₂|/q, the error is
≪ (max w/avg w)·(1/(Q√(ρ₁ρ₂)) + Q/(X√(ρ₁ρ₂)))·(log X)³.
- For children Q ≥ B = L[3/8, a] and 1/ρ = L[3/8, c_ρ] with c_ρ ∝ 1/a. So the error is L[3/8]^{−(a − c_ρ)}: super-
  polylogarithmically small once a² exceeds an absolute constant. **This is the only place where a large factor-base
  constant is essential.**
- The weights w_q need only polylogarithmic domination (max w/avg w ≤ (log N)^{O(1)}); under DW this is explicit, §6.3.
- The corollary extends to A_i restricted to a class mod Δ, with smooth weights and additive twists, via 4.1(c). The
  bilinear form then runs over (q, ψ, a mod Δ).

### 4.3 The principal part and its covariance

Fix a pair type and the lattice Ω ∋ (n₁, n₂) of index K, smooth-weighted by ω. Decompose by classes mod K:
  Σ_{(n₁,n₂)∈Ω} ω 1_S(n₁)1_S(n₂) = ∫∫ ω̂(ξ) Σ_{(a₁,a₂)∈Ω̄} S(ξ₁; K, a₁) S(ξ₂; K, a₂) dξ.
Write S(ξ; K, a) = M(ξ; K, δ(a)) + E(ξ; K, a), where M is the average of S over the gcd class δ(a) = gcd(a, K).
- **Cross terms vanish exactly.** For each a₂, Σ_{a₁ : (a₁,a₂)∈Ω̄} M(ξ₁; K, δ(a₁)) depends only on δ(a₂), by §3's
  scaling invariance, and Σ_{a₂ ∈ class δ} E = 0.
- **The E₁E₂ term** is ≤ (fibre bound)·(Σ_{a₁}|E₁|²)^{1/2}(Σ_{a₂}|E₂|²)^{1/2}, i.e. BDH-type.
- **The principal part** is Σ_δ π(δ) ∫∫ ω M₁M₂, which is the same as for two *independent* smooth variables in the
  gcd class δ.

**Lemma PCV (second-order principal covariance) [outline].** The principal part minus the product of the marginals
consists of three pieces.
- (i) **Local.** A covariance over gcd classes, Σ_δ π(δ)(τ₁(δ) − 1)(τ₂(δ) − 1)·(…), with
  τ(δ) := Ψ(x/δ; coprime to K/δ)/(δ^{−1}Ψ(x; coprime to K)). By Hildebrand–Tenenbaum Theorem 3 (δ ≤ y) and the coprime
  local factors Π_{p|m}(1 − p^{−α}) (de la Bretèche–Tenenbaum 2005; UNVERIFIED citation), τ(δ) − 1 =
  O((1 − α)log δ + 1/u)·δ^{1−α}. For p > y, τ = 0 on p | n, with weight 1/p. Hence (i) ≪ ((1−α) + 1/u)²
  Σ_{p|K}(1 + log p)²/p + Σ_{p|K, p>y} 1/p.
- (ii) **Primitivity of g** (Möbius over e). This contributes −Σ_{p∤K} p^{−2}(τ_p − 1)² + …, which is
  ≪ ((1−α) + 1/u)². E1 shows it at toy scale.
- (iii) **Archimedean.** With E(t) := Ψ(t) − (Ψ(x)/x)t, the archimedean covariance is exactly the second-order
  expression ∫∫ω dE dE − (∫ω₁ dE)(∫ω₂ dE)/∫∫ω, since the first-order terms cancel identically. Integrating by parts
  twice and using |E(t)| ≪ ((1−α)log(x/t) + 1/u)Ψ(t) for t ≥ x/y (HT Theorem 3) gives
  ≪ ((1−α) + 1/u)²·κ, where κ = |P′||P″|/|det(P′, P″)| is the pair's condition number.
  - Over a box the density of κ is ≍ κ^{−2}, so pairs with κ ≤ K₀ contribute ≪ ((1−α)+1/u)² log K₀.
  - Pairs with κ > K₀ have weight ≍ 1/K₀. For them we use only an upper bound joint/product ≤ C(log N)^C. That needs
    an upper bound for smooth numbers in intervals of relative length ≥ 1/κ: Hildebrand (Quart. J. Math. 1985):
    Ψ(t + h, y) − Ψ(t, y) ≍ hΨ(t, y)/t for h ≥ t/y^{5/12} (UNVERIFIED statement), plus the trivial bound beyond.
  - With K₀ = D^{3/8+C+γ} the tail is negligible. The thinnest pairs (κ ≥ y^{5/12}) have weight ≤ y^{−5/12}, and this
    beats the trivial loss 1/ρ once a is large.

Averaged over the pair family (E_Δ Σ_{p|Δ}(log p)²/p = O(1)), the principal covariance is
**O(((1−α) + 1/u)² log D) = O(D^{−3/4} log D)**. Note that 1/u ∝ a·D^{−3/8} and (1−α) ∝ D^{−3/8}/a.

*Fourier tails.* ω̂ is not compactly supported. Truncating at |ξ| ≤ D^A/x costs D^{−Ak} times the *smooth* pair count
over a slightly wider region, so a crude upper bound with polylog loss is needed. For product boxes it follows from the
exact class decomposition with BDH and LSQ applied to sharp intervals, with no Fourier step at all. A Whitney cover then
gives it for general regions [outline].

## 5. The rational two-point block

**Proposition RB2 [outline; inputs: BDH-S, LSQ, PCV, smooth g-law].** Use the design with g drawn from a smooth
density supported in [M/2, M]² (randomised thresholds for the sharp filters, §7.1). Let C(P′, P″) be the relative
covariance of 1_S(G(P′)) and 1_S(G(P″)), weighted by smooth size weights h₁(|G(P′)|)h₂(|G(P″)|); these weights come from
the algebraic block, §6. Then:
- (a) **child pairs:** Σ_𝔮 w_𝔮 E_{P′≠P″ ∈ Box(𝔮)} |C(P′, P″)| ≪ W·(D^{−3/4} log D + ε);
- (b) **parent pairs of a label**, with the rational parent condition G(P) ∈ J·S_B: the analogous average of the pair
  covariance over P₁ ≠ P₂ ∈ Λ_𝔮 ∩ 𝒫 is ≪ D^{−3/4} log D + ε.

Here ε := (log N)^{O(1)}(ε_B(u) + L[3/8]^{−(a−c_ρ)}).

*Outline.*
- Lemma LAT puts each pair on a lattice of index qΔ. Characters mod q split off: n_i are units mod q.
- ψ₀ (mod q): the mod-Δ problem. Decompose by classes (§4.3); cross terms vanish; the E₁E₂ terms are bounded by
  BDH-S (4.1(a–c)), summed over the pair family. The family of Δ's is dense:
  #{pairs with det/q = Δ} ≪ #pairs·τ(Δ)log/Δ_max, a standard lattice-point count. The principal part is covered by
  Lemma PCV.
- ψ ≠ ψ₀: bilinear in (n₁, n₂) after Fourier separation of the weight. Corollary LSQ over q with the w-weights (max/avg
  polylogarithmic).
- (b): the parent values n = ℓs with ℓ ∈ J and s B-smooth with u_s ≍ D^{3/8}. Character sums factor as
  Σ_ℓ χ(ℓ)Ψ(x/ℓ, B; χ), so BDH-S applies to each x/ℓ. Δ₁₂ ≤ 2V₁²/q ≪ Ψ, so the level is fine. LSQ needs
  Q ≫ 1/ρ(J·S) = L[3/8]^{O(1/a)}, as before.
- Möbius over primitive g, and the q²|D pairs (weight 1/q, bounded crudely), are lower order. ∎(outline)

*Why this is not Hyp RB.* RB had three forms G(P)/ℓ, G(P′), G(P″) of a common g with Plücker coefficients up to
x^{O(D^{−1/8})}; it is a genuinely ternary problem whose minor arcs are scrambled by the dilations (multidescent §3.4;
codex blind §2, which also shows that the principal arc does not reproduce the conditional marginals exactly). After DC
no configuration contains a parent together with a child, and every statement is binary.

## 6. The algebraic block

### 6.1 The obstruction in the original design

In Theorem U's design, w is a random prime, F̂_{g,Nw} depends non-linearly on w, and the forced residue
F(P) ≡ Nw b^d g₁^{−d} (mod G(P)) ([pp] Lemma homog) is spread only up to a Brun–Titchmarsh constant (Lemma W-U).
- A two-point statement then needs BV/BDH over the moduli G(P′), G(P″). These are the rational values themselves,
  weighted by the smoothness of the partner, an endogenous and exponentially sparse family.
- The best-known bounds are absolute, not density-relative. Dropping the partner's indicator loses 1/ρ ≫ ε^{−1}.
- So a Brun–Titchmarsh bound for smooth numbers in APs to moduli ≈ x^θ would be needed. We know of none.

This is the "density-relative BV" obstruction of multidescent §3.4 and codex blind §5. We do not solve it. We remove it
by a design change.

### 6.2 Design change DW and Lemma DW

**Design DW.**
- F_A := w·F̂_{g,N} + G·R_c, where F̂_{g,N} is Lemma GB's polynomial for N (so F_A(g₀, g₁) = wN; the resultant is ±wN as
  before).
- w is uniform on the integers in (W, 2W] coprime to g₁, which keeps the leading coefficient a unit mod g₁ as in
  Lemma GB.
- W := ⌈N^{2μ+4ν+2η₀}⌉ ≥ (2MV)²·2V²·N^{η₀}, and the height satisfies H ≥ W·N^{η₀}·dV. Then |wF̂_N(P)| ≤ X_F N^{−η₀/2},
  and h ≥ 2μ + 4ν + 3η₀ replaces Theorem U's floor h ≥ 1 − dμ + ν + η₀ (now F̂_N's constant coefficient is ≤ M).
- θ changes by a constant factor (≈ 3/(d+1) → ≈ 4/(d+1) at the optimum), which affects only the constant of L[3/8, c].

**Lemma DW [proved; TV bookkeeping as in [pp] Lemma smooth1/kern].** Fix g, and points P₁, P₂ with |P_i| ≤ V,
t_i := |G(P_i)| ≥ 1, gcd(N, t₁t₂) = 1 and gcd(b_i, g₁) = 1.
- (a) **One point.** F(P₁) ≡ w·F̂_N(P₁) ≡ w·N b₁^d g₁^{−d} (mod t₁), and w mod t₁ is uniform up to
  O(t₁ τ(g₁)/W). So the law of F(P₁) over (w, c) is within TV O(N^{−η₀/4} + d²V/H) of
  ℒ₁ := (uniform residue mod t₁) ⊕ t₁·ϑ. Here ϑ is the smooth unimodal law of ⟨c, v(P₁)⟩ ([pp] Lemma smooth1). In
  particular Pr_{w,c}[F(P₁) ∈ 𝒮] equals the density of 𝒮 against a smooth unimodal weight at scale X_F, **pointwise in
  g**, with no equidistribution of smooth numbers in progressions.
- (b) **Two points.** (F(P₁), F(P₂)) = w·v + (t₁z₁, t₂z₂) with v := (F̂_N(P₁), F̂_N(P₂)) and
  (z₁, z₂) = (⟨c, v(P₁)⟩, ⟨c, v(P₂)⟩) ∈ I. The image lattice I has index |det(P₁, P₂)|·γ (Kern-U), with γ a gcd factor
  bounded by the common factors of the coordinates. Hence (F(P₁), F(P₂)) lies in the lattice Ω := Zv + (t₁, t₂)·I, and
  its law is within TV O(N^{−η₀/4} + d²V/H + K τ(g₁)/W) of the smooth-weighted uniform law on Ω. Ω has index
  K ≤ t₁t₂·|det|·γ ≤ W N^{−η₀}, and Ω ⊂ {F₁ ≡ λ^d F₂ (mod det)} by Lemma LAT(c).
- (c) **Root event.** Under R_𝔮 the same holds with c restricted to an affine condition mod q, and q | F(P_i) for
  P_i ∈ Λ_𝔮. For the quotients F(P_i)/q the q-part imposes no further ratio condition, since det has q to the first
  power generically.

*Proof.*
- (a) F is linear in (w, c). Reduction mod t₁ kills G·R_c, and w ↦ w·F̂_N(P₁) permutes Z/t₁ because F̂_N(P₁) is a
  unit mod t₁ ([pp] Lemma homog). The archimedean part: t₁⟨c, v(P₁)⟩ has the smoothed law of [pp] Lemma smooth1, and
  the shift w·F̂_N(P₁) has size ≤ X_F N^{−η₀/2}, so it moves the law by TV ≤ N^{−η₀/4}.
- (b) Linearity again. The pair is an affine image of (w, c), so its support is the lattice Ω. The equidistribution on
  Ω-cosets comes from two sources: w uniform modulo the order of v in Z²/((t₁, t₂)I), which divides K ≤ W N^{−η₀},
  and the two-point smoothing of c on I ([pp] Lemma kern, Kern-U). ∎

*Consequence.* For fixed g, the algebraic pair is a smooth-pair count on an explicit lattice Ω of index K. The
scaling invariance of §3 applies to Ω, so §4.3's decomposition applies verbatim. The moduli are det-related
(≤ 2V²·(local)), not the rational values: t₁, t₂ enter only through the archimedean scale of F and through local
factors at primes dividing gcd(t_i, det). **The algebraic block is no longer coupled to an endogenous modulus
family.**

### 6.3 Proposition AB2 [outline]

Under DW, for child pairs (resp. parent pairs of a label), the root-conditioned algebraic pair covariance given g,
averaged over the pair family, is O(D^{−3/4}log D + ε) relative, uniformly in g. The values are B-smooth for children;
for parents they are B₂-smooth with the acceptance shape (distinct large primes, which is a subtraction of density
≲ B^{−1/2}).
- *Outline.* Lemma DW(b, c), then §4.3 with K = index(Ω). BDH over the pair family: Lemme 5 at y = B (children) or
  y = B₂ (parents), with u ≍ D^{3/8} in both cases, so ε is super-polylogarithmic. No mod-q ratio appears. The
  principal part is Lemma PCV with α replaced by the algebraic saddle point α_F (resp. β₂).
- The w-weights of the labels are explicit under DW: E Z_𝔮 ∝ Σ (1/q)ρ_𝒜(X_F/q)·(…), a smooth function of q. So the
  domination needed by LSQ (max/avg polylog) is a computation, not a hypothesis (Lemma CL's algebraic analogue).
- **Combination with the rational block.** For fixed g the algebraic pair probability is π₂(g) = h₂(|t₁|, |t₂|; type)
  ·(1 + e(type, scales)), with h₂ explicit and smooth. So the full pair probability is a rational pair count with the
  smooth weight h₂ (RB2), plus Σ_types ē(type)·(rational pair probability), where ē := sup over scales of |e|. The second
  sum is controlled by BDH over the types, with the rational pair probabilities bounded on average over q by LSQ.

### 6.4 Side consequence for Theorem U and Hyp U′ (not our files; for `uprime`)

Under DW, Lemma DW(a) makes the algebraic factor of Theorem U's first moment pointwise in (g, P). FM-U then needs
neither def:good nor Lemma 8.5′ (Drappeau Lemme 5) for its algebraic half. The rational half (fibres, Lemma FIB) is
unchanged.
- LD*-U's algebraic pair becomes the lattice problem of Lemma DW(b) with det = ℓ·q′, so the ratio condition is mod ℓ
  and the large sieve averages over ℓ ∈ J, as in [pp] Theorem L2.
- Codex adversarial #1 (multidescent) located U′'s hybrid-BV difficulty in Lemma 8.5′ with B₂-smooth values and
  B-smooth modulus cofactors; DW bypasses that lemma in FM.
- **Unchecked:**
  - Lemma IR, LV6-U and POSTPROCESS with w composite. Lemma IR's proof uses only c; character primes are tested against
    D_A by gcd.
  - Lemma GB with w·F̂_N in place of F̂_{Nw}.
  - The new height floor.

  We flag DW for re-check in the internal chain.

## 7. Assembly, exponents, and what remains

### 7.1 Theorem PC₃′ and Theorem 2L′ [outline; conditional]

**Design 2L′** = design 2L with three changes:
- DW;
- g drawn from a smooth density on [M/2, M]² (primitive);
- the sharp filters (box radii, b-window, |G| ≥ MV/32) applied with thresholds drawn uniformly from ranges of
  relative width 1/2. The analysis then sees smooth weights; the search cost changes by O(1).

The factor-base constant a is large enough that a² > C_LSQ and 5a/12 > c_ρ, and the child box is large enough that
m_𝔮 ≥ D^{3/4}. This costs a polylogarithmic factor in search.

**Theorem PC₃′ [outline].** In design 2L′, Corollary DC's two conditions hold for N ≥ N₀:
- by RB2(a) + AB2, ϕ_𝔮 ≤ 1/m_𝔮 + O(D^{−3/4}log D) on w-average;
- by RB2(b) + AB2, σ_𝔮² ≤ 1/z_𝔮 + O(D^{−3/4}log D).

Hence Prop UB's conclusion E N₂ ≥ μ₁/2 holds, without Hyp PC₃ or Hyp RB.

**Theorem 2L′ [conditional on Hyp U′ for design 2L′, on Theorem PC₃′ (outline), and on the internal chain].** FACTOR_2L
in design 2L′ factors every N in expected time L_N[3/8, c].

**Internal lemmas used:**
- multidescent: Prop UB, Theorem 2L's proof;
- lv5unif: §5 (E3, E6, E8, LS′), Lemma U1, Lemma SM1/Kern-U, Lemma FIB;
- [pp]: Lemmas GB, homog, smooth1, kern, IR (to be re-checked under DW);
- this note: Prop DC, Lemmas 0, LAT, LSQ, DW, PCV, Props RB2, AB2.

**Published inputs:**
- Harper 1208.5992 Thm 2 (BDH; preprint), or Drappeau 2015 Lemme 5 + the multiplicative large sieve (IK 7.13);
- Hildebrand–Tenenbaum 1986 Thm 3; de la Bretèche–Tenenbaum 2005 (coprime local factors);
- Hildebrand 1985 (short intervals);
- Rademacher–Menshov.

### 7.2 Margins

| quantity | needed (Cor. DC) | supplied (outline) | heuristic | margin |
|---|---|---|---|---|
| ϕ (child box, root-conditioned) | ≤ η/4 ≍ D^{−3/8} | 1/m + O(D^{−3/4}log D) + ε | ≍ D^{−3/4} | D^{3/8}/log D |
| σ²ϕ | ≤ η²/16 ≍ D^{−3/4} | O(D^{−3/2}log²D) | ≍ D^{−3/2} | D^{3/4}/log² D |
| LSQ error | ≪ D^{−3/4} | L[3/8]^{−(a−c_ρ)} | — | super-polylog |
| BDH error | ≪ D^{−3/4} | e^{−cD^{3/8}/log²D} | — | super-polylog |

The binding inputs are the constants (a large) and the unverified assembly. Where the outlines use the trivial bound,
the precision obstruction of multidescent §3.4 (secondary terms O(1/u) ≍ D^{−3/8}) is avoided in two ways: the
covariances are written in second-order form (§4.3(iii)), and they are needed only to first-order precision.

### 7.3 Deeper descents, and α ∈ (3/8, 2/5)

- **j ≥ 3 levels [model].** At depth i the union bound needs failure ≤ η_i ≍ D^{−iα}. DC then needs ϕ_i ≤ η_i and
  σ_i²ϕ_i ≤ η_i². With heuristic ϕ ≍ σ² ≍ D^{−2α}, the second condition reads D^{−4α} ≤ D^{−2iα}, i.e. i ≤ 2. At depth 2
  of a 3-level chain (α = 4/11) it fails by D^{−α}: σ_i would need to be ≤ D^{−2α}, below its natural size D^{−α}.
  multidescent's direct route had margin 0 there (§3.5). **So DC helps only for j = 2**, and 4/11 and below remain
  open by both routes.
- **3/8 < α < 2/5 [model; codex blind §4].** A two-level chain with fewer children, m ≍ D^{3−7α}, needs precision
  D^{−(3−7α)} > D^{−α}. Then secondary-term accuracy suffices without covariance form, and DC's margins grow. This does
  not change what must be proved, only the room.

### 7.4 What is not proved (the honest gap)

1. Props RB2, AB2 and Lemma PCV are **outlines**. The likeliest failure points are:
   - (a) the archimedean treatment of nearly parallel pairs (the Hildebrand short-interval input and the κ-tail);
   - (b) the Fourier-tail upper bound;
   - (c) the local structure at primes dividing gcd(t_i, det) in Lemma DW(b), which gives a lattice of slightly
     larger index;
   - (d) the dependence of BDH-S on Harper's preprint.
2. Design DW is unchecked against Lemma IR, LV6-U, POSTPROCESS and Lemma GB's sizes.
3. Hyp U′ (level-1 moments at 3/8) is untouched by us beyond the remark in §6.4.
4. If DW is rejected, the algebraic two-point block in the original design remains open. Its exact obstruction is a
   density-relative BV, or a Brun–Titchmarsh bound, for smooth numbers to the endogenous moduli G(P′) at level θ < 1/2
   (§6.1). We know no route to it.
5. Hyp RB itself (the ternary statement) is neither proved nor needed. The dilated ternary circle method of the task
   (major arcs at effective denominators ≤ R, but θ-denominators up to qHR) is not carried out (codex blind §2).

## 8. Answers to the task, item by item

1. **Averaging.** Yes. Averaged over the pair family and the child primes, BDH and the large sieve give the needed
   two-point statements at super-polylogarithmic precision, for the non-principal parts. The principal parts are
   second-order. The essential averaging is over pairs inside one child box (dense Δ's) and over children q.
2. **Circle method.** For the two-point problems the "principal arc" is the gcd-class/archimedean principal part. It
   does not reproduce the product of marginals exactly; the residual is second-order (Lemma PCV, codex blind §2). For
   the ternary RB we did not carry out the dilated circle method; DC makes it unnecessary.
3. **Weaker variants.** None needed for 3/8 in design 2L′, subject to the outlines. Below 3/8, DC does not help (§7.3).
4. **Algebraic block.** Clean reduction under DW (Lemma DW, Prop AB2). In the original design the obstruction is named
   in §7.4.4.
5. **Obstruction.** If the outlines fail, they fail at the items of §7.4.1. None of them is a known open problem in
   the sense of RB (smooth ternary equations with growing coefficients).

## 9. Experiments (pre-registered, `notes/pc3_prereg.md`, commit `1025bda`)

Toy, rational block only, exact smoothness. y = 150, M = 2000, S = 20, labels q ∈ (500, 1000], 20 labels per
experiment. Child values ≈ 10⁶ (u ≈ 2.9); parent values up to 8·10⁷.

### 9.1 E1: two-point structure of the child box (20 labels, 511 points each, 3·10⁵ g per label)

| statistic | value | 95% CI (label bootstrap) | band | verdict |
|---|---|---|---|---|
| P1a mean residual ⟨C_emp − C_mod⟩ | −0.0002 | [−0.0008, 0.0005] | ±0.01, CI ⊂ ±0.02 | held |
| P1b slope cov(C_emp, C_mod)/var(C_mod) | 1.18 | [1.16, 1.21] | [0.8, 1.3] | held |
| P1c excess dispersion var(C_emp) − noise − var(C_mod) | 0.0070 (var(C_mod) = 0.0150) | — | ≤ 0.5·var(C_mod) | held (barely) |
| P1d box mean ⟨C_emp⟩ | 0.0005 | [−0.0007, 0.0018] | ±0.01 | held |

Reading:
- Individual pair covariances are large at toy scale: sd 0.15 across pairs, bins from −0.28 to +0.23. They are
  explained mostly by the principal proxy (local factors including the primitivity of g, times the size density).
- The proxy under-disperses by ≈ 18%: the slope is 1.18 and the excess dispersion is about half of var(C_mod). Likely
  causes are the non-multiplicativity of local factors at p ≥ 17 and the local × archimedean interaction. Neither is
  modelled.
- Over the whole box they cancel to 5·10⁻⁴: the box is "complete mod p", so the number of box points with p | G is
  nearly independent of g. This is a toy illustration that ϕ_𝔮 ≈ 1/m. It is the DC1 premise, not evidence for any
  asymptotic claim.

### 9.2 E2: conditioning on the parent vs the unconditional box (20 labels × K ∈ {20, 80, 320})

| statistic (pooled over labels; label bootstrap CI) | K = 20 | K = 80 | K = 320 | band | verdict |
|---|---|---|---|---|---|
| mean m, z (parents per g) | 40.5, 80.5 | 40.5, 62.6 | 40.5, 39.5 | — | — |
| ϕ_ex = ϕ − 1/m (unconditional) | −0.0043 | −0.0043 | −0.0043 | — | — |
| P2a Φ_ex − ϕ_ex | −0.0013 [−0.0015, −0.0012] | −0.0014 [−0.0015, −0.0012] | −0.0013 [−0.0015, −0.0012] | ≤ 0.01 in abs. | held |
| P2c σ²_ex = σ² − 1/z | −0.0018 [−0.0021, −0.0014] | 0.0024 [0.0020, 0.0029] | 0.0029 [0.0023, 0.0036] | [−0.02, 0.05] | held |
| P2d shift ⟨E[N∣E_P]/m − 1⟩ | 0.0043 | 0.0045 | 0.0046 | ≤ 0.03 | held |
| P2e dc − ϕ (dc = E[Z(N−m)²]/(E Z·m²)) | −0.0012 | −0.0012 | −0.0012 | ≤ 0.01 | held |
| P2b Φ_ex(320) − Φ_ex(20) | −3·10⁻⁶ [−3·10⁻⁵, 3·10⁻⁵] | | | ≤ 0.01 | held |

Reading (fixed in advance):
- At toy scale the parent-conditioned relative variance of the child box, which is PC₃'s object, equals the
  unconditional one to 0.0014. The parent-weighted variance, Prop DC's object, does too.
- The coefficient scale K of the parent does not matter: a 16-fold change in |(i₀, j₀)| moves Φ_ex by 3·10⁻⁶.
- The small systematic offset (−0.0013, conditioning slightly *reduces* the variance) and the shift +0.0045 are
  third- and second-order effects of the shared g.
- The child box is sub-Poisson (ϕ_ex < 0). This is the completeness-mod-p cancellation seen in E1.
- The label's parent count is nearly Poisson: σ² − 1/z ≤ 0.003, while 1/z ≈ 0.01–0.025.
- Nothing here is asymptotic; the algebraic block is absent.

## 10. Literature

- **BDH/BV for smooth numbers.** Harper, arXiv 1208.5992, Theorems 1–2 (preprint; Theorem 2 is the BDH form quoted
  in 4.1). Drappeau, Compositio 151 (2015) 828–862, Lemme 5 (small conductors). Fouvry–Tenenbaum 1996; Granville 1993.
- **Large sieve.** Iwaniec–Kowalski, Analytic Number Theory, Theorem 7.13; maximal version by Rademacher–Menshov.
- **Local behaviour.** Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Theorem 3. de la Bretèche–Tenenbaum, "Propriétés
  statistiques des entiers friables", Ramanujan J. 9 (2005) (coprime local factors; UNVERIFIED theorem number).
  Hildebrand, Quart. J. Math. Oxford 36 (1985) (short intervals; UNVERIFIED exact range).
- **Ternary smooth equations.** Harper, Compositio 152 (2016), Cor. 1 and Theorems 1–2 (minor arcs, restriction).
  Lagarias–Soundararajan, Proc. LMS 104 (2012). Drappeau, MPCPS 154 (2013), Thm 1.3. All as checked in multidescent
  §7 and codex blind.
- **Codex correction.** Matthiesen–Wang Theorem 1.3 assumes Ψ(X, y) > X/log₈X with log₈ the eighth iterated
  logarithm, not (log X)^8. multidescent §7 should be read accordingly; still far outside our density regime.
- **Decoupling.** No prior use found by us of "parent count of a label" decoupling in descent analyses. This is a
  search result, not a novelty claim.

## 11. Reviews

[to fill: codex adversarial pc3_2; Opus proof critic; Opus literature critic]

## 12. Reproduction

```
.venv/bin/python -m pytest -q tests/test_pc3.py
.venv/bin/python -m src.pc3.experiments E1 experiments/raw/pc3_e1.jsonl 20
.venv/bin/python -m src.pc3.experiments E2 experiments/raw/pc3_e2.jsonl 20
.venv/bin/python -m src.pc3.analyze E1 experiments/raw/pc3_e1.jsonl
.venv/bin/python -m src.pc3.analyze E2 experiments/raw/pc3_e2.jsonl
```

## 13. Candidate card

- **Mechanism.** Two-level descent at 3/8 with child boxes closed by a union bound (multidescent Prop UB). The
  parent–child conditioning is removed by decoupling through the label's parent count (Prop DC), and the field
  multiplier is made linear (DW).
- **Result.**
  - Proved: Prop DC; Lemmas 0, LAT, LSQ, DW.
  - Outline: PCV, RB2, AB2 ⇒ Theorem PC₃′ (Prop UB's conclusion without Hyp PC₃/RB) ⇒ Theorem 2L′: L[3/8], conditional
    on U′.
- **Failure scope.** Not proved: the outlines (§7.4.1), DW's compatibility with the internal chain, U′. DC fails at
  depth ≥ 2 of deeper chains.
- **Next step.** Check §7.4.1(a)–(c) line by line. Re-check DW in [pp] Lemmas GB/IR/LV6. Hand DW to `uprime`.
