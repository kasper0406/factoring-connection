# c12alg: the algebraic two-point blocks of C₁/C₂ under design DW — results

Agent `c12alg`, phase 27, 2026-10-24. **v2.1** (v2 + §3.6 Theorem AB₂). v2: the codex adversarial round (MAJOR REVISION) and the Opus proof critic (Lemma η MAJOR; the rest ACCEPT or ACCEPT WITH FIXES) are applied (§10). The Opus literature critic is recorded in §10 v2.1.

Files:
- this note; pre-registration `notes/c12alg_prereg.md` (commit `20c882b`, code frozen there);
- code: `src/c12alg/dw.py` (DW objects, Lemma TP construction, exact laws), `e1.py` (E1), `dw_pipeline.py` and
  `harness_e2.py` (E2: the thmA toy pipeline with DW fields), `analyze.py`, `params.py` (leading-order margins of the
  window design, 25 checks);
- tests: `tests/test_c12alg.py` (13 tests);
- data: `experiments/raw/c12alg_e1a.jsonl`, `c12alg_e1b.jsonl`, `c12alg_e2.jsonl` and `_summary.json` files;
- codex: `codex/c12alg_1.{md,out}` (blind round, launched before any analysis was written).

Labels as in pc3: **proved** (from the project's internal unrefereed chain plus cited published results; anything
resting on Harper arXiv 1208.5992 rests on an unrefereed preprint, and is marked [H]); **sketch** (every estimate
named and checked at leading order, constants not written to the end); **heuristic**, **open**, **experimentally
supported** (toy), **UNVERIFIED** (literature statement not checked by us).

Notation: λ := log log N, D := log N/λ, L[α, c] = e^{(c+o(1))D^αλ}; "polylog" = (log N)^{O(1)} = e^{O(λ)}. The
Hildebrand–Tenenbaum saddle point is written β (α is the L-exponent).

---

## 0. Verdict

**Plain language.** On the algebraic side, design DW works as hoped. It does **not** by itself give a rigorous
algorithm below L[2/5]: the result is conditional on two unreviewed sketches, ours (Theorem AB/AB₂) and sister agent
c12rat's rational blocks.

What we have:
- DW is compatible with the whole chain (proved).
- Its two-point law is an explicit lattice-parallelogram mixture (proved).
- The child boxes are never empty on average (proved).
- The algebraic two-point blocks of C₁/C₂ hold to the precision the descent needs (sketch):
  - at **first order** on the window 3/8 < α < 2/5, where first order is enough (Theorem AB);
  - at **second order** on [3/8, 2/5), including α = 3/8 (Theorem AB₂, added in v2.1), by a cell-covariance device
    (the covariance of two one-point functions is at most the product of their first-order standard deviations).
    The device was suggested in c12rat's uncommitted v1 note.

If both sketches survive, the two-level algorithm runs in L[3/8]. If only the first-order versions survive, it runs in
L[α] for every α > 3/8.

Precisely:
1. **Window precision budget [proved].** In the design 2L_DW(α) (§1), 3/8 < α < 2/5, Prop DC_e's dead mass is
   ≤ μ₁/2 as soon as
   - **A ≤ D^{−α}λ^{C}** (first order at the child scale, u_child ≍ D^α), and
   - **B ≤ D^{7α−3}λ^{C}** (first order at the parent scale, u₂ ≍ D^{3−7α}).

   This uses Lemma η, proved here in v2 by a pointwise ratio with Hildebrand's theorem at B₂: W/μ₁ ≤ 3u₂, with **no
   polylog loss**. v1's proof lost (log N)^{2C}, which is fatal in this window (codex FATAL #1, proof critic MAJOR).
   Both conditions are weaker than "A, B ≤ η/4". Codex reached the same budget blind.
2. **DW compatibility (task 2) [proved].** With w uniform among *all* integers in (W_m, 2W_m], with no primality and
   no coprimality to g₁, every step of the chain survives after explicit repairs:
   - GB with the multiplier outside;
   - Lemma CU: F̂_N(P) is a unit mod |G(P)| whenever gcd(G(P), N) = 1, with no condition on b, g₁ or w;
   - Lemma CT: the content divides w, and "content > 1" is rejected, with probability ≤ 2^{2−d};
   - Lemma IR-DW: certificate primes dividing w are lost, but averaging over w costs only 2^{1−|𝒫_IR|};
   - Lemma W-DW: LD*'s Lemma W holds with C_W = 2;
   - DW(a) makes FM's algebraic factor pointwise, so Lemmas 8.5′/8.5H are not needed at level 1;
   - POSTPROCESS/LV6-U, the row identity, the square root in Z[β] and the factor exits are unchanged.

   One claim of codex's blind round is rebutted: it said LD* is "not verbatim valid", but LD* Step 4 goes through with
   the all-residue sum.
   E2 ran the thmA toy pipeline with composite w on 32 instances and found no algebraic violation.
3. **Two-point smoothing DW(b) (task 1).**
   - **Lemma TP / TP-W / TP-q [proved].** For filtered g, the law of (F(P₁), F(P₂)) over (w, c) is within TV
     3N^{−η₀/2} of an explicit mixture of uniform measures on the DW-L lattice intersected with parallelograms. Under
     the root event q | F, the bound is q·3N^{−η₀/2}, on the sublattice x₁ ≡ λ^d x₂ (mod Δ/q). The TV bound does not
     degrade for nearly parallel pairs. That rests on two elementary lemmas: a maximal-minor bound (Lemma MM) and a
     Bezoutian bound (Lemma BZ).
   - **Theorem AB [sketch; no conceptual gap found by either reviewer; several bookkeeping errors of v1 corrected].** In the window, the algebraic parts of C₁ and C₂ (law of total variance with
     respect to g) satisfy
     - Σ(w/W)|ϕ^alg| ≤ D^{−α}λ^C(1 + A^rat);
     - Σ(w/W)|σ^alg| ≤ D^{7α−3}λ^C(1 + B^rat).

     The inputs are:
     - Lemma TP;
     - the BV theorem for smooth numbers [H, Thm 1] on fibres, averaged over pairs with uniform multiplicities;
     - Hildebrand's short-interval asymptotic;
     - Hildebrand–Tenenbaum Thm 3 and de la Bretèche–Tenenbaum Thm 2.4;
     - [pp] Lemma DS (generalised to long progressions);
     - a weight-domination lemma (Lemma WD, proved here under DW);
     - an explicit treatment of thin pairs and skewed labels.
   - **Skewed labels, a by-product [proved, given the pair upper bounds of Lemma PU].** Labels of skewness > T have
     total weight·ϕ ≤ D^{−α}. This holds for the whole ϕ, rational part included, so it disposes of pc3 §7 item 2 for C₁.
4. **One-point child lower bound (task 3) [proved].** m_𝔮 ≥ m₂(q)(log N)^{−C} ≥ e^{λ²}/(log N)^C ≫ 8/η, for **every**
   label. This includes skewed labels and fields with q | G(r, 1).
   - It requires rescanning S₂ with DW's height (codex).
   - It uses DW(a) at the child level: the root coset conditioning costs a factor q in TV, which N^{η₀} absorbs.
   - The rational side uses a per-point fibre count.
   - uprime's "routine" is confirmed, with two corrections: the child filter |G(P′)| ≥ M|P′|N^{−η₀/4} must be added to
     Good(P′), and the scan must use DW's height.
5. **Achieved α.** **No rigorous α < 2/5.** v2.1: conditionally on Theorem AB₂ (sketch, §3.6) and Hyp RAT₂ (the
   rational blocks at second order, which c12rat's unreviewed v1 claims), the same assembly gives **L[3/8]**
   (Theorem 2L-DW₂, §5).
   - **Theorem 2L-DW(α) [proved conditional; §5].** For every fixed α ∈ (3/8, 2/5), FACTOR_2L in design 2L_DW(α)
     factors every N in expected time L_N[α, c(α)]. It is conditional on:
     - (i) the internal chain;
     - (ii) **Hyp RAT(α)**: the rational parts of C₁ at precision D^{−α}λ^{C} and of C₂ at D^{7α−3}λ^{C}, for labels
       of skewness ≤ T, both with the explicit smooth algebraic weights;
     - (iii) Theorem AB (sketch).
   - Hyp RAT(α) is c12rat's object, pc3 §7 items 1 and 4, at first order only.
   - If c12rat proves RAT(α) and Theorem AB survives review, a rigorous (internal) L[α] for every α > 3/8 follows.
6. **Toward 3/8 [analysis; superseded in v2.1 by §3.6].** The window argument needs α > 3/8 in exactly two places:
   - the first-order budget D^{−α} ≪ η;
   - √(AB) ≪ η.

   At α = 3/8, Theorem AB's local factor (Lemma LF) and its BV/BDH errors are already good enough: they are of second
   order or super-polylogarithmically small. The archimedean comparison (Lemma AR combined with the first-order
   short-interval asymptotics) is not: it must be done to second order. That is pc3 §7 item 4 (Saias + Lemma DS), the
   same obstacle as on the rational side.
7. **Experiments (pre-registered; toy).** 10 of 12 predictions held. P1c failed: the registered band was miscalibrated,
   because the TV bound is 20–50× pessimistic. P2a failed in E2b: one instance hit the field cap, and post hoc it was
   factored correctly with a larger cap.
   - E1a: TV ≤ bound in 72/72 exact cases.
   - E1b: the DW pair law matches "lattice gcd classes × archimedean sizes". Σ J/Σ J_full = 1.0001 [0.9995, 1.0006],
     with slope 0.9999.
   - E2: 31/32 factorisations, 0 violations of the row identity, the congruence, Lemma CU (5 911 atoms) or the forced
     residue. Content rejections are 5 392 vs 5 398.6 predicted (z = −0.09; v1's 5 399.3 came from a truncated Möbius
     sum). DD fractions are 0.46 (ω = 2) and 0.72 (ω = 3).

**Probabilities** (ours; reviewers in brackets):
- Theorem AB as sketched can be written out without a new idea: 55–65% [codex adversarial 40%; proof critic 65%].
- DW is compatible with the chain: 90% [codex 93%].
- Theorem 2L-DW(α) is a correct conditional theorem after the v2 repairs: 65% [codex 70%; proof critic 65–70%].
- A rigorous fixed α < 2/5, including the rational blocks: 25–35% [codex blind 30%].

---

## 1. The window design 2L_DW(α) and the precision budget

### 1.1 Parameters

Fix rationals α ∈ (3/8, 2/5) and a, κ, d₀, b₂ > 0.
- **Degrees and special-q.** d := least odd ≥ max(3, d₀D^{1−2α}); k := max(⌈κD^α⌉, d + 1); ν := 1/(4k); μ := (1 − ν)/d
  (μ_opt, since k < (d² + 1)/4); Y, J and ℓ ≈ N^{1/(2k)} = L[1 − α] as in lv5unif.
- **Smoothness bounds.** B := ⌊e^{aD^αλ}⌋; children are the primes in (B, B₂], B₂ := ⌊e^{b₂D^{4α−1}λ}⌋ (codex pc3_1 §4: the
  upper end of the feasible interval 2 − 4α ≤ ζ₁ ≤ 4α − 1, which opens exactly at 3/8).
- **Height margin.** η₀ := 1/⌈D^{3/2−3α}⌉, i.e. η₀ log N ≍ D^{3α−1/2}λ. This lies strictly between max(α, 4α − 1) and
  2α, so N^{η₀} beats every L[α]^{O(1)} and every B₂^{O(1)}, and η₀ = o(μ). (Codex proposed η₀ ≍ C D^{4α−2}; either
  works. uprime's D^{−1/2} does not: L[1/2] < B₂ in the window.)
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
  σ₂ = (α/a)(4/d₀ + d₀b₂). This agrees with codex.
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
D^{8α−3} margins (codex #1 FATAL, proof critic MAJOR). v2 bounds the ratio pointwise, with **no polylog loss**.

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
  Ψ(J/q, B₂) = (|J|/q)ρ(u_x − v_q)(1 + O(1/u₂)) and Ψ(J, B₂) = |J|ρ(u_x)(1 + O(1/u₂)), where v_q = log q/log B₂ ≤ 1. The
  hypotheses hold: B₂ ≥ exp((log log x)^{5/3+ε}), and relative length 1/(64dK₁) ≥ B₂^{−5/12}.
- **Dickman.** −ρ′/ρ is increasing, and equals ξ(u)(1 + O(1/u)) at u (de Bruijn). So
  ρ(u − v)/ρ(u) ≤ e^{vξ(u)(1+O(1/u))} for v ≤ 1. Hence
  Σ_q(1/q)ρ(u_x − v_q)/ρ(u_x) ≤ (1 + o(1))∫_{log B}^{log B₂}e^{(ξ/log B₂)t}dt/t + O(1) ≤ (1 + o(1))e^{ξ}/ξ + O(1)
  ≤ (1 + o(1))u_x, using e^ξ = 1 + uξ. And u_x ≤ u₂.
- **The part |n| < X_F/K₁.** Its Ψ-mass relative to Ψ(I) is ≤ 64dK₁·K₁^{1−β₂}/K₁·C ≤ CdK₁^{−3/4}. On it m ≤ m̄ ≍ D^{2−4α}, so
  it contributes ≤ CdD^{−3}D^{2−4α} = o(1).
- **Total.** R(I) ≤ (1 + o(1))u₂ + o(1) ≤ 3u₂. ∎

The deterministic cap m̄ ≍ D^{2−4α} would not suffice: D^{−α} ≪ D^{−(2−4α)} fails for α < 2/5 (codex). The input
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
- Codex's example (pc3_1 out §1) shows that w coprime to g₁ restricts F(P) mod p for p | gcd(t, g₁). The law is then not
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

This replaces [pp] lem:coprime, which used w prime > t. No hypothesis on b or g₁ is needed; codex reached the same
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
(g, 𝔮) contribute a dead environment of conditional probability ≈ 1/q ≤ 1/B, which is negligible against η (codex, same).

### 2.4 Lemma CT (content) and Lemma IR-DW [proved]

**Lemma CT.** For a prime p ∤ N, p divides the content of F iff p | w and p | c_i for every i. In particular the content
e (coprime to N) divides w, and
Pr_{w,c}[content > 1, coprime to N] ≤ Σ_{p≤2W_m}(1/p + 1/W_m)(1/p + 1/(2H + 1))^d ≤ 2^{2−d}. The sum is over p | w, so
p ≤ 2W_m. For p > 2H + 1, "p | c" means c = 0, of total probability ≤ (2H + 1)^{−d}.
*Proof.*
- If F ≡ 0 (mod p), then wN = F(g₀, g₁) ≡ 0, so p | w (as p ∤ N).
- If p | w: F ≡ G·R_c, and G ≢ 0 (mod p), so F ≡ 0 iff R_c ≡ 0, i.e. iff p | c for every coordinate. ∎

SAMPLE_FIELD-DW: gcd-check the content e against N (a proper divisor is a factor exit, since e ≤ |f_d| < N). If e > 1,
reject. This replaces "content ∈ {1, w}; reject if w". E2 observed 5 392 content rejections against 5 398.6 predicted (z = −0.09). *This probability is not N^{−Ω(η₀)}* (codex). Rejection is therefore used
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

A worst-case w divisible by every certificate prime exists (codex), so averaging over w is necessary. Codex's sketch is
the same. E2: 0 violations of claim 2.

**Lemma reject-DW.** Under the raw DW law, the total rejection probability is
Pr[ℛ] ≤ exp(−log H/(26d log log H)) + 3/(2H + 1) + 2^{2−d}.
It collects: no certificate; N | f′(m) (unchanged: coefficient g₁ of c_{d−1}); degree drop; content > 1. In the window
this is e^{−Θ(D^{4α−1})} + e^{−Θ(D^{1−2α})} ≤ γ₂/4.

### 2.5 Lemma W-DW (the second moment's multiplier step) [proved]

**Lemma W-DW.** For 1 ≤ t ≤ W_m and every Φ₀: Z/t → [0, ∞), E_w Φ₀(wu mod t) ≤ (2/t)Σ_{c mod t}Φ₀(c), for every u with
gcd(u, t) = 1. (v1 said "every u", which is false for u ≡ 0; proof critic. In LD*, u = F̂_N(P) is a unit by Lemma CU.)
*Proof.* w ↦ wu permutes Z/t. And w is uniform on W_m consecutive integers, so each class mod t has probability
≤ (⌊W_m/t⌋ + 1)/W_m ≤ 2/t. ∎

*Use in LD*-U′ (answering codex OPEN (f)).* The steps of [pp] Thm LD / lv5unif LD*-U are:
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
- **Range restatement (proof critic).** uprime's U1₂, DS-U′ and CS′ are stated for log z/log B₂ ≥ D^{1/4}. In the window
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
G(P), F(P) and the child primes) and transcript invariance (T) are unchanged (codex agrees). E2 ran the thmA pipeline
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

E1a checked this exactly in 72 cases: TV ≤ bound always. The actual TV is 2–4% of the bound, so the constant
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
   - Hence |σ|_∞ ≤ 8d⁴2^dW_mN^{η₀/2}, **independently of the angle between P₁ and P₂**. This is codex's (1), obtained
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
- q² ∤ Δ in the child box (Δ = qΔ′ with Δ′ < q), so no mod-q ratio remains (codex). ∎

So the **DW(b) TV statement is proved**. The proof critic's M4 and pc3 §7 item 3 asked for this, and it holds uniformly,
nearly parallel pairs included.

### 3.3 Theorem AB (the algebraic blocks in the window) [sketch, no gap known]

**Theorem AB.** In design 2L_DW(α), 3/8 < α < 2/5, there are constants C and N₀ such that, cell by cell:
- (i) Σ_{𝔮: s(𝔮)≤T}(w_𝔮/W)|ϕ^alg_𝔮| ≤ D^{−α}λ^C(1 + Σ_{𝔮: s(𝔮)≤T}(w_𝔮/W)ϕ^rat_𝔮);
- (ii) Σ_𝔮(w_𝔮/W)|σ^alg_𝔮| ≤ D^{7α−3}λ^C(1 + Σ_𝔮(w_𝔮/W)σ^rat_𝔮);
- (iii) Σ_{𝔮: s(𝔮) > T}(w_𝔮/W)ϕ_𝔮 ≤ D^{−α}, where s(𝔮) = |e₂|/|e₁| is the skewness of the reduced basis and
  T := Ξ³(log N)^{C}D^α.

Here Ξ := exp(C₁(D^{3α−1} + D^{1−2α} + D^{2α−1/2})λ²) collects every L^{o(1)} loss used below:
- the large-modulus loss of Lemma PU, e^{O(D^{3α−1}λ)};
- the spread of ϱ over the physical sizes |P′| in a box, e^{O(D^{1−2α}λ)}: d·log(S₂s) ≍ D^{1−α}λ, times
  1 − β ≍ D^{−α};
- **(v2, proof critic)** the spread over |G(P′)| down to the filter threshold M|P′|N^{−η₀/4}: e^{O((1−β)η₀ log N)} =
  e^{O(D^{2α−1/2}λ)}. This is the dominant term in the window.

Since max(3α − 1, 1 − 2α, 2α − 1/2) < α for 1/3 < α < 1/2, Ξ^{O(1)} ≤ B^{1/3} and Ξ^{O(1)} ≪ e^{cD^α/λ²}
(`params.py` checks this).

The proof assembles the lemmas below. Every input is a published theorem, an internal lemma, or a lemma proved in this
section. It is a sketch because the grid/step bookkeeping of Lemmas FB and PP is not written to the last constant.

**Lemma FB (fibre BV) [proved [H]].** Let ν := Unif(L′ ∩ Π) with L′ = {x₁ ≡ κx₂ (mod Δ′)}, κ a unit, and Π a
parallelogram with x₂-fibres I(x₁). Write E(t; Δ′, a) := Ψ(t; Δ′, a) − Ψ(t; gcd(a, Δ′))/φ(Δ′/gcd(a, Δ′)), where
Ψ(t; δ) counts B-smooth n ≤ t with gcd(n, Δ′) = δ. Put
- 𝔐(Π) := (1/#)Σ_{x₁∈𝒮}m(I(x₁); gcd(x₁, Δ′)), where m(I; δ) := Ψ(I; δ)/φ(Δ′/δ);
- ℰ(Π) := (2/#)Σ_{x₁∈𝒮} max_{t∈𝒳_K}max_a|E(t; Δ′, a)|.

Then:
- |ν(x₁ ∈ 𝒮, x₂ ∈ 𝒮) − 𝔐(Π)| ≤ ℰ(Π) + (grid error).
- The grid error is ≤ 4·(1/K)·𝔐-scale·(log N)^C, by [pp] Lemma DS on the grid cells. The grid is
  𝒳_K := {x₀(1 + 1/K)^i}, and modulus Δ′ ≤ B^{C₀} and cell length x/K ≥ xB^{−C₀} are needed.
- For a family of pairs with moduli Δ′ ≤ Q ≤ √Ψ and multiplicity ≤ μ_max per value of Δ′, Harper's Theorem 1
  (BV for smooth numbers, max over residues) at each grid point gives
  Σ_{pairs}#·ℰ/(#{x₁ ∈ 𝒮}) ≤ 2μ_max|𝒳_K|[Ψ(x)(e^{−cu/log²(u+1)} + B^{−c}) + Q√Ψ(x)(log x)^{7/2}].
- Non-reduced classes reduce to modulus Δ′/δ at endpoint t/δ, with the multiplicity in Δ″ = Δ′/δ inherited from Δ′.
  Summing Harper's Theorem 1 at (t/δ, Δ″ ≤ Q/δ) over B-smooth δ ≤ Q, with Ψ(t/δ) ≪ δ^{−β}Ψ(t) (dBT 2.4(i)) and
  Σ_{δ≤Q}δ^{−β} ≤ Q^{1−β}(1 + log Q) = e^{O(λ)} (as (1 − β)log Q = O(λ) for Q = L[α]^{O(1)}), costs only a polylog.
  The BV saving absorbs it. v1's "Σ_δ δ^{−β} < ∞" was false for β < 1 (codex #3, proof critic). For parents
  Q = L[1 − α] and Q^{1−β₂} ≤ Ξ₂.

*Proof.* Expand ν(both) = (1/#)Σ_{x₁∈𝒮}Ψ(I(x₁); Δ′, κ^{−1}x₁), and bound |E| at the endpoints of I(x₁) by the nearest
grid points plus a DS term. ∎

The averaging is over **pairs**, as pc3 found for the rational block. The algebraic modulus Δ′ = |ij′ − ji′| depends only
on the coordinates, and its multiplicity in a dyadic coordinate class (|(i, j)| ≍ 2^{c₁} ≤ |(i′, j′)| ≍ 2^{c₂}) is
≤ 2^{c₁+c₂+6}. This is ≪ #pairs/Q_c, uniformly (pc3 §3, refined per dyadic class). So the relative error per class is
≤ K·K₀·(log N)^C(e^{−cD^α/λ²} + Q/√Ψ), where K₀ bounds the relative fibre length of non-thin pairs.

**Lemma LF (local factor) [proved].** Let π*(δ) := Pr(gcd(n, Δ′) = δ) in the multiplicative model p | n with probability
p^{−β}, restricted to p ≤ B. Let π₀(δ) := φ(Δ′/δ)/Δ′. Then:
- LF := Σ_δ π*₁(δ)π*₂(δ)/π₀(δ) = ∏_{p|Δ′, p≤B}(1 + (p^{1−β₁} − 1)(p^{1−β₂} − 1)/(p − 1)) for squarefree Δ′, with the
  analogous geometric factors otherwise;
- primes p | Δ′ with p > B contribute the factor p/(p − 1), since smooth values avoid p (codex). Their product is
  ≤ e^{ω(Δ′)/B};
- 1 ≤ LF ≤ 1 + D^{−2α}λ^C, using Δ′ ≤ L[α]^{O(1)}: every p | Δ′ has (1 − β)log p = O(λ), and p ≤ (log N)^C has
  p^{1−β} − 1 = O(D^{−α}λ²).

*Proof.* Multiplicativity, then the per-prime identity τ²/p + (1 − τ/p)²/(1 − 1/p) = 1 + (τ − 1)²/(p − 1) (and its
two-variable form). ∎ (`test_local_factor_product_formula`.)

**Lemma AR (archimedean factor) [proved].** Let ν be a mixture of the measures of Lemma TP-W (or Cor TP-q). Then:
- each marginal satisfies ν(|X_i| ≤ e^{−s}S_i) ≤ 8de^{−s}, where S_i := max|X_i|. This is Lemma MM's column claim:
  each component's marginal density is ≤ 2/((2H + 1)|G_i||P_i|^{d−1}).
- Put Z_i := ϱ(|X_i|), ϱ(t) := Ψ(t, B)/t. Then |E[Z₁Z₂]/(EZ₁·EZ₂) − 1| ≤ C((1 − β)²(log d)² + 1/u²). This is second
  order: the 1/u errors of HT enter only through Var(Z_i), not linearly (proof critic). The atom at X_i = 0 has mass
  ≤ 2/((2H + 1)|G_i||P_i|^{d−1}) and is removed first (codex).

*Proof.*
- HT Theorem 3 gives ϱ(S/c) = c^{1−β}ϱ(S)(1 + O(1/u)) for 1 ≤ c ≤ B; iterate for c ≤ B^{C}, and use
  de la Bretèche–Tenenbaum Thm 2.4(i) beyond.
- So Z_i = ϱ(S_i)Y_i(1 + O(1/u)) with Y_i := (S_i/|X_i|)^{1−β}.
- Cov(Z₁, Z₂) ≤ √(VarZ₁·VarZ₂), and VarZ_i/(EZ_i)² ≤ C(VarY_i + 1/u²). Moreover VarY ≤ E[(Y − 1)²] ≤ C(1 − β)²E[log²(S/|X|)e^{2(1−β)log(S/|X|)}] ≤ C(1 − β)²(log 8d)²,
  by the exponential tail. ∎

**Lemma PP (principal part ≈ product of marginals, first order) [sketch].** For a non-thin pair (sin∠(P₁, P₂) ≥ 8d/K₀,
K₀ ≤ B^{1/3}) and every filtered g:
- 𝔐 = A₁A₂⁰·LF·AR·(1 + O(λ^C/u)), where A₂⁰ is the principal part of A₂. Then |A₂ − A₂⁰| is again an ℰ-type error, so
  𝔐 = A₁A₂(1 + O(λ^C/u)) + O(ℰ).

*Proof sketch.*
- Cut the projection of Π into intervals J of relative length 1/K. On J, the fibre I(x₁) moves by relative O(1/K).
- On each J, and on fibres of relative length ≥ 1/K₀ ≥ B^{−1/3} ≥ B^{−5/12}, Hildebrand's short-interval theorem
  (Ψ(x + z) − Ψ(x) = zρ(u)(1 + O(log(u + 1)/log B)) for x B^{−5/12} ≤ z ≤ x) applies. Combined with Möbius over e | m,
  e ≤ L[α]^{O(1)}, and HT Theorem 3, it gives, in every gcd class:
  Ψ_m(I/δ) = (|I|/δ)ϱ(x/δ)∏_{p|m,p≤B}(1 − p^{−β})(1 + O(λ^C/u)).
- The Möbius error sum is ≤ (λ^C/u) × the main term, since ∏_{p|m}(1 + p^{−β})/(1 − p^{−β}) ≤ λ^C for m ≤ L[α]^{O(1)}.
- Summing over x₁ ∈ 𝒮 ∩ J in each gcd class gives the same structure at the x₁ scale.
- The result factorises as (Σ_δ π*₁π*₂/π₀) × (E_ν[ϱϱ]/(E_νϱ E_νϱ)) × (1 + O(λ^C/u)).
- **Short fibres and corners (v2; codex #4).** Fibres of a parallelogram shrink to 0 at its vertices. Let
  ε := D^{−α}/(Ξ(log N)^C). The x₁ whose fibre is shorter than ε × (maximal fibre) carry ν-mass ≤ 4ε. Their
  smooth-weighted contribution is bounded by [pp] Lemma DS on fibres enlarged to length max(|I|, XB^{−C₀}) (the
  enlargement costs ≤ κB^{−C₀}) and on the x₁-interval. This gives ≤ (4ε + κB^{−C₀})(log N)^CΞ × A₁A₂ ≤ D^{−α}A₁A₂.
  Every remaining fibre has relative length ≥ ε/κ ≥ 1/(K₀D^α(log N)^C) ≥ B^{−5/12}, and moves by relative ≤ K₀/(εK)
  across a grid cell; K is chosen accordingly.
- **Parent version (v2; codex #5, proof critic).** At y = B₂ the moduli reach L[1 − α] ≫ B₂^{O(1)}, so HT cannot be
  iterated over all divisors.
  - Truncate gcd classes δ and Möbius divisors e at E := B₂^{1/2}, i.e. v ≤ 1/2, where
    ρ(u − v)/ρ(u) = e^{vξ}(1 + O(1/u₂)).
  - The tail weight is ≤ Σ_{δ|Δ, δ>E}δ^{−β₂} ≤ E^{−1/4}∏_{p|Δ}(1 + 2p^{−1/2}) ≤ B₂^{−1/8}e^{O(D^{(1−α)/2}λ^{1/2})}, which is
    negligible since (1 − α)/2 < 4α − 1.
  - The same bound handles e: Möbius over e | m, e > E, has weight ≤ Σ e^{−β₂}e^{O(v_e²/u₂)}, and
    log e ≤ log Δ ≤ D^{1−α}λ ≪ u₂log²B₂.
  - The quadratic error v²/u₂ that codex flagged arises only beyond E, where it is killed by e^{−β₂log e}. ∎

*Why first order suffices and the κ-amplification disappears.* Codex (c12alg_1 §1) warned that with Saias-type
endpoint differences the archimedean error grows like 𝒦²/u, with a 1/𝒦 thin tail, which balances only at u^{−1/3}.
Hildebrand's theorem is a *direct* short-interval asymptotic with relative error O(log(u + 1)/log y) for every relative
length ≥ y^{−5/12}. So the per-pair relative error is O(λ^C/u) **uniformly in the thinness** up to K₀ ≤ B^{1/3}. Only the
BV error terms are amplified by K₀ (Lemma FB), and those carry the super-polylog saving. This is codex's "sketch:
Hildebrand could replace Saias", carried out.

**Lemma PU (pair upper bounds with loss Ξ) [proved; v2 repairs].** For every label, every filtered g and all P₁ ≠ P₂
in the box:
- (alg) A₁₂(g) ≤ Ξ(log N)^C ϱ(S₁(g)/q)ϱ(S₂(g)/q), with the g-dependent scales S_i(g) := H|G_g(P_i)||P_i|^{d−1};
- (rat) Pr_g[Rat₁ ∧ Rat₂] ≤ Ξ(log N)^C ϱ(M|P₁|)ϱ(M|P₂|).

For parents the same holds with Ξ₂Ξ_r := e^{O(D^{2−5α}λ + D^{1−2α}λ²)} and ϱ at B₂ and B.

*Proof.* [pp] Lemma DS, generalised to large moduli and long progressions: count n = d₁m with d₁ ∈ (D/B, D] smooth, where
D := min(#terms, X/B^{C₀}).
- Reduce first by e = gcd(residue, modulus), as in DS. Then n/e has the lower endpoint X/(eB^{C₀}), and D/B ≥
  X/(B₂^{C}L[α]^{C}) stays in the range of Lemma U1.
- This gives (2A/Δ)·8log₂(2B)·ϱ(D/B), and ϱ(D/B) ≤ ϱ(X)(XB/D)^{1−β}(log N)^C by HT Theorem 3 iterated. The number of
  steps is ≤ log(X/D)/log B, each costing 1 + O(1/u).
- (alg) On a fibre of the Cor TP-q model with modulus Δ′:
  - every fibre shorter than X·B₂^{−3} is first enlarged to that length. The total enlargement is ≤ κB₂^{−3} ≤ 1 relative
    to the area, since κ ≤ 16dqS₂⁴ ≤ B₂²;
  - then (XB/D)^{1−β} ≤ (BB₂³S₂²)^{1−β} ≤ Ξ;
  - the x₁-count uses DS on the projection, an interval of length ≥ X₁/(4d);
  - max fibre ≤ 2 × average fibre.
- (rat) On the g-line through a fixed G(P₁): modulus qΔ′, M/|P₁| terms. Here (XB/D)^{1−β} ≤ (BV₂²)^{1−β} ≤ Ξ.
- Here Ξ ≥ exp((1 − β)log(B₂^{O(1)}L[α]^{O(1)})), and (1 − β)log B₂ ≍ (λ/(aD^αλ))b₂D^{4α−1}λ.
- *Scales (proof critic).* The algebraic scale is S_i(g), not the envelope. Statements (iii) and the assembly use a
  dyadic split in |G_i|/(M|P_i|) ≥ N^{−η₀/4}. The dyadic sum Σ_τ τ^{4(1−β)−1} converges, and the max/typical ratio is in
  Ξ's third term. ∎

**Lemma TH (thin pairs).** In a box of skewness s, the fraction of pairs with sin∠ < 8d/K₀ is ≤ 30·8d·s/K₀ + 3/S₂².
*Proof.*
- |P′||P″| ≤ (8/√3)sqS₂², so sin∠ = q|Δ′|/(|P′||P″|) ≥ √3|Δ′|/(8sS₂²).
- The number of coordinate pairs with |Δ′| = k is ≤ Σ_{(i,j)}2(2S₂/|(i, j)|_∞ + 2) ≤ 48S₂² for each k.
- Divide by ≥ S₂⁴ pairs. ∎

**Lemma WD (weight domination) [proved].** For each prime q ∈ (B, B₂] and all r, s mod q: w_{(q,r)} ≤ (log N)^C w_{(q,s)}.
In particular the residues r with skewness s(q, r) > T carry a weight fraction ≤ C(log N)^C/T of their window.
*Proof.*
- E[Z_{(q,r)}] = Σ_ℓΣ_{P∈𝒫, a≡rb (q)}Pr[Acc_ℓ(P) ∧ q child].
- Algebraic factor: by DW(a)⁺ with the root coset (upper bound by [pp] DS at B₂; lower bound by HT Theorem 3 and the
  conditional SQF subtraction ≤ (1/q)(Σ_p p^{−2} + 2e(d − 1)/q)), it is ≍_{polylog} (1/q)ϱ₂(X_F/q), uniformly in P ∈ 𝒫.
  The b-window keeps X_F(P) within a factor 2^{O(d)} of X_F, which moves ϱ₂ by 1 + o(1).
- Rational factor: by FIB and DS it is ≍_{polylog} ϱ(X_r)/ℓ.
- Count: #{P ∈ 𝒫 : a ≡ rb (q)} = #𝒫/q·(1 + O(qd log V₁/V₁)), uniformly in r (V₁ ≫ q, b prime ≫ q).
- The residues with λ₁(Λ_{(q,r)})² < (2/√3)q/T number ≤ (6/π)q/T + 1.
- *Per cell (v2; proof critic, codex #9).* Cells are functions of g alone. Within a cell e, the parent filter holds for
  P in a union R_e of angular sectors of size ≍ V₁. Instead of a uniform pointwise comparison over P (false for P nearly
  parallel to the cell's direction), sum over P ∈ R_e ∩ 𝒫 with a ≡ rb (q). This count is #(R_e ∩ 𝒫)/q·(1 + O(q d log V₁/V₁)),
  uniformly in r, because R_e is a union of sectors of radius V₁ ≫ q. So WD holds cell by cell. ∎

Codex flagged weight domination as OPEN and asked for it per prime window; DW(a)⁺ supplies the pointwise algebraic upper
bound it lacked.

### 3.4 Proof of Theorem AB (assembly) [sketch]

Fix α and choose the auxiliary parameters:
- T := Ξ³(log N)^{C}D^α (skewness threshold);
- K₀ := 240d·Ξ⁴T(log N)^{C}D^α ≤ B^{1/3} (thin threshold; v2: v1's choice left a factor d(log N)^C, codex #2, proof
  critic);
- K := K₀D^α(log N)^{C_DS}/ε (grid; ε as in Lemma PP).

Within a box, the largest and the typical values of A₁A₂ (and of ϱϱ) differ by ≤ Ξ (Ξ's second term). Every comparison
below is made after paying this factor.

**(iii) Skewed labels.** For every label, Lemma PU gives
E[N_𝔮² | R] ≤ m + Σ_{pairs}Ξ²(log N)^Cϱϱϱϱ ≤ m + Ξ³(log N)^{C}m²,
using m ≥ Σϱϱ/(log N)^C (§4) and the max/typical factor. So ϕ_𝔮 ≤ Ξ³(log N)^C. Lemma WD (per cell) then gives
Σ_{s>T}(w/W)ϕ ≤ Ξ³(log N)^{2C}/T ≤ D^{−α}.

**(i) Labels with s ≤ T.** For each g, split the pairs:
- *Thin pairs (Lemma TH).* Fraction ≤ 240dT/K₀ + 3/S₂². With Lemma PU and the max/typical factor, their contribution to
  |ϕ^alg| is ≤ (240dT/K₀)Ξ³(log N)^C ≤ D^{−α}/Ξ.
- *Non-thin pairs.* A₁₂ − A₁A₂ = (𝔐 − A₁A₂) + (A₁₂ − 𝔐), up to TV errors 7qN^{−η₀/2}.
  - Lemma PP: |𝔐 − A₁A₂| ≤ (λ^C/u)A₁A₂ + O(ℰ).
  - Summed against Rat₁Rat₂ and divided by m², the first term gives ≤ (λ^C/u)E[m(g)²]/m² = (λ^C/u)(1 + ϕ^rat_𝔮).
  - The ℰ-terms: Σ_{pairs}E_g[Rat₁Rat₂ℰ] ≤ max_{pair}Pr_g(Rat₁Rat₂)·Σ_{pairs}sup_gℰ.
    - sup_g ℰ ≤ 2A₁(Δ′/ℓ_fib)ℰ*(Δ′; window) (Lemma FB) depends on g only through the scale window. The grid 𝒳_K is
      common to all g and all pairs.
    - Summing over the (coordinate class × scale window) subclasses, Σ_c Q_cμ_max(c) ≤ 2⁹#pairs, gives
      Σ_{pairs}sup_gℰ/(A₁A₂) ≤ 2^{10}K₀|𝒳_K|(e^{−cD^α/λ²} + Q(log N)^{4}/√Ψ)·#pairs.
    - With Lemma PU (rat) and the max/typical factor, the ℰ-terms are ≤ Ξ²(log N)^C K K₀e^{−cD^α/λ²}·m² ≤
      e^{−c′D^α/λ²}m².
  - TV errors: ≤ #pairs·7qN^{−η₀/2} ≤ N^{−η₀/3}m², since m₂ ≥ S₂²ϱϱ and N^{η₀} beats B₂ and L[α]^{O(1)}.
- Multiply by w/W and sum. 1 + ϕ^rat_𝔮 enters linearly, so its weighted average is all that is needed.

**(ii) Parents.** The same steps, with these differences:
- the pairs are (ℓ₁, P₁), (ℓ₂, P₂), P_i ∈ 𝒫 ∩ Λ_𝔮, and P₁ ≠ P₂ automatically (G(P)/ℓ₁ smooth excludes a second ℓ₂ > B);
- **same-ℓ pairs (v2; codex #7).** If ℓ₁ = ℓ₂ = ℓ, then ℓ | det(P₁, P₂), and the rational weight is ≍ 1/ℓ rather than
  1/ℓ². These pairs are not averaged by BV. Their total contribution to E[Z²] is bounded by Lemma PU: per ℓ, there are
  ≤ #(𝒫 ∩ Λ_𝔮)²/ℓ·(1 + o(1)) pairs of weight ≤ Ξ₂Ξ_r(log N)^C ϱ⁴/ℓ. Summed over ℓ ∈ J this is ≤ Ξ₂Ξ_r(log N)^C z²/π_J.
  That is negligible, since π_J = L[1 − α]^{1−o(1)};
- for ℓ₁ ≠ ℓ₂, a pair with ℓ₁ | det has P₂ ∈ Λ_{ℓ₁}(g), hence ℓ₁ | G(P₂) = ℓ₂·(B-smooth), which is impossible. So
  the rational weight of the different-ℓ class is ≍ 1/(ℓ₁ℓ₂) uniformly, and BV is applied within that class;
- the modulus is Δ₁₂/q ≤ 2V₁²/q = L[1 − α] ≤ √Ψ(X_F/q, B₂), with multiplicity ≤ 8#(𝒫 ∩ Λ_𝔮) per value (b-window: ≤ 4
  partners per line);
- smoothness is over B₂ with u₂ ≍ D^{3−7α};
- the first-order error is λ^C/u₂ = D^{7α−3}λ^C;
- the losses are Ξ₂Ξ_r ≤ e^{O(D^{1−2α}λ)}, below the BV saving e^{−cD^{3−7α}/λ²} (as 1 − 2α < 3 − 7α) and inside
  Hildebrand's range at B₂ (as 1 − 2α < 4α − 1);
- there are no skewed labels for parents (the parent set is not a reduced-basis box), and the thin parent pairs have
  fraction ≤ C/K₀.
- Acceptability (squarefree large part, children ∤ bf_d disc(f)N) is removed by a two-point version of Lemma SQF with
  loss Ξ₂: the removed mass is ≤ (log N)²Ξ₂/B ≪ D^{7α−3}. ∎

### 3.5 What does not extend to 3/8

At α = 3/8:
- **Fine at 3/8.** The precision needed is η ≍ D^{−3/8} = 1/u. Lemmas TP, FB (super-polylog saving), LF (second order,
  ≤ D^{−3/4}λ^C), PU, TH, WD and part (iii) still work.
- **Lemma PP fails.** Its first-order error O(λ^C/u) is ≍ η·λ^C, so it does not suffice: the archimedean factor must be
  computed to relative error o(η). That requires either:
  - (a) a second-order short-interval asymptotic for smooth numbers in gcd classes, i.e. Saias's Λ(x, y) with derivative
    bounds and the κ-uniformity Hildebrand provides only at first order; or
  - (b) a cancellation between the first-order terms of 𝔐 and A₁A₂. Both are expanded around the same saddle point, so
    the leading 1/u terms of Hildebrand's error may cancel in the ratio. We did not prove this. It is pc3 §4.3(ii) and
    §7 item 4.
- **√(AB) ≪ η fails** with the present first-order bounds. At 3/8, Prop DC needs A + √(AB) ≤ η/2, a statement
  about constants.

So **the present first-order bounds do not close the endpoint** (codex #12). Possible routes are small constants,
cancellation between the first-order terms of 𝔐 and A₁A₂, or second-order (Saias-type) archimedean asymptotics.
Lemma AR is itself second order (v2), so the obstruction on the algebraic side is the first-order error of the
Hildebrand/HT counting step in Lemma PP. That is the same kind of estimate as pc3 §7 item 4. We have not shown it is
the only obstacle at 3/8. (SQ) is also borderline there, as in uprime.

### 3.6 Theorem AB₂: the algebraic blocks at second order, via cell covariances [sketch; added in v2.1]

*Source of the device.* c12rat's v1 note (`notes/c12rat_results.md` §8.1 Prop PD, read in the working tree, not yet
committed or reviewed when this was written) bounds the rational pair covariance with a cell/Cauchy–Schwarz identity.
In its §9.3 it suggests that the same device applies to the algebraic pair. We carry this out here. Under DW the
algebraic pair is simpler than the rational one: it is exactly a mixture of uniform measures on a lattice intersected
with parallelograms (Lemma TP-W), and there is no primitivity constraint. The point is that **no second-order asymptotic
is needed**:
- the covariance of two functions of *separate* coordinates is bounded by the product of their one-point standard
  deviations;
- each standard deviation is first order;
- so the covariance is second order.

**Setting.** One mixture component ν_y = Unif(L′ ∩ Π_y), L′ = {x₁ ≡ κx₂ (mod Δ′)}, of Cor TP-q. Cells are
c = (I₁ × I₂, δ), where I_i run over a common grid of intervals of relative length 1/K and δ | Δ′ is the gcd class
gcd(x₁, Δ′) = gcd(x₂, Δ′). This class is common to both coordinates, since κ is a unit. Define the **exact** one-point
cell densities
  d_i(I, δ) := Ψ(I; gcd(n, Δ′) = δ)/#{n ∈ I : gcd(n, Δ′) = δ}  (B-smooth n; for parents B₂).

**Step 1 (exact factorisation in a cell).** For a cell inside Π_y:
- #{x ∈ L′ ∩ c : x₁, x₂ ∈ 𝒮} = N(c)d₁(I₁, δ)d₂(I₂, δ) + Err(c), where N(c) := #(L′ ∩ c), |Err(c)| ≤ √(V(I₁)V(I₂)),
  V(I) := Σ_{a mod Δ′, gcd(a,Δ′)=δ}|E(I; a)|². This is the pair-sum identity: cross terms vanish exactly
  (`test_pair_sum_identity_PS`). The class-δ principal term is φ(Δ′/δ)m₁m₂ = N(c)d₁d₂(1 + O(Δ′/|I|)), with
  m_i := Ψ(I_i; δ)/φ(Δ′/δ).
- The singles are exact, with **no** character error: #{x ∈ L′ ∩ c : x₁ ∈ 𝒮} = N(c)d₁(I₁, δ)(1 + O(Δ′/|I₂|)). Each
  x₁ has |I₂|/Δ′ + O(1) partners.

**Step 2 (covariance identity).**
- Let ν̂ be the probability on cells with ν̂(c) = N(c)/#(L′ ∩ Π_y) for interior cells, plus an atom ∗ (boundary cells)
  where d_i := 0.
- Then A₁₂(y) − A₁(y)A₂(y) = Cov_ν̂(d₁, d₂) + (interior Err) + (boundary terms) + tiny.
- Under ν̂, δ is π₀-distributed independently of the rectangle, up to O(Δ′/|I|): the residue pair of a uniform point
  of L′ is uniform on {(κa, a)}.
- Cauchy–Schwarz gives |Cov_ν̂(d₁, d₂)| ≤ (Var d₁·Var d₂)^{1/2}.
- The boundary terms are handled exactly as in c12rat's Prop PD, by Cauchy–Schwarz against fourth moments of d_i over
  the one-point marginal. That gives ≤ C h^{1/4}A₁A₂, where h ≍ K₀/K is the boundary-cell mass. With the short-fibre
  truncation of Lemma PP: K := K₀⁴D^{8α}(log N)^C.
- The mixture over y adds Cov_y(A₁(y), A₂(y)) ≤ sd_y(A₁(y))sd_y(A₂(y)), again a product of one-point deviations.

**Step 3 (one-point deviations are first order) [proved from first-order inputs].** Under the one-point marginal
(and over y),
  Var(d_i)/(E d_i)² ≤ C(D^{−2α}λ^C + λ^C/u²)  (child; for parents, u₂ and B₂).
- d_i(I, δ) = ϱ(x_I)τ*(δ)(1 + e(I, δ)) with |e| ≤ λ^C/u. This uses Hildebrand's short-interval theorem and Möbius
  over e | Δ′ (Lemma PP), or for parents the truncation at B₂^{1/2}.
- log ϱ(x_I) varies by (1 − β)·log(scale ratio). Its variance under the marginal is ≤ C(1 − β)²log²d, by the
  exponential tail of Lemma AR.
- Var_{π₀}τ* ≤ LF − 1 ≤ D^{−2α}λ^C (Lemma LF).
- Var(e) ≤ λ^C/u².
- So Var(d_i)/(Ed_i)² ≤ C(D^{−2α}λ^C + λ^C/u²).

The first-order errors e enter only through variances. That is why no second-order asymptotic is needed.

**Step 4 (character errors).** Σ_c|Err(c)| over all cells of all non-thin pairs of a dyadic coordinate class is
controlled by Harper's Theorem 2 (BDH), or by Drappeau's Lemme 5 since Δ′ ≤ L[α]^{O(1)} = x^{o(1)}. This is summed over
the K² grid cells with the multiplicities of Lemma FB. The relative error is ≤ K^{O(1)}e^{−cD^α/λ²}, as before.

**Theorem AB₂.** In design 2L_DW(α) for every fixed α ∈ [3/8, 2/5), and with the same thin-pair, skew and TV
treatment as Theorem AB:
- (i) Σ_{s≤T}(w/W)|ϕ^alg| ≤ D^{−2α}λ^C(1 + Σ_{s≤T}(w/W)ϕ^rat);
- (ii) Σ(w/W)|σ^alg| ≤ (D^{−2α} + D^{−2(3−7α)})λ^C(1 + Σ(w/W)σ^rat);
- (iii) as in Theorem AB, with T and K₀ enlarged by D^{α}.

At α = 3/8 the design is 2L_DW(3/8). This is uprime's 2L′ with DW, η₀ ≍ D^{5/8−1}, and uprime's (SQ) constraint
ab₂ > 3d₀/(32κ), which is now needed, since (SQ) is borderline-automatic only for α > 3/8.

*Sketch status.*
- Steps 1, 2 (interior), 3 and 4 use only lemmas proved here or published inputs.
- Step 2's boundary bookkeeping is c12rat's Prop PD argument transposed. It is not re-derived line by line here, and
  c12rat's note is itself unreviewed at the time of writing.
- The margins at 3/8: max(3α − 1, 1 − 2α, 2α − 1/2) = 1/4 < 3/8 = α. (1 − α)/2 = 5/16 < 4α − 1 = 1/2. 2 − 5α = 1/8 and
  1 − 2α = 1/4 < 3 − 7α = 3/8. `params.margins2` encodes them.

**Consequence for the budget.** With A and B both second order:
- A ≤ D^{−2α}λ^C ≪ η ≍ D^{−(3−7α)}, which holds iff α > 1/3;
- √(AB) ≤ D^{−α−min(α, 3−7α)}λ^C ≪ η, which holds for all α ≥ 3/8.

So Prop PB's conclusion holds on the closed range [3/8, 2/5). §3.5's diagnosis is superseded: the first-order obstacle
at 3/8 is removed by the covariance device, not by second-order asymptotics.

---

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
   - (v2; codex #8, proof critic) **No angular discard.** v1 discarded points within angle D^{−A₀} of the cell's
     direction. That set can be almost the whole box for skewed labels (codex's example r = −1).
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
- (a) The child filter must be part of Good(P′), otherwise DW(a) fails at points with tiny |G(P′)| (codex pc3 #2). Its
  pass fraction is 1 − O(N^{−η₀/4}), so it introduces no O(1) environment term (unlike the parent filter, pc3 M1).
- (b) S₂ must be rescanned with DW's larger height (codex). Keeping uprime's S₂ would lose an L[α]-scale factor.

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
    the convolution changes its smoothness average at first order (1 − β)log d (codex #10).
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
     w_𝔮/z_𝔮 = Pr(R_𝔮)/#J (v2, proof critic): Σ_𝔮(w_𝔮/W)/z_𝔮 = Σ_𝔮Pr(R_𝔮)/(#J·W) ≤ dπ(B₂)/(#J·μ₁) =
     L[4α − 1 − (1 − α) + o(1)] → 0. This uses Σ_r Pr(R_{(q,r)}) = E[#roots of F mod q] ≤ d. No lower bound on z_𝔮 is
     needed.
   - Theorem AB(i), (iii) and RAT(R1) give A ≤ D^{−α}λ^{C}. RAT(R1) is restricted to s ≤ T, matching AB(i).
   - Theorem AB(ii) and RAT(R2) give B ≤ D^{7α−3}λ^{C}.
2. Prop PB: E N₂ ≥ μ₁/2.
3. U′-DW(α) (§2.7) gives E[N₁²] ≤ (log N)^C(μ₁² + μ₁) and μ₁ ≥ (log N)^{−C}, so γ₂ ≥ (log N)^{−C} by Paley–Zygmund
   (Prop UB).
4. Theorem K-U and Prop select-U with Lemma reject-DW give the split cost (γ₂/(16 log N))^{2k}/(32k²) = L[α, −O(κ)].
5. Rows (nested Lemma RW), POSTPROCESS-DW (§2.8), Theorem DD with Lemma T_multi, and Lemma F_α at α finish the proof, as
   in multidescent Theorem 2L. ∎

**Theorem 2L-DW₂ [proved, conditional; v2.1].** For every fixed α ∈ [3/8, 2/5), FACTOR_2L in design 2L_DW(α) (at 3/8
with uprime's (SQ) constraint ab₂ > 3d₀/(32κ)) factors every N in expected time L_N[α, c(α)]. It assumes:
- (i) the internal chain;
- (ii) **Hyp RAT₂(α)**: (R1), (R2) with D^{−2α}λ^C on both right-hand sides (c12rat's Theorem RAT, unreviewed v1);
- (iii) **Theorem AB₂** (§3.6, sketch).

*Proof.* As Theorem 2L-DW(α), with A, B ≤ D^{−2α}λ^C. The budget (§3.6, "Consequence") holds on [3/8, 2/5), and the
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
  2L-DW₂, also Theorem AB₂ (§3.6), whose boundary bookkeeping transposes c12rat's Prop PD (unreviewed v1).
- *Proved here:* Lemma η (v2 rewrite), Prop PB, Lemma VD, Lemmas GB-DW, CU, RE, CT, IR-DW, reject-DW, W-DW, DW(a)⁺, MM, BZ, TP,
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
  - Hildebrand, J. Number Theory 22 (1986) 289–307 (Ψ = xρ(u)(1 + O(log(u + 1)/log y)) and the short-interval form for
    x y^{−5/12} ≤ z ≤ x; Jain, QJM, Lemma 3.2(ii), attributes it to Hildebrand's Theorem 3, per codex). It is used at
    y = B (Lemma PP, child) and at y = B₂ (Lemma η, parent PP);
  - de Bruijn: the monotonicity of −ρ′/ρ and ρ(u − 1)/ρ(u) = uξ(u)(1 + O(1/u)), for Lemma η;
  - Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Theorems 1, 2, 3;
  - de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), Thm 2.4(i);
  - Canfield–Erdős–Pomerance;
  - Tenenbaum's Ψ_q theorem;
  - Iwaniec–Kowalski Thm 7.13;
  - Montgomery–Vaughan;
  - Lagarias–Odlyzko, Stark, Landau (LV6);
  - AKS, Wiedemann/Kaltofen–Saunders, Cantor–Zassenhaus, Bernstein;
  - Trudgian 2016 (CS′).

  Not needed: Drappeau's Lemme 5 and Lemma 8.5′/8.5H, which are replaced by DW(a). Theorem U itself still uses them, for
  the recursive value factorisations (FACTOR_U).

**Achieved α:**
- none unconditionally;
- conditionally on RAT(α) and Theorem AB: every α > 3/8;
- conditionally on RAT₂ and Theorem AB₂: α = 3/8.

---

## 6. Experiments (pre-registered, `notes/c12alg_prereg.md`, commit `20c882b`)

### 6.1 E1a: exact TV of Lemma TP (72 pairs, 6 configurations)

| (d, V, H) | median TV | max TV | bound | median TV(2H + 1)/(dV²) | max TV/bound |
|---|---|---|---|---|---|
| (3, 2, 40) | 0.0103 | 0.0133 | 0.630 | 0.069 | 0.021 |
| (3, 2, 80) | 0.0042 | 0.0067 | 0.317 | 0.056 | 0.021 |
| (3, 2, 120) | 0.0035 | 0.0039 | 0.212 | 0.069 | 0.018 |
| (3, 3, 80) | 0.0083 | 0.0157 | 0.689 | **0.0495** | 0.023 |
| (3, 3, 120) | 0.0077 | 0.0186 | 0.461 | 0.069 | 0.040 |
| (4, 2, 40) | 0.0102 | 0.0176 | 0.840 | 0.051 | 0.021 |

- **P1a held** (72/72, exact rationals).
- **P1b held** (‖x(z)‖_∞ ≤ 1 + 4V², 72/72).
- **P1c failed as registered.** The monotonicity part held, but the band [0.05, 4] was missed at (3, 3, 80) by 0.0005.
  The registered lower edge was too high: the true constant is ≈ 0.05–0.07 × dV², i.e. Lemma TP's d(1 + 4V²) is ≈ 50×
  pessimistic at toy scale. This affects no claim.

### 6.2 E1b: joint smoothness under DW (6 fields × 2·10⁶ samples, 602 pairs; d = 3, y = 512, u ≈ 3.2)

| statistic | value | 95% CI (field bootstrap) | band | verdict |
|---|---|---|---|---|
| P1d Σ J/Σ J_gcd | 1.0172 | [1.0160, 1.0184] | [0.95, 1.08] | held |
| P1e Σ J/Σ J_full | 1.0001 | [0.9995, 1.0006] | [0.97, 1.03], CI ⊂ [0.95, 1.05] | held |
| P1f slope of (J/M₁M₂ − 1) on (J_full/M₁M₂ − 1) | 0.9999 | [0.994, 1.005] | [0.8, 1.2] | held |
| P1g mean (J − J_full)/(M₁M₂) | 1.4·10⁻⁴ | [−5·10⁻⁴, 7·10⁻⁴] | ±0.03 | held |

Descriptive:
- the per-pair excess J/(M₁M₂) − 1 ranges over [−0.037, 0.203], pooled 1.075;
- the archimedean proxy ranges over [−0.029, 0.117];
- the gcd-only slope is 1.095.

Reading, fixed in advance: at toy scale the DW pair law shows no correlation beyond "lattice gcd classes ×
archimedean sizes", the structure of Lemmas LF, AR and PP. J_gcd uses empirical gcd-class marginals from the same
samples. So E1b does not independently validate LF's saddle-point model (codex). Both factors are large at u ≈ 3: (1 − β)² is not small at
toy scale. Asymptotically they are second order. The toy does not test precision.

### 6.3 E2: the thmA toy pipeline with DW fields (composite w, W = 64M)

- **P2a failed as registered.** E2 gave 24/24, but E2b gave 7/8: the instance semi 44 b seed 7024 hit the
  field-selection cap 20 000. Post hoc, not registered: with cap 4·10⁵ it was factored correctly after 21 568 fields.
  Cause: at toy height H = 4 the term wF̂_N dominates the values, so DW values are ≈ 64× larger than thmA's, smoothness
  is rarer and field selection is slower. This is a toy-scale effect: in the design, H ≥ W_mN^{η₀} makes the w-part
  negligible.
- **P2b held.** 0 row-identity, 0 congruence, 0 CU (5 911 atoms checked), 0 forced-residue, 0 IR-DW violations.
- **P2c held.** 120/120 tested kernel vectors were squares.
- **P2d held.** 5 392 content rejections among 54 915 fields, against 5 398.6 predicted (z = −0.09). The exact
  computation replaced a truncated Möbius sum after codex #13; the verdict is unchanged.
- **P2e held.** DD fruitfulness ω = 2: 26/56 = 0.46 (z = −0.53); ω = 3: 23/32 = 0.72 (z = −0.41).

---

## 7. Answers to the task

1. **DW(b) (item 3).**
   - The TV lemma is proved (Lemma TP, TP-W, Cor TP-q), uniformly in the angle.
   - The two-point smoothing for smooth numbers on the lattice cosets is reduced (Theorem AB, sketch) to Harper's BV
     (on average over pairs), Hildebrand's short intervals and the HT/dBT local estimates, at first order.
   - The relative error achieved is λ^C/u, enough in the window. At 3/8 the exact failing estimate is the second-order
     archimedean comparison (§3.5).
   - Lee–Venkatesan's Lemma 8.5 and Lemma 8.5H are not needed: DW makes the moduli the determinants, not G(P′). LV's
     Lemma 8.5 assumes X^ε ≤ Q (codex), so it would not apply anyway.
2. **DW compatibility (item 5).** Compatible, with explicit repairs (§2). The alternative designs examined do not keep
   DW-L (§2.9).
3. **Child lower bound (item 6).** Proved for every label (Prop CLB), with the corrections in §4.

## 8. Literature

Inputs as in §5's dependency list. Specific notes:
- **Harper 1208.5992.** Thm 1 is BV with max over residues, Q ≤ √Ψ, error Ψ(e^{−cu/log²(u+1)} + y^{−c}) +
  Q√Ψ(log x)^{7/2}. Thm 2 is BDH. Neither has an inner maximum over x (hence the grid). Statement forms as in pc3 §4.1
  and codex; the Thm 1 form is UNVERIFIED by us against the preprint.
- **Hildebrand 1986 (JNT 22).** The direct short-interval asymptotic is what removes the κ-amplification (§3.3). The
  theorem number is UNVERIFIED.
- **Rigorous-NFS precedent for an integer multiplier.** None found by pc3 or by us; this is a search result. The
  Knuth–Schroeppel multiplier kN is the heuristic analogue.
- **Pascadi** (Compositio 2025): smooth numbers in APs to large moduli, for a fixed residue (codex). It does not replace
  Harper's maximum over residues, which Lemma FB needs.
- **Lee–Venkatesan.** Their randomisation f̂ + G·R_c (JNT 187 (2018)) is the source of the c-part. Their stochastic
  deepening (Lemma 5.6) avoids second moments in a different way.

## 9. Reproduction

```
.venv/bin/python -m pytest -q tests/test_c12alg.py
.venv/bin/python -m src.c12alg.e1 E1a experiments/raw/c12alg_e1a.jsonl
.venv/bin/python -m src.c12alg.e1 E1b experiments/raw/c12alg_e1b.jsonl
.venv/bin/python -m src.c12alg.harness_e2 experiments/configs/c12alg_e2.json experiments/raw/c12alg_e2.jsonl
.venv/bin/python -m src.c12alg.analyze E1a experiments/raw/c12alg_e1a.jsonl   # likewise E1b, E2
.venv/bin/python -c "from src.c12alg import params as P; from fractions import Fraction as F; print(P.margins(F(7,18)))"
```

## 10. Reviews

**Codex blind `c12alg_1`** (xhigh, launched before any analysis here). Verdict OPEN. Points it found independently:
- the TV lemma with the interpolation/determinant cancellation;
- the window budget A ≪ D^{−α}, B ≪ D^{−s};
- η₀ ≍ D^{4α−2} and σ₂ = (α/a)(4/d₀ + d₀b);
- unrestricted w as the clean design;
- the IR averaging;
- the warning that first-order slack does not by itself remove Saias (𝒦²/u amplification).

It listed weight domination and LD* as OPEN. Probabilities: TV lemma 95%, DW compatible 85%, α < 2/5 30%.

**Codex adversarial `c12alg_2`** (xhigh, on v1): **MAJOR REVISION.**

| # | sev. | finding | action |
|---|---|---|---|
| 1 | FATAL | Lemma η loses (log N)^C, not λ^C | accepted; Lemma η rewritten (pointwise ratio, Hildebrand at B₂, W/μ₁ ≤ 3u₂) |
| 2 | MAJOR | thin-pair parameters leave d(log N)^C | accepted; T, K₀, K re-chosen |
| 3 | MAJOR | non-reduced classes: Σδ^{−β} diverges; weighted grid inequality | accepted; δ ≤ Q sum = polylog; per-subclass statement in §3.4 |
| 4 | MAJOR | fibres shrink at vertices | accepted; short-fibre truncation with DS enlargement (Lemma PP) |
| 5 | MAJOR | parent divisors up to L[1−α]: v²/u₂ error | accepted; truncation at B₂^{1/2} with divisor-moment tail |
| 6 | MAJOR | PU's uniform fibre length false | accepted; fibre enlargement and gcd reduction (Lemma PU) |
| 7 | MAJOR | same-ℓ parent pairs (weight 1/ℓ) | accepted; separated, contribution ≤ z²/π_J |
| 8 | MAJOR | CLB angular discard false for skewed labels | accepted; no discard, per-point interior intervals |
| 9 | MAJOR | cells, WD per cell, parent mean lower bound | accepted; WD per cell; exact identity w/z = Pr(R)/#J replaces the z lower bound |
| 10 | MAJOR | RAT weight is not a single ϱ value | accepted; idealised weight defined |
| 11 | MINOR | W-DW needs (u,t) = 1; content cutoff; q² \| Δ possible (S₂² > B for a = 3/5); LF p > B; AR atom | all accepted; constraint (A1) stated |
| 12 | MINOR | 3/8 diagnosis overstated | accepted; §3.5 softened |
| 13 | MINOR | 10/12 not 9/11; p_content truncated (5 398.6, z = −0.09) | accepted; analysis fixed |

Codex withdrew its blind OPEN on LD* (the all-residue rebuttal is correct). Probabilities: AB as sketched 40%, DW 93%,
conditional theorem 70%.

**Opus proof critic** (read-only, on v1).
- §2 lemmas, MM, BZ, TP, TP-W, Cor TP-q and Lemma VD: ACCEPT.
- Lemma η: MAJOR (same as codex #1; its suggested fix is the one adopted).
- W-DW, DW(a)⁺ edge case ab = 0, §2.7 range restatement (uprime's D^{1/4} hypotheses fail for α > 11/28), PU
  (g-dependent scales), WD (per cell), CLB (skew remark): ACCEPT WITH FIXES, all applied.
- Ξ was missing the dominant filter-threshold spread e^{O(D^{2α−1/2}λ)}: added to Ξ and to `params.py`.
- Lemma AR is second order: restated.
- Theorem AB(i)'s right side must be restricted to s ≤ T: done.
- No hidden first-order term and no circularity.
- Probabilities: AB writable 65%; conditional theorem 15% as v1, 65–70% after repairs.

**Opus literature critic:** see v2.1 below.

## 11. Candidate card

- **Mechanism.** Two-level descent at α ∈ (3/8, 2/5). The integer field multiplier (DW) makes the algebraic pair law an
  explicit lattice-parallelogram mixture (Lemma TP). Prop DC's two-point statements are then needed only to first order.
- **Result.**
  - Proved: DW compatibility; Lemma TP; the child lower bound; the window budget; weight domination; skewed labels for
    C₁.
  - Sketch: the algebraic blocks of C₁ and C₂ in the window (Theorem AB).
  - Conditional: L[α] for every α > 3/8, given Hyp RAT(α).
- **Failure scope.** No rigorous α < 2/5 until Theorem AB/AB₂ are written out and reviewed and the rational blocks
  (c12rat) are reviewed. Below 3/8 the two-level design is infeasible (multidescent §3.1). Three levels are borderline
  (pc3 §7.5).
- **Next step.** Write Lemmas FB/PP and §3.6's boundary terms to the last constant; an independent review of AB₂
  together with c12rat's Prop PD.
