# c12rat: the rational two-point blocks of C₁/C₂ (pc3 §7 items 1, 2, 4) — results

Agent `c12rat`, phase 27, 2026-10-24. **v2.1** (literature critic applied, §12.1). **v2**: codex adversarial round (MAJOR REVISION) and Opus proof critic (W: ACCEPT
WITH FIXES) applied, §12. The v2 repairs have not been re-reviewed.

Files:
- this note; pre-registration `notes/c12rat_prereg.md` (commit `b9c2b8d`, code frozen there);
- code: `src/c12rat/geom.py` (Lemmas LS, MULT, TF), `bts.py` (Lemma BTS), `e3.py` (E3: the cell decomposition),
  `params.py` (leading-order exponent margins of §1.3), `experiments.py`, `analyze.py`;
- tests: `tests/test_c12rat.py`;
- data: `experiments/raw/c12rat_e{1,2,3}.jsonl` and `_summary.json`;
- codex: `codex/c12rat_1.{md,out}` (blind round, launched before any analysis was written).

Labels as in pc3/c12alg: **proved** (from the project's internal, unrefereed chain plus the cited published results);
**sketch**; **heuristic**; **open**; **experimentally supported** (toy); **UNVERIFIED** (literature statement not
checked by us). Notation: λ := log log N, D := log N/λ, L[α, c] = e^{(c+o(1))D^αλ}, "polylog" = e^{O(λ)}. The
Hildebrand–Tenenbaum saddle point is β = β(x, y) and γ := 1 − β. The determinant of a pair is written **K**
(pc3 writes D, which here is log N/λ), and the unit of pc3's Lemma LAT is **ι**.

---

## 0. Verdict

**Plain language.** The rational two-point blocks of C₁ and C₂ are **proved** (internal, unrefereed; v2 after one codex
adversarial round and an Opus proof critic, whose repairs are applied but not re-reviewed). These are the parts of the
two pair-variance hypotheses that concern only the smoothness of the linear forms G_g(P′) = g₁a′ − g₀b′. The bound is
**second order**: the incidence-weighted relative pair covariance of child points is ≤ D^{−2α}(log D)^C, and for parents
≤ (D^{−2α} + D^{−2(3−7α)})(log D)^C. The window 3/8 < α < 2/5 needs only first order, and the bound also suffices at
α = 3/8.

How it works:
- **Cells.** For a pair of child points, cut the joint range of (G′, G″) into cells: a small rectangle of values times a
  local "type" (valuations at small primes and at the primes dividing the pair's determinant).
- **Exact factorisation.** Inside a cell, the pair count factorises exactly into the product of two one-point cell
  densities, plus an error made of non-principal characters.
- **Cauchy–Schwarz.** The covariance of the two cell densities is at most the product of their one-point standard
  deviations.
- **First-order inputs.** One-point densities are known to first order (Hildebrand's short-interval theorem,
  Hildebrand–Tenenbaum). So first-order one-point asymptotics give a second-order two-point bound, with no Fourier tails,
  no pair main term ever computed, and no transversality assumption.
- **Character errors.** These are controlled on average over the whole incidence-weighted family by Drappeau's published
  Lemme 5 alone (which Drappeau states without proof, as a consequence of §3.3 of Harper's unrefereed preprint
  arXiv:1208.5992), through a weighted application theorem whose only inputs are a pointwise upper bound for the incidence
  weights and lattice multiplicities.
- **Skew and thin pairs.** Skewed labels are negligible by weight. Nearly parallel pairs are bounded by an elementary
  smooth-numbers-in-progressions estimate, whose loss e^{O(D^θλ)} (θ < α) is beaten by the BDH-type saving e^{−cD^α/λ²}.

**This does not give an unconditional algorithm below L[2/5].** C₁ and C₂ also have algebraic blocks, which are c12alg's.
With c12alg's theorems (internal chain + the rational blocks + their algebraic Theorem AB/AB₂, which are **sketches**):
- L_N[α] for every fixed α ∈ (3/8, 2/5): c12alg's Theorem 2L-DW(α), conditional on AB (first order).
- L_N[3/8]: c12alg's Theorem 2L-DW₂, conditional on AB₂ (second order, with four listed open bookkeeping obligations).
  AB₂ uses the cell device of this note on the algebraic side.

So the honest status is: **conditional on c12alg's algebraic sketches and on the internal chain**, and on this note's
v2 repairs surviving re-review.

Precisely:
1. **Theorem RAT [proved].** Assume c12alg's window design 2L_DW(α), for any fixed α ∈ [3/8, 2/5). (pc3's 2L′ at 3/8 is
   covered only under a proviso, §1.1.) Work cell by cell, with environment cells of relative size D^{−A₀}, A₀ ≥ 1, and
   with the smooth algebraic weights of DW(a)⁺. Then:
   - (RAT-1) Σ_{𝔮,e: s(𝔮)≤T₀}w_{𝔮,e}Φ_{𝔮,e} ≤ C D^{−2α}λ^C W + e^{−cD^α/λ²}μ₁;
   - (RAT-2) Σ_{𝔮,e}w_{𝔮,e}Φ^par_{𝔮,e} ≤ C(D^{−2α} + D^{−2(3−7α)})λ^C W + e^{−cD^α/λ²}μ₁;
   - (RAT-0) the labels of skewness s(𝔮) > T₀ = e^{C_TD^θλ²} carry weight ≤ T₀^{−1/2}μ₁.

   Here Φ is the off-diagonal part of the g-variance in c12alg's Lemma VD (the law of total variance; codex's (15)), the
   one-point terms are c12alg's, and θ = 2α − 1/2 < α.
2. **Corollaries [proved].**
   - Hyp RAT(α) of c12alg holds for every α ∈ (3/8, 2/5).
   - Hyp RAT₂(α), in the form c12alg's Theorem 2L-DW₂ uses (A ≪ D^{−2α}λ^C, B ≪ (D^{−2α} + D^{−2(3−7α)})λ^C), holds on
     [3/8, 2/5).
3. **The tools (task items 1, 2, 4).**
   - *Weighted application theorem* (Theorem WA, after codex's transfer lemma), with Lemmas FAM, DOM and SPR. The only
     published input is Drappeau's Lemme 5. The large sieve over child primes is not needed: every modulus is a small
     conductor.
   - *Skew* (Prop SK, Lemmas LS, TF, PUB). Skewed labels have weight ≤ 3/T₀ of their window. Inside non-skewed labels,
     thin pairs are a fraction ≤ 60T₀/κ₀ and are bounded with loss Ξ.
   - *Archimedean second order* (Prop PD, Lemmas FO, OPV, ADM). Saias is not needed. The rate O(log(u + 1)/log y) in
     Hildebrand's short-interval theorem is essential (with o(1) only, the bound would not be second order).
4. **Achieved α.** No unconditional α < 2/5. Conditionally, as above: every α ∈ (3/8, 2/5), and 3/8 itself. The
   rational side imposes no constraint on [3/8, 2/5) (§9.2).
5. **Experiments (pre-registered; toy): 12 of 13 predictions held** (§10). P2b failed: it was a mis-specified sanity
   statistic, explained post hoc. E3 shows the cell decomposition capturing toy pair covariances, with a residual that
   cancels over the box.

**Probabilities (ours, after review).**
- That Theorem RAT in design (W) is correct after the v2 repairs: 70%. (Proof critic: 80%; codex: 60%.)
- That the conditional L[α] (α > 3/8) is correct given c12alg's Theorem AB: 50%. (Codex: 55%.)

---

## 1. Setting

### 1.1 Designs

Two designs are covered; both use the DW field law of pc3 §6.2 as repaired by c12alg §2 (w uniform on all integers of
(W_m, 2W_m], c uniform on C_H, F = wF̂_N + G·R_c, g uniform on primitive points of [M/2, M]²).
- **(W)** c12alg's window design 2L_DW(α), α ∈ (3/8, 2/5) (c12alg §1.1): d ≍ D^{1−2α}, k ≍ D^α, B = e^{aD^αλ},
  B₂ = e^{b₂D^{4α−1}λ}, η₀ log N ≍ D^{3α−1/2}λ, special-q ℓ ∈ J = (Y, 2Y], Y = L[1 − α], level-1 box radius V₁ with
  V₁² ≍ YS₁², child boxes of radius S₂ (S₂² = L[α, σ₂]) in Gauss-reduced coordinates, child primes q ∈ (B, B₂].
- **(P)** pc3/uprime's 2L′ at α = 3/8 with DW: B₂ = L[1/2, b₂], η₀ = D^{−1/2}, otherwise as uprime §1. **Proviso**
  (proof critic M2): here N^{η₀/2} = L[1/2, 1/2] and q ≤ L[1/2, b₂], so the TV terms qN^{−η₀/2} of DW(a)⁺ are small only
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

**Lemma VD (c12alg §1.3; codex c12rat_1 (15)).** Conditional on R_𝔮 ∩ e,
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
θ = 2α − 1/2 (proof critic m13). **On [3/8, 2/5), θ < α**,
so Ξ^{O(1)} = ε^{−o(1)}, and also Ξ^{O(1)} ≤ y^{1/12}. Further, with P₀ := ⌈D^{7/10}⌉:
- every modulus is a small conductor: log K ≤ (η_D/4) log x;
- ω(K) = o(P₀) for every child and parent modulus, so Σ_{p | K, p > P₀} p^{−β} → 0;
- P₀# = e^{(1+o(1))P₀} ≤ x^{o(1)}, far below every cell length.
`params.margins` checks these (and the window/3/8 comparisons of §1.2) for α = 3/8 (both designs) and on a grid of
(3/8, 2/5); at α = 2/5 only "B₂ ≪ V₁" becomes borderline.

### 1.4 The algebraic weights (regularity used)

We use only the weights a_P of DW(a)⁺ (children, with the root coset) and the parent weights a_P of §1.2 — not an
abstract class (codex c12rat_2 #7: a pointwise regularity alone says nothing about relative amplitudes across points;
the cross-point spread is proved for these weights in Lemma SPR, §5.7). Their regularity:

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
expansion over a ratio ≤ y² is O(log²(y²)/(log x log y)) = O(D^{−α}), inside E_a (codex c12rat_2 #1). Global: Fact dBT1.
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
cancellation (Lemma SI): no cross terms, and the error is a product of two one-coordinate L² discrepancies. Codex's
blind round found the same coupling principle in L² form (c12rat_1 (13)).

**Corollary CPW (slowly varying weights).** If moreover a(m) = ā·u(m)·(1 + ϵ(m)), b(m) = b̄·v(m)(1 + ϵ′(m)) with
u, v ≥ 0, |ϵ|, |ϵ′| ≤ ε′ ≤ 1/4, then with Ẽ_u := Σ_{χ≠χ₀}|û(χ)|²/(Σu)²,
  |Σ a b 1_C − (|C|/φ²)AB| ≤ (|C|/φ²)AB·[3(Ẽ_u + Ẽ_v) + 4ε′²].
*Proof.* Class sums: a(r) = ā u(r) + e(r) with |e(r)| ≤ ε′āu(r) (positivity of u). So ‖α‖ ≤ ā(1 + ε′)‖u − Σu/φ‖
+ ε′āΣu/√φ, and similarly for β; A ≥ (1 − ε′)āΣu. Insert in Lemma CP and use xy ≤ (x² + y²)/2. ∎
*The point:* the weights' oscillation enters **squared**, and only multiplied by the main term. (A naive treatment of
the weights inside the character sums would cost 1/ρ; codex's warning "slowly varying does not supply an inner
supremum for free" is met by this positivity argument.)

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
- **Fact H1, H2** (Hildebrand, J. Number Theory 22 (1986) 289–307, Theorem 1 (H1) and Theorem 3 (H2); restated with
  these ranges and this error term in Jain, arXiv:2502.10530, Lemma 3.2(i)–(ii), and in the Hildebrand–Tenenbaum survey,
  JTNB 5 (1993) §5 (Thm 5.1); the original paper was not opened by us or our critics). For y ≥ exp((log log x)^{5/3+ε}):
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
is why the weighted application (§6) needs the moduli to be **spread** (multiplicity domination), as codex warned.

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
Rankin (proof critic m7, codex #3: the small and mid primes must be treated separately). Write e = e_s·e_m with
e_s | P₀# and e_m | 𝔪₁. Since s ≤ S₀ = y^{o(1)}, se > y² forces e_s > y^{0.9} or e_m > y. With σ := 1/log P₀ (p^σ ≤ e for
p ≤ P₀): Σ_{e_s>y^{0.9}}e_s^{−β} ≤ y^{−0.9σ}∏_{p≤P₀}(1 + e p^{−β}) ≤ e^{−0.9 log y/log P₀}(log P₀)^{C}, and
Σ_{e_m}e_m^{−β} ≤ exp(Σ_{p|𝔪₁}p^{−β}) ≤ e; with exponent 1/2 instead: Σ_{e_m>y}e_m^{−β} ≤ y^{−1/2}∏_{p|𝔪₁}(1 + p^{−1/2+γ})
≤ y^{−1/2}e^{ω(K)P₀^{−0.49}} ≤ y^{−1/2}e^{D^{0.3}}, and Σ_{e_s}e_s^{−β} ≤ (log P₀)^{C}. Both tails are smaller than every power of D (log y/log P₀ ≍ D^α and aD^αλ/2 ≫ D^{0.3}). So the count is
(z/x)s^{−β}Ψ(x)[g_𝔪(β) + O((E_H + 1/u)∏_{p|𝔪}(1 + p^{−β}))], and ∏(1 + p^{−β})/∏(1 − p^{−β}) ≤ exp(3Σ_{p|𝔪}p^{−β})
≤ C(log P₀)³. Finally (z/x)s^{−β}Ψ(x) = (z/s)s^{γ}ρ̂(x). ∎

### 5.4 Lemma BTS (smooth numbers in a progression in a short interval) [proved; elementary]

Let k ≥ 1, a ∈ Z, I = (X₀, X₀ + Z] with W := Z/k ≥ 2y, and suppose gcd(a, k)·W ≤ X₀ (proof critic m8; in every
application gcd(a, k) ≤ K ≪ X₀/W). Then
  #{n ∈ I : n ≡ a (mod k), P⁺(n) ≤ y} ≤ 4W log(2y)·max_{W/y≤w≤2W}ρ̂(w),
and if log W ≥ (1/2)log(X₀ + Z), this is ≤ e^{O(λ)}·W·((X₀ + Z)/W)^{γ}ρ̂(X₀ + Z) (Facts dBT1, HT(α2)).
*Proof.* Put δ := gcd(a, k); n = δn′ with n′ ≡ a/δ (mod k/δ) a unit, n′ ∈ I/δ (if δ is not y-smooth there is nothing
to count). Each y-smooth n′ > W has a divisor m₁ ∈ (W/y, W] (take the largest product of its smallest prime factors
that is ≤ W). For fixed m₁ (smooth, coprime to k/δ), m₂ := n′/m₁ lies in an interval of length Z/(δm₁) in a fixed
class mod k/δ: ≤ W/m₁ + 1 values. Sum over m₁: W Σ_{W/y<m₁≤W}1_S(m₁)/m₁ + Ψ(W), and dyadic blocks give the bound.
(The exact inequality, before the dyadic step, is `bts.rhs`; E2 checks it on 2000 instances.) ∎
This is codex's (11) (found blind) and the large-modulus case of [pp] Lemma DS; c12alg's Lemma PU uses the same device.

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
mid-prime term is ≤ D^{−2α}λ²·o(1). (We centre at the explicit c, not at S′; codex c12rat_2 #4.)
*Proof.*
- *Ordinary blocks* (relative width ≤ h^{1/2}; Lemma CELL(iii) puts all others — width ≤ 1/8, at |n₁| ≤ 128z₁h^{−1/2} —
  into a set of mass ≤ Ch^{1/2}) and types with s_t ≤ S₀: Lemma FO with 𝔪 = P₀#·rad(K/δ)_mid, Lemma FIB for #U, and
  Lemma ADM (local) across the block give ψ′(I, t) = c·X(ξ)·Y(t)·Z(t)·(1 + e), |e| ≤ C(E_FO + E_a + γ_ah^{1/2}), where
  X(ξ) := a′(ξ)ρ̂(ξ)/(a′(x₁)ρ̂(x₁)) (ξ the block's end nearer 0), Y(t) := s_t^{γ}, and
  Z(t) := ∏_{p mid, j_p = k_p}(1 − 1/p)/(1 − p^{−β}): for a mid prime with j_p = k_p the cofactor is no longer coprime
  to p (p ∤ M_t), so the factor (1 − p^{−β})/(1 − p^{−1}) present in c is cancelled. Thus Z ≥ 1 and
  Z − 1 ≤ 2Σ_{p mid, j_p=k_p}(p^{−β} − p^{−1}) ≤ 2Σ_{p mid, p^{k_p} | n₁}(p^γ − 1)/p (codex c12rat_2 #2: v1 had the
  direction reversed). (Replacing β(ξ, y) by β(x₁, y) inside s^γ costs O(log s·log(x₁/ξ)/(log x·log y)), negligible.)
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

(Codex c12rat_2 #6–7, proof critic m10–m11: Lemma FAM needs **pointwise** comparisons across box points.)
- (children) For a label with s(𝔮) ≤ T₀, a cell e and every P′ ∈ Box(𝔮), let c_{P′} be OPV's c. Then
  max_{P′}c_{P′}/min_{P′}c_{P′} ≤ Ξ, and ψ′ ≤ Ξc_{P′} on every cell of every pair; hence S′S″/m_𝔮² ≤ 4Ξ³/n² and
  ψ′ψ″/m_𝔮² ≤ Ξ⁴/n².
- (parents) For every label and cell, call P ∈ 𝒫 ∩ Λ_𝔮 *active* if |G_g(P)| ≥ MV₁/16 on all of Ω_e. Then at least half
  of the points are active, every active P has S_P ≥ S_typ/polylog, and every P has ψ_P ≤ S_typ·polylog on every cell,
  where S_typ := ϱ̂_J(MV₁)·ρ₂(X_F)-scale. Hence m^par_𝔮 ≥ n·S_typ/polylog and ψψ/(m^par)² ≤ polylog/n²; **inactive P can
  have S_P = 0** (codex's example a = q⌊b/q⌋, b′ ≈ a′), so no lower bound for individual S_P is claimed or used.
*Proof.* (children) x₁(P′) ∈ [cδ_eM|P′|, 2√2M|P′|] (the cell is a square of side δ_eM, so max|G′| over it is
≥ δ_eM|P′|/4), and |P′| varies over the box by a factor ≤ 4S₂s(𝔮) ≤ 4S₂T₀ (Lemma TF's proof). So ρ̂(x₁) varies by
≤ C(16S₂T₀/δ_e)^{γ} (Fact dBT1) = e^{o(D^θλ)}, and a′(x₁) — whose argument is x₁L_{P′} with L_{P′} ∝ |P′|^{d−1} — by
≤ C((16S₂T₀/δ_e)^{d})^{γ_a} = e^{O(dγ_a(log S₂ + log T₀))} = e^{O(D^{1−2α}λ + D^{1−3α+θ}λ²)} ≤ Ξ (proof critic m10: this
is the term D^{1−2α}λ of θ). The local factors in c differ by ≤ e^{O(ω(K)/P₀)} = 1 + o(1). The pointwise bound
ψ′ ≤ C·X·Y·Z·c with X ≤ C(x₁/X⁻)^{γ+γ_a} ≤ Ξ^{1/2}, Y ≤ S₀^γ = 1 + o(1), Z ≤ 2. m ≥ n·min S′ ≥ n min c/2.
(parents) For g ∈ Ω_e, |G_g(P)| = |P||g||sin(θ_P − θ_g)|, |P| ≥ 2V₁/3, |g| ≥ M/√2, so P is active when
|sin(θ_P − θ_e)| ≥ 1/4 + O(δ_e); the angles θ_P of 𝒫 ∩ Λ_𝔮 are equidistributed on [π/4, 3π/4] (for each prime b in the
window, a ≡ rb (mod q) runs over an AP of difference q ≪ V₁), and θ_e ∈ [0.46, 1.11], so at least half the points are
active. For active P, the values lie in [MV₁/16, 2MV₁] and F(P) at scale [X_F/(64d), X_F] (c12alg Lemma η), so
ρ̂_J and ρ₂ vary by (64d)^{γ+γ₂} = 1 + o(1) (proof critic m11: the variation of |P|^{d−1} is ≈ 2^d, absorbed into
this range, not "≤ e"); S_P ≥ (pass fraction 1)·S_typ/polylog by Lemma FO. Upper bounds hold for all P by BTS/dBT1. ∎

---

## 6. The weighted application theorem (pc3 §7 item 1)

pc3 §5 (codex pc3_2 #1) asked for "one theorem stating the exact measure, its multiplicity and domination bounds".
Given Lemma CP, the needed statement is linear in the error energies, and it is elementary.

**Theorem WA (weighted application; after codex c12rat_1 §1, transfer lemma (2)–(3)).** Let 𝔏 be a finite *library* of
finitely supported sequences u_j: Z → [0, 1], independent of the moduli, with normalisers U_j > 0 such that
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ≠χ₀ (k)}|û_j(χ)|² ≤ ε U_j²   for every j ∈ 𝔏.
Let (c_t, j(t), k(t))_{t∈𝒯} be a finite family of terms, c_t ≥ 0, j(t) ∈ 𝔏, k(t) ≤ Q, and H_j(k) := Σ_{t: j(t)=j,
k(t)=k}c_t. Then
  Σ_t c_t·Σ_{χ≠χ₀ (k(t))}|û_{j(t)}(χ)|²/U_{j(t)}²  ≤  ε·𝔇,   𝔇 := Σ_{j∈𝔏}max_{k≤Q}φ(k)H_j(k).
*Proof.* Σ_jΣ_k[φ(k)H_j(k)]·φ(k)^{−1}Σ_χ|û_j(χ)|²/U_j² ≤ Σ_j max_k[φH_j]·ε. ∎

The content is in bounding 𝔇 for the true family (Lemma FAM), which needs only an **upper** bound for the incidence
weights (Lemma DOM) and the lattice multiplicities (Lemma MULT); errors are measured against the dominating mass
Σŵ ≤ (log N)^Cμ₁, not against W. Codex's concern that "small label means remain dangerous" is met because every
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
(k ≤ 4Q₂S₂², q > Q₂ ≥ B) and ≤ log(2V₁²)/log B ≍ D^{1−2α} for parents (polylog; proof critic m10). The pair's K is
one of the ≤ τ(s)2^{ω(s)} numbers kδ/r_δ with δ | s, r_δ | rad(δ) (M_τ = (K/δ)r_δ and δ = gcd(s, K)); each gives
Δ = K/q, with ≤ 10S² + 2S ordered pairs per label (Lemma MULT). For one pair and one dilation, the cells with endpoint ξ
in coordinate i and s_i = s have ν-mass ≤ 2Pr[s | n_i] ≤ 4/s (Lemma OPV's local law). Summing over labels (all r, e)
and dilations (Σϖ_{e′} ≤ 2): H_j(k) ≤ O(1)·max_qŴ_q·Ξ⁴n^{−2}(10S² + 2S)·τ(s)2^{ω(s)}(8/s)·h_*^{−2}(log P₀)⁴, with
Ŵ_q := Σ_{r,e}ŵ_{q,e}.
*The max over k.* φ(k) ≤ k ≤ 4Q₂S₂² (M_τ ≤ K). So max_kφ(k)H_j(k) ≤ CQ₂·max_qŴ_q·Ξ⁴h_*^{−2}(log P₀)⁴τ(s)2^{ω(s)}/s.
*The sum over j.* The grid points in the relevant range number ≤ log(2x_max/X⁻_min)/h_* ≤ CD^{2α}λ/h_* (proof critic
m10), and Σ_{s≤S₀}τ(s)2^{ω(s)}/s ≤ (log S₀)⁴ = D^{2(θ+α)}. So 𝔇 ≤ CD^{2α+2θ+2α}λQ₂max_qŴ_qΞ⁴h_*^{−3}(log P₀)⁴.
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
e^{O(D^θλ)} points, each a separate library element (codex: "a grid of size exp(o(u/log²u)) is affordable"); (d) LSQ
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
decomposition, whose only κ-dependence is the factor h_*^{−3} ∝ κ₀³ in Lemma FAM. (Codex's example in uprime §7.3 item 3 —
a label with 1 transverse point in 719 — is a skewed label in this sense.) c12alg's Theorem AB(iii) gives the same
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
  covariance moves by O(S′S″/q) (proof critic m9). ∎

**Lemma EXC (exceptional types) [proved].** Let ℰ be the n ∈ L ∩ Π with (a) gcd(n₁n₂, K_large) > 1, or (b) (parents) ℓ | n₁n₂
for one of the ≤ 2 primes ℓ ∈ J dividing K, or (c) s₁ > S₀ or s₂ > S₀, or (d) n_i in the single grid interval containing
a filter threshold (X⁻_{P′} for children; MV₁/32 and 2MV₁ for parents; proof critic m2: these intervals are not
library intervals). Then Z^{−1}Σ_{n∈ℰ}f′(n₁)f″(n₂) ≤ (e^{−λ} + κ₀^{−1/2})c′c″ and Z^{−1}Σ_{n∈ℰ}f′(n₁) ≤ e^{−λ}c′, and
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
together. No pair main term is ever computed or compared with the product of marginals. This is codex's "(13)" and
"common normalized model", with the model chosen per pair: since S′ = E_νψ′ holds exactly for *each* pair's own cell
structure, different partners may use different cells (codex's caution "a different one-point approximation for every
partner does not suffice" concerned approximations of the marginals). For interior cells the marginal identity is exact;
boundary and exceptional cells change it by the terms B′, P_∂, P_ℰ above, which are bounded separately (codex c12rat_2 #4).

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
primes contribute ≤ e^{−1/(4γ)}; and Σ_{p>P₀}log²p/p² ≍ log P₀/P₀, Σ_{p>P₀}p^{−2} ≍ 1/(P₀ log P₀) (proof critic M3, codex
#8: v1 replaced every log p by log y, which loses). So the sum is ≤ C(γ²log²P₀ + E_FO²)/(P₀ log P₀)·S′S″ ≤
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
and the **hard filter** |G(P)| ∈ [MV₁/32, 2MV₁], which needs care (codex c12rat_2 #6, proof critic m11):
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

## 9. Composition and the achieved exponent

### 9.1 Theorem 2L-RAT(α) [proved conditional on c12alg's Theorem AB (sketch) and the internal chain]

**Theorem 2L-RAT(α).** Let α ∈ (3/8, 2/5) be fixed and take c12alg's design 2L_DW(α) and c12alg's threshold T = Ξ_alg³(log N)^CD^α
(Ξ_alg = e^{O(D^θλ²)}); choose C_T so that T₀ ≥ T, so that c12alg's (R1), a sum over s(𝔮) ≤ T, is a sub-sum of (RAT-1).
Assume:
- (i) the project's internal chain, as in Theorem U and multidescent Theorem 2L;
- (ii) c12alg's Theorem AB (§3.3 there; **sketch**, "no gap known") and the lemmas c12alg labels proved.
Then FACTOR_2L in design 2L_DW(α) is Las Vegas and factors every N in expected time L_N[α, c(α)].
*Proof.* c12alg's Theorem 2L-DW(α) assumes (i), Hyp RAT(α) and Theorem AB. Corollary RAT(α) proves Hyp RAT(α) (with
room: D^{−2α}λ^C for A, (D^{−2α} + D^{14α−6})λ^C for B). ∎

**Theorem 2L-RAT₂ (α = 3/8) [proved conditional on c12alg's Theorem AB₂ (sketch, obligations O1–O4) and the internal
chain].** FACTOR_2L in design 2L_DW(3/8) (with c12alg's constraint (O4)) factors every N in expected time L_N[3/8, c].
*Proof.* c12alg's Theorem 2L-DW₂ assumes the internal chain, Hyp RAT₂ and Theorem AB₂; Theorem RAT supplies Hyp RAT₂ in
the form used there (A ≪ D^{−2α}λ^C, B ≪ (D^{−2α} + D^{−2(3−7α)})λ^C, which at 3/8 are both D^{−3/4}λ^C). ∎
(c12alg's §0 cites "c12rat's Theorem RAT" for RAT₂ with D^{−2α} on both sides on all of [3/8, 2/5); our B-bound is
D^{−2α} only at 3/8, but c12alg's proof of 2L-DW₂ already uses the B-bound with D^{−2(3−7α)}, so nothing breaks.)

**Complete dependency list.**
- *This note (proved here, internal):* Lemmas LC, LS, MULT, TF, CELL, TYPE, CP, CPW, FIB, EN, FO, BTS, PUB, OPV, ADM,
  DOM, FAM, PRIM; Theorem WA, Cor WA′; Props PD, SK; Theorem RAT; Corollary RAT(α); `params.margins`.
- *c12alg (internal):* Lemma VD (law of total variance), Lemma RE (root event uniform in g), Lemma DW(a)⁺ and Cor TP-q
  (one-point law with root coset; the form of the weights a_P), Lemma η, Prop CLB (one-point lower bound m_𝔮), Prop PB
  (precision budget), U′-DW(α), Lemmas GB-DW, CU, CT, IR-DW, reject-DW, W-DW, POSTPROCESS-DW, and **Theorem AB
  (sketch)**.
- *pc3 (internal):* Lemma 0, Prop DC, Prop DC_e, Corollary DC, Lemma LAT, Lemma DW-L.
- *uprime, multidescent, nested, escapeR2, Theorem U's chain:* as listed in c12alg §5 (FM/LD*/D-U′, CS′, SQF, child
  scan, recursion and IR at γ₂; Prop UB, Theorem 2L's proof; Lemmas RW, T_multi; §7.3; [pp] §§3, 5–9, Kern-U, L2-U,
  Lemma S″, Lemmas OS, W, F, T_rec, Theorem DD, LV6-U).
- *Published, used here:* Drappeau, Compositio 151 (2015), Lemme 5 (Fact D; stated there without proof, as a
  consequence of §3.3 of Harper's unrefereed preprint arXiv:1208.5992, so its justification rests on that preprint) and
  Lemme 1(i); de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), Thm 2.4(i); Hildebrand–Tenenbaum,
  Trans. AMS 296 (1986), Thms 2 and 3; Hildebrand, J. Number Theory 22 (1986), Thms 1 and 3 (Ψ = xρ(1 + O(E_H)) and the
  short-interval form, uniform for xy^{−5/12} ≤ z ≤ x; restated in Jain, arXiv:2502.10530, Lemma 3.2); elementary: Brunn–Minkowski in dimension 2, Gauss reduction, inclusion–exclusion.
  Not used: Harper arXiv 1208.5992's BV/BDH theorems, Saias, the large sieve (IK 7.13), de la Bretèche–Drappeau.
- *Published, used by the rest of the chain:* as Theorem U (lv5unif §10) and c12alg §5.

**Achieved α.** No unconditional α < 2/5. **Conditionally on c12alg's Theorem AB (a sketch) and the internal chain:
L_N[α] for every fixed α ∈ (3/8, 2/5); conditionally on Theorem AB₂ (a sketch with open obligations): L_N[3/8].** The
rational side imposes no further condition, but its v2 repairs are themselves unreviewed (codex: "calling only AB
unfinished is premature").

### 9.2 How much precision the rational blocks leave

Theorem RAT gives A_rat ≤ D^{−2α}λ^C and B_rat ≤ (D^{−2α} + D^{−2(3−7α)})λ^C (codex c12rat_2 #10: the second term
dominates for α > 3/8; at α = 7/18 it is D^{−5/9}). The window's budget is D^{−α}λ^C for A and D^{7α−3}λ^C for B
(c12alg Prop PB). So the margins are D^{α} for C₁ and D^{min(9α−3, 3−7α)} for C₂ — positive on (3/8, 2/5), and at
α = 3/8, against pc3's η ≍ D^{−3/8}, D^{3/8} for both. In the language of the task: the rational blocks do not
constrain α on [3/8, 2/5); c12alg's Theorem 2L-DW₂ uses B ≪ (D^{−2α} + D^{−2(3−7α)})λ^C, exactly what RAT-2 gives.

### 9.3 Toward 3/8

At α = 3/8 the rational blocks hold at second order (Theorem RAT in design (P), or (W) at α = 3/8). What remains is
c12alg §3.5: (a) the algebraic block's archimedean comparison must be second order (their Lemma PP has first-order error
λ^C/u ≍ η·λ^C), and (b) √(AB) ≪ η, which needs both A and B at precision D^{−3/8−ε}. Both would follow from a second-order
algebraic block. **Heuristic (not attempted; c12alg's file):** Prop PD should apply to the algebraic pair essentially
verbatim. Given g, Lemma TP-W/Cor TP-q describe (F(P₁), F(P₂))/q as a mixture over y of uniform measures on
L′ ∩ (parallelogram_y) with L′ = {x₁ ≡ λ^dx₂ (mod Δ/q)}; each parallelogram is convex, so Prop PD bounds the
y-conditional covariance by one-point variances plus character errors (moduli Δ/q ≤ 2S₂², averaged over pairs as in
Lemma FAM), and the mixture covariance Cov_y(S′_y, S″_y) is again bounded by one-point variances over y (Cauchy–Schwarz).
Then no first-order term survives, as here. The design (P) at 3/8 additionally needs c12alg's DW repairs re-checked
there ((SQ) is borderline in (W) at 3/8, satisfiable in (P)).

---

## 10. Experiments (pre-registered, `notes/c12rat_prereg.md`, commit `b9c2b8d`)

Toy checks only: nothing here tests asymptotic smallness. 13 predictions; 12 held, 1 failed (P2b, mis-specified).

### 10.1 E1: label shapes and thin pairs (exact; q = 10007 and 100003 all r; 1000 labels at q = 10007, S = 10)

| prediction | result | verdict |
|---|---|---|
| P1a [det] F(σ) ≤ π(√(σq) + 1)²/(2q), all q, σ | 0 violations (8 cells) | held |
| P1b F(σ)/σ ∈ [0.85, 1.05] at q = 100003 | 0.9542, 0.9547, 0.9543 (σ = 1/64, 1/16, 1/4); 3/π = 0.9549 | held |
| P1c [det] MULT ≤ 5S² + S and Lemma TF bound, every label | 0 violations (1000 labels × 3 K) | held |
| P1d E f(8)/E f(2), E f(32)/E f(8) ∈ [0.15, 0.45] | 0.348, 0.303 (E f = 0.504, 0.175, 0.053) | held |
| P1e skewed (σ < 1/16) share of E f(32) ∈ [0.3, 0.9] | 0.437 (60 skewed labels of 1000) | held |

Reading: the label-shape law is exactly the one-short-vector count (3/π), as Lemma LS uses. Thin pairs decay like
1/K: for regular labels (σ ≥ 1/4) K·f(K) ≈ 0.80, 0.77, 0.67 for K = 2, 8, 32 (the box's finite size truncates large K).
6% of labels are skewed at σ < 1/16 and carry 44% of the thin pairs at K = 32, as Prop SK expects.

### 10.2 E2: Lemma BTS (2000 random instances, y = 50)

| prediction | result | verdict |
|---|---|---|
| P2a [det] the explicit BTS inequality | 0 violations (also 0 for the dyadic form) | held |
| P2b median count/((Z/φ(D))·density) ∈ [0.8, 1.25] | 0.153 | **failed** |

*Disclosed deviation.* The frozen E2 code drew infeasible instances (Z up to 10⁸ against a 5·10⁷ table) and crashed after
5 instances. After the crash I changed `run_E2` to redraw (D, Z/D) until the instance fits (rejection sampling of the
registered distribution conditional on feasibility), and reran. The predictions were not changed.
*P2b post hoc (not pre-registered; `src/c12rat/posthoc.py`).* The registered "expected count" omitted the factor
Ψ_D/Ψ ≈ ∏_{p|D}(1 − p^{−β}) for smooth numbers coprime to D, which is large at y = 50. With the correct expectation
Ψ_D(I)/φ(D), the median ratio is 0.81 (10%–90%: 0, 2.03; the counts are small). The statistic was a sanity check of
equidistribution, not of the lemma. The lemma's loss (bound/count) has median 656 (10%–90%: 313–1527): the bound is
lossy at toy scale, as expected (it loses y^{O(γ)}log y, with γ ≈ 0.4 here).

### 10.3 E3: cell (Cauchy–Schwarz) decomposition of pair covariances (16 labels × 2 laws; y = 150, M = 1200, S = 8)

| prediction | result | verdict |
|---|---|---|
| P3a [det] \|C_ν\| ≤ B₁ for every label and law | 32/32 | held |
| P3b all-g law: \|R_lab\| ≤ 0.015 for every label | max 0.0133 | held (narrowly) |
| P3c prim law: \|R_lab − R_prim\| ≤ 0.02 for every label | max 0.0130 | held |
| P3d all-g law: median sd_pairs(R/S′S″)/\|R_lab\| ≥ 3 | 3.60 | held (narrowly) |
| P3e all-g law: mean ϕ_ex(skew) − mean ϕ_ex(uniform) ≥ −0.02 | +0.075 | held |

Means over the 16 labels: all-g law ϕ_ex = 0.167, C_ν = 0.160, R = 0.008, B₁ = 0.74; primitive law ϕ_ex = −0.021,
C_ν = 0.068, R = −0.089 (model R_prim within 0.013 on every label).
Reading:
- The cell decomposition captures the pair covariance: under the all-g law the principal (cell) part C_ν accounts for
  ϕ_ex up to 0.008 on average, and the within-cell residual cancels over the box (per-pair sd ≥ 3.6× the box mean).
- Under the primitive law, the residual is the primitivity correlation at the untyped primes 7 < p ≤ 150, predicted by
  the explicit local formula to ±0.013. In the theory these primes are typed (p ≤ P₀) or handled by Lemma PRIM
  (p > P₀); at toy scale γ ≈ 0.3 makes the term large (−0.09), asymptotically it is (γ²/P₀)·polylog.
- The Cauchy–Schwarz bound B₁ is ≈ 4.4× the cell part C_ν at toy scale: it is a bound, not an estimate.
- Skewed boxes have larger ϕ_ex (+0.075): their child values share the factor G(e₂), so their archimedean sizes are
  more correlated. This is a real second-order effect (bounded by the same one-point variances), not a failure of the
  method; the proof's difficulty with skew was technical (the BDH error's κ-dependence), and Prop SK removes it by
  weight.

---

## 11. Literature

- **Inputs:** §5.1 and §9.1, checked by the Opus literature critic (OK WITH FIXES; §12.1). Fact D: Drappeau's Lemme 5
  has the 1/φ(q) weight, Q ≤ x, conductors ≤ x^η and effectivity for A < 1, and is stated without proof as a consequence
  of Harper's preprint. Hildebrand's two theorems, with the rate O(log(u + 1)/log y) uniform for xy^{−5/12} ≤ z ≤ x: Jain,
  arXiv:2502.10530, Lemma 3.2(i)–(ii), and the HT survey §5. dBT Thm 2.4(i) via Drappeau Lemme 1(i); HT Thm 3 as in
  paper2/uprime.
- **Coupled pairs / covariance by one-point norms.** The inequality |∫(τ₁ − 1)(τ₂ − 1)dγ| ≤ ‖τ₁ − 1‖‖τ₂ − 1‖ for any
  coupling γ is Cauchy–Schwarz; codex found it blind (c12rat_1 (13)). Its use to get a *second-order* covariance from
  *first-order* one-point asymptotics, with exact cell factorisation through diagonal-invariant coupling sets (Lemma CP),
  is not something we found in the smooth-number literature (search result, not a novelty claim). The nearest relatives
  give first-order counts: de la Bretèche–Drappeau, JEMS 22 (2020), Cor. 4.2 (a first-order pair upper bound, small
  coefficients; Jain Lemma 3.3); Harper, Compositio 2016 (linear equations in smooth numbers, unbounded u);
  Balog–Blomer–Dartyge–Tenenbaum 2012, Lachand 2017 and Matthiesen–Wang 2025 (linear forms). None gives covariance-level
  control on families of large-modulus lattices.
- **BDH-type averages.** The L¹ → L² use of Drappeau's Lemme 5 is pc3 §4.1's; Harper's BDH (1208.5992 Thm 2) would give
  the same with an extra QΨ term and an unrefereed source. The parent-weighted family is handled here by pointwise
  domination (Lemma DOM; c12alg's Lemma WD is the two-sided version).
- **Short progressions of large modulus.** Lemma BTS is the standard divisor-switching argument (elementary; proof given),
  the large-modulus case of [pp] Lemma DS; codex (11) and c12alg's Lemma PU have the same statement.
- **Skewed lattices.** Lemma LS is the elementary count of short vectors; the label-shape law 3/π is visible in E1.

---

## 12. Reviews

**Codex blind `c12rat_1`** (xhigh, launched before any analysis here). Verdict OPEN. It independently found: Drappeau's
Lemme 5 suffices and LSQ is optional; a transfer lemma (our Theorem WA) with the domination quantity 𝔇 and the warning
that weights must be dominated pointwise; the BTS bound (its (11)) and the thin-pair losses e^{O(D^{1/8}log D)} /
e^{O(D^{1/4}log D)} at 3/8; the shape count; the coupling inequality (its (13)); and the idea that a common normalised
model would give second order from first-order inputs, which it left open. Lemmas TYPE/CP supply that step. Its
probabilities: 40% (3/8), 70% (some α < 2/5).

**Codex adversarial `c12rat_2`** (xhigh, on v1): **MAJOR REVISION.** "I found no counterexample to the intended
DW-specific RAT inequalities"; TYPE/CP/CPW checked on 59 624 congruence cases and 60 000 class-sum cases without
violation; WA and EN valid; library independence and the Σŵ comparison legitimate.

| # | sev. | finding | action (v2) |
|---|---|---|---|
| 1 | MAJOR | ADM's single exponent over the whole filter range ignores saddle-point drift | ADM split into local (ratio ≤ y²) and global one-sided parts (§1.4) |
| 2 | MAJOR | OPV: mid-prime factor Z has the wrong direction; fourth-moment line invalid | Z ≥ 1 corrected; moments recomputed (§5.6) |
| 3 | MAJOR | FO/EXC Rankin applies p^σ ≤ e to mid primes | small/mid split (§5.3) |
| 4 | MAJOR | S′ = E_νψ′ not exact with the boundary atom | OPV centred at c; B′, P_∂, P_ℰ carried (§§5.6, 8.1) |
| 5 | MAJOR | CELL counts use target lengths; union of rectangles not convex | actual lengths; convex Π⁺ (§3) |
| 6 | MAJOR | parent densities can vanish (hard filter): S_min argument fails | Lemma SPR (active points), absolute PUB, π_P-centring (§§5.5, 5.7, 8.3) |
| 7 | MAJOR | FAM needs pointwise spread bounds; "(Adm) weights" too general | Lemma SPR; theorem restricted to the DW weights |
| 8 | MAJOR | PRIM's (log y)²/P₀ is not second order | weighted prime sum: (γ²log²P₀ + E²)/(P₀ log P₀) (§8.2) |
| 9 | MAJOR | T₀ = e^{CD^θλ} cannot dominate c12alg's T (λ²) | T₀ := e^{C_TD^θλ²} |
| 10 | MAJOR | "common D^{−2α}" contradicts RAT-2's D^{−2(3−7α)} | margins recomputed (§§0, 9.2, 14) |
| 11 | MINOR | 𝒱_mid's extra D^{−0.2} wrong (−0.075 for parents at 3/8) | stated as o(1) |
| 12 | MINOR | Lemma LS's simplified constant false | fixed |

Codex's probabilities: 60% (Theorem RAT after repairs, with RAT-2's rate), 55% (Theorem 2L-RAT given AB).

**Opus proof critic** (read-only, on v1): design (W) **ACCEPT WITH FIXES**; the (P) clause and the statement of
Theorem 2L-RAT **MAJOR REVISION**. It checked the core mechanism line by line (TYPE's bijection; CP's exact split; ψ′
depending only on (I₁, t₁); PD's covariance identity with the atom; Cauchy–Schwarz against one-point variances; CPW's
positivity), found no place where Ξ multiplies a polylog-small quantity, and estimated 80% for Theorem RAT in (W).

| # | sev. | finding | action (v2) |
|---|---|---|---|
| M1 | MAJOR | threshold mismatch with c12alg (λ²) | as codex #9 |
| M2 | MAJOR | design (P): DOM(ii) needs b₂ < 1/4; DW lemmas unchecked at η₀ = D^{−1/2} | (P) only under a stated proviso (§1.1) |
| M3 | MAJOR | PRIM's bound too weak | as codex #8 |
| m1 | minor | PRIM: PD on Ω/pp′ (side δ_eM/y²); FO for p ≤ y via HT3; primes > y for F = 1 | fixed |
| m2 | minor | filter thresholds are not library endpoints | threshold intervals into ℰ (Lemma EXC(d)) |
| m3 | minor | small blocks may violate CPW's ε′ ≤ 1/4 | blocks of relative width ≤ min(h^{1/2}, 1/8) |
| m4 | minor | LS constant | fixed |
| m5 | minor | 𝒱_mid exponent; S′ = Eψ′; local law for mid powers; Z's index set | fixed |
| m6 | minor | ADM only for t′/t ≤ y² | as codex #1 |
| m7 | minor | FO Rankin for mid primes | as codex #3 |
| m8 | minor | BTS missing case δW > X₀ | hypothesis added |
| m9 | minor | PD bookkeeping (ℰ terms, ν(∗), Π̃, root-event sublattice, constants) | fixed |
| m10 | minor | FAM grid count, typo /q, C₂ window primes, |P′|^{d−1} spread | fixed; spread in θ and SPR |
| m11 | minor | parent spread justification ("≤ e" false) | Lemma SPR |
| m12 | minor | headline vs RAT-2; SK(ii) conditional; Hildebrand's rate essential | fixed; flagged for the literature critic |
| m13 | — | exponent bookkeeping checks out (θ = 2α − 1/2) | — |

### 12.1 Opus literature critic (read-only, on v2): OK WITH FIXES

Every published input of §5.1/§9.1 holds in the form used, **including Hildebrand's rate O(log(u + 1)/log y), uniform
for xy^{−5/12} ≤ z ≤ x with a constant depending only on ε** (Hildebrand JNT 1986 Thms 1, 3, via Jain arXiv:2502.10530
Lemma 3.2 and the HT survey §5). No range violation was found; **Theorem RAT survives; no verdict changed.** Fixes
applied in v2.1: Drappeau states Lemme 5 without proof, as a consequence of §3.3 of Harper's preprint (caveat added in
§§0, 5.1, 9.1); effectivity for A < 1 (A = 1/2); Pascadi citation flagged; Hildebrand theorem numbers and uniformity;
dBT via Drappeau Lemme 1(i); HT(α2) labelled a consequence of HT Thm 2; nearest relatives rewritten (dlB–D Cor. 4.2,
Harper 2016; first-order only); unreferenced "Friedlander", "Harper 2016 SNR3" and "Lee–Venkatesan Lemma 8.5" removed
from §11; Jain cited directly instead of codex for Hildebrand.

---

## 13. Reproduction

```
.venv/bin/python -m pytest -q tests/test_c12rat.py
.venv/bin/python -m src.c12rat.experiments E1 experiments/raw/c12rat_e1.jsonl
.venv/bin/python -m src.c12rat.experiments E2 experiments/raw/c12rat_e2.jsonl
.venv/bin/python -m src.c12rat.experiments E3 experiments/raw/c12rat_e3.jsonl
for E in e1 e2 e3; do .venv/bin/python -m src.c12rat.analyze ${E^^} experiments/raw/c12rat_$E.jsonl experiments/raw/c12rat_${E}_summary.json; done
.venv/bin/python -m src.c12rat.posthoc experiments/raw/c12rat_e2.jsonl experiments/raw/c12rat_e2_posthoc.json
.venv/bin/python -c "from fractions import Fraction as F; from src.c12rat import params as P; print(P.all_hold(F(3,8),'P'), P.all_hold(F(7,18),'W'))"
```

## 14. Candidate card

- **Mechanism.** Rational two-point blocks of C₁/C₂ via cells (value rectangles × local types), with exact factorisation
  through diagonal-invariant coupling sets (Lemma CP). Cauchy–Schwarz against one-point cell variances turns first-order
  inputs into a second-order output. Character errors are handled by Drappeau's Lemme 5 under pointwise weight
  domination, thin pairs by divisor switching, skewed labels by weight.
- **Result.** Theorem RAT (proved, internal; v2 repairs unreviewed): A_rat ≤ D^{−2α}λ^C and
  B_rat ≤ (D^{−2α} + D^{−2(3−7α)})λ^C on [3/8, 2/5), in design (W). Consequences:
  - Hyp RAT(α) and Hyp RAT₂ of c12alg follow.
  - Conditional on c12alg's sketches: L_N[α] for every α ∈ (3/8, 2/5) (via AB), and L_N[3/8] (via AB₂).
- **Failure scope.** No unconditional α < 2/5, because the algebraic blocks are sketches. Design (P) needs b₂ < 1/4 and
  a re-check of DW.
- **Next step.** A re-review of v2. Then c12alg's AB/AB₂ written to the last constant.
