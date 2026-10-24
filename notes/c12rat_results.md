# c12rat: the rational two-point blocks of C₁/C₂ (pc3 §7 items 1, 2, 4) — results

Agent `c12rat`, phase 27, 2026-10-24. **v1** (before the adversarial round and the critics; §12 records them).

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

**Plain language.** The rational two-point blocks of C₁ and C₂ — the parts of the two pair-variance hypotheses that
concern only the smoothness of the linear forms G_g(P′) = g₁a′ − g₀b′ — are **proved** (subject to this note's review),
and they are proved at **second order**: the incidence-weighted relative pair covariance is ≤ D^{−2α}(log D)^C, not
merely ≤ D^{−α}. That is far more than the window 3/8 < α < 2/5 needs (first order), and it is also enough at α = 3/8
itself. The idea is simple. For a pair of child points, cut the joint range of (G′, G″) into small cells (a rectangle in
value space times a local "type": valuations at small primes and at the primes of the pair's determinant). Inside a
cell the pair count factorises *exactly* into the product of two one-point cell densities, up to an error made of
non-principal characters. The covariance of two functions of the two coordinates is bounded by the product of their
one-point standard deviations (Cauchy–Schwarz), and a one-point density is known to first order by Hildebrand's
short-interval theorem and de la Bretèche–Tenenbaum. So first-order one-point asymptotics give a second-order
two-point bound, without Fourier tails, without computing any principal term, and uniformly in how nearly parallel the
two points are. The character errors are controlled on average over the pair family by Drappeau's published Lemme 5
(no Harper preprint at the level of statements), with a weighted application theorem whose only inputs are a pointwise
upper bound for the incidence weights and lattice multiplicities. Skewed labels and nearly parallel pairs are handled by
an elementary smooth-numbers-in-progressions upper bound whose loss e^{O(D^{θ}λ)}, θ < α, is beaten by the BDH saving.

**This does not yet give a rigorous algorithm below L[2/5].** C₁ and C₂ also have algebraic blocks (c12alg). With
c12alg's Theorem 2L-DW(α), whose hypotheses are the internal chain, **Hyp RAT(α)** (proved here) and its Theorem AB
(**a sketch**), the two-level algorithm would run in L_N[α] for every fixed α ∈ (3/8, 2/5). So the honest status is:
**L[α] for every α ∈ (3/8, 2/5), conditional on c12alg's Theorem AB being written out** (and on the internal chain). At
α = 3/8 the rational side is no longer the obstacle: what remains is the second-order algebraic block, for which the
cell/Cauchy–Schwarz device of this note is the natural tool (§9.3, not attempted: c12alg's file).

Precisely:
1. **Theorem RAT [proved].** In c12alg's window design 2L_DW(α), every fixed α ∈ [3/8, 2/5), and in pc3's design 2L′
   with DW at α = 3/8, cell by cell (environment cells of relative size D^{−A₀}, A₀ ≥ 1), with the smooth algebraic
   weights of DW(a)⁺ (or any weights satisfying the regularity (Adm) of §1.4):
   - (RAT-1) Σ_{𝔮,e : s(𝔮) ≤ T₀} w_{𝔮,e} ϕ^rat_{𝔮,e} ≤ (Σ w_{𝔮,e}/m_{𝔮,e}) + C D^{−2α}λ^{C} W + e^{−cD^α/λ²} μ₁;
   - (RAT-2) Σ_{𝔮,e} w_{𝔮,e} σ^rat_{𝔮,e} ≤ (Σ w_{𝔮,e}/z_{𝔮,e}) + C (D^{−2α} + D^{−2(3−7α)})λ^{C} W + e^{−cD^α/λ²} μ₁
     (the second term is the square of the first-order error 1/u₂ ≍ D^{7α−3} of the parents' B₂-smooth algebraic weight;
     at α = 3/8 both are D^{−3/4});
   - (RAT-0, skew) labels with skewness s(𝔮) = |e₂|/|e₁| > T₀ = e^{C D^{θ}λ} carry incidence weight
     ≤ e^{−CD^θλ}(log N)^{C}μ₁.

   Here ϕ^rat, σ^rat are the g-variances of c12alg's Lemma VD (law of total variance; codex's (15)), and
   θ = θ(α) < α is the loss exponent of §1.3 (θ = 1/4 at α = 3/8).
2. **Corollary [proved].** Hyp RAT(α) of c12alg §5 holds for every α ∈ (3/8, 2/5), with D^{−α} replaced by D^{−2α}λ^C
   in (R1) and D^{7α−3} by (D^{−2α} + D^{14α−6})λ^C in (R2). Its second-order form (D^{−3/4}λ^C) also holds at α = 3/8.
3. **The tools (task items 1, 2, 4).**
   - *Weighted application theorem* (Theorem WA, after codex's transfer lemma) with the family accounting Lemma FAM and
     the weight domination Lemma DOM: the non-principal errors of all cells of all pairs of all labels, weighted by the
     true incidence measure and each label's own normalisation, are ≤ e^{−cD^α/λ²}μ₁. Published input: Drappeau's
     Lemme 5 only (Fact D). The large sieve over child primes (pc3 Lemma LSQ) is not needed: every modulus is a small
     conductor.
   - *Skew* (Prop SK, Lemmas LS, TF, PUB): skewed labels have weight ≤ C/T₀ of their window (pointwise domination),
     their dead mass is bounded through C₂ (whose parent boxes are Euclidean, so never skewed), and inside non-skewed
     labels the thin pairs are a fraction ≤ CT₀/κ₀ bounded by Lemma PUB with loss Ξ = e^{O(D^θλ)}.
   - *Archimedean second order* (Prop PD, Lemmas FO, OPV): no Fourier tails and no smooth-pair upper bound inside the
     main argument; the boundary of the value region is charged by Cauchy–Schwarz to the cells' fourth moments.
     Saias is not needed; Hildebrand's short-interval theorem and dBT Thm 2.4 suffice.
4. **The achieved α.** None unconditionally. **L[α] for every fixed α ∈ (3/8, 2/5), conditional on c12alg's Theorem AB
   (sketch) and on the internal chain** (Theorem 2L-RAT(α), §9, with its complete dependency list). The rational blocks
   allow any precision loss up to D^{9α−3−ε} (relative to η ≍ D^{−(3−7α)}), i.e. D^{3/8−ε} at α = 3/8.
5. **Experiments (pre-registered; toy).** See §10.

**Probabilities (ours).** That Theorem RAT survives review as stated (possibly with repaired constants): 70%. That the
whole two-level algorithm in the window becomes rigorous once c12alg's Theorem AB is written: 55% (codex's blind
numbers: 70% for the rational blocks at some α < 2/5, 40% at 3/8).

---

## 1. Setting

### 1.1 Designs

Two designs are covered; both use the DW field law of pc3 §6.2 as repaired by c12alg §2 (w uniform on all integers of
(W_m, 2W_m], c uniform on C_H, F = wF̂_N + G·R_c, g uniform on primitive points of [M/2, M]²).
- **(W)** c12alg's window design 2L_DW(α), α ∈ (3/8, 2/5) (c12alg §1.1): d ≍ D^{1−2α}, k ≍ D^α, B = e^{aD^αλ},
  B₂ = e^{b₂D^{4α−1}λ}, η₀ log N ≍ D^{3α−1/2}λ, special-q ℓ ∈ J = (Y, 2Y], Y = L[1 − α], level-1 box radius V₁ with
  V₁² ≍ YS₁², child boxes of radius S₂ (S₂² = L[α, σ₂]) in Gauss-reduced coordinates, child primes q ∈ (B, B₂].
- **(P)** pc3/uprime's 2L′ at α = 3/8 with DW: B₂ = L[1/2, b₂], η₀ = D^{−1/2}, otherwise as uprime §1.

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
it is the exponent of every "loss" Ξ := e^{O(D^θλ)} below ((|P′||P″|)^γ, (V₁²)^γ, (X/X⁻)^γ). **On [3/8, 2/5), θ < α**,
so Ξ^{O(1)} = ε^{−o(1)}, and also Ξ^{O(1)} ≤ y^{1/12}. Further, with P₀ := ⌈D^{7/10}⌉:
- every modulus is a small conductor: log K ≤ (η_D/4) log x;
- ω(K) = o(P₀) for every child and parent modulus, so Σ_{p | K, p > P₀} p^{−β} → 0;
- P₀# = e^{(1+o(1))P₀} ≤ x^{o(1)}, far below every cell length.
`params.margins` checks these (and the window/3/8 comparisons of §1.2) for α = 3/8 (both designs) and on a grid of
(3/8, 2/5); at α = 2/5 only "B₂ ≪ V₁" becomes borderline.

### 1.4 Admissible weights

A family of weights a_P: [X⁻_P, ∞) → [0, 1] is **admissible with exponent γ_a and error E_a** if for all t ≤ t′ in the
range of |G_g(P)|, g ∈ Ω_e:
  (Adm)  a_P(t′)/a_P(t) = (t/t′)^{γ_P}(1 + O(E_a + log(t′/t + 1)/log t)),  with 0 ≤ γ_P ≤ γ_a.

**Lemma ADM [proved].** The child weights a_{P′} of DW(a)⁺ (with the root coset) and the parent weights a_P are
admissible: for children γ_a ≤ C/(aD^α) and E_a ≤ Cλ³D^{−α}; for parents γ_a ≤ C/(b₂D^{4α−1}) and
E_a ≤ Cλ³(D^{7α−3} + D^{1−4α}) (the first-order errors 1/u + E_H of the respective target sets).
*Proof.* By DW(a)⁺/Cor TP-q (c12alg §§2.6, 3.2), given g with |G_g(P)| = t, F(P) (resp. F(P)/q) is within TV
ε_TV ≤ 7qN^{−η₀/2} of the law of t·(Z + U), where Z has a law independent of t (the remaining coordinates of the
smoothing decomposition, plus the w-shift, which is ≤ N^{−η₀/4}·scale and moves the law by ≤ ε_TV), and U is uniform on
an interval [0, ℓ_P/t) of fixed relative length ≥ 1/(2d) of the support. So a_P(t) = E_Z[Ψ_𝒯(t(Z + [0, L)))/(tL)] + O(ε_TV)
for the target set 𝒯 (𝒮_B for children, the acceptable B₂-smooth set for parents). For each Z the interval t(Z + [0, L))
has relative length ≥ 1/(2d) ≥ y^{−5/12}; by Hildebrand's two theorems (Facts H1, H2 below) its 𝒯-density is
ρ̂_𝒯(t|Z|_+)(1 + O(E_H)), and by Fact HT3 (iterated ≤ 3 times; Fact dBT1 for ratios beyond y²) ρ̂_𝒯(t′x)/ρ̂_𝒯(tx) =
(t/t′)^{γ}(1 + O(1/u)) for t′/t ≤ y², while the ratios t′/t > y² occur on a set of mass ≤ 40y^{−2} (Lemma LC) where
only the upper bound is used. (For 𝒯 = qZ ∩ 𝒮_B the factor q is fixed, so 𝒯-densities scale as 𝒮_B-densities at x/q.)
The layer |Z| ≤ N^{−η₀/8} has ϑ-mass ≤ dN^{−η₀/8} (small layers, [pp] smooth1) and is absorbed into E_a. For the
acceptable set, the squarefree/discriminant subtraction is a relative O((log N)²/B) (uprime Lemma SQF). ∎

The DW(a)⁺ weights depend on g only through t (c12alg). Arithmetic factors (e.g. a design with gcd(b′, g₁) = 1, which
makes the weight depend on gcd(G′, b′)) are also allowed if they are constant on the P₀-types of §4 and within
1 ± ε_h, ε_h ≤ 1/P₀, otherwise (§8.1, remark).

---

## 2. Main statements

For a label 𝔮, a cell e and child points P′ ≠ P″ write S′ := E_{𝒢_e} f′(G′), Pair := E_{𝒢_e} f′(G′)f″(G″),
κ := |P′||P″|/K (K = |det(P′, P″)|; κ ≥ 1, large for nearly parallel pairs).

**Theorem RAT.** In designs (W) (α ∈ [3/8, 2/5)) and (P), with admissible weights, there are constants C, c > 0 and
parameters T₀ = κ₀^{1/3} = e^{C₀D^θλ} (§7) such that:
- (RAT-1) Σ_{𝔮,e: s(𝔮)≤T₀} w_{𝔮,e}Φ_{𝔮,e} ≤ C D^{−2α}λ^{C}·W + e^{−cD^α/λ²}μ₁;
- (RAT-2) Σ_{𝔮,e} w_{𝔮,e}Φ^par_{𝔮,e} ≤ C (D^{−2α} + D^{−2(3−7α)})λ^{C}·W + e^{−cD^α/λ²}μ₁;
- (RAT-0) Σ_{𝔮,e: s(𝔮)>T₀} w_{𝔮,e} ≤ e^{−C₀D^θλ/2}μ₁.

Here W = Σ w_{𝔮,e} is the expected number of parent–child incidences, μ₁ = E N₁, and W ≥ (1 − o(1))μ₁ (Lemma DOM(iii)).

**Corollary RAT(α).** For every fixed α ∈ (3/8, 2/5), Hyp RAT(α) of c12alg §5 holds: (R1) Σ_{s(𝔮)≤T}(w/W)ϕ^rat
≤ D^{−2α}λ^C + Σ(w/W)/m and (R2) Σ(w/W)σ^rat ≤ (D^{−2α} + D^{14α−6})λ^C + Σ(w/W)/z, with T := T₀; note
D^{14α−6} ≪ D^{7α−3} for α < 3/7. [The one-point terms are c12alg's
Prop CLB and z_𝔮 ≍ L[1 − α]·m₁/polylog; c12alg's (R1)/(R2) contain them implicitly.]
*Proof.* Divide Theorem RAT by W ≥ μ₁/2. ∎

The proof occupies §§3–8. Its skeleton, for one non-thin pair (κ ≤ κ₀) and one cell, is **Prop PD** (§8.1):
  |Pair − S′S″| ≤ √(V′V″) + Σ_c |Err(c)| + O(h^{1/4}S′S″) + tiny,
where V′, V″ are one-point variances of cell densities (≤ CD^{−2α}λ⁶S′², pointwise in the pair: Lemma OPV) and Err(c)
are character errors (Theorem WA); primitivity at primes > P₀ adds O((γ²/P₀)λ² + P₀^{−2})S′S″ (Lemma PRIM). Thin pairs (κ > κ₀) use Lemma PUB, skewed labels Prop SK.

---

## 3. Geometry [all proved; elementary]

**Lemma LC (log-concave marginals).** Let Ω be a closed square of side δM ≥ 8|P| in R² and P = (a, b) primitive, and
let 𝒢 ⊂ Ω ∩ Z² contain at least a fraction 1/2 of Ω ∩ Z² (true for the P₀-primitive points of Ω). Put
n(g) := det(P, g) = G_g(P), x_max := max_Ω |n| and H := length of n(Ω) ≥ δM|P|. Then for every interval I,
  Pr_{g∈𝒢}[n(g) ∈ I] ≤ 8(|I| + 1)/H,
and for 0 < t ≤ 1/2, Pr_{g∈𝒢}[|n(g)| ≤ t·x_max] ≤ 40t + 16/H. Consequently, for k ≤ 8 and 0 ≤ σ ≤ 1/2,
E_{g∈𝒢}[log^k(x_max/|n|)(x_max/|n|)^σ·1(|n| ≥ X⁻)] ≤ C_k + (16/H)(log(x_max/X⁻))^k(x_max/X⁻)^σ.
*Proof.* The g ∈ Z² with n(g) = n form a coset of ZP; its intersection with Ω has ≤ c(n)/|P| + 1 points, where c(n) is
the length of the chord {det(P, ·) = n} ∩ Ω. Lines with consecutive n are at distance 1/|P|, so |Ω| = |P|^{−1}∫c. The
function c is concave on n(Ω) (Brunn–Minkowski in dimension 2), so ∫c ≥ c_max·H/2. Hence
#{g ∈ Ω∩Z² : n(g) ∈ I} ≤ (|I| + 1)(c_max/|P| + 1) ≤ 4(|I| + 1)|Ω|(1 + |P|/c_max)/H, and c_max ≥ δM ≥ 8|P|,
#(Ω∩Z²) ≥ |Ω|/2·… give the first bound. For the second: if 0 ∈ n(Ω), H ≥ x_max; if not, the event is empty unless the
nearest endpoint is ≤ tx_max, in which case H ≥ (1 − t)x_max. The moment bound follows by integrating the tail. ∎

**Lemma LS (label shapes).** For a prime q and T ≥ 1, #{r mod q : s(q, r) > T} ≤ (π/2)(√(2q/(√3T)) + 1)² ≤ 2q/T + 3.
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

**Lemma CELL (tilings).** Fix a non-thin pair (κ ≤ κ₀) and a cell e, and let Π := A(Ω_e), A(g) := (det(P′, g),
det(P″, g)), a parallelogram of area K|Ω_e|. Let ℓ_r (ℓ_c) be the longest chord of Π parallel to the n₁-axis (n₂-axis),
H₁, H₂ the widths of Π in the two coordinates, x′_max := max_Π|n₁|, and put z₁ := hℓ_r, z₂ := hℓ_c with h := D^{−A₁}.
Let 𝒳 := {±X_*(1 + h_*)^k : k ∈ Z} ∪ {0} with h_* := hδ_e/(4κ₀) and X_* := 1 (the **universal grid**, the same for all
labels, cells and pairs). Then:
- (i) z₁ ≥ hδ_e x′_max/(√2κ) ≥ 2h_*x′_max, and z_i ≫ K; so each coordinate axis (on |n_i| ≤ 2x_max) can be partitioned
  into **blocks**: unions of consecutive grid intervals of total length in [z_i, 2z_i] (blocks inside [−z_i, z_i] are
  allowed to be shorter, "small blocks").
- (ii) The rectangles R = I₁ × I₂ (I_i blocks) meeting Π but not contained in Π ("boundary rectangles") cover an area
  ≤ 20h|Π|, and each rectangle contains (z₁z₂/K)(1 + O(K/z_min)) points of L := A(Z²) = {n₁ ≡ ιn₂ (mod K)}.
- (iii) The set of n with |n₁| ≤ z₁/h^{1/2} ("small blocks" and their neighbours) has μ-mass ≤ 40h^{1/2} + tiny, μ the
  law of A(g) for g ∈ 𝒢_e (by Lemma LC, since H₁ ≥ ℓ_r).
*Proof.* (i) A chord of Π parallel to the n₁-axis is the image of a chord of Ω_e in direction P″, along which n₁ moves
by K per unit of the parameter t in g = g₀ + tP″; the central chord has length ≥ δ_eM, so ℓ_r ≥ Kδ_eM/|P″|, while
x′_max ≤ √2M|P′|. So z₁/x′_max ≥ hδ_e/(√2κ). The grid spacing at |n| ≤ 2x_max is ≤ 2h_*x_max ≤ z₁/2. And
z_i ≥ hδ_eM/κ₀ ≫ K since M = N^μ = L^{ω(1)}. (ii) Convexity gives |Π| ≥ ℓ_rH₂/2 and |Π| ≥ ℓ_cH₁/2, and ℓ_rℓ_c ≤ ℓ_rH₂.
The union of boundary rectangles lies within box-distance (2z₁, 2z₂) of ∂Π, of area ≤ 4z₁H₂ + 4z₂H₁ + 16z₁z₂ ≤
(16h + 32h²)|Π|. L contains KZ², and each K×K square contains exactly K points of L. (iii) Lemma LC with I = [−z₁h^{−1/2},
z₁h^{−1/2}], H₁ ≥ ℓ_r = z₁/h. ∎

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
  |Σ a b 1_C − (|C|/φ²)AB| ≤ (|C|/φ²)AB·[2(Ẽ_u + Ẽ_v) + 4ε′²].
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

- **Fact D** (Drappeau, Compositio 151 (2015), Lemme 5; Pascadi, Compositio 161 (2025), Lemma 3.2; = paper2 Fact `fact:D`).
  There are η_D, c_D, δ_D > 0 such that for (log x)^{c_D} ≤ y ≤ x and Q ≤ x,
  Σ_{k≤Q}φ(k)^{−1}Σ_{χ mod k, 1<cond χ≤x^{η_D}}|Ψ(x, y; χ)| ≪ Ψ(x, y){H(u)^{−δ_D} + y^{−δ_D}} (with A = 1/2 and the
  factor (log x)^{−A} discarded), H(u) = e^{u/log²(u+1)}. Its published proof refers to §3.3 of Harper's unpublished
  arXiv 1208.5992 (paper2 discusses this). We use **only** this statement: no BDH theorem of Harper is used.
- **Fact dBT1** (de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), Thm 2.4(i); Drappeau Lemme 1(i)). For 2 ≤ y ≤ x and
  1 ≤ d ≤ x, Ψ(x/d, y) ≤ C Ψ(x, y)d^{−β(x,y)}; equivalently ρ̂(x/d) ≤ C d^{γ(x)}ρ̂(x), ρ̂(x) := Ψ(x, y)/x.
- **Fact HT3** (Hildebrand–Tenenbaum, Trans. AMS 296 (1986), Thm 3). Ψ(cx, y) = c^{β(x,y)}Ψ(x, y)(1 + O(1/u + log y/y))
  for x ≥ y ≥ 2, 1 ≤ c ≤ y. Iterating at most k times (with HT Thm 2's control of the drift of β, as in lv5unif Lemma
  U1): Ψ(x/d, y) = d^{−β(x,y)}Ψ(x, y)(1 + O_k(1/u)) for 1 ≤ d ≤ y^k, x/d ≥ x^{1/2}.
- **Fact HT(α2)** (HT Thm 2, as paper2 (α2)): (1 − β(z, y)) log y ≤ log(u_z log y) + C_α log log y; in our ranges
  γ = 1 − β ≤ C₁/(aD^α), and y^{γ} = e^{O(λ)}.
- **Fact H1, H2** (Hildebrand, J. Number Theory 22 (1986) 289–307; restated as Theorem 5.1 of the Hildebrand–Tenenbaum
  survey, JTNB 5 (1993); theorem numbers in the original UNVERIFIED). For y ≥ exp((log log x)^{5/3+ε}):
  Ψ(x, y) = xρ(u)(1 + O_ε(E_H)) and, for xy^{−5/12} ≤ z ≤ x, Ψ(x + z, y) − Ψ(x, y) = zρ(u)(1 + O_ε(E_H)),
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
Rankin with σ := 1/log P₀ (so p^σ ≤ e for p ≤ P₀),
Σ_{e|𝔪, se>y²}(se)^{−β} ≤ s^{−β}(y²/s)^{−σ}∏_{p|𝔪}(1 + e·p^{−β}) ≤ s^{−β}·e^{−log y/log P₀}(log P₀)^{C};
since log y/log P₀ ≍ D^α, this tail is smaller than every power of D. So the count is
(z/x)s^{−β}Ψ(x)[g_𝔪(β) + O((E_H + 1/u)∏_{p|𝔪}(1 + p^{−β}))], and ∏(1 + p^{−β})/∏(1 − p^{−β}) ≤ exp(3Σ_{p|𝔪}p^{−β})
≤ C(log P₀)³. Finally (z/x)s^{−β}Ψ(x) = (z/s)s^{γ}ρ̂(x). ∎

### 5.4 Lemma BTS (smooth numbers in a progression in a short interval) [proved; elementary]

Let k ≥ 1, a ∈ Z, X₀ ≥ W ≥ 2y, and I = (X₀, X₀ + Z] with W := Z/k. Then
  #{n ∈ I : n ≡ a (mod k), P⁺(n) ≤ y} ≤ 4W log(2y)·max_{W/y≤w≤2W}ρ̂(w),
and if log W ≥ (1/2)log(X₀ + Z), this is ≤ e^{O(λ)}·W·((X₀ + Z)/W)^{γ}ρ̂(X₀ + Z) (Facts dBT1, HT(α2)).
*Proof.* Put δ := gcd(a, k); n = δn′ with n′ ≡ a/δ (mod k/δ) a unit, n′ ∈ I/δ (if δ is not y-smooth there is nothing
to count). Each y-smooth n′ > W has a divisor m₁ ∈ (W/y, W] (take the largest product of its smallest prime factors
that is ≤ W). For fixed m₁ (smooth, coprime to k/δ), m₂ := n′/m₁ lies in an interval of length Z/(δm₁) in a fixed
class mod k/δ: ≤ W/m₁ + 1 values. Sum over m₁: W Σ_{W/y<m₁≤W}1_S(m₁)/m₁ + Ψ(W), and dyadic blocks give the bound.
(The exact inequality, before the dyadic step, is `bts.rhs`; E2 checks it on 2000 instances.) ∎
This is codex's (11) (found blind) and the large-modulus case of [pp] Lemma DS; c12alg's Lemma PU uses the same device.

### 5.5 Lemma PUB (universal pair upper bound) [proved]

For every label (or parent set), cell Ω (square of side ≥ δ_eM/y), and every pair P′ ≠ P″ (thin or not):
  Pair ≤ Ξ·S′S″,  Ξ := e^{O(λ)}(|P′||P″|/δ_e)^{γ}·(a″_max/a″_min)·y^{γ} ≤ e^{C_ΞD^θλ},
with a″_max/a″_min the ratio of the weight over the range of |G″| on Ω (≤ (X″/X″⁻)^{γ_a}(1 + O(E_a)) by (Adm)).
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
c12alg's Lemma PU (alg and rat parts) is the same bound with the same exponent; we need only its rational part.

### 5.6 Lemma OPV (one-point variance of cell densities) [proved]

Fix a non-thin pair, its cell structure (Lemma CELL, types relative to K), and the first coordinate. For a block I
and a type t let ψ′(I, t) := Σ_{m∈U(I,t)}f′(s_t m)/#U(I, t), U(I, t) := {m ∈ I/s_t : P₀-free, (m, M_t) = 1}. Let ν₁ be
the law of (block(G′), t(G′)) for g ∈ 𝒢 (P₀-primitive points of Ω). Then
  E_{ν₁}[(ψ′ − S′)²] ≤ (𝒱₀ + 𝒱_mid(K))·S′²,
  𝒱₀ := C[(γ + γ_a)²(log P₀)² + E_FO² + E_a² + γ_a²h^{1/2} + h^{1/4}],  𝒱_mid(K) := Cγ²(log P₀)²ω(K)/P₀.
In our ranges 𝒱₀ ≤ CD^{−2α}λ^{6} (with A₁ ≥ 8α + 1, i.e. h^{1/4} ≤ D^{−2α}), and **pointwise in K**,
𝒱_mid(K) ≤ CD^{−2α}λ²·ω(K)/P₀ ≤ CD^{−2α−0.2}λ² (children: ω(K) ≤ D^{4α−1+o(1)} resp. D^{1/2+o(1)}; parents: ω(K) ≤
D^{1−α+o(1)}; P₀ = D^{7/10}; §1.3). Moreover E_{ν₁}[(ψ′/S′)⁸] ≤ C (bounded relative moments).
*Proof.* Let x₁ := max_Ω|G′| and c := a′(x₁)ρ̂(x₁)·g_{P₀#}(β)∏_{p≤P₀}(1 − 1/p)^{−1}·∏_{p|K, mid}((1 − p^{−β})/(1 − 1/p))
·∏_{p|K, large}(1 − 1/p)^{−1}. Since Var ≤ E[(ψ′ − c)²] and S′ = E_{ν₁}ψ′, it suffices to bound E[(ψ′/c − 1)²] by
(𝒱₀ + 𝒱_mid)/4 (then S′ ≥ c/2 and the claim follows).
- *Regular blocks* (left end ξ ≥ z₁h^{−1/2}, so relative width ≤ h^{1/2}) and types with s_t ≤ S₀: Lemma FO with
  𝔪 = P₀#·rad(K/δ)_mid (Σ_{p|K mid}p^{−β} ≤ 2ω(K)/P₀ → 0 by §1.3), Lemma FIB for #U, and (Adm) for the weight across
  the block give ψ′(I, t) = c·X(ξ)·Y(t)·Z(δ)·(1 + O(E_FO + E_a + γ_ah^{1/2})), where
  X(ξ) := a′(ξ)ρ̂(ξ)/(a′(x₁)ρ̂(x₁)), Y(t) := s_t^{γ}, Z(δ) := ∏_{p|δ, mid}(1 − 1/p)/(1 − p^{−β}) ∈ [1 − 4Σ_{p|δ}p^{−β}, 1].
  (We replaced β(ξ, y) by β(x₁, y) in s^{γ}: the drift is O(log s·log(x₁/ξ)/(log x·log y)), negligible.)
- *Archimedean factor.* For ξ ≥ x₁y^{−2}: X = (x₁/ξ)^{γ+γ_{P′}}(1 + O(1/u + E_a)) (Fact HT3 for ρ̂, (Adm) for a′). For all
  ξ: X ≤ C(x₁/ξ)^{γ+γ_a} (Fact dBT1). Lemma LC then gives E[(X − 1)⁴] ≤ C(γ + γ_a)⁴ + C(1/u + E_a)⁴ (the part
  ξ < x₁y^{−2} has mass ≤ 40y^{−2} and fourth moment ≤ Σ_{k≥2log y}40e^{−k(1−4γ−4γ_a)} ≤ y^{−1}), and E[X⁸] ≤ C.
- *Local factor.* Under ν₁, Pr[s_t(G′) ≡ 0 (mod s)] ≤ 2/s for every P₀-smooth s ≤ S₀·y (G′ is equidistributed mod s
  on Ω up to a factor 1 + O(s|P′|/(δ_eM)), and P₀-primitivity changes local densities by ≤ 1 + 1/p). Hence
  E[(Y − 1)⁴] ≤ Cγ⁴(log P₀)⁴ + CΣ_{p|K mid}(2/p)(p^{γ} − 1)⁴ (the small-prime part is a finite Euler product,
  Σ_{s P₀-smooth}s^{−1+4γ}log⁴s ≤ (log P₀)^{C}, since P₀^{4γ} = O(1)). For a mid prime, (p^γ − 1)^k/p ≤
  max(γ^k log^k p/p, p^{kγ−1}) ≤ γ^k(log P₀)^k/P₀ on (P₀, e^{1/γ}] and ≤ e^{1−1/(2γ)} beyond (kγ ≤ 1/2), so the mid
  part is ≤ Cγ⁴(log P₀)⁴ω(K)/P₀ ≤ 𝒱_mid(K)². Further Pr[s_t > S₀] ≤ S₀^{−σ}E[s^σ] ≤ e^{−D^{(θ+α)/2}/(2log P₀)}, and
  E[(Z − 1)²] ≤ CΣ_{p|K mid}p^{−1}(p^{−β} − p^{−1})²p² ≤ 𝒱_mid(K).
- *Small blocks* (|G′| ≤ z₁h^{−1/2}) carry mass ≤ 40h^{1/2} (Lemma CELL(iii)) and ψ′/c ≤ C·X·Y·(1 + O(1)) there, so they
  contribute ≤ (40h^{1/2})^{1/2}(E[(ψ′/c)⁴] + 1)^{1/2} ≤ Ch^{1/4}.
- Combine with (u − 1) = (X − 1)YZ(1+e) + (Y − 1)Z(1+e) + (Z − 1)(1+e) + e and Cauchy–Schwarz. ∎

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
s(𝔮) ≤ T₀, cell e, ordered non-thin pair, cell c of the pair, coordinate i, endpoint), library j = (endpoint ξ ∈ 𝒳,
s ≤ S₀) with u_j = 1_{𝒜_{ξ,s}}, U_j = Ψ(ξ/s, y), modulus k(t) = M_c, and weights
  c_t := ŵ_{q,e}·m_{𝔮,e}^{−2}·(N(c)/Z)ψ′(c)ψ″(c)·C(ξ/z_i)²(log P₀)⁴
— one has 𝔇 ≤ Ξ^{C}·(κ₀/(hδ_e))³·Σ_{q∈win,r,e}ŵ_{q,e}. The same holds for the C₂ family (parent pairs, library
𝒜^J_{ξ,s}).
*Proof.* *Size of a term.* N(c)/Z ≤ ν-mass of the cell; ψ′ψ″ ≤ Ξ²S_typ² and m_𝔮 ≥ n·S_min ≥ nS_typ/Ξ (one-point
lower and upper bounds from Lemmas FO, OPV and Fact dBT1 for every box point; the spread over a box of skewness ≤ T₀ is
≤ Ξ by §1.3), so c_t ≤ ŵ_{q,e}Ξ⁴n^{−2}·ν(c)·(ξ/z_i)²(log P₀)⁴.
*H_j(k) for fixed j = (ξ, s) and k.* The prime q of the window dividing k: ≤ 1 + log(4S₂²)/log B = O(1) choices
(k ≤ 4Q₂S₂², q > Q₂ ≥ B). The pair's K is one of the ≤ τ(s)2^{ω(s)} numbers kδ/r_δ with δ | s, r_δ | rad(δ)
(M_τ = (K/δ)r_δ and δ = gcd(s, K)); each gives Δ = K/q, which has ≤ 10S² + 2S ordered pairs per label (Lemma MULT).
For one pair, the cells with an endpoint ξ in coordinate i and s_i = s have total ν-mass ≤ 2Pr[s_i(G) = s] ≤ 4/s
(Lemma OPV's local law), and (ξ/z_i)² ≤ 2κ₀²/(hδ_e)² for regular blocks (≤ h^{−1} for small blocks). Summing over the
labels (all r and e): H_j(k) ≤ O(1)·max_qŴ_q·Ξ⁴n^{−2}(10S² + 2S)·τ(s)2^{ω(s)}(4/s)·(2κ₀²/(hδ_e)²)(log P₀)⁴, with
Ŵ_q := Σ_{r,e}ŵ_{q,e}.
*The max over k.* φ(k) ≤ k ≤ 4Q₂S₂² (M_τ ≤ K). So max_kφ(k)H_j(k) ≤ C Q₂·max_qŴ_q·Ξ⁴(κ₀/(hδ_e))²(log P₀)⁴·τ(s)2^{ω(s)}/s.
*The sum over j.* #(𝒳 ∩ relevant range) ≤ log(x_max/x_min)/h_* ≤ D·4κ₀/(hδ_e), and Σ_{s≤S₀}τ(s)2^{ω(s)}/s ≤
(log S₀)⁴ = D^{2(θ+α)}. So 𝔇 ≤ C D^{1+2(θ+α)}Q₂ max_qŴ_q Ξ⁴(κ₀/(hδ_e))³(log P₀)⁴.
*Comparison.* Σ_{q∈win}Ŵ_q ≥ (Q₂/(2 log Q₂))min_qŴ_q, and max/min of Ŵ_q over a dyadic window is ≤ C (ŵ is an explicit
smooth function of q: q^{−1}ρ₂(X_F/q)/q up to (1 + O(q/V₁)), and ρ₂(X_F/q) varies by ≤ 2^{γ₂}·(1 + o(1)) on the window).
Hence 𝔇 ≤ C D^{3}log Q₂·Ξ⁴(κ₀/(hδ_e))³(log P₀)⁴·Σ_{q∈win}Ŵ_q. For C₂ the same computation uses the parent
multiplicity ≤ 5#(𝒫 ∩ Λ_𝔮) per K, K ≤ 2V₁², #(𝒫 ∩ Λ_𝔮)² pairs and #(𝒫∩Λ_𝔮) ≍ #𝒫/q: the ratio
Q·5#(𝒫∩Λ)/(#(𝒫∩Λ)²·#q_win) ≤ C d log V₁ log Q₂ is again polylogarithmic. ∎

**Corollary WA′.** For each window, Σ over the family of ŵ_{q,e}m^{−2}·Σ_c|Err(c)| (the character part of Cor CPW) is
≤ ε_lib·Ξ^{C}(κ₀/(hδ_e))³·Σ_{win}ŵ; summed over the ≤ log₂B₂ windows and with Lemma DOM(ii), ≤ e^{−cD^α/λ²}μ₁.
*Proof.* Cor CPW bounds |Err(c)| by N(c)ψ′ψ″·2(Ẽ′ + Ẽ″) with Ẽ for the interval sequence 1_{𝒜}·1_{I/s} mod M_c. An
interval sum is a difference of two endpoint sums, so |û_J|² ≤ 2|û_{ξ₁,s}|² + 2|û_{ξ₂,s}|², and the principal sum is
≥ (z/ξ₂)Ψ(ξ₂/s)(log P₀)^{−2}/C (Lemma FO with 𝔪 ⊃ rad(M_c)_y; Σ_{p|M_c mid}p^{−β} → 0). This gives the weights c_t
above; apply Theorem WA with ε = ε_lib (Lemma EN, Q = 4Q₂S₂² ≤ x^{η_D/2}) and Lemma FAM. Ξ^{C}(κ₀/(hδ_e))³ = e^{O(D^θλ)}
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

Put T₀ := e^{C_TD^θλ} and κ₀ := T₀³ (C_T large, but C_T D^θλ ≤ (c_ε/10)D^α/λ²).

**Prop SK [proved].**
- (i) Σ_{q,r,e: s(q,r)>T₀}w_{(q,r),e} ≤ (2/T₀ + 3/B)(log N)^{C}μ₁ (Lemmas LS, DOM).
- (ii) The dead-incidence mass (Prop DC_e) of the skewed labels is ≤ W_sk + (W_sk·Σ_{all}wσ²)^{1/2}, W_sk the mass in (i);
  with Σ_{all}wσ² ≤ (log N)^{C}μ₁·(B_rat + B_alg + one-point) (C₂ for all labels: RAT-2 and c12alg's Theorem AB(ii)),
  this is ≤ T₀^{−1/3}μ₁.
- (iii) In a label with s(𝔮) ≤ T₀, the thin pairs (κ > κ₀) are a fraction ≤ 60T₀/κ₀ = 60T₀^{−2} of all pairs, and
  contribute ≤ 120T₀^{−2}Ξ³ ≤ T₀^{−1} to Φ_𝔮 (Lemma PUB and the spread ratio).
- (iv) Parent sets 𝒫 ∩ Λ_𝔮 lie in a Euclidean square; their thin pairs are a fraction ≤ 3/κ₀ + O(q/V₁) for **every**
  label (skewed child lattices included), contributing ≤ 4Ξ_p/κ₀ to Φ^par_𝔮.
*Proof.* (i) By DOM(i) the weight of a residue r is ≤ ŵ_{q,e}, independent of r, and by LS at most 2q/T₀ + 3 residues
are skewed; sum over e and q and use DOM(ii) (Σ_{r}ŵ = qŵ). (ii) Prop DC_e: the dead mass of a set 𝔖 of (label, cell)
is ≤ Σ_𝔖w(p + σ√p) ≤ Σ_𝔖w + (Σ_𝔖w)^{1/2}(Σ_𝔖wσ²)^{1/2} by p ≤ 1 and Cauchy–Schwarz. (iii) Lemma TF with s(𝔮) ≤ T₀;
each thin pair has |Cov| ≤ Pair + S′S″ ≤ (Ξ + 1)S′S″ ≤ (Ξ + 1)Ξ²m²/n² (spread ratio). (iv) Lemma TF (parents) and PUB. ∎

So skewness is handled without any shape-dependent analytic estimate: a skewed child lattice costs nothing beyond its
(pointwise dominated) weight, because the companion statement C₂ is insensitive to the child lattice's shape, and
transversality is never used — nearly parallel pairs are either thin (upper bound PUB, loss Ξ) or handled by the cell
decomposition, whose only κ-dependence is the factor (κ₀/(hδ_e))³ in Lemma FAM. (Codex's example in uprime §7.3 item 3 —
a label with 1 transverse point in 719 — is a skewed label in this sense.) c12alg's Theorem AB(iii) gives the same
conclusion for the whole ϕ via its Lemmas PU and WD; the two arguments are independent.

---

## 8. Assembly

### 8.1 Prop PD (one non-thin pair) [proved]

Let Ω be a square of side ≥ δ_eM/y, 𝒢 its P₀-primitive points, P′ ≠ P″ child (or parent) points with κ ≤ κ₀, f′, f″
admissible, S′ := E_𝒢f′(G′), Pair := E_𝒢f′(G′)f″(G″). Then
  |Pair − S′S″| ≤ (V′V″)^{1/2} + 𝔈 + C(h^{1/4} + L^{−c})S′S″,
where V′ := E_{ν₁}[(ψ′ − S′)²] + C h S′² ≤ (𝒱₀ + 𝒱_mid(K))S′² (Lemma OPV) and likewise V″, and
𝔈 := Z^{−1}Σ_c N(c)ψ′(c)ψ″(c)·2(Ẽ′(c) + Ẽ″(c)) (sum over interior and boundary rectangles c, Cor CPW).
*Proof.* Let Z := #𝒢 and use the cells c = (R, τ) of Lemma CELL × Lemma TYPE (admissible types; the n with
gcd(n₁n₂, K_large) > 1, with ℓ | n for a prime ℓ ∈ J dividing K (parents), or with s_i > S₀, form an exceptional set
ℰ, handled by Lemma EXC below).
- *Interior cells.* By Lemma TYPE(i) the n ∈ c ∩ 𝒢 are the (s₁m₁, s₂m₂) with (m̄₁, m̄₂) ∈ C_τ, m_i ∈ U(I_i, t_i). Lemma
  CP/CPW with a(m) := f′(s₁m)1_{U₁}(m), b likewise, and Lemma FIB give Σ_{n∈c}f′f″ = N(c)ψ′(c)ψ″(c)(1 + O(tiny))
  + Err(c), |Err(c)| ≤ N(c)ψ′ψ″[2(Ẽ′ + Ẽ″) + 4ε′²], with ε′ ≤ C(γ_ah^{1/2} + E_a) the weight oscillation over the
  block (for small blocks ε′ ≤ C, but their mass is ≤ 40h^{1/2}); and Σ_{n∈c}f′(n₁) = N(c)ψ′(c)(1 + O(tiny)) (the fibre
  over each m₁ has the same size, Lemma FIB). Here ψ′(c) depends on c only through (I₁, t₁) (it is the ψ′ of Lemma OPV).
- *Boundary rectangles.* 0 ≤ Σ_{n∈c∩𝒢∩Π}f′f″ ≤ Σ_{n∈c∩𝒢}f′f″, and the full rectangle obeys the interior formula; the
  same for the singles.
- *Covariance identity.* Let ν be the probability measure on {interior cells} ∪ {∗} with ν(c) = N(c)/Z and the atom ∗
  carrying the rest (boundary and ℰ), and put ψ′(∗) = ψ″(∗) := 0. Write S′ = E_νψ′ + B′ + O(tiny) with
  0 ≤ B′ ≤ b′ := Z^{−1}Σ_{∂}N(c)ψ′(c) (boundary part), and similarly S″. Then
  Pair − S′S″ = Cov_ν(ψ′, ψ″) + Z^{−1}Σ_{int}Err(c) + P_∂ − (E_νψ′)B″ − B′(E_νψ″) − B′B″ + tiny,
  with 0 ≤ P_∂ ≤ Z^{−1}Σ_∂(N(c)ψ′ψ″ + |Err(c)|). And |Cov_ν(ψ′, ψ″)| ≤ (Var_νψ′·Var_νψ″)^{1/2}, Var_νψ′ ≤ E_ν[(ψ′ − S′)²]
  ≤ E_{ν₁}[(ψ′ − S′)²] + ν(∗)S′² (ν's first marginal is dominated by the one-point law ν₁ of Lemma OPV), ν(∗) ≤ 20h + e^{−λ²}.
- *Boundary terms.* By Cauchy–Schwarz, b′ ≤ (Z^{−1}Σ_∂N)^{1/2}(Z^{−1}Σ_∂Nψ′²)^{1/2} and Z^{−1}Σ_∂Nψ′ψ″ ≤
  (Z^{−1}Σ_∂Nψ′²)^{1/2}(Z^{−1}Σ_∂Nψ″²)^{1/2}; and Z^{−1}Σ_∂Nψ′² ≤ (Z^{−1}Σ_∂N)^{1/2}(Z^{−1}Σ_{R⊂Π̃}Nψ′⁴)^{1/2}, where Π̃ ⊃ Π
  is the union of all rectangles meeting Π. The last sum is a one-point fourth moment over a region of area ≤ (1 + 20h)|Π|
  whose n₁-marginal is again the marginal of a convex set, so ≤ CS′⁴ by the bounded relative moments of Lemma OPV
  (moments, not maxima: the spread ratio Ξ never enters). With Z^{−1}Σ_∂N ≤ 20h (Lemma CELL(ii)) all boundary terms are
  ≤ Ch^{1/4}S′S″.
- *Weights.* 4ε′² ≤ C(γ_a²h + E_a²) + (small-block mass)·C ≤ 𝒱₀; this part of Err(c) sums to ≤ 𝒱₀S′S″·… and is put into
  V′ (it is ≤ (V′V″)^{1/2} after adjusting C). ∎

**Lemma EXC (exceptional types) [proved].** Let ℰ be the n ∈ L ∩ Π with (a) gcd(n₁n₂, K_large) > 1, or (b) (parents) ℓ | n₁n₂
for one of the ≤ 2 primes ℓ ∈ J dividing K, or (c) s₁ > S₀ or s₂ > S₀. Then Z^{−1}Σ_{n∈ℰ}f′(n₁)f″(n₂) ≤ e^{−λ}S′S″ and
Z^{−1}Σ_{n∈ℰ}f′(n₁) ≤ e^{−λ}S′, and ν(ℰ) ≤ e^{−λ}.
*Proof.* (a) contributes 0 to every f-sum (a y-smooth n has no prime > y). For (b) and (c), the n with σ | n₁ for a fixed
integer σ (σ = ℓ, or σ a P₀-smooth s > S₀ dividing the P₀-part) come from the g on the lines g* + ZP′ with σ | G′(g);
running Lemma PUB's proof on these lines only gives Z^{−1}Σ_{σ|n₁}f′f″ ≤ Ξ·S′_σS″, where S′_σ := Z^{−1}Σ_{σ|G′}f′ ≤
e^{O(λ)}σ^{−1+γ}S′ (Lemma FO applied to the class σ | n₁; for σ = ℓ ∈ J the density of J·𝒮_B among multiples of ℓ is
≤ (log Y)·(its global density)). Summing, (b) gives ≤ Ξ(log Y)²Y^{−1}S′S″ and (c) gives ≤ ΞΣ_{s>S₀}s^{−1+γ}S′S″ ≤ Ξ·S₀^{−1/(2log P₀)}(log P₀)^C S′S″
(Rankin, exponent 1/log P₀); both are ≤ e^{−λ}S′S″ because Y and S₀^{1/log P₀} = e^{D^{(θ+α)/2}/log P₀} exceed every
power of Ξ = e^{O(D^θλ)}. The singles are the
same computation without the partner. ∎

**Remark (why first order suffices).** Every first-order quantity (E_H, 1/u, E_a, the archimedean factor's deviation
γ log(x₁/|n|), the local factor's deviation γ log s) appears only inside one-point L² norms, which are multiplied
together. No pair main term is ever computed or compared with the product of marginals. This is codex's "(13)" and
"common normalized model", with the model chosen per pair: since S′ = E_νψ′ holds exactly for *each* pair's own cell
structure, different partners may use different cells (codex's caution "a different one-point approximation for every
partner does not suffice" concerned approximations of the marginals; here the marginals are exact).

### 8.2 Lemma PRIM (primitivity at primes > P₀) [proved]

Let Ω, P′, P″, f′, f″ be as in Prop PD and let 𝒢^prim ⊂ 𝒢 be the primitive points. Then
  |Cov_{𝒢^prim}(f′, f″) − Cov_𝒢(f′, f″)| ≤ C[(γ² + E_FO²)(log y)²/P₀ + P₀^{−2}]S′S″ + Σ_{1<e≤y², P⁻(e)>P₀}e^{−2+2γ}|R_e|,
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
≤ CΣ_p p^{−2}(γ log p + E_FO)²p^{2γ}a₁b₁ + C(Σ_pp^{−2}…)² ≤ C(γ² + E_FO²)(log y)²P₀^{−1}S′S″. ∎
(At toy scale, where γ ≈ 0.3, this term is visible: E3's R ≈ −0.08 under the primitive law is exactly the
untyped-prime primitivity covariance, §10.3.)

### 8.3 Proof of Theorem RAT

Choose h := D^{−A₁} with A₁ := 8α + 1, δ_e := D^{−A₀} with A₀ ≥ 1, P₀ := ⌈D^{7/10}⌉, T₀ := e^{C_TD^θλ} and κ₀ := T₀³,
with C_T so large that T₀ ≥ Ξ^{10} and so small that κ₀^{10} ≤ ε_lib^{−1/10} (possible since θ < α).

*(RAT-1).* Fix a window (Q₂, 2Q₂], a label with s(𝔮) ≤ T₀ and a cell e; all box points have one-point densities S′
within a factor Ξ of each other (§1.3), and S′ = E ψ′ exactly for each pair's cell structure. Split the ordered pairs:
- **Thin pairs** (κ > κ₀): by Prop SK(iii) their contribution to Φ_𝔮 is ≤ 120T₀^{−2}Ξ³ ≤ T₀^{−1}.
- **Non-thin pairs.** Prop PD and Lemma PRIM give, for each pair,
  |Cov(f′, f″)| ≤ (V′V″)^{1/2} + 𝔈 + Σ_{e}e^{−2+2γ}(…)_{Ω/e} + C[h^{1/4} + (γ² + E_FO²)λ²/P₀ + P₀^{−2} + L^{−c}]S′S″.
  Here V′ ≤ (𝒱₀ + 𝒱_mid(K))S′² and V″ likewise, with the **same** K, so (V′V″)^{1/2} ≤ (𝒱₀ + 𝒱_mid(K))S′S″ ≤
  C D^{−2α}λ^{6}S′S″ (Lemma OPV, pointwise in K). Summing over pairs and dividing by m² = (Σ S′)² gives at most
  CD^{−2α}λ⁶ + C(h^{1/4} + λ²D^{−2α}/P₀ + P₀^{−2}), i.e. ≤ C D^{−2α}λ⁶ by the choice of A₁ (h^{1/4} ≤ D^{−2α−1/4}) and
  P₀^{−2} = D^{−7/5} ≤ D^{−2α} on α ≤ 7/10.
- **Character errors.** The 𝔈-terms of all non-thin pairs of all labels with s(𝔮) ≤ T₀ and all cells (including the
  Ω/e terms of Lemma PRIM, with weights e^{−2+2γ}), multiplied by w_{𝔮,e}/m_{𝔮,e}² ≤ ŵ_{q,e}/m², form the family of
  Lemma FAM; Cor WA′ bounds their total by e^{−cD^α/λ²}μ₁.
Multiplying the per-label bounds by w_{𝔮,e} and summing: Σ_{s(𝔮)≤T₀}w_{𝔮,e}Φ_{𝔮,e} ≤ (CD^{−2α}λ⁶ + T₀^{−1})W +
e^{−cD^α/λ²}μ₁, which is (RAT-1). ∎

*(RAT-2).* Identical, with the parent library 𝒜^J_{ξ,s}, the parent multiplicities, no skew (Prop SK(iv)), and the
filter |G(P)| ∈ [MV₁/32, 2MV₁]. The one-point densities of parents are within a **polylogarithmic** factor of each
other (|P| ∈ [2V₁/3, √2V₁] by the b-window, so the algebraic weight's |P|^{d−1} varies by ≤ e, and the filter confines
|G(P)| to a range of ratio 64), so S-weighted averages over P cost only polylogs. The filter's one-point relative
variance r_P is O(1/pass fraction) for the P whose threshold MV₁/32 falls inside the range of |G(P)| on Ω_e
("straddling" P) and ≤ Cδ_e otherwise; straddling P are a fraction ≤ Cδ_e of 𝒫 for each cell, and S_P·r_P = O(S_typ),
so (V′V″)^{1/2} ≤ S′S″(r′ + r″)/2 sums to ≤ m·Σ_PS_Pr_P ≤ Cδ_eλ^{C}m² (pc3's O(1) environment term M1 is exactly what
the cells remove; A₀ ≥ 1 > 2α makes it ≤ D^{−2α}), and the J-structure (the ≤ 2 primes ℓ ∈ J dividing K are in ℰ: by PUB their
contribution is ≤ Ξ_p(log Y)/Y ≪ L^{−1}). The one-point FO for J·𝒮_B: Σ_{ℓ∈J}#{m ∈ J/(sℓ) : …}, Lemma FO for each ℓ. The
mid-prime term uses Lemma MULT (parent): Pr_{P′}[p | det(P, P′)] ≤ (1/p)(1 + O(pq/V₁)) — exact up to 1 + o(1), so no
spread ratio appears. Choose δ_e = D^{−A₀} with A₀ ≥ 1. ∎

*(RAT-0).* Prop SK(i). ∎

---

## 9. Composition and the achieved exponent

### 9.1 Theorem 2L-RAT(α) [proved conditional on c12alg's Theorem AB (sketch) and the internal chain]

**Theorem 2L-RAT(α).** Let α ∈ (3/8, 2/5) be fixed and take c12alg's design 2L_DW(α) with T := T₀ (any T₀ = e^{C_TD^θλ}
with C_T large enough to dominate c12alg's T = Ξ_{c12alg}³D^αλ^C works for both notes). Assume:
- (i) the project's internal chain, as in Theorem U and multidescent Theorem 2L;
- (ii) c12alg's Theorem AB (§3.3 there; **sketch**, "no gap known") and the lemmas c12alg labels proved.
Then FACTOR_2L in design 2L_DW(α) is Las Vegas and factors every N in expected time L_N[α, c(α)].
*Proof.* c12alg's Theorem 2L-DW(α) assumes (i), Hyp RAT(α) and Theorem AB. Corollary RAT(α) proves Hyp RAT(α) (with
room: D^{−2α}λ^C in place of D^{−α}λ^C and D^{7α−3}λ^C). ∎

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
- *Published, used here:* Drappeau, Compositio 151 (2015), Lemme 5 (Fact D; its proof refers to Harper's unpublished
  computations) and Lemme 1(i); de la Bretèche–Tenenbaum, Ramanujan J. 9 (2005), Thm 2.4(i); Hildebrand–Tenenbaum,
  Trans. AMS 296 (1986), Thms 2 and 3; Hildebrand, J. Number Theory 22 (1986) (Ψ = xρ(1 + O(E_H)) and the short-interval
  form for z ≥ xy^{−5/12}); elementary: Brunn–Minkowski in dimension 2, Gauss reduction, inclusion–exclusion.
  Not used: Harper arXiv 1208.5992's BV/BDH theorems, Saias, the large sieve (IK 7.13), de la Bretèche–Drappeau.
- *Published, used by the rest of the chain:* as Theorem U (lv5unif §10) and c12alg §5.

**Achieved α.** No unconditional α < 2/5. **Conditionally on c12alg's Theorem AB (a sketch) and the internal chain:
L_N[α] for every fixed α ∈ (3/8, 2/5).** The rational side imposes no further condition.

### 9.2 How much precision the rational blocks leave

Theorem RAT gives A_rat, B_rat ≤ D^{−2α}λ^C. The window's budget is D^{−α}λ^C for A and D^{7α−3}λ^C for B (c12alg
Prop PB), and pc3's raw requirement is η ≍ D^{−(3−7α)}. So the rational side tolerates a precision loss up to
D^{9α−3−ε} relative to η (D^{3/8−ε} at α = 3/8, D^{3/5−ε} at α = 2/5). In the language of the task: the rational blocks
do not constrain α at all on [3/8, 2/5).

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

- **Inputs:** §5.1 and §9.1. Fact D's statement was checked by codex (c12rat_1) against Drappeau's paper (the 1/φ(q)
  weight, Q ≤ x, effectivity for A < 1); Hildebrand's short-interval theorem is confirmed through the HT survey
  Theorem 5.1 (codex), original numbering UNVERIFIED; dBT Thm 2.4(i) and HT Thm 3 as in paper2/uprime.
- **Coupled pairs / covariance by one-point norms.** The inequality |∫(τ₁ − 1)(τ₂ − 1)dγ| ≤ ‖τ₁ − 1‖‖τ₂ − 1‖ for any
  coupling γ is Cauchy–Schwarz; codex found it blind (c12rat_1 (13)). Its use to get a *second-order* covariance from
  *first-order* one-point asymptotics, with exact cell factorisation through diagonal-invariant coupling sets (Lemma CP),
  is not something we found in the smooth-number literature (search result, not a novelty claim). The nearest relatives
  are the local-global factorisations in linear-forms problems (Balog–Blomer–Dartyge–Tenenbaum, Lachand,
  Matthiesen–Wang) and the W-trick of additive combinatorics, which need bounded u or bounded coefficients.
- **BDH-type averages.** The L¹ → L² use of Drappeau's Lemme 5 is pc3 §4.1's; Harper's BDH (1208.5992 Thm 2) would give
  the same with an extra QΨ term and an unrefereed source. Density-relative BV/BT bounds exist (Lee–Venkatesan Lemma
  8.5; Harper 2016 SNR3); the parent-weighted family is handled here by pointwise domination (Lemma DOM; c12alg's
  Lemma WD is the two-sided version).
- **Short progressions of large modulus.** Lemma BTS is the divisor-switching argument (Friedlander; [pp] Lemma DS;
  Harper 2016's SNR3 has the same (Xq/Z)^{1−α} shape); codex (11) and c12alg's Lemma PU have the same statement.
- **Skewed lattices.** Lemma LS is the elementary count of short vectors; the label-shape law 3/π is visible in E1.

---

## 12. Reviews

*(To be completed: codex adversarial round `codex/c12rat_2`; Opus proof critic; Opus literature critic.)*

**Codex blind `c12rat_1`** (xhigh, launched before any analysis here). Verdict OPEN. It independently found: Drappeau's
Lemme 5 suffices and LSQ is optional; a transfer lemma (our Theorem WA) with the domination quantity 𝔇 and the warning
that weights must be dominated pointwise (its (6)–(7); our Lemma DOM); the BTS bound (its (11)) and the thin-pair loss
e^{O(D^{1/8}log D)} (children) / e^{O(D^{1/4}log D)} (parents) at 3/8; the shape count (its (8)); the coupling
inequality (13) and the idea that a common normalised model would give second order from first-order inputs, which it
left open ("the exact missing step: the actual principal term and both actual marginals must be represented by the same
normalized model"). Our Lemma CP/TYPE supplies exactly that step (the marginals are exact for each pair's cells). Its
probabilities: 40% (3/8), 70% (some α < 2/5).

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

- **Mechanism.** Rational two-point blocks of C₁/C₂ by cells (value rectangles × local types) with exact factorisation
  through diagonal-invariant coupling sets (Lemma CP), Cauchy–Schwarz against one-point cell variances (first-order
  inputs ⇒ second-order output), character errors by Drappeau's Lemme 5 under pointwise weight domination, thin pairs by
  divisor switching, skewed labels by weight.
- **Result.** Theorem RAT (proved, internal): A_rat, B_rat ≤ D^{−2α}λ^C on [3/8, 2/5). Corollary: Hyp RAT(α) of c12alg.
  Conditional on c12alg's Theorem AB (sketch): rigorous internal L_N[α] for every α ∈ (3/8, 2/5).
- **Failure scope.** No unconditional α < 2/5 (algebraic block). At 3/8 the algebraic block is needed at second order.
- **Next step.** Write c12alg's Theorem AB to the last constant; apply Prop PD to the algebraic pair for 3/8 (§9.3).
