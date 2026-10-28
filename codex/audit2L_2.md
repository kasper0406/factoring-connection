You are an adversarial referee (analytic number theory, smooth numbers). Read-only repository at
/Volumes/git/factoring-connection; scratch under your current directory. This is a focused follow-up to your own blind
report codex/audit2L_1.out (read it), on notes/audit2L_blind_claims.md.

Your weakest link was RAT-2's parent normalization (blind excerpt B §8.3, the hard-filter bullet, and Lemma SPR
(parents) in B §5.7): the proof bounds Σ_{P,P′} c_P c_{P′} r_P by (m^par)²·polylog·(δ_e + 𝒱₀), with
polylog = e^{O(λ)} = D^{O(1)}, but claims λ^C. You also noted that S_typ is built on ρ₂(X_F) while the parent weight
lives at the root-conditioned scale X_F/q, which differs by q^{1−β₂} ≍ u₂ log u₂ (a power of D).

The auditor proposes the following repair. Decide whether it is correct and complete (closing RAT-2 at the claimed
rate C(D^{−2α} + D^{−2(3−7α)})λ^C on [3/8, 2/5)), and whether anything else in RAT-2 still loses a power of D.

**Lemma SPR-P′ (within-label comparison of parent OPV constants).** Fix a label 𝔮 = (q, r) and an environment
cell e. Let 𝒫_e(𝔮) := {P ∈ 𝒫 ∩ Λ_𝔮 : π_P > 0} (P passes the hard filter for some g ∈ Ω_e). Then there is an
absolute C with c_P ≤ C·c_{P′} for all P, P′ ∈ 𝒫_e(𝔮). Moreover, for every P ∈ 𝒫 ∩ Λ_𝔮 with π_P = 0, f_P ≡ 0 on
Ω_e.

*Proof.* The OPV constant is c_P = a_P(x₁)·ρ̂_J(x₁)·Λ_P(K), where x₁ = max_{Ω_e}|G_g(P)| and Λ_P(K) is the product
of local factors (small, mid and large primes of the pair's K). For P ∈ 𝒫_e(𝔮):
- (i) x₁ ∈ [MV₁/32, 2MV₁], because |G_g(P)| ≤ |g||P| ≤ 2MV₁ and P passes somewhere. So ρ̂_J(x₁)/ρ̂_J(x₁′) ∈ [64^{−γ}, 64^{γ}]·(1 + O(1/u)), by Fact HT3/dBT1 at y = B. That is 1 + o(1).
- (ii) a_P(t) is the B₂-smoothness probability of the DW(a)⁺ model law at |G| = t, restricted to qZ and divided by
  q, at scale tL_P/q with L_P = (2H + 1)|P|_∞^{d−1}. All P in one label share the **same q**. Their |P|_∞ lie in
  ((1 − 1/d)V₁, V₁] (b-window), so L_P/L_{P′} ∈ [e^{−1}, e]. By Lemma ADM (global one-sided bound plus local two-sided
  bound), a_P(t)/a_{P′}(t′) ∈ [C^{−1}, C] for t, t′ ∈ [MV₁/32, 2MV₁]. The scales differ by a factor ≤ 64e(d + 1),
  costing (64e(d + 1))^{γ₂}(1 + O(E_a)) = 1 + o(1). Here γ₂ ≍ D^{1−4α}, so γ₂ log d → 0. The model-law average
  over Z contributes the same bounded factor for every P, because the law of Z is that of a sum of d independent
  segments whose shape depends on P only through the normalised coefficient vector; its |Z|^{−γ₂}-moment lies in
  [(Cd)^{−γ₂}, C] = [1 − o(1), C].
- (iii) Λ_P(K) = exp(O(ω(K)/P₀ + ω(K)/y)) = 1 + o(1).
Hence c_P/c_{P′} ≤ C. If π_P = 0 then the filter indicator in f_P vanishes on Ω_e. ∎

**Consequence (RAT-2 hard-filter step, corrected).**
- By the angle argument of SPR, at least half of 𝒫 ∩ Λ_𝔮 is *active* (π_P = 1), and every active P has
  S_P ≥ c_P/2 (Lemma OPV).
- So m^par ≥ (n/2)·min_{𝒫_e} c_P/2 ≥ n·max_{𝒫_e} c_P/(4C).
- With r_P := π_P(1 − π_P) + 𝒱₀^par, and r_P irrelevant for π_P = 0:
  Σ_{P≠P′}|Cov(f_P, f_{P′})| ≤ Σ_{P,P′∈𝒫_e} c_Pc_{P′}(r_P + r_{P′})/2 + (𝔈 and smaller terms)
  ≤ n·max c·Σ_{P∈𝒫_e} c_P r_P ≤ (n max c)²(Cδ_e + 𝒱₀^par) ≤ 16C²(m^par)²(Cδ_e + 𝒱₀^par).
  This uses #{P : 0 < π_P < 1} ≤ Cδ_e·n per cell, the SPR angle argument for straddlers.
- So the only losses are absolute constants. Every remaining polylog loss sits inside 𝔇 (character errors) or inside
  thin-pair terms, where exponentially small savings absorb it.
- No common comparator S_typ is needed, so the X_F vs X_F/q issue does not arise: all comparisons are within one
  label, at one q.

Questions:
1. Is Lemma SPR-P′ true as stated? In particular, is (ii) justified by Lemma ADM and the DW(a)⁺ model law (blind
   excerpt B §1.4, F §2.6), uniformly in the label? Does the Z-law average really give a P-uniform bounded factor?
2. Is the straddler count ≤ Cδ_e·n per (label, cell) correct, including for skewed child lattices? Parent sets are
   Euclidean squares, but the label's lattice Λ_𝔮 restricts P.
3. Does the corrected chain give RAT-2 at C(D^{−2α} + D^{−2(3−7α)})λ^C, given RAT's other lemmas? Does any other
   step of RAT-2, or of AB §§9–11, multiply 𝒱₀ (or another merely D-polynomially small quantity) by a power of D?
4. Your MINOR items (EXC(c) mid-prime part, PRIM parent J-content, the W-P(P3) citation): are they now the only
   remaining items, and is each repairable without changing an exponent?

Output: a one-line verdict on the repair (CORRECT / CORRECT WITH FIXES / INCORRECT), then a short issue table, then
your updated verdicts for (i) α ∈ (3/8, 2/5) and (ii) α = 3/8 (ACCEPT / ACCEPT WITH FIXES / MAJOR REVISION / REJECT).
Be concrete.
