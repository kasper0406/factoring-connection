# audit2L: the two blind Opus critics (verbatim)

Both critics were given only `notes/audit2L_blind_claims.md` plus the underlying source notes. They were told not to
read any Reviews/Verdict section or probability estimate, any audit note, or any codex file before forming their
verdict. Both ran read-only and one at a time. The reports below are verbatim. The only change is formatting: the
harness's indentation has been removed.

---

## Opus critic 1: the analytic core (RAT, AB/AB₂, the cell device, Drappeau Lemme 5 / Theorem WA, ranges)

Agent id a73a82bb1c2026998. Run time about 15 min.

**Verdict (my scope only: the analytic core, RAT + AB/AB₂):**
- (i) α ∈ (3/8, 2/5): **ACCEPT WITH FIXES.**
- (ii) α = 3/8 with a > 2c₂: **ACCEPT WITH FIXES.** The fixes are the same as for (i). The a > 2c₂ step in §11 is correct.

Both verdicts treat Drappeau's Lemme 5 as a black box. I found no fatal gap. One step in RAT-2 does not prove what it claims as written; it can be repaired without changing any exponent.

## Issue table

| # | Sev. | Location | Issue | Exponent or constant? | Fix |
|---|---|---|---|---|---|
| 1 | **MAJOR** (gap as written, fixable) | Excerpt B §8.3 (RAT-2, hard-filter bullet); Lemma SPR (parents) | The text derives Σ_{P,P′} c_P c_{P′} r_P ≤ (m^par)²·**polylog**·(δ_e + 𝒱₀) and then concludes C(D^{−2α} + D^{−2(3−7α)})λ^C. In these notes "polylog" means e^{O(λ)} = D^{O(1)}λ^{O(1)}. A D^C factor of unknown size times 𝒱₀ ≍ D^{−2α}λ⁶ fails both budgets. The B-budget margin is only D^{9α−3} (D^{3/8} at 3/8), and at 3/8 the η-margin is only D^{3/8}. This is exactly the fatal pattern your item 3 asks about. Source: SPR(parents) states only S_P ≥ S_typ/polylog and c_P ≤ c_typ·polylog. | Neither, once fixed. As written, an unknown D-power. | The polylog is not real. Every P with π_P > 0 has x₁(P) ∈ [MV₁/32, 2MV₁] and an algebraic scale within a factor 64ed. So ρ̂_J and ρ̂_{B₂} vary by (64ed)^{γ+γ₂} = 1 + o(1). The local factors are e^{O(ω(K)/P₀)} = 1 + o(1). By OPV, S_P ≥ c_P/2 for active P. So c_P ≤ C·c_typ for all P with π_P > 0, and m^par ≥ n·c_typ/C. This gives Σ c c r ≤ C(m^par)²(Cδ_e + 𝒱₀). Points with π_P = 0 have f_P ≡ 0 and drop out. State this O(1) comparison as a lemma and delete "polylog" from SPR's lower bound. Polylog losses then remain only inside 𝔇, where they are harmless. |
| 2 | MINOR (external risk) | Fact D, Lemma EN, EN-A | I checked Drappeau arXiv 1307.7554v3, Lemme 5, verbatim. Hypotheses: (log x)^c ≤ y ≤ x and Q ≤ x. The sum is Σ_{q≤Q} φ(q)^{−1} Σ over **all** χ mod q with 1 < cond χ ≤ x^η, imprimitive included. The bound is Ψ(x,y){H(u)^{−δ}(log x)^{−A} + y^{−δ}}, with η, c, δ independent of A and effective for A < 1. This matches every use here: Q = 2S₂² resp. 2V₁²/q ≤ x^{η_D}, and the energies over all χ ≠ χ₀ mod 𝔨 are of that form. The lemma has **no proof** in the paper ("découle des calculs de la section 3.3 de [Har12b]"). | — | Nothing to fix internally. Disclose it as the one unrefereed published input, as the notes already do. |
| 3 | MINOR | AB §10(iii) and FAM-A step 3: "Prop CLB per point" | Excerpt G proves CLB with a global envelope X_{F2} and says "In all cases \|P′\| ≤ V₂". That is false for skewed labels: \|e₂\| can be ≍ q, so \|P′\| can reach ≍ S₂q ≫ V₂ ≍ S₂√q. (iii) needs the per-point bound m ≥ Σ_P ϖ_P(log N)^{−C} at each point's own scale. That bound is plausible: fibres M/(8\|P′\|) ≥ N^{η₀/2} still hold, and the DW(a)⁺ intervals have relative length ≥ 1/(2d) at the point's own scale. But it is asserted, not written. | No | Write CLB per point with X_{F2}(P′) := (2d+1)MH\|P′\|^d/q, and drop the "≤ V₂" claim. |
| 4 | MINOR | EXC-A large classes; FO tails | The Rankin bound uses P₀^{−1/2+γ} as if p^γ were bounded for all p \| Δ′, p ≤ y. In fact p^γ ≤ y^γ = e^{O(λ)}. | No | Split at p = e^{1/γ}: Σ p^{−1/2+γ} ≤ eωP₀^{−1/2} + ω·e^{O(λ)−cD^α}. |
| 5 | MINOR | FAM-A step 1 | It cites "OPV-A step 6" for the pointwise bound Ψ_i ≤ e^{O(λ)}δ^γc_i. That bound is step 5. | No | Correct the reference. |
| 6 | MINOR | AB §10 / Lemma TF | TF uses Euclidean κ = \|P′\|\|P″\|/K, while AB's ϑ uses \|·\|_∞. | Constant only | Note the factor 2. |
| 7 | MINOR | LF-A | The bound "≤ 16γ log Δ′ = O(λ + D^{2−5α}λ)" is loose. Harmless, since each such prime contributes e^{−cD^α} (children) or e^{−cD^{4α−1}} (parents). | No | — |
| 8 | MINOR | RAT-2 δ_e term | The bound is only valid after fix 1 (A₀ ≥ 1 > 2α). With a genuine polylog in front, it would need A₀ larger than an unknown constant. | No | Covered by fix 1. |

## Single weakest link in scope

**RAT-2's hard-filter π_P-centring (excerpt B §8.3, resting on SPR(parents)).** It is the only place where a factor e^{O(λ)} multiplies the D-polynomially small one-point variance 𝒱₀. As written, the conclusion does not follow. The repair (an O(1) comparison of OPV constants across the parents of a label) is routine but must be stated.

The external weakest link is Fact D: Drappeau's Lemme 5 is stated without proof and rests on §3.3 of Harper's unrefereed preprint.

## Items 1–5

**1. Theorem RAT.** All but the last two lemmas check out:
- LC, CELL, TYPE, CP, CPW, FIB, FO, BTS, PUB, ADM, SPR(children), EXC, PRIM, WA, DOM, FAM, Cor WA′, Prop PD and SK are correct as structured.
- CP is an exact Cauchy–Schwarz plus Parseval identity on a diagonal-invariant set.
- CPW's constant 3(Ẽ_u + Ẽ_v) + 4ε′² checks out for ε′ ≤ 1/4.
- WA is a one-line Hölder step that is valid for any assignment j(t), k(t).
- In FAM, the q-part of the modulus is spread over ≍ Q₂/log Q₂ window primes, and DOM(i) is r-uniform, so 𝔇 ≲ e^{O(D^θλ²)}Σŵ.
- PRIM's signed-weight identity is correct, and its second-order term is (γ²log²P₀ + E_FO²)/(P₀ log P₀).
- RAT-1 and RAT-0 hold.
- RAT-2 holds modulo issue 1.
- Not re-derived: TYPE's exhaustive casework (I trust it plus the tests) and PUB's Hölder line argument (plausible).

**2. Theorem AB/AB₂.**
- ML, PG, BL, CL, CF, LF-A, EXC-A, PU-A, PD-A and its corollary, EN-A, FAM-A, W-P, §10 and §11 are correct. In LF-A, E τ_p = 1 telescopes and Var τ_p ≍ (p^γ − 1)²/p.
- OPV-A step 6: p₁ ≥ c₁(1 − Ch − e^{−λ} − 𝒱^{1/2})(1 − tiny). This uses only that ν̂ on interior cells is dominated by ν, and it is not circular.
- Step 1's β-drift: |Δβ| ≪ C_X/log X₁.
  - Children: drift × log S₀ = D^{ζ−2α}/λ = D^{−(α+1/2)/2}/λ, whose square D^{−α−1/2} is ≤ D^{−2α}.
  - Parents: D^{4α−2}, whose square D^{8α−4} is ≤ D^{−2α} for α ≤ 2/5.
  - It enters only through e, so it appears squared.
- EXC-A's dyadic BTS covering is correct. λ_min ≤ X₂⁻ ≤ t, so J fits in [t, 2t]. W ≥ 2y, log W = (1 − o(1))log X₂, and the total covering length is ℓ + O(log X₂)λ_min.
- §11's subtraction is correct. The perturbation is |Δcov| ≤ 6D^C/B, taken relative to a₁a₂ ≥ s₁s₂/4 ≥ ρ̂_{B₂}(X_F)²D^{−C}. That gives exp(−aD^αλ + 2c₂D^{3−7α}λ) = L^{−a+o(1)} on the open window. At α = 3/8 it gives L[3/8, −(a − 2c₂)], so a > 2c₂ is exactly what is needed, and it implies (SQ).
- Parents use their own thresholds ϑ₀^p and h_*^p, which is necessary: the child h_* would fail for α ≥ 7/18.

**3. The Cauchy–Schwarz cell device.**
- On an interior cell (R, δ) the count is exactly Σ a₁(m₁)a₂(m₂)1[m₁ ≡ ιm₂ (mod Δ_δ)]. This uses CL(i): both coordinates have the same gcd with Δ′. So the count is A₁A₂/φ = n·d₁d₂·(1 + O(2^ω Δ′/|J|)) plus an error of size at most √(E₁E₂)/φ.
- The singles carry no character error, because every m₁ has |J₂|/Δ_δ + O(1) partners.
- The decomposition p₁₂ − p₁p₂ = Cov_ν̂ + (boundary, exceptional and character terms) is exact.
- No polylog multiplies 𝒱 anywhere in AB. The bound is pointwise relative to A₁A₂, and Σ Rat₁Rat₂A₁A₂ ≤ m(g)². In RAT-1, (V′V″)^{1/2} ≤ 4𝒱₀S′S″. The single exception is issue 1 in RAT-2.
- Every first-order quantity enters only through 𝒱 or ε′², that is, squared: γ log d, γλ, E_FO = E_H + 1/u, E_a and the β-drift. These are all ≍ D^{−α}λ^C for children and D^{1−4α} or D^{7α−3} for parents.
- The terms that enter linearly are all ≪ D^{−2α}: h^{1/2} = D^{−5}, h^{1/4}, e^{−λ} = (Dλ)^{−1}, P₀^{−2} = D^{−7/5}, κ₀^{−1/2}, T₀^{−1}, D^C/B, and the TV terms N^{−Ω(η₀)}. The exceptional bounds are in fact super-polynomially small.

**4. Drappeau Lemme 5 and Theorem WA.**
- Hypotheses verified. Children: y = B, x ≥ X⁻/S₀ with log x ≍ D^{2α}λ, Q = 2S₂² = x^{o(1)}. Parents: x = X_F/q with log x ≍ D^{2−3α}λ, Q = 2V₁²/q with log Q ≍ D^{1−α}λ = o(log x) since α < 1/2.
- The library {(ξ, δ): ξ ∈ universal grid 𝒳, δ ≤ S₀} is independent of the moduli. Block endpoints are snapped to 𝒳, which is essential.
- Multiplicity: MULT gives φ(𝔨)·#pairs/n² ≤ 24/δ for children and ≤ 80d·log V₁/δ for parents.
- Children: 𝔇 = Ξ^C h_*^{−3}D^C = e^{O(D^θλ²)} with θ = 2α − 1/2 < α, so it is beaten by e^{−cD^α/λ²}.
- Parents: 𝔇_p = e^{O(D^{θ_p}λ²)}. On the window, θ_p = max(1 − 2α, 2 − 5α) < 3 − 7α (at 3/8: 1/4 < 3/8; equality only at 2/5). So it is beaten by e^{−cD^{3−7α}/λ²}.
- Circularity is resolved: C_T is chosen after C_Ξ. ϑ₀ never enters Ξ except through ϑ₀^γ = e^{o(1)}.

**5. Ranges.**
- Hildebrand (y ≥ exp((log log x)^{5/3+ε}), relative length ≥ y^{−5/12}, error log(u+1)/log y):
  - Children: log B = aD^αλ, and h_* ≥ B^{−5/12} because θ < α.
  - Parents: log B₂ = b₂D^{4α−1}λ (D^{1/2}λ at 3/8), and h_*^p ≥ B₂^{−5/12} because θ_p < 4α − 1.
- HT Thm 3 (c ≤ y, error 1/u + log y/y): iterated ≤ C_X + 1 times, with drift controlled by φ₂ ≍ log x log y.
- dBT 2.4(i): confirmed verbatim as Drappeau Lemme 1(i), Ψ(x/d, y) ≪ d^{−α}Ψ(x, y) for 1 ≤ d ≤ x.
- FO's Rankin tails at B₂ need 1 − α − 0.343 < 4α − 1, which holds (α > 0.332).
- ω(Δ′) = o(P₀) needs α > 0.3, which holds.

## Probability

I put **0.6** on the analytic core (RAT + AB₂) being correct as a proof once fix 1 is made, taking Fact D as given. Most of the remaining doubt is in internal lemmas I did not re-derive (Cor TP-q, DW(a)⁺, TYPE's casework, PUB) and in unwritten details such as per-point CLB. If you also discount for Fact D resting on an unproved, unrefereed source, I would put it at about 0.5.
