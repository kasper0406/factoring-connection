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

---

## Opus critic 2: the assembly and the parameters (DC, U′, multidescent, cost ledger, window, a > 2c₂, gaps, c(α))

Agent id ad4f5b675b9dbeb75. Run time about 19 min. Code: `src/audit2L/critic2_params.py`, written from scratch; run it
with `.venv/bin/python -m src.audit2L.critic2_params`.

**Verdict (my scope only, the assembly and the parameters): (i) α ∈ (3/8, 2/5): ACCEPT WITH FIXES. (ii) α = 3/8: ACCEPT WITH FIXES.** I found no fatal or exponent-changing gap in the reduction chain or the cost ledger. The one MAJOR item is a write-up gap: Hyp U′ is carried over to design DW by assertion, not by proof. My leading-order audit of that transfer passes.

I formed this verdict without reading any Reviews, Verdict or probability sections. One disclosure: a `grep -n polylog notes/c12rat_results.md` printed a single line from c12rat's Reviews section (line 945). I did not use it.

## Issue table

| # | Sev | Location | Issue | Changes exponent or constant? | Fix |
|---|---|---|---|---|---|
| 1 | MAJOR (write-up) | c12alg §2.7 ("only these replacements", "range restatement … not re-proved line by line"); thmAB §13, the 3/8 bullet ("a leading-order re-check, not a line-by-line re-proof of uprime §3 at the new η₀") | Hyp U′ (μ₁, E N₁², hence γ₂) is the base of the whole chain. Its proofs were written for design 2L′: d ≍ D^{1/4}, η₀ = D^{−1/2}, height floor h = μ + 4ν + η₀, and w coprime to t. For DW they are transferred by assertion: for α ∈ (11/28, 2/5), where u₂ < D^{1/4} and U1₂/DS-U′/CS′ are out of their stated range, and at 3/8 with η₀ = D^{−3/8}. I re-checked every η₀- and range-sensitive inequality I could identify and all hold strictly (see item 1). | Neither, if true | Write U′-DW(α) as a lemma that lists each changed inequality. |
| 2 | MINOR | §13: c(α), the "degenerate" paragraph, the 5.39 illustration | In the open window c(α) carries no information. For a rational α′ ∈ (3/8, α), L[α′, c′] = o(L[α, ε]) for every ε > 0, so the window theorem already gives time L_N[α, o(1)] for every α ∈ (3/8, 2/5). The infimum 0 is a consequence of monotonicity in α, not of C_γ behaviour. Also, the C_γ = 0 illustrations are lower bounds on the formula (C_γ > 0), not achieved values. | Presentation | State "L_N[α, o(1)] for every α > 3/8"; keep an explicit c only at 3/8. |
| 3 | MINOR | thmAB §11, §14 (endpoint constraint O4) | a > 2c₂ is an artifact of the crude absolute subtraction \|(a₁₂−a₁a₂)−(s₁₂−s₁s₂)\| ≤ 6D^C/B measured against s₁s₂ ≍ ρ̂_{B₂}². It is sufficient: relative error exp(−(a−2c₂)D^{3/8}λ), far below the D^{−3/4} budget. It is satisfiable, since nothing else constrains (a, κ, d₀, b₂) at 3/8 now that (A1) is dropped. Every other 3/8 use (FM-U′, Lemma η's denominator, the weight bridge ε_p) needs only (SQ) a > c₂. The text's "iff" should be "if". | Constant only. On my rational grid with C_γ = 0, the best c(3/8) is 4.7062 under a > 2c₂ (exactly the claimed (4/5, 1/8, 1, 2)) versus 4.5625 under a > c₂. | Prove the mixed bound Pr(bad₁ ∧ smooth₂) ≪ D^C B^{−1}s₂, or keep the constraint. |
| 4 | MINOR | §13 step 5; multidescent (vi) | The (a+s)/(1−2^{−α}) term assumes the shape test of DESCEND_ONCE calls FACTOR_2L itself. "FACTOR_U with … replaced", read literally, keeps FACTOR_U's self-call, which would cost L_{√N}[2/5] ≫ L[α]. The calls on parent cofactors may use FACTOR_U (this needs α > 4/11). | Exponent, if misread | Say explicitly which routine each call uses. |
| 5 | MINOR | Theorem statement | Says "every fixed α", but the design (c12alg §1.1) is defined for rational α. | None (item 2 covers irrational α) | Say "rational", or cite monotonicity. |
| 6 | MINOR | Lemma η | "Hildebrand's error ≪ 1/u₂ because u₂ log u₂ ≪ log B₂ (α > 4/11)" is more than needed. With slack factor 3, an o(1) relative error suffices. The real requirements are: R(I) = u(1+O(1/log u)) + O(λ), so u₂ ≫ λ (true, since 3−7α ≥ 1/5), and the small-|n| part d·D^{−3}·D^{2−4α} = D^{−6α} → 0. | None | Optional. |
| 7 | MINOR | uprime §4.1 (CS′) | The certified quadrature that computes R_{B₂} (needed by the algorithm's S₁ scan) is "routine, not written out". | None | Write it out. |
| 8 | MINOR | uprime §6 vs DW | uprime's σ₂ = (3/(8a))(3/d₀ + d₀b₂) and "c ≤ 4.23" belong to 2L′. Under DW, h ≈ 2μ + 2ν + 2η₀ adds 1/d₀ to the child envelope, giving σ₂ = (α/a)(4/d₀ + d₀b₂). §13's 4.71 correctly uses 4/d₀. | Constant | Do not cite 4.23 for 2L-final(3/8). |
| — | Outside my scope (flag only) | c12rat RAT-2, hard-filter step: "(m^par)²·polylog·(δ_e + 𝒱₀)", where "polylog" is defined as e^{O(λ)} = (log N)^{O(1)} | The result is then reported with only a λ^C loss. A genuine (log N)^C = D^Cλ^C factor would wipe out B's budget at every α in [3/8, 2/5). Lemma SPR's own proof gives an active-point spread of (64d)^{γ+γ₂} = 1 + o(1), so this is almost certainly a slip (it should read O(1)). It matters because the assembly tolerates no (log N)^{O(1)} loss in A or B; this is why Lemma η is loss-free. | Exponent, if real | Fix the wording in c12rat. |

## Single weakest link in my scope
Hyp U′ under DW, i.e. issue 1. It is the only input to γ₂ that is transferred by assertion rather than by proof, at two places: the window (11/28, 2/5) and the 3/8 endpoint at η₀ = D^{−3/8}. Structurally, the most fragile point is Prop PB's margin, D^{−3/8} at 3/8 (D^{3−8α} in the window with the first-order AB forms). It has no room for any (log N)^{O(1)} loss anywhere in A, B or W/μ₁.

## 1. Reduction chain

- **Lemma 0 / Prop DC / DC_e / Cor DC.** Correct.
  - Lemma 0: a primitive P′ ∈ Λ_𝔮 has q ∤ b′, so F(P′) ≡ b′^dF(r,1) mod q. Both Good and child-hood force q ∥ F.
  - Prop DC: E[Z1{N=0}] ≤ zp + sd(Z)√p by Cauchy–Schwarz, and p ≤ ϕ by Chebyshev.
  - The second Cauchy–Schwarz over (𝔮, e) gives W(A + √(AB)).
  - Cells are functions of g only, so DC_e's hypotheses hold. W = Σ_{𝔮,e} w_{𝔮,e} is the same W as in Lemma η.
- **Prop UB.** E N₂ ≥ μ₁ − dead mass. N₂ ≤ N₁ gives E N₂² ≤ E N₁². With μ₁ ≥ (log N)^{−C}, both branches (μ₁ ≥ 1 and μ₁ < 1) give γ₂ ≥ (log N)^{−C_γ}. Matches.
- **Lemma VD.** Law of total variance in g. E[N | g, R] = m(g) exactly, and the diagonal Σ Rat·A(1−A) ≤ m(g) gives the 1/m term. Parents likewise. The P₁ = P₂ ⇒ ℓ₁ = ℓ₂ reduction is fine.
- **Root-event corrections.** ϕ_full = ϕ_good/t + (1−t)/t is correct (I checked the mixture identity), with t = 1 − O(1/q).
- **One-point terms.**
  - w/z = Pr(R)/#J is exact. Σ_r Pr(R_{(q,r)}) ≤ d + 2.
  - The parent one-point term is ≤ B₂(log N)^C/#J = exp(b₂D^{4α−1}λ − D^{1−α}λ/(2κ)) → 0. This needs 4α−1 < 1−α, i.e. α < 2/5 (a tie at 2/5).
  - The child one-point term is ≤ e^{−λ²}(log N)^C via Prop CLB.
  - CLB's conservative envelope V₂ ≈ S₂q covers skewed labels; its fibre condition needs V₂ ≪ M, i.e. 4α−1 < 2α.
- **Lemma η.**
  - The ratio-of-sums reduction is valid: the rational part of Acc is a function of g.
  - The mixture-over-intervals ratio is ≤ the maximum over components. The denominator subtraction needs (SQ).
  - Dickman sum: ∫_{v₀}^1 e^{ξv}dv/v = e^ξ/ξ(1+O(1/ξ)) + O(λ), with v₀ = log B/log B₂ → 0 since α > 1/3. This gives W/μ₁ ≤ 3u₂ with no (log N)^C loss.
- **Prop PB.** With the AB forms: A/η = D^{3−8α}λ^C and √(AB)/η = D^{(3−8α)/2}λ^C, both → 0 iff α > 3/8. With the AB₂ forms: max(D^{3−9α}, D^{−α})λ^C → 0 on [1/3, 2/5). At 3/8 both equal D^{−3/8}λ^C. Correct.
- **Interfaces.**
  - The AB₂(i)/(iii) and RAT thresholds are unified at T₀, with the same θ = 2α − 1/2.
  - The RAT-to-true weight bridge holds: (1+5ε)ϕ̃ + 5ε, with ε = N^{−η₀/3} for children. For parents ε_p = L^{−(a−c₂)} at 3/8, which needs only (SQ).
  - RAT is divided by W ≥ (1 − L^{−c})μ₁ (DOM(iii), since u_B ≍ D^{2−4α} ≫ u₂).
- **Field rejections (reject-DW).** e^{−Θ(D^{4α−1})} + 2^{2−d} = e^{−Θ(D^{1−2α})} ≪ (log N)^{−C}. Rejection is applied globally, so select-U goes through (𝒢* has measure ≥ γ₂/4).
- **Split acceptance.** (γ₂/(16 log N))^{2k}/(32k²) = L[α, −2κ(1+C_γ)]. (P6) holds: log Y ≍ D^{1−α}λ/(2κ) ≫ ω(N) ≤ D^{1−α}/a.
- **Rows / POSTPROCESS-DW / LV6.** e^{O(d log d)} = e^{O(D^{1−2α}λ)} = L^{o(1)}. d(μ + h) = O(1), since dν → 0 and dη₀ = D^{α−1/2} → 0. ‖f‖ ≤ 3MH < N.
- **Fruitfulness and recursion.** Fruitfulness rests on the deterministic search plus Lemma OS. The parent-value recursion costs L_{X_F}[2/5] = exp(O(D^{(4−6α)/5}λ)), which is L^{o(1)} iff α > 4/11.
- **Assumed rather than proved.** Only the items in rows 1 and 7 of the table.

## 2. Cost ledger
From the design definitions alone I re-derived:

- log M ~ D^{2α}λ/d₀
- log H ~ 2D^{2α}λ/d₀
- log X_F ~ (d₀/(4κ))D^{2−3α}λ, from the d·log V₁ term; this needs 2−3α > 2α, i.e. α < 2/5
- log X_r ~ D^{2α}λ/d₀
- log X_{F2} ~ (3/d₀ + d₀b₂)D^{2α}λ and log X_{r2} ~ D^{2α}λ/d₀ (conservative box)
- u₂ log u₂ = (d₀(3−7α)/(4κb₂))D^{3−7α}λ, which equals c₂ = 3d₀/(32κb₂) at 3/8

These give σ₁ = α/(ad₀) (+ c₂ at 3/8) and σ₂ = (α/a)(4/d₀ + d₀b₂); sympy confirms both match the claimed forms exactly. The cost terms follow from lv5unif (10.1):

- M_rows = L^a
- 1/p = L^s, with s = 2κ(1 + C_γ)
- T_SQ = W₁(1 + m̄W₂)·L^{o(1)} = L^{max(a, σ₁+σ₂)}
- T_post = L^{2a}

Lemma F_α then gives c(α) as claimed.

- **Illustrative values.** Reproduced exactly. α = 7/18 gives c = 2909/540 ≈ 5.3870, with the special-q term binding; the recursion term is 3.5974. α = 3/8 gives c = 753/160 ≈ 4.7063, with the recursion term at 4.5873 and a > 2c₂ satisfied (0.8 > 0.75).
- **Finiteness and C_γ.** c is finite for every admissible choice. Leaving C_γ unquantified is legitimate for an "L_N[α, O(1)]" claim; it is moot in the window (row 2). At 3/8, c(3/8) is a real but unknown constant, and 4.71 is only its C_γ = 0 floor.

## 3. Window
I computed every leading-order comparison in the excerpt as an exact linear inequality; there are 61, all strict on the open window.

- **Lower end, 3/8.** Two ties occur there:
  - 3 − 7α = α. This is (SQ), the level-1 algebraic B₂-smoothness cost entering σ₁, and the §11 subtraction. It is the j = 2 size-model wall (2 − 4α ≤ ζ₁ = 4α − 1).
  - The first-order PB budget (8α − 3). It is not binding, because AB₂ moves PB's threshold down to 1/3.
- **Upper end, 2/5.** Ties for: θ_p vs 3−7α (gap 2−5α), B₂ vs V₁, the one-point parent term, d log V₁ vs log MH, the μ_opt branch, and μ − (d−2)ν < 0. Beyond 2/5 Theorem U applies instead.
- **The two named gaps.** The child gap is ½ − α > 0 and the parent gap is 2 − 5α > 0, with the other parent loss component 2−5α below 3−7α by 1−2α.
- **Other comparisons, all strict on [3/8, 2/5):**
  - h_* ≥ B^{−5/12} and h_*^p ≥ B₂^{−5/12}
  - ζ, ζ_p and the Rankin tail
  - ω(Δ′) = o(D^{7/10}), using 1−α < 7/10 and 4α−1 < 7/10
  - Δ′ = x^{o(1)}
  - N^{η₀} ≫ L[α] and ≫ B₂, using 3α − ½ > max(α, 4α − 1)
  - u₂ log u₂ ≪ log B₂ and the recursion, both needing α > 4/11
  - LV6, and the rejection probabilities vs γ₂
- **The range restatement** is genuinely needed only for α > 11/28.

**Command:** `cd /Volumes/git/factoring-connection && .venv/bin/python -m src.audit2L.critic2_params`

**Key output:**
```
match claimed sigma1, sigma2: True True        (alpha = 3/8, 7/18, 79/200)
u2 log u2 = 3*d0/(32*b2*kappa) * D^3/8 lambda  (alpha = 3/8)
open-window failures (should be none): []
3/8 TIE: PB-with-AB (A/eta, sqrt(AB)/eta), (SQ)/subtraction/sigma1 ; everything else ok
2/5 TIE: theta_p vs 3-7a, B2<<V1, one-point parent, d log V1 vs log MH, mu_opt, h_* branch
alpha=7/18 (3/5,1/8,1,2): c = 2909/540 = 5.38704 ; alpha=3/8 (4/5,1/8,1,2): c = 753/160 = 4.70625
1/(1-2^(-3/8)) = 4.368823 ; 3/8 grid C_gamma=0: a>2c2 -> 4.7062 ; a>c2 -> 4.5625
degenerate sequence at 7/18 (C_gamma=0): c = 5.29, 2.65, 1.32, 0.66, 0.33
```

## Probability
Given the analytic core (RAT, AB/AB₂ as stated) and Theorem U's chain, I put the probability that the assembly is correct at **0.80**. The main residual risks are:

- the unwritten U′-DW transfer: about 0.08
- an interface mismatch not visible in the excerpt: about 0.07
- the DW design lemmas (CT, IR-DW, DW(a)⁺) used at level 1: about 0.05

The file I created is `/Volumes/git/factoring-connection/src/audit2L/critic2_params.py`; nothing else was edited or committed.
