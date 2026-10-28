# audit2L: is Theorem 2L-final (rigorous Las Vegas factoring in L_N[α], 3/8 ≤ α < 2/5) true?

Agent `audit2L`, phase 31, 2026-10-28. I audited Theorem 2L-final(α) for α ∈ (3/8, 2/5) and Theorem 2L-final(3/8)
(notes/thmAB_results.md §13). The chain behind them:
- multidescent (Prop UB, Theorem 2L);
- uprime (U′);
- pc3 (Prop DC/DC_e);
- c12rat (Theorem RAT, v2.1);
- c12alg (design DW, Lemma TP, Lemma η, Prop PB, Prop CLB);
- thmAB (Theorem AB/AB₂, v3).

All of it rests on Theorem U's one-level chain (paper2/main.tex, lv5unif). I followed auditU's pattern.

Files:
- This note.
- `notes/audit2L_blind_claims.md`: the blind excerpt every reviewer received. It holds the statements and proofs of
  the chain, with review sections, verdicts, probabilities and attributions removed. It was built by
  `src/audit2L/blind/make_blind2.py` and `patch_blind.py` (copied from scratch; they read paths under `scratchpad/audit2L/`), plus two hand edits.
- `notes/audit2L_own_verdict.md`: my verdict, committed (6fdec2b) before I read any reviewer output.
- Reviewer records:
  - `notes/audit2L_codex1_report.md`: the GPT-6.1-sol blind referee (prompt `codex/audit2L_1.md`, output
    `codex/audit2L_1.out`).
  - `notes/audit2L_codex2_report.md`: the GPT-6.1-sol follow-up on my repair of the RAT-2 gap (`codex/audit2L_2.*`).
  - `notes/audit2L_critics.md`: Opus critic 1 (analytic core) and critic 2 (assembly and parameters), verbatim.
- Code, all exact and written from scratch:
  - `src/audit2L/exponents.py`: mine, 46 exponent comparisons and c(α). Run
    `.venv/bin/python -m src.audit2L.exponents experiments/raw/audit2L_exponents.json`.
  - `src/audit2L/critic2_params.py`: critic 2's, 61 comparisons and the cost ledger. Run
    `.venv/bin/python -m src.audit2L.critic2_params`; the output is in
    `experiments/raw/audit2L_critic2_params.txt`.
  - `src/audit2L/codex_scratch/referee_checks_audit2L_1.py`: codex's finite counterexamples, copied from its scratch
    directory and re-run by me.

## Verdict

**(i) α ∈ (3/8, 2/5): ACCEPT WITH FIXES. (ii) α = 3/8 (with a > 2c₂): ACCEPT WITH FIXES, at lower confidence.**
Both verdicts are relative to Theorem U's internal chain and to Drappeau's Lemme 5. The fixes include one MAJOR gap
in the written proof, which this audit repairs (§2).

In plain terms:
- **As written, the proof had a gap that would have been fatal.** RAT-2 (c12rat §8.3, the parent hard-filter step)
  bounds the parent rational covariance by (m^par)²·**polylog**·(δ_e + 𝒱₀). In these notes polylog means
  e^{O(λ)} = (log N)^{O(1)}. The proof then reports only a λ^C loss.
  - One genuine factor of log N there gives √(AB)/η ≫ D^{1/2−α} → ∞, so no α < 2/5 would follow (codex's
    computation; I confirmed it).
  - All four reviewers found this independently: me (pre-registered, rated MINOR), Opus critic 1 (MAJOR, fixable),
    Opus critic 2 (flagged out of scope) and codex (MAJOR, hence MAJOR REVISION).
  - Codex added a second point: SPR's comparator S_typ is built on ρ₂(X_F), while the parent weights live at the
    root-conditioned scale X_F/q. These differ by q^{1−β₂} ≍ u₂ log u₂, a power of D.
- **The gap is repairable, and the repair is written and checked.**
  - Lemma SPR-P′ (§2) compares the OPV constants of the parents **within one label**, at one q, and finds they agree
    up to an absolute constant.
  - With it, the RAT-2 step costs only C(δ_e + 𝒱₀^par), with no power of D.
  - Codex's focused follow-up (`audit2L_2`) rated the repair CORRECT WITH FIXES. Its three corrections are applied
    in §2.
  - With the repair, codex's verdict moves to ACCEPT WITH FIXES for both (i) and (ii). So all four reviewers now
    stand at ACCEPT WITH FIXES.
- **Everything else holds up.**
  - The cell-covariance device, the character-error accounting through Drappeau's Lemme 5, thin pairs and skewed
    labels, the parent thresholds, the endpoint constraint, and the reduction chain UB → DC → VD → η → PB all pass.
  - The exponent bookkeeping passes in two independent exact scripts (46 and 61 comparisons).
  - The four thmAB round-3 repairs are correct. The c12rat v2 repairs are correct, apart from the RAT-2 step
    above and two MINOR missing cases found by codex (§1.2).
- **The window and the endpoint.**
  - The window is open at both ends.
  - At 3/8 the first-order budget ties (8α − 3 = 0); second order (AB₂ + RAT-2) clears it with margin D^{3/8}.
  - The level-1 B₂-smoothness and acceptability subtraction also tie at 3/8 (3 − 7α = α). a > 2c₂ closes that tie.
    It is sufficient but probably not necessary: (SQ) a > c₂ would do if the mixed bound were proved.
  - At 2/5 the parent loss θ_p meets the saving 3 − 7α (gap 2 − 5α → 0), and B₂ meets V₁.
  - Gaps: children ½ − α (1/8 at 3/8), parents 2 − 5α (1/8 at 3/8).
- **What the claim still rests on.**
  - Theorem U's chain: [pp]'s unrefereed second moment.
  - Drappeau's Lemme 5. I checked it verbatim in arXiv 1307.7554v3, and so did critic 1 and codex. It is stated
    without proof, as following from §3.3 of Harper's unrefereed preprint 1208.5992.
  - About 4000 lines of internal, unrefereed notes that only AI has reviewed.
  - One assertion not yet written as a proof: the transfer of Hyp U′ to design DW, in the window and at 3/8 with
    η₀ = D^{−3/8} (issue 2). Everyone who checked it at leading order found it holds.
- **Why (ii) gets lower confidence.** At 3/8 every block must be genuinely second order. The acceptability tie is
  closed only by a constant constraint. And the U′-DW transfer at η₀ = D^{−3/8} has only a leading-order re-check.
- **Is c(α) meaningful?** The constant c(α) is finite for every admissible parameter choice but not quantified,
  since C_γ is unknown. In the open window it carries no information: the window theorem already gives
  L_N[α, o(1)] for every α ∈ (3/8, 2/5), by monotonicity in α (critic 2). Only c(3/8) is a meaningful constant;
  4.706 is its C_γ = 0 floor.

**No weaker α is needed.** No issue is FATAL once the §2 repair is adopted. Without it, nothing below 2/5 follows,
because the defect is a power-of-D loss in the parent budget at every α in the window. Theorem U's L[2/5] is unaffected.

### Reviewer verdicts (all blind; my own verdict was committed before I read any reviewer)

| Reviewer | Scope | (i) window | (ii) 3/8 | Weakest link | Prob. |
|---|---|---|---|---|---|
| audit2L (me) | whole chain | ACCEPT WITH FIXES | ACCEPT WITH FIXES | U′-DW glue; RAT-2 polylog (rated MINOR before review) | 70% / 55% |
| GPT-6.1-sol `audit2L_1` (xhigh) | whole chain, Theorem U given | MAJOR REVISION | MAJOR REVISION | RAT-2 parent normalization | — |
| GPT-6.1-sol `audit2L_2` (xhigh) | the RAT-2 repair | ACCEPT WITH FIXES | ACCEPT WITH FIXES | (repair CORRECT WITH FIXES) | — |
| Opus critic 1 | analytic core (RAT, AB₂, Lemme 5/WA, ranges) | ACCEPT WITH FIXES | ACCEPT WITH FIXES | RAT-2 hard-filter π_P-centring | 60% (≈50% with Lemme 5 discounted) |
| Opus critic 2 | assembly, parameters, cost | ACCEPT WITH FIXES | ACCEPT WITH FIXES | U′ under DW (transfer by assertion) | 80% (assembly, given the core) |

Reconciliation. Codex's MAJOR REVISION and critic 1's MAJOR-fixable are the same finding, and I agree it is a real
gap in the written proof. I had seen it but under-rated it as MINOR. I now rate it MAJOR as written, because the
literal proof does not give the theorem. With Lemma SPR-P′ the gap is closed at the claimed rate. That lemma uses only
ingredients already proved in the chain (OPV, ADM/DW(a)⁺, HT3, the SPR lattice count), and its one re-review accepted
it. So the overall verdict is ACCEPT WITH FIXES, with F1 as the mandatory fix.

My updated probabilities, correct relative to Theorem U's chain and Lemme 5:
- (i): 70%.
- (ii): 55%.
Unconditionally, the internal chain and Harper/Lemme 5 must also be discounted.

---

## Issue log

"Exp/const" says whether the issue changes the exponent α or only the constant c. "Status" is its state after this
audit. Fixes recorded here are for the notes' owners to apply; I edited no other agent's note.

| # | Sev. | Location | Issue | Exp/const | Status | Found by |
|---|---|---|---|---|---|---|
| 1 | **MAJOR** (as written) | c12rat §8.3 RAT-2 hard-filter bullet; §5.7 Lemma SPR (parents) | Σ_{P,P′}c_Pc_{P′}r_P ≤ (m^par)²·**polylog**·(δ_e + 𝒱₀), with polylog = (log N)^{O(1)}, is reported as λ^C. The comparator S_typ uses ρ₂(X_F), but the weights live at X_F/q (a factor q^{1−β₂} ≍ u₂ log u₂). | **Exponent if unrepaired** (an extra log N gives √(AB)/η ≫ D^{1/2−α}) | **Repaired here** (§2, Lemma SPR-P′; codex audit2L_2: CORRECT WITH FIXES, fixes applied) | me (as MINOR), critic 1, critic 2, codex |
| 2 | MAJOR (write-up) | c12alg §2.7 ("only these replacements", range restatement "not re-proved line by line"); thmAB §13 3/8 bullet (leading-order re-check at η₀ = D^{−3/8}) | Hyp U′ (μ₁, E N₁², hence γ₂) is transferred from design 2L′ (d ≍ D^{1/4}, η₀ = D^{−1/2}, h = μ + 4ν + η₀) to DW by assertion. This affects α ∈ (11/28, 2/5), where u₂ < D^{1/4}, and the endpoint 3/8. | Neither, if true. Every leading-order check passes (me, critic 2, codex). | Open: write U′-DW(α) as a lemma listing each changed inequality | critic 2, me |
| 3 | MINOR | c12rat Lemma EXC(c) | s_i includes mid-prime factors; the Rankin argument covers only a large P₀-smooth part. Codex example: P₀=5, S₀=7, y=17, q=19, P=(1,1),(1,210), g=(1,23), s=(22,11). I re-ran it. | No | Fix: add the mid-prime Euler product ∏(1 − p^{−1+γ+σ})^{−1} = 1 + o(1) | codex |
| 4 | MINOR | c12rat Lemma PRIM (parents) | "f′f″ = 0 unless p ≤ y" fails for a common content prime ℓ ∈ J not dividing K. Codex example: y=7, J=(20,40], q=13, P=(4,5),(3,7), g=(23,23). I re-ran it. | No (≤ Ξ_pD^C log Y/Y) | Fix: add the J-content term | codex |
| 5 | MINOR | thmAB §13 step 1 (parent weight bridge) | Cites W-P(P3), which concerns active points and an aggregate, for a pointwise lower bound a_P^{id} ≥ ρ_a(q)D^{−C}. | No | Fix: cite the DW(a)⁺ pointwise lower bound for every filtered P, then subtract | codex |
| 6 | MINOR | thmAB §10(iii), FAM-A step 3: "Prop CLB per point" | c12alg Prop CLB is written with a global envelope. The per-point bound m ≥ Σϖ_P(log N)^{−C} is asserted, not written out. (Critic 1 also called "\|P′\| ≤ V₂" false for skewed labels. That part is wrong: V₂ = 4Q₂S₂ ≥ 2qS₂ ≥ \|P′\|, as critic 2 notes.) | No | Fix: write CLB per point with X_{F2}(P′) | critic 1 |
| 7 | MINOR | thmAB EXC-A, FO tails | The Rankin step uses P₀^{−1/2+γ} as if p^γ were bounded for all p ≤ y. | No | Fix: split at p = e^{1/γ} | critic 1 |
| 8 | MINOR | thmAB FAM-A step 1 | Cites "OPV-A step 6" for the pointwise bound; it is step 5. | No | Fix the reference | critic 1 |
| 9 | MINOR | thmAB §10 vs c12rat TF | The Euclidean κ and AB's ∞-norm ϑ differ by a factor ≤ 2. | Constant | Note it | critic 1 |
| 10 | MINOR | thmAB §13, the c(α) paragraph | In the open window c(α) carries no information: monotonicity in α already gives L_N[α, o(1)]. The C_γ = 0 values 5.39 and 4.71 are floors, not achieved values. | Presentation | Fix: say so; keep c only at 3/8 | critic 2 |
| 11 | MINOR | thmAB §11, §14 | a > 2c₂ is sufficient, not necessary: the text's "iff" should be "if". The other 3/8 uses need only (SQ) a > c₂. Best grid c(3/8) at C_γ = 0: 4.7062 under a > 2c₂, 4.5625 under a > c₂. | Constant | Fix the wording. Optional: the mixed bound | critic 2 |
| 12 | MINOR | thmAB §13 step 5; multidescent (vi) | Name the routine for each recursive call: FACTOR_2L self-call on the halves (giving (a+s)/(1−2^{−α})) and FACTOR_U on parent cofactors (needs α > 4/11). | Exponent if misread | Fix: state it | critic 2 |
| 13 | MINOR | thmAB §13 statement | The design is defined for rational α. | None (by monotonicity) | Fix: say "rational", or cite monotonicity | critic 2, me |
| 14 | MINOR | uprime §4.1 (CS′) | Certified quadrature for R_{B₂} is not written out. | No | Write it out | critic 2 |
| 15 | MINOR | uprime §6 | Its σ₂ = (3/(8a))(3/d₀ + d₀b₂) and "c ≤ 4.23" belong to 2L′. Under DW, σ₂ = (α/a)(4/d₀ + d₀b₂), which thmAB uses. | Constant | Do not cite 4.23 for 2L-final(3/8) | critic 2 |
| 16 | MINOR | c12rat §1.3 vs SPR | Ξ is defined as e^{O(D^θλ)}, but SPR needs e^{O(D^θλ²)}. Harmless, since T₀ = e^{C_TD^θλ²}. | No | Unify | me |
| 17 | MINOR | c12rat Lemma BTS uses | An interval crossing 0 needs a split, or an extra Ψ(W, y) term. | No | Fix: add the split | codex |
| 18 | MINOR | thmAB LF-A | "16γ log Δ′ = O(λ + D^{2−5α}λ)" is loose. Harmless. | No | — | critic 1 |
| 19 | scope | Fact D | Drappeau's Lemme 5 is the single essential published input; its proof rests on Harper 1208.5992 §3.3 (unrefereed). Checked verbatim by me, critic 1 and codex: 1/φ(q) weight, all χ mod q with 1 < cond ≤ x^η, Q ≤ x, (log x)^c ≤ y, effective for A < 1. | — | Not repairable here; phase-31 agent `mf2L` is trying to remove it | all |

No issue is FATAL after fix F1 (= issue 1, §2). Issues 3–18 are text or proof-detail corrections for the owners.
Issue 2 is the main unwritten piece.

---

## 1. Targeted re-review (task 1)

### 1.1 thmAB round-3 repairs (v3; never re-reviewed before this audit)

| Repair | Verdict | Evidence |
|---|---|---|
| **Marginal lower bound** (OPV-A step 6: p₁ ≥ c₁(1 − Ch − e^{−λ} − 𝒱^{1/2})) | **Correct** | On interior cells, CF(ii) gives the exact single count n(R,δ)d₁(1 + tiny). So p₁ ≥ (1 − tiny)E_ν̂Ψ₁ by positivity. Then E_ν̂[Ψ₁/c₁] ≥ ν̂(int) − ν̂(int)^{1/2}E_ν̂[(Ψ₁/c₁ − 1)²1_int]^{1/2}. The second factor is ≤ 𝒱, because ν̂ restricted to interior cells is dominated by the law of (I(x₁), δ(x₁)) under ν on non-exceptional points: the mass of an interior cell is the ν-mass of its points. 𝒱 itself (steps 1–4) uses PG(d) and PG(e) and never uses p₁, so there is no loop. Critic 1 and codex agree. |
| **Parent RAT bridge** (value-only a_P^{id}; ε_p ≤ D^C/(Bρ_a(q)) + N^{−η₀/3}) | **Correct**, with a citation fix (issue 5) | We have z_true(g) = z_ideal(g)(1 + O(ε_p)) pointwise, since every filtered P has a^{id} ≥ ρ_a(q)D^{−C} by DW(a)⁺. Perturbing a relative variance gives σ_true ≤ (1 + O(ε))σ_ideal + O(ε). At 3/8, ε_p = L[3/8, −(a − c₂)] needs only a > c₂. |
| **BTS placement** (EXC-A large classes: cover each dyadic piece by J ⊂ [t, 2t] with \|J\| = max(\|piece\|, λ_min)) | **Correct** | λ_min = X₂y^{−C_X}/ϑ ≤ X₂⁻ ≤ t, so J fits. gcd(a, Δ′)·W = δ\|J\|/Δ′ ≤ \|J\| ≤ t ≤ inf J. W ≥ λ_min/Δ′ ≥ 2y. log W = (1 − o(1))log X₂. The loss (4ϑΔ′y^{C_X})^γ is e^{O(λ)} for non-thin child pairs. The total covering length is ℓ(x₁) + O(log X₂)λ_min. |
| **β-drift** (D^{ζ−2α}/λ children, D^{4α−2}/λ parents) | **Correct** | \|Δβ\| ≪ log(X₁/\|ξ\|)/(log X₁ log y) ≤ C_X/log X₁, because φ₂ ≍ log x log y. Times log S₀ = D^ζ this gives D^{ζ−2α}/λ. It enters only through e, so it appears squared: D^{−α−1/2} ≤ D^{−2α}, and for parents D^{8α−4} ≤ D^{−2α} for α ≤ 2/5. |

### 1.2 c12rat v2 repairs (never re-reviewed before this audit)

| Item | Verdict | Evidence |
|---|---|---|
| **Prop DC / RAT bridge** | **Correct** | DC_e conditions on R_𝔮 ∩ e. RAT works on 𝒢_e (q ∤ G_g(r,1)). On the complement both counts vanish, and ϕ_full = ϕ_good/t + (1 − t)/t with t = 1 − O(1/q). I checked the mixture identity, and so did critic 2. The idealised child weight differs from the true A_{P′}(g) by N^{−η₀/3} relatively, which costs (1 + N^{−η₀/4})ϕ + N^{−η₀/4}. W ≥ (1 − L^{−c})μ₁ (DOM(iii)). |
| **Lemma TYPE** | **Correct** | The valuation casework is consistent. Admissibility (b) excludes p \| g. The extra p in M_τ encodes the primitivity non-congruence. Every condition is r₁ ≡/≢ c·r₂, so the set is diagonal-invariant. Critic 1 did not re-derive the casework; it trusts the exhaustive test. |
| **Lemma CP / CPW** | **Correct** | Cauchy–Schwarz on a diagonal-invariant C plus Parseval: \|C\|/φ²·(Σ\|â\|²)^{1/2}(Σ\|b̂\|²)^{1/2}. CPW's positivity argument gives 3(Ẽ_u + Ẽ_v) + 4ε′². |
| **Lemma DOM** | **Correct** | (i): w = #J^{−1}Σ_{ℓ,P∈𝒫∩Λ_𝔮}Pr[e, Acc_ℓ, q child] ≤ Pr[Ω_e]·α̂(q)·ϱ̂·#(𝒫∩Λ_𝔮)/#J, and this bound does not depend on r. (ii): q residues × Σ_q α̂(q) ≈ ρ₂(X_F)·u₂(log N)^C ≤ (log N)^Cμ₁. The TV term Σ_q 7q²N^{−η₀/2} is negligible. (iii): parents with no child need B-smooth F(P), and u_B ≍ D^{2−4α} ≫ u₂. |
| **Lemma FAM / Theorem WA / Cor WA′** | **Correct** | WA is a one-line Hölder step. FAM: for a fixed modulus k, only O(1) window primes divide k (children; ≤ D^{1−2α} for parents). φ(k)·(pairs with that K)/n² ≤ 40Q₂. Since ŵ is r-uniform, Q₂·max_qŴ_q ≤ C log Q₂·Σ_{win}Ŵ_q. So 𝔇 = e^{O(D^θλ²)}Σŵ, which ε_lib = e^{−cD^α/λ²} beats. Critic 1 and codex agree. |
| **Lemma SPR** | children **correct**; parents **as written insufficient** (issue 1) | Children: γ log(S₂T₀) = O(λ), and dγ_a log(S₂T₀) = O(D^{1−2α}λ² + D^{1/2−α}λ³) ≤ D^θλ², with equality of exponents at 3/8 (1 − 2α = θ = 1/4). Parents: see §2. |
| **Thin pairs** (Lemma TF; Prop SK(iii)/(iv)) | **Correct** | Children: κ ≤ 2√3 sS²/Δ and MULT give a thin fraction ≤ 60s/κ₀, each costing (Ξ+1)·4Ξ³m²/n². Parents: ≤ 3/κ₀ + O(q/V₁). Codex prefers a per-b count here (MINOR wording). |
| **Skewed labels** (Lemma LS + DOM; thmAB §10(iii)) | **Correct** | #{r : s > T} ≤ 3q/T + 2, and DOM(i) is r-uniform, so the skewed weight is ≤ (3/T₀ + 2/B)(log N)^Cμ₁. thmAB's per-point bound E[N(N−1)] ≤ Ξ^C(Σϖ_P)² against m ≥ Σϖ_P(log N)^{−C} never compares across the box. That needs per-point CLB (issue 6). |
| Also checked: **PRIM** (Bonferroni with signed weights; (γ²log²P₀ + E²)/(P₀ log P₀)), **OPV** (Z ≥ 1, moments), **FO** (small/mid Rankin split), **ADM**, **PUB** | **Correct**, except EXC(c) mid primes (issue 3) and PRIM parent J-content (issue 4) | Codex's two finite counterexamples are genuine; I re-ran them. Each is a missing case whose contribution is negligible. |

---

## 2. Fix F1: Lemma SPR-P′ (closing the RAT-2 parent normalization)

**The defect.** c12rat §8.3 writes m^par ≥ n·c_typ/polylog and c_P ≤ c_typ·polylog, both against a common comparator
S_typ that is built on ρ₂(X_F). It then concludes Σ_{P,P′}c_Pc_{P′}r_P ≤ (m^par)²·polylog·(δ_e + 𝒱₀). A factor
D^K on 𝒱₀^par ≍ (D^{−2α} + D^{−2(3−7α)})λ⁶ breaks the budget: √(AB)/η ≍ D^{K/2−α}. This is exactly the
precision-class violation the chain cannot afford; Lemma η was rewritten loss-free for the same reason.

**Lemma SPR-P′.** Fix a label 𝔮 = (q, r) and a cell e. Let 𝒫_e(𝔮) := {P ∈ 𝒫 ∩ Λ_𝔮 : π_P > 0}, the parents that
pass the hard filter for some g ∈ Ω_e. Then:
- c_P ≤ C·c_{P′} for all P, P′ ∈ 𝒫_e(𝔮) and all pair moduli K, with C absolute;
- if π_P = 0, then f_P ≡ 0 on the averaging support.

*Proof (with codex audit2L_2's corrections applied).* c_P = a_P(x₁)·ρ̂_J(x₁)·Λ_P(K), where x₁ = max_{Ω_e}|G_g(P)|.
- (i) x₁ ∈ [MV₁/32, 2MV₁], since |G_g(P)| ≤ |g||P| ≤ 2MV₁ and P passes somewhere. So ρ̂_J varies by
  64^{γ}(1 + O(1/u)) = 1 + o(1) (HT3/dBT1 at y = B).
- (ii) Put A_P(t) := t(2H + 1)|P|_∞^{d−1}/q, the root-conditioned scale.
  - DW(a)⁺ (c12alg §2.6) writes the law of F(P)/q as a mixture of translates of one interval of length A_P(t), inside
    [−CdA_P(t), CdA_P(t)].
  - Uniformly over the translates (split at 0 when needed), Hildebrand's short-interval theorem at B₂ gives each
    interval a B₂-smooth density comparable, up to absolute constants, to ρ̂_{B₂}(A_P(t)). The relative length is
    ≳ 1/d ≫ B₂^{−5/12}.
  - So a_P(t) ≍ ρ̂_{B₂}(A_P(t)) with absolute constants, **whatever the shape of the other segments**. This avoids
    ADM's fixed-P comparison and any negative moment of the shift Z_P alone, which can concentrate near 0.
  - All P in the label share **the same q**, and |P|_∞ ∈ ((1 − 1/d)V₁, V₁]. Hence A_P(t)/A_{P′}(t′) ∈ [(64e)^{−1}, 64e],
    and HT3 at B₂ costs (64e)^{γ₂} = 1 + o(1).
- (iii) Λ_P(K) = exp O(γ log P₀ + ω(K)/P₀ + ω(K)/y) = 1 + o(1), uniformly in K.
- If π_P = 0, the filter indicator in f_P vanishes on Ω_e. ∎

**Straddlers and active points.** Fix a prime b in the window. Straddling, i.e. 0 < π_P < 1, confines a to at most two
intervals of total length O(δ_eV₁). Intersected with a ≡ rb (mod q), that leaves O(δ_eV₁/q + 1) points, against
2V₁/q + O(1) parents for this b. Hence #straddlers ≤ C(δ_e + q/V₁)n = O(δ_e n), uniformly in the label, skewed child
lattices included. The same count shows that at least half of 𝒫 ∩ Λ_𝔮 is active (π_P = 1), and active P have
S_P ≥ c_P/2 by OPV.

**Corrected RAT-2 step.**
- m^par ≥ (n/2)·min_{𝒫_e}c_P/2 ≥ n·max_{𝒫_e}c_P/(4C).
- Var_{ν₁}(ψ_P) ≤ 2c_P²(π_P(1−π_P) + 𝒱₀^par).
- Hence Σ_{P≠P′}|Cov(f_P, f_{P′})| ≤ (Σ_{P∈𝒫_e}Var(ψ_P)^{1/2})² + (𝔈 and smaller terms)
  ≤ 4(n max c)²(Cδ_e² + 𝒱₀^par) ≤ C′(m^par)²(δ_e² + 𝒱₀^par) + ….
- So RAT-2 ≤ C(D^{−2α} + D^{−2(3−7α)})λ^C + δ_e², with A₀ ≥ 1 > α. The polylog losses that remain sit only inside 𝔇
  (character errors) and inside the thin/straddling remainders, where exponentially small savings absorb them.
- No common comparator is needed, so the X_F vs X_F/q discrepancy does not arise.
- Codex audit2L_2 also checked AB §§9–11: there the principal covariance term is C𝒱_p*·E[z(g)²]/z² =
  C𝒱_p*(1 + σ^rat), with no normalization loss. Coarse D^C factors there enter only character, thin-pair and
  acceptability errors, whose margins absorb them.

Status: a repair proposed by the auditor and verified once (codex audit2L_2: CORRECT WITH FIXES, corrections applied
above). The owner of c12rat should replace the parent bullet of Lemma SPR and the hard-filter bullet of §8.3 with this
lemma.

---

## 3. Blind excerpts and reviewers (tasks 2–3)

- **The excerpt.** `notes/audit2L_blind_claims.md` has 2738 lines. It contains thmAB §§1–14, c12rat §§1–8 and §9.2,
  pc3 §§1–3 and §6.2, multidescent §§3.1–3.3, uprime §§1–6, c12alg §§1–3.2 and §§4–5, and a neutral header stating
  the claim and its logical structure.
  - Attributions were removed mechanically: balanced parenthetical groups that mention codex, a critic, a round or a
    version tag were dropped. A patch list then handled the rest.
  - Every line containing a percentage was dropped.
  - A residual grep for codex/critic/review/round/% is clean, apart from the header's own instruction.
- **Who saw what.**
  - Codex round 1 and both Opus critics were told not to read the source notes' Reviews/Verdict sections, any audit
    note, or any codex file.
  - Critic 2 disclosed one accidental grep line from c12rat's Reviews; it did not use it. Codex disclosed one
    escapeR2 sentence; it did not use it.
  - The Opus critics ran one at a time. Codex round 1 overlapped critic 1.
  - Codex round 2 was **not** blind. It is a focused verification of my repair and saw its own round-1 report.

## 4. Exact checks (task 3c)

- **Mine** (`src/audit2L/exponents.py`): 46 leading-order comparisons, with every D-exponent exact.
  - All hold on (3/8, 2/5): grid 3/8 + k/400, k = 1…9, and 7/18.
  - At 3/8 exactly three fail, all expected: the two first-order budget checks (A/η and √(AB)/η, which second order
    replaces and which pass in their second-order form), and the tie 3 − 7α = α, which a > 2c₂ closes.
  - At 2/5 there are three ties: 2 − 3α vs 2α, B₂ vs V₁, and θ_p vs 3 − 7α.
  - Gaps: children ½ − α, parents 2 − 5α. At 7/18: 1/9 and 1/18.
  - c(7/18; 3/5, 1/8, 1, 2; C_γ = 0) = 5.38704 and c(3/8; 4/5, 1/8, 1, 2) = 4.70625. Both reproduce thmAB §13.
- **Critic 2** (`src/audit2L/critic2_params.py`, independent): 61 comparisons, all strict on the open window, with the
  same ties at 3/8 and 2/5.
  - It re-derives σ₁ and σ₂ symbolically from the design sizes: c = 2909/540 at 7/18 and 753/160 at 3/8.
  - Best grid c(3/8) at C_γ = 0: 4.7062 under a > 2c₂, 4.5625 under a > c₂.
  - A degenerate parameter sequence at 7/18 drives c → 0, which is consistent with L[α, o(1)].
- **Codex** (`src/audit2L/codex_scratch/referee_checks_audit2L_1.py`): the two finite counterexamples (issues 3 and
  4), and the effect of one extra log N on √(AB)/η: D^{1/8} at 3/8 and D^{1/9} at 7/18. That is why issue 1 is
  exponent-relevant if left unrepaired.

## 5. Bottom line

- Theorem 2L-final(α), α ∈ (3/8, 2/5), and Theorem 2L-final(3/8) are true **if** all of the following hold:
  - Theorem U's chain;
  - Drappeau's Lemme 5;
  - Lemma SPR-P′ (F1), which replaces c12rat's parent normalization;
  - the U′-DW transfer (issue 2), still to be written out.
- No reviewer found a counterexample or an exponent failure in the repaired chain. The analytic novelty, the
  cell-covariance device turning first-order one-point asymptotics into second-order pair covariances, survived four
  independent checks.
- The headline is now as reliable as two things: [pp]'s unrefereed second moment, and roughly 4000 lines of internal
  notes that only AI has reviewed. One of those notes had a precision-class slip, which would have been fatal, that
  two prior review rounds missed.
- Recommended next steps:
  1. The owners apply F1 and issues 3–18.
  2. Write U′-DW(α) out as a lemma (issue 2).
  3. Remove Lemme 5 (agent `mf2L`).
  4. Get an expert human referee for c12rat §§4–8 and thmAB §§5–9.

