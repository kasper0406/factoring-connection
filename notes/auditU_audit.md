# auditU: is Theorem U (rigorous Las Vegas factoring in L_N[2/5, c]) true?

Agent `auditU`, phase 12, 2026-10-20. This note audits `notes/lv5unif_results.md` (Theorem U and the τ = 0 claim). It
rests on the fixed-degree chain: [pp] = `paper2/main.tex`, `notes/thmA_proof.md`, `notes/lem85_results.md` and
`notes/recurse_results.md`.

Files:
- This note: `notes/auditU_audit.md`.
- The blind excerpt given to every reviewer: `notes/auditU_blind_claims.md`. It is lv5unif_results.md with §15 (Reviews),
  the codex-verdict bullet and the reviewer attributions removed.
- Reviewer records:
  - `notes/auditU_codex1_report.md`, the verbatim GPT-6.1-sol report (prompt `codex/auditU_1.md`, raw output `codex/auditU_1.out`);
  - `notes/auditU_critics.md`, the two Opus critics.
- Exact checks: `src/auditU/constants_check.py`, written from scratch; it imports lv5unif only to compare at the end.
  - Command: `.venv/bin/python -m src.auditU.constants_check experiments/raw/auditU_constants.json` (about 95 s).
  - Output: `experiments/raw/auditU_constants.json`.

## Verdict

**ACCEPT WITH FIXES, conditional on the internal fixed-degree chain.** I reached this verdict myself, and so did all
three blind reviewers: GPT-6.1-sol and two Opus critics. **Nobody found a FATAL or MAJOR issue.** Every issue logged
below is MINOR. None changes α = 2/5, the cost formula, or the finiteness of c.

In plain terms:
- **The growing-degree step is sound.** Theorem U is a correct consequence of [pp]'s fixed-degree relation-supply chain
  plus Lemma 8.5′ and the recursion.
  - The feared d^{O(d)} and 2^{O(d)} losses really do come from four write-up choices, and the four repairs really do cost only poly(d).
  - I checked each repair against the actual [pp] proofs, not against lv5unif's summaries.
  - I also checked that the computable S makes m = S²ρ_rρ_F polylogarithmic. This matters because γ is raised to the power 2k ≍ D^{2/5}.
- **The step from 1/2 to 2/5 is mostly bookkeeping.** This is the most useful structural observation (critic 2; I confirmed it).
  - Relative to the algebraic value size x = X_F = N^{θ}, with θ ≍ D^{−1/5}, the smoothness bound is still B = L_x[1/2, Θ(1)].
  - So every analytic input (Hildebrand–Tenenbaum, CEP, Tenenbaum's Ψ_q, Drappeau's Lemme 5, the large sieve) is used in
    the same regime, relative to x, as at fixed degree.
  - What is new is only explicit constants, the four repairs, Lemma CS, Lemma S″ and the factorial-size character primes.
- **So the whole claim is only as strong as [pp]'s unrefereed second-moment machinery.**
  - That means LD*, the multiplier w, L2_fam and divisor switching, together with Lemma 8.5′.
  - It also rests on Drappeau's Lemme 5, whose published proof refers to Harper's unpublished preprint 1208.5992, §3.3.
  - If that chain is right, rigorous L[2/5] follows. If it has a hole, Theorem U and Theorem A fall together.
  - This is a dramatic claim, and it still needs an outside expert referee of [pp]. All reviews in this project are AI reviews.
- **Weakest link, as named by the reviewers.**
  - Lemma 8.5′(R′) / Drappeau Lemme 5: codex and Opus critic 1.
  - Theorem D-U (γ ≥ (log N)^{−C_γ}): Opus critic 2.
  - My own answer combines them: the relation-supply estimate γ ≥ (log N)^{−C_γ}, i.e. LD* together with Lemma 8.5′. It is inherited from [pp], and Theorem U adds no new risk to it beyond the regime transfer of §5.
- **The asymptotics start astronomically late.**
  - Every N^{η₀} margin beats L[2/5]-scale quantities only by a factor D^{1/10}. For example, N^{η₀/4} ≥ B^{10} needs D = log N/log log N ≳ 6·10^{13} (A8).
  - This is not a proof gap. The algorithm never uses N₀, because it alternates SPLIT with trial division. But c is a purely asymptotic constant, and it is not effective.

### Reviewer verdicts (all blind to every prior verdict)

| Reviewer | Scope | Verdict | Weakest link |
|---|---|---|---|
| auditU (me) | whole theorem | **ACCEPT WITH FIXES** (conditional on the internal chain) | relation supply γ (LD* + Lemma 8.5′(R′)), inherited from [pp] |
| GPT-6.1-sol `codex/auditU_1` (xhigh) | whole theorem, fixed-d chain stipulated | **ACCEPT WITH FIXES** | Lemma 8.5′(R′): the small-conductor → smooth-modulus transfer |
| Opus critic 1 (analytic, uniformity) | relation-supply chain | **ACCEPT WITH FIXES** | Lemma 8.5′(R′) (Drappeau Lemme 5; D^{1/10}-thin level gap) |
| Opus critic 2 (algebraic, algorithmic, cost) | LV6, split, recursion, cost | **ACCEPT WITH FIXES** (in scope; "not an ACCEPT of the whole theorem") | Theorem D-U, i.e. the analytic chain |

### Issue table

| # | Sev. | Location (lv5unif_results.md) | Issue | Fix | Found by |
|---|---|---|---|---|---|
| 1 | MINOR | §13 (R3) | "exp((1 − α)θ log N·log log N/log B)": at θ ≍ D^{2α−1}, u ≍ D^α, so log u ∼ α log log N; the coefficient is α. The 2/5 threshold is unchanged, and σ = 2θ₀/(5a) in §10 already uses the right value. | Replace (1 − α) by α. | codex |
| 2 | MINOR | §4 Lemma FIB | The Möbius error 2^{ω(ℓ\|n\|)} is compared, as worded, with N^{η₀/2}. But 2^ω can be ≫ N^{η₀/2}. The correct comparison is with M/V = N^{μ−ν−o(1)}, since log 2^ω = O(μ log N/λ) = o((μ−ν)log N). | Reword. The conclusion holds. | Opus 1, codex |
| 3 | MINOR | Lemma P, (P2), (P10), (R′4), Lemma W-U | Each N^{η₀}-type margin beats Λ only by D^{1/10}, so N₀ is astronomical (A8). | State this, and that N₀ is not used by the algorithm. | Opus 1, me |
| 4 | MINOR | §1 "Common V" | "ℓ ≤ 2Y < V²/2" fails at S = 1. The inequality actually needed, ℓ < 2(1−1/d)V², holds for all S. S_ℓ ∈ [S(1−o(1)), √2S] (the stated range is weaker but still fine). | Reword. | Opus 1, Opus 2 |
| 5 | MINOR | §6 LD*-U Step 5 | U = 4M/V·ρ_r(log N)^C, not 3M/V (two signs). | Constant only. | Opus 1 |
| 6 | MINOR | §5 Lemma 8.5′(R′) | Given as a list of substitutions into [lem85], not a rewritten proof. Every substitution checks: I re-checked E7, SM(a), SM(b), Part 1, Part 2 and (i) against lem85 §§2–6, and so did Opus 1 and codex. | Write it out in full. | Opus 1 |
| 7 | MINOR | §8 Theorem K-U | "2k dyadic windows": k suffice. | Cosmetic. | Opus 2 |
| 8 | MINOR | §8 Prop select-U(ii) | The displayed (b/p)·2(…) chain is garbled. The argument (size-biased law, E[1/ς_{A*}] ≤ 1/b) is right. | Rewrite. | Opus 2 |
| 9 | MINOR | §10 cost table | 3.517, 4.047, 5.038, 6.498 are valid upper bounds but not grid optima. My grid gives 3.5095, 4.0395, 5.0350, 6.4722, and Opus 2's gives 3.5103, 4.028, 5.031, 6.493. | Write "≤". | Opus 2, me (A6) |
| 10 | MINOR | §9, [pp] Fact LO | The Lagarias–Odlyzko theorem number is unpinned (also noted in lv5unif §14). | Pin it. | Opus 2 |
| 11 | MINOR (docs, not lv5unif) | notes/nested_results.md Verdict | "nor any fixed exponent below 1/2" and "best α under ERH: 1/2" read literally as conflicting with Theorem U. They refer to nested's own nested-descent route. nested's Prop R2 itself puts designs without algebraic children at the one-level 2/5. | Owners should qualify the wording ("this route"). No mathematical conflict. | me (orchestrator request) |
| 12 | scope | §14, lem85 | The single external analytic input of substance is Drappeau Lemme 5, whose proof points to unpublished Harper §3.3. It is used only for conductors ≤ B^K, in Harper's range (log x)^K ≤ y. | Not repairable here. Inherited from Theorem A. | all |

No issue is FATAL or MAJOR. Items 1–10 are text or constant corrections inside lv5unif_results.md, whose owner should apply
them. I did not edit it.

---

## 1. The four τ = 0 repairs (task 1)

**(a) b-window ((1−1/d)V, V] with a shared V. True, loss poly(d).**
- **Effect on Kern.** (V²/bb′)^{d−2} ≤ (1−1/d)^{−2(d−2)} ≤ e^{8/3}, so the old 4^{d−2} disappears.
- **Effect on #𝒫.**
  - #𝒫 ≍ V²/(d log V). Primes in a window of relative width 1/d: PNT with de la Vallée Poussin's error, since log d ≪ √log V.
  - The upper bound #𝒫 ≤ 5V²/(d log V) is Brun–Titchmarsh.
  - The 1/d appears in FM as a loss and in LD* as a gain. The net is (2d+1)²/d ≤ 9d in γ's denominator.
- **Effect on FM.** The narrower window does not hurt FM-U.
  - The fibre-chord argument of [pp] Lemma `lem:fibres` needs only V′ ≥ V/2.
  - One-point smoothing needs H ≥ 3dV′.
  - The good-modulus layer argument loses V′^{d−1} ≥ V^{d−1}/e. Even the old 2^{−d} would only have shifted the argument of ϱ, which costs 2^{d(1−α)} = 1 + o(1).
- **Effect on Theorem K.** The window does not enter Theorem K-U at all. K-U depends on the field only through δ_A = |𝒟_A|/π_J.
- **Common V.**
  - Removes the (V_max/V₋)^{2d} ≤ 2^d of old Step 7.
  - Simplifies the bad-modulus ranges to one dyadic family covering [MV/32, 2MV].
  - Keeps Lemma `lem:reg` at ≤ 179V²/ℓ ≤ 358S² ≤ 400S² points and keeps the pair bound q ≤ 4S². I re-derived both; see issue 4.

**(b) Envelope X_F = (2d+1)MHV^d. True.**
- In [pp] Lemma `lem:valuesizes`, |F̂(P)| ≤ 2^{2d+5}MHV^dN^{−η₀/2} ≤ MHV^dN^{−η₀/4}, since 2d + 5 ≪ η₀ log N/4 = (log N log log N)^{1/2}/4.
- Also |G(P)z| ≤ 2dMHV^d for |z| ≤ dHV^{d−1}.
- The exponential base-expansion constant is absorbed before any smoothness density is taken. Only (2d+1)² reaches Ψ(X_F)².

**(c) Kern-U. True, with constant e uniformly in d.**
- **Proof check.**
  - Gauss's lemma gives the kernel basis e_m = Q·x^{d−3−m}y^m, with entries bb′, −(ab′+a′b), aa′.
  - So coordinate i of a fibre point depends on λ_i through bb′λ_i plus earlier λ's. That gives the triangular count (⌊2H/bb′⌋ + 1)^{d−2}.
  - (1 + 1/d)^{d−2} ≤ e when bb′ ≤ (2H+1)/d.
- **Brute force (A3).** 75 random cases (d = 3–5), enumerated over all of C_H independently of lv5unif's code: 0 violations; maximum normalised constant 1.21 < e.
- **The (3/2)^d family.** At H = b² the normalised constant grows (1.26, 1.60, 1.78 for d = 3, 4, 5), and those cases violate Kern-U's hypothesis, as they must.
- **Image index (A4, Smith normal form, 100 pairs, d = 3–7).** The index of Φ(Z^d) in Z² is exactly |Δ|.
  - So Kern-U's point bound exceeds the uniform-on-image density by about d²V²/|Δ| = d²S²/q.
  - The q-sum pays for this: there are ≤ 10 partners per q, and Σ_{q≤4S²} S²/q ≍ S² log S. This is a log factor, consistent with Step 7. Opus 2 made the same probe.
- **Step 7 power count (A2, sympy).** (#𝒫/|𝒢|)·kern·U·Ψ(X′₊)Ψ(X_F)²/(Y(MV)²) = (2d+1)²V⁴/(dY²) = O(d)·4S⁴.
  - The powers of M and H cancel exactly, and no 2^{O(d)} remains.
  - All three reviewers re-derived the same.

**(d) Lemma 8.5′ at Y = L[3/5] in regime (R′). True as applied.**
- I read lem85 §§2–6. Regime (R) is used only through the consequences lv5unif §5 lists:
  - (α1)/(α2);
  - E7 (ϱ ≥ y^{−1−K₀});
  - Lemme 5's hypotheses (log z)^{c_D} ≤ y and y^K ≤ z^{η_D};
  - the large-sieve level conditions Y ≥ y^K and Q ≤ √x·e^{−Γ};
  - the count of 𝒯_Q.
- Each holds in the design:
  - u ≍ D^{2/5} ≫ (log log N)⁴;
  - u log(4u log B)/log B → (4/5)θ_F₀/a² < K₀;
  - log Y ≍ D^{3/5}λ ≫ K log B;
  - Γ ≥ (η₀/2) log N − O(Λ), which exceeds O(log B) by a factor ≍ D^{1/10}.
- The decisive comparison: Drappeau's saving exp(−δ_D u/(2 log²u)) = exp(−Θ(D^{2/5}/λ²)) beats the transfer loss (log x)^{K(1+C_α)} = e^{O(λ)}, with K fixed.
- lem85's remark that "this works exactly at y = L[1/2]" refers to y = L_x[1/2, Θ(1)] relative to x. That is exactly the situation here (B = L_x[1/2, a/√θ_F₀]).

**Side effects of the repairs.**
- **H ≥ dV² does not move θ.**
  - h_* already contains the kernel floor 2ν (H ≥ V² at fixed d).
  - The extra factor d·N^{η₀/2} is absorbed by the existing slack η₀ = D^{−1/2} = o(θ).
  - In the binding regime κ < d₀²/4 (where the cost optimum lies) the kernel floor is the active one: h = 2ν + η₀.
  - Kern-U's hypothesis then needs N^{η₀} ≫ dS², i.e. (log N log log N)^{1/2} ≫ σΛ. True, with a D^{1/10} margin.
- **The level floor is automatic.** h ≥ μ − (d−2)ν + η₀ gives log(√x/Q) ≥ (η₀/2)log N − O(Λ) in both regimes.

## 2. Uniformity at growing (d, k) (task 2)

I went through every "N ≥ N₀", every C_d and every o(1) in [pp] §§2–9, 11 (via grep for O_d, C_d, L^{o(1)}, h₂, ω(N)) and
in lem85. Result:

| Input | Status at d ≍ D^{1/5}, k ≍ D^{2/5} |
|---|---|
| HT Thms 1, 2(i), 3 | uniform for x ≥ y ≥ 2; error 1/u + log y/y → 0. Shifts by B^C cost (u log u·…)^C = (log N)^{O(C)} with C fixed (C₀ ≤ 2σ/a + 3, C₁ = ⌈σ/a⌉ + 3); shifts by 2^{O(d)} cost 1 + o(1). |
| CEP / de Bruijn | u ≍ D^{2/5} ≪ log z/log log z ≍ D^{4/5}; used only for sizes (σ). de Bruijn's second-order term would be fatal if dropped (it is ≫ log log N); Lemma CS never uses it. |
| Lemma CS (computable S) | correct. R*/ϱ = α√(2πφ₂)(1+O(1/u)) ≤ log N (HT Thm 1); the bisection error is O(1) by convexity; existence is proved before the scan; m ∈ [(log N)^{−2C_R}/4, 64], since doubling S changes X_F by 2^d ≤ B. |
| Tenenbaum Ψ_q | needs ω(t) ≤ y^{c/log(1+u)} = e^{Θ(D^{2/5})}; ω(t) ≤ log N. Error O(λ log log N/D^{2/5}) → 0. Ψ_t/φ(t) vs Ψ/t: ∏_{p\|t}(1−p^{−α})/(1−1/p) = 1 + o(1). |
| Drappeau Lemme 5, Harper | used for conductors ≤ B^K only, z ≥ x^{1/2}; (log z)^{c_D} ≤ B. Range fine; constants absolute. |
| Large sieve (MV), Brun–Titchmarsh | exact inequalities. |
| Lemma W | C_W = 24/η₀ = O((log N)^{1/2}); needs η₀ log N ≫ log S (D^{1/10} margin). |
| L2_fam | the cushion B (not L[1/2]) keeps C₀ fixed; classes (B), (C) save N^{−h−(d+1)ν}B^{O(1)}d^{O(1)}, N^{−2ν}B^{O(1)}, N^{(5−d)ν−h}B^{O(1)}. |
| FM-U, LD*-U, D-U | losses d, d log V, log log N, η₀^{−1}, (log N)^{O(1)}: C_γ = 2C_FM + C_LD + 2C_R + 5 depends only on (a, κ, d₀) and absolute constants. |
| Theorem R | used only by Prop MF (not by Theorem U). Lemma `lem:reg` with common V gives ≤ 400S² partners. |
| Lemma S″ | verified line by line. Coprime same-class pairs ≥ \|S\|²/(2·2^ω) (Cauchy–Schwarz + shared-prime count, needs \|𝒟_A\| ≥ 4k³2^ω, (P6), margin λ = log log N); injective into Q_N; \|Q_N\| = φ(N)/2^ω cancels. No 2^ω loss and no C(2k, k). |
| Split probability | (δ_A/(4 log N))^{2k}/(32k²) re-derived exactly. With γ ≥ (log N)^{−C_γ}: log(1/p) = 2κ(1+C_γ)Λ(1+o(1)), i.e. L[2/5, s]. **Any** 2^{d} or (log N)^{d} loss in γ would be fatal here (k·d ≍ D^{3/5} ≫ Λ; A7), which is why τ = 0 is the whole question. |
| Rankin box radius | V = ⌊√(2Y)S⌋, S² ≤ e^{(σ+o(1))Λ}. S^d = N^{o(θ)} (d log S ≍ D^{3/5}λ ≪ θ log N ≍ D^{4/5}λ). |
| η₀ = D^{−1/2} | every N^{η₀} margin = L_N[1/2, Θ(1)] beats L[2/5]^{O(1)}, 2^{O(d)}, S^{O(1)} (D^{1/10}); η₀ = o(θ) (θD^{1/5} → θ₀; A1 shows θ·D^{1/5} → θ₀ with the actual integer d, k). |
| Lemma IR, rejections | the irreducible count r^{d−1}/(2d) holds for r ≥ 16, d ≥ 3; CRT distortion (1+1/d)^d ≤ e (Q_IR = H^{1/2+o(1)}); Pr ≤ exp(−Θ(D^{2/5})) (kernel regime) = o(γ). |

## 3. LV6 at growing degree (task 3)

All correct, and the cost is L^{o(1)} on top of L^{2a}.
- **Stark's threshold.**
  - Stark's Lemma 8 is applied with the general constant 4·n! = 4(2d)! for L = K(√h).
  - d is odd, so K has no quadratic subfield and L has at most one (4 ∤ 2d).
  - Kummer injectivity plus Landau (conductors ≤ Δ̄^{1/d}, ν₁ ≤ c_L/(2 log Δ̄)) then give ≤ 1 exceptional class for every field and every d.
- **Character primes of s = ⌈(2d)!(log N + log Δ̄ + d)⁴⌉ bits.**
  - They satisfy ν₁s ≫ log s and Lagarias–Odlyzko's s ≥ 20d(log Δ̄)², with log Δ̄ = O(d log N).
  - They also control the error c₁x·exp(−c₂√(s/2d)).
- **Cost.** poly(s) = exp(O(d log d)).
- **Characters.** t_χ = O(d² log N).
- **Square roots.** Lemma SQ costs 2^{d−1} sign vectors × Newton lifts of L^{a+o(1)}-bit numbers = e^{O(D^{1/5})}L^{a+o(1)}.
- **Linear algebra.** L^{2a+o(1)}.
- **The orchestrator's question (nested's correction).** nested shows that at α = 1/3 the unconditional Stark-threshold repair costs L[1/3, Θ(δ)]. This does **not** carry over to Theorem U.
  - Exactly (A9): (2d log d)/Λ_α ≍ d₀δ·D^{δ−α}. This is a positive constant iff δ = α (nested: d = δD^{1/3} at α = 1/3) and → 0 iff δ < α.
  - Theorem U has δ = 1/5 < α = 2/5, so the factor is D^{−1/5}: exp(O(D^{1/5}λ)) = L[2/5]^{o(1)}.
  - Concretely, 7·2d ln 2d/Λ = 0.19, 0.028, 2.8·10^{−4} at D = 10⁶, 10^{10}, 10^{20} (d₀ = 1).
  - There is no hidden δ-dependence: d₀ only multiplies a vanishing ratio.
  - (LV Remark 6.7 states the same L(1/3, 4δ/3) for their own degree; checked in LV's text.)
- **nested's Props SI and R2 are consistent with Theorem U.**
  - Theorem U is exactly a design without algebraic children: every algebraic value is B-smooth, and the special-q targets ℓ ∈ J are exogenous.
  - Its relation supply uses route (ii) of Prop SI (an average over the field for a field-independent target set J). Its split uses the counting-plus-uniform-group principle (Lemma S″), which Prop SI names as legitimate.
  - Prop R2 says such designs cannot go below 2/5 in the size model (its "window" is non-empty iff α < 2/5). Theorem U sits at the boundary 2/5, as §13 says.
  - nested's "α = 1/2 under ERH" concerns nested descent only (issue 11).
- **Fruitfulness ((T)).**
  - Every quantity before X depends on the squares R_i and on fresh coins.
  - The recursive FACTOR_U calls on |u₀| and |v₀| use fresh coins and canonical outputs (Lemma OS/T_rec).
  - So X/Y is uniform on μ_N, and the round succeeds with probability 1 − 2^{1−ω} ≥ 1/2, uniformly in d and k.

## 4. Recursion at Λ_{2/5} (task 4)

- **Lemma F_α.**
  - Λ_α(√(2x)) ≤ 2^{−α}Λ_α(x) + O((log x)^{α−1}λ^{1−α}).
  - The induction with η′ = η(1−2^{−α})/4 closes.
  - Finiteness of F holds by the trial-division alternation.
- **The shape test on the halves.**
  - |u₀|, |v₀| ≤ λ₁(Λ_R) ≤ (4/3)^{1/4}√N.
  - Its cost is F(√(2N)) = L_N[2/5, 2^{−2/5}c].
  - The adaptive sub-calls are handled by Wald, not by independence.
- **Recursion inequality (10.1).** F(x) ≤ 2log₂x[P(x) + D(x)F(√(2x))] with D = 2(2M_rows + 1)/p = L^{a+s}. All three reviewers re-derived it.

## 5. Cost formula and exponent (task 5): exact script results (`experiments/raw/auditU_constants.json`)

- **A1 (θ).**
  - Exhaustive breakpoint minimisation of θ(μ) at the actual integers d = odd ≥ d₀t, k = max(⌈κt²⌉, d+1) (t = D^{1/5} = 10, …, 10⁴) gives θ·t → θ₀ = 3/d₀ or 2/d₀ + d₀/(4κ).
  - For example, (d₀, κ) = (1, 1/8): 4.09, 4.020, 4.002, 4.0002 → 4.
- **σ = 2θ₀/(5a).** u_r + u_F = (θ₀/a)D^{2/5} and log u ∼ (2/5)λ.
- **A5 (α floor).**
  - The exponent LP is: minimise α subject to κ′ + τδ ≤ α (split), −δ ≤ 2α − 1 and δ − κ′ ≤ 2α − 1 (smoothness).
  - Its optimum is α = 2/5 at δ = 1/5, κ′ = 2/5 (τ = 0). The grid minima for τ = 1, 2 approach 3/7 and 4/9 (= (2+τ)/(5+2τ)).
- **A6 (the constant).**
  - 1/(1 − 2^{−2/5}) ∈ [329877/79877, 13195079/3195079], decided exactly via x⁵ vs 4.
  - c = max(2a, a + s + max(a, σ), (a + s)/(1 − 2^{−2/5})), s = 2κ(1 + C_γ). My grid optima (a, d₀, κ) are below; lv5unif's own points reproduce their values exactly.

    | C_γ | c ≤ | (a, d₀, κ) |
    |---|---|---|
    | 0 | 3.5095 | (23/40, 17/16, 35/256) |
    | 1 | 4.0395 | |
    | 5 | 5.0350 | |
    | 20 | 6.4722 | |
    | 100 | 8.855 | |

  - At the optimum the special-q term and the recursion term balance; linear algebra (2a) never binds.
- **Is c finite independently of the uncomputed C_γ?**
  - Yes: c is a finite piecewise-affine function of C_γ, and C_γ is a finite constant depending only on (a, κ, d₀) and on absolute analytic constants.
  - The algorithm never uses C_γ: S is found by scanning, and N₀ is avoided by the alternation.
  - c is **not effective**, because C_γ contains HT's and Drappeau's unspecified constants. The numbers above are illustrations.
- **A7 (margins, leading-term exact arithmetic).** All 18 of my comparisons hold. Two of them are negative controls that show why τ = 0 is necessary:
  - a 2^d loss in γ would cost k·d ≍ D^{3/5} ≫ Λ;
  - 2^{ω(N)} is not L^{o(1)}, which is why Lemma S″ is needed.
- **A8 (where the asymptotics start).**
  - N^{η₀/4} ≥ B^C needs D ≥ (4Ca)^{10}: 3.7·10⁸ for C = 3 and 6·10^{13} for C = 10.
  - The C_W margin needs D ≳ 10⁹.

## 6. Sanity and literature (task 6)

**Why has this not been done before?**
- Lee–Venkatesan (arXiv 1805.08873, JNT 2018) prove, unconditionally and at growing degree, that a randomised NFS finds congruences of squares in L[1/3]. Their argument is a first moment with L^{o(1)} precision.
- They do not get fruitfulness. Their Thm 2.3 covers only p ≡ q ≡ 3 (mod 4) and assumes their Conjecture 7.1. §7 says that for general n "it is unclear how (even notionally) one might show that the linear algebraic step may produce non-trivial congruences." (LV text checked.)
- The project's route buys fruitfulness from Dixon's hidden square roots. The price has two parts:
  - **Geometric.** The residues must be split into L[3/5]-size primes and descended in one level, and that geometry gives exactly 2/5 (§13; nested Prop R2).
  - **Analytic.** It needs a per-special-q descendability statement: γ = Pr_{ℓ,A}[ℓ descends] ≥ (log N)^{−C}, with polylog precision, because γ is raised to 2k.
- That needs a second moment over random fields: the multiplier w, Kern, L2_fam and the restricted BV over moduli ℓ·(smooth). LV have no second moment, and LV's L^{o(1)}-precision first moment would not survive the power 2k.
- So the novelty, and the risk, sits in [pp]'s second-moment chain. The growing-degree bookkeeping is genuine but routine once the four write-up artefacts are seen.
- Heuristically, a one-level Dixon descent at L[2/5] is the natural exponent for this family. I checked this myself:
  - split u = log√N/log Y, so log Y ≳ (log N)^{1−α};
  - norms ≈ exp(√(log N log Y)) at the optimal d;
  - smoothness at L[α] then needs 1 − 3α/2 ≤ α, i.e. α ≥ 2/5.

**Known obstruction?**
- None found. The folklore L[1/2] barrier of random-squares methods comes from asking √N-size residues to be L[α]-smooth. Descent removes that requirement, as the recursion audit (auditRA) noted.
- nested's Prop R2 is consistent with Theorem U: it bounds such designs below by 2/5, not 1/2.

**Gordon 1993.**
- I could not obtain the text. Codex read it (dmgordon.org) and reports that its individual-logarithm step is heuristic and stated at L[1/3]. A search found no source that attributes a one-stage L[2/5] to Gordon (lv5unif §14 also leaves this unverified).
- No Gordon-type descent supplies a rigorous relation-supply estimate or a factoring fruitfulness theorem.

**Literature search** (WebSearch, extended; the codex web search; lv5unif's literature critic).
- No rigorous L[α < 1/2] for factoring arbitrary N was found. Lenstra–Pomerance 1992, L[1/2, 1 + o(1)], is still described as the best proven bound.
- No published L[2/5] Dixon-plus-one-level-descent algorithm, heuristic or rigorous, was found. This is a search result, not a novelty judgement.

## 7. Bottom line

- Theorem U is true **if** [pp]'s fixed-degree chain, Lemma 8.5′ and Drappeau's Lemme 5 are true. Given them, I found no remaining gap at growing degree, and neither did three blind AI reviewers.
- The extension adds no analytically new input: relative to x, the regime is the fixed-degree one. So the dramatic claim (rigorous L[2/5], beating three decades of L[1/2]) is exactly as reliable as [pp]'s unrefereed second moment.
- The next step that would most change confidence is an expert human referee of [pp] §§5–6 (FM, LD*, L2_fam) and of Appendix A, not further AI review of the growing-degree layer.
