# thmA audit: is Theorem A true?

Agent `thmA`, phase 10, 2026-10-19.
- Proof roadmap: `notes/thmA_proof.md`.
- Pre-registration: `notes/thmA_prereg.md`, commit `9c67cc6`.
- Code: `src/thmA/`. Tests: `tests/test_thmA.py`.
- Data: `experiments/raw/thmA_e1.jsonl` and `_summary.json`.
- Codex reports: `codex/thmA_{1,2}.out`; full texts copied to `notes/thmA_codex{1,2}_report.md`.

## Verdict

**No FATAL or MAJOR error was found, by me or by any of the five blind referees (2 GPT-6.1-sol, 3 Opus). All five returned
ACCEPT WITH FIXES, and every fix is local and leaves the exponent unchanged.**

So, as far as this audit can tell, the internal chain does prove the following:
- for each fixed (d, k, μ, η₀) with θ < 1, a Las Vegas algorithm factors every N in expected time
  L_N[1/2, (1 + θ)/√2 + o(1)];
- at d = 9 this is below 1 from k = 13 on (c ≈ 0.9972 + η₀/√2);
- the constant tends to 0.9192 as k → ∞ and η₀ → 0, and to 1/√2 as d grows.

Four caveats apply:
1. **0.919 is an infimum over fixed parameters, not a single algorithm's "+ o(1)".** The claim as phrased in
   the brief ("L[1/2, (1 + 3/(d+1))/√2 + o(1)]") is an overstatement unless a growing-k diagonalisation is
   carried out, and nobody has done that.
2. **The weakest link is external-plus-transplant.** FM's restricted Bombieri–Vinogradov (Lemma 8.5′) imports
   the absolute-character majorant from the *proof* of Lee–Venkatesan's Lemma 8.5. It relies on Harper's
   Proposition 2 (arXiv 1208.5992, a preprint) and on LV's loosely written Page/Siegel absorption, used at
   u ≍ (log N/log log N)^{1/2}. Two analytic referees independently named this the weakest link. Both checked
   it against the sources and found it holds in this regime. It deserves a fully written proof before any
   public claim.
3. **The whole chain is internal and unrefereed** (lemmaM, hypD, pcest, lv6, ldbound, and this audit). Blind
   AI referees are not peer review.
4. **The literature search found no rigorous c < 1 result and no stated obstruction at c = 1** (§5). A result of
   this kind would be significant, which is a reason for more scrutiny, not less.

The toy run is consistent with every checkable claim. It is a correctness smoke test only.
- 48/48 complete factorisations.
- 0 row-identity and 0 congruence violations.
- 352/352 tested kernel vectors were squares.
- Theorem DD's fruitfulness rate matched exactly: 0.476 for ω = 2 (z = −0.82, n = 288) and 0.766 for ω = 3
  (z = +0.29, n = 64).

## 1. Issue table

Severity: FATAL (breaks the theorem), MAJOR (a gap needing a new argument), MINOR (local fix). Status:
FIXED-HERE means fixed in `notes/thmA_proof.md`; the source notes belong to other agents and were not edited.

| # | location | issue | severity | found by | status |
|---|---|---|---|---|---|
| 1 | statement (brief; ld §6 "Consequence") | "L[1/2, 0.919 + o(1)]" conflates an infimum over fixed k with o(1). Constants depend on k through (log N)^{O(k)}, C(2k, k) and N₀(k). At d = 9 the best fixed k gives 0.9972. | MINOR (wording); MAJOR if the o(1) form is claimed | me; thmA_1 §18; thmA_2 §18 | FIXED-HERE (thmA_proof §2 S1). Diagonalisation open |
| 2 | thmA_proof §6 U4, Step 5 | γ ≥ 1/(log N)^{O(1)} is false for the computable S. de Bruijn's second-order term makes m = S²ρ_rρ_F = L^{o(1)} but smaller than every inverse power of log N. Only γ ≥ L^{−o(1)} holds, and it suffices for fixed k. | MINOR (my error) | thmA_1 §15; Opus analytic F1 | FIXED-HERE |
| 3 | ld §5 hypothesis; pcest Thm 4.3 | LD* and PZ are stated for S² = (log N)^{C₃}/(ρ_rρ_F), which is not computable. The proofs hold for any S = L^{σ+o(1)}. PZ gives γ ≥ 1/(P³ + P/m) without E N ≥ 1. | MINOR | me (U5); thmA_1 §14; Opus analytic F2 | FIXED-HERE |
| 4 | pcest §2.1 (T), FM Step 2(d)/3 | Tenenbaum's hypothesis is misquoted (y^{c/log(1+u)}, not y^c/log(1+u)). It is also applied to a non-smooth modulus t = ℓ|n|; apply it to |n| instead. The range is ample. | MINOR | thmA_1 §2; Opus analytic F3 | FIXED-HERE (Step 3 note) |
| 5 | thmA_proof Step 3 | FM is proved averaged over ℓ, not for each ℓ. Only the average is consumed. | MINOR | thmA_1 §5; thmA_2 §8 | FIXED-HERE |
| 6 | ld §5 Step 0 | "=" should be "≤" (filters, F ≠ 0) | MINOR | Opus analytic F4 | noted |
| 7 | pcest §2.2 Part 1 ← LV Lemma 8.5 proof, Harper Prop. 2 | Imports the W-majorant from LV's *proof*, not its statement. Relies on Harper's preprint and LV's loosely written Page/Siegel step at u ≍ (log N/log log N)^{1/2}. Checked against the sources: consistent, ranges hold with fixed power margins, u/log²u ≫ log log x. | MINOR as found; **weakest link** | thmA_1 (weakest link); Opus analytic (weakest link, F5) | open: needs a self-contained write-up |
| 8 | thmA_proof Steps 6–7 (lv6 §4, ld §1) | Factor exits: E[g_A 1_{𝓡ᶜ}] ≥ γ/2 under the raw law does not survive conditioning on non-terminating fields. A search may also exit with a factor on an earlier candidate. Repaired by the absorbing-outcome bookkeeping, p = δ + b, cost C·C_trial/p. | MINOR (bookkeeping; thmA_2 calls it its weakest link) | thmA_2 §5 | FIXED-HERE |
| 9 | hypD K(ii)/(iii) | K(ii)'s display uses the fourth-power Theorem S, while the roadmap uses the S′ form. The shape count should use T_J = ∪_{j≤k}(j-products), not S_J. | MINOR | thmA_2 §4 | FIXED-HERE (Step 7 note) |
| 10 | lv6 (1.1) | The row identity omits sign(u₀v₀). Example: R = −1, split (1, −1), q = −1. | MINOR | me (U6); Opus algebraic F1; thmA_2 §9 | FIXED-HERE; the toy implements it, P2 holds |
| 11 | lv6 §2 hypothesis | "r_i uniform on (Z/N)^×" is false for accepted rows. What holds, and what DD needs, is uniformity on the fibre r_i·μ₂ given the certificates. | MINOR | thmA_2 §13 | FIXED-HERE |
| 12 | lemmaM §1 step 1 | Drawing r from [2, N−2] breaks exactness (N = 15, R = 1). The final algorithm uses [1, N−1]. | MINOR | thmA_2 §1 | correct in thmA_proof |
| 13 | thmA_proof Step 10 | A field-selection attempt costs T_shape + 2k·T_SQ, not T_shape + polylog. Same exponent. | MINOR | thmA_2 §17 | FIXED-HERE |
| 14 | lv6 §3.3 sampler | x must be the *least* admissible integer, so log x = O((log N)²). Stark's interval is narrower than LO's; the residual zero is handled by the same decay estimate. | MINOR | thmA_2 §12 | FIXED-HERE |
| 15 | thmA_proof §3 SQ | The Hensel prime must satisfy p ∤ D, p odd. | MINOR | Opus algebraic F2 | FIXED-HERE |
| 16 | lv6 §3.3 Stark constant | 4·n! checked only through restatements (LV Fact 6.22). Any O_d(1) constant works. | MINOR (citation) | Opus algebraic F3; thmA_2 §11 | noted |
| 17 | lv6 §1, §3.8 | Still written with M = N^{1/(d+1)}. The facts actually used are ‖f‖ ≪ MH < N, 0 < \|f_d\| < N and 2H+1 < N, all from μ + h < θ ≤ 1. | MINOR | Opus algebraic F4 | noted (ld §6.1 is correct) |
| 18 | lv6 §3.4 | Uniform kernel sampling is not needed. Two independent kernel vectors chosen independently of the roots suffice. "One fixed coordinate" means the coefficient of the basis element 1. | MINOR (presentation) | Opus algebraic F5; thmA_2 §13 | noted |
| 19 | constants | Re-derived exactly by me (`src/thmA/params.py`, independent of ldbound's module, agreeing on all (d, k) ≤ 200), by thmA_1 (d ≤ 51, k ≤ 500) and by thmA_2. k = 13 is minimal at d = 9: θ(12) = 23/54 > √2 − 1, θ(13) = 16/39, (1 + θ)² = 3025/1521 < 2; requires η₀ < 0.0039572. | — (verified) | me; thmA_1; thmA_2; complexity critic | — |
| 20 | Theorem R (pcest §3) | Not on the critical path: LD* re-derives the partner's rational weight in its Step 5 | — (verified) | me; thmA_1 §17 | — |
| 22 | thmA_proof §3 shape test | LP's output primes should be certified by AKS, so that the algorithm is Las Vegas | MINOR | complexity critic C1 | FIXED-HERE |
| 23 | thmA_proof §3 N₀ branch | N₀ must be computable. It is effective; alternatively dovetail SPLIT with trial division. "FACTOR(m) repeated e times" wastes polylog. | MINOR | complexity critic E2 | noted |
| 21 | toy design | The toy uses the relaxed split shape (B-smooth part × ≤ 2 primes in (B, 4B]) and w ≈ 2M. Theorem A's pure k-product shape is infeasibly rare at 30–50 bits. The analytic devices (Lemma W, the floors) are not tested. | disclosed limitation | me (prereg) | — |

Interface checks the brief requested, with results:
- **S/S′/K under F = Nw and floating M.** Valid verbatim. They treat G_A abstractly, provided G_A is the set
  recognised by the actual filtered search (thmA_2 §§3–4; complexity critic).
- **(T) with public w, g₀, g₁ and descent.** No exposed quantity depends on r beyond R = r². The roots first
  enter at X_S. Field retention biases the field law but not the fibre signs. Confirmed by me, thmA_2 §2, the
  Opus algebraic critic §11 and the complexity critic.
- **PZ / Theorem 4.3, filtered.** Correct when applied to the filtered count. The filter and prime-b are part of
  G_A's definition (thmA_1 §16; thmA_2 §16).
- **FM, R, Kern, L2_fam, Lemma W, LD* in (a‴_k; μ).** I re-derived the power counts (thmA_proof §§4–5).
  - FM gives (V²/ℓ)ρ_rρ_F.
  - LD* gives (H²V^{2d−4})^{−1}M^{−2}V²(M/V)(MV/Y)(MHV^d)²/(Y(MV)²) = V⁴/Y² ≍ S⁴, with exact cancellation of
    M and H.
  - Margins: μ − ν, h − 2ν, (ξ/2 − q₀) = (h − μ + (d−2)ν)/2, η₀/2, h + (d−5)ν, all ≥ fixed positive constants.
  - The Opus analytic critic checked the d = 9, k = 13 exponents: q₀ = 20/156 < ξ/2 = 25/156.
- **lv6 with content e and w.** e | Nw and e < N, so e ∈ {1, w} after factor exits. e = w rejects a field that
  has zero counted relations. w | f_d or w | disc f only enlarges S_D. Nothing needs (w, f_d·disc f·index) = 1.
- **Uniformity in N.**
  - Trial division to B gives 2^ω = L^{o(1)}.
  - Prime powers and perfect powers are removed first.
  - Non-squarefree N is fine (|μ₂| = 2^ω).
  - The cofactor uses its own parameters.
  - N₀ depends only on (d, k, μ, η₀).
  - No exceptional set of N appears. The effective Stark/Landau/LO constants make N₀ effective but huge.
- **Hidden L[1/2] loss.** None found by any reviewer.
  - Rows (d+1)π(B); per row (T_shape + T_SQ)L^{o(1)}.
  - Batch Õ(B + S²) per box with B ≤ S² for θ ≤ 1.
  - LV6, the Hensel sieve and characters are each ≤ L^{2a+o(1)}.
- **Expected time vs probability.** All loops are Las Vegas with verified outputs.
  - E[1/q_A | kept] ≤ 1/μ_sel.
  - Fresh-row rounds are geometric with mean ≤ 2.
  - Factor exits are handled by item 8.
  - LP is used through its expected-time bound for complete factorisation of every n (checked in LP 1992,
    thmA_2 §6).

## 2. Review verdicts (blind: no access to earlier verdicts or review sections)

| reviewer | scope | whole-theorem verdict | weakest link named |
|---|---|---|---|
| GPT-6.1-sol thmA_1 (xhigh) | analytic: FM, Lemma 8.5′, W, L2_fam, LD*, Kern, DS, PZ, constants | **ACCEPT WITH FIXES** (analytic half) | pcest §2.2 Part 1: the absolute-character majorant imported from LV Lemma 8.5's proof (induced and exceptional characters) |
| GPT-6.1-sol thmA_2 (xhigh) | algebraic/algorithmic: DD + (T), S/S′/K, LV6 transfer, IR, BS, cost, uniformity | **ACCEPT WITH FIXES** | Steps 6–7: raw-law supply → field selection with factor exits (repaired, item 8) |
| Opus critic, analytic number theory | FM, 8.5′, W, L2_fam, LD*, PZ; LV/Harper statements checked against the arXiv sources | **ACCEPT WITH FIXES** (analytic half only) | Lemma 8.5′ Part 1 rests on LV Lemma 8.5's Page/Siegel absorption at u ≍ (log N/log log N)^{1/2}; no published fallback |
| Opus critic, algebraic NT / NFS | LV6 (C, H, E, X, SQ, congruence), w/content, (T) | **ACCEPT WITH FIXES** (post-processing only) | Theorem LD* (novel second moment; any L^ε loss there breaks D_k). It did not audit LD* |
| Opus critic, complexity / probability / descent | DD, S/S′/K, field selection, LP shape test, cost, uniformity, constants | **ACCEPT WITH FIXES** | LV §6 adaptation (lv6 Lemmas H, X, E): a new unpublished fixed-degree specialisation that it did not itself re-derive. The Opus algebraic critic and thmA_2 did re-derive it and found only minor fixes. |
| Opus literature reviewer | prior art and barriers | — (search) | see §5 |

Two reviewers per half returned ACCEPT WITH FIXES with no MAJOR finding. The two analytic referees, working
independently, named the same weakest link (item 7).

## 3. Toy end-to-end run (E1/E1b; pre-registered `9c67cc6`; all predictions hold)

Pipeline at d = 3, μ = 1/3:
- random (g₀, g₁), random prime w ∈ (2M, 4M], so F(g₀, g₁) = Nw;
- uniform r, R = r², Gauss split, shape test, special-q descent of each large prime with batch smoothness per
  box;
- LV6 columns (sign, rational primes, one parity column, first-degree (r, s) with r ∤ f_d·disc f) and 40
  characters;
- GF(2) kernel, Hensel square root in Z[β], X² ≡ Y², gcd.

The deviations are listed in the prereg.

| prediction | result |
|---|---|
| P1 correctness (E1) | **48/48** complete factorisations (20 semiprimes, 12 three-prime, 4 unbalanced, 4 small × semi, 4 p²q, 4 powers/prime); E1b 12/12 |
| P2 row identity | **0** violations |
| P3 congruence X² ≡ Y² | **0** violations |
| P4 square fraction | **352/352** tested kernel vectors were squares in K |
| P5 DD, ω = 2 | **137/288 = 0.476**, predicted 1/2, z = −0.82 |
| P5 DD, ω = 3 | **49/64 = 0.766**, predicted 3/4, z = +0.29 |
| P6 accidental exits (E1) | 32 observed ≤ threshold 118 (null λ = 89; the null overcounts, since split halves are not uniform) |
| P7 preprocessing | 12/12 |

Step rates (descriptive, E1, per split):
- shape acceptance 0.11–0.26;
- descent success given shape 0.05–0.11;
- special-q success γ̂ = 0.06–0.13;
- 35–260 fields tried per kept field;
- rows 83–743, kernel dimension 37–216;
- 0.15–3.1 s per split.

Factor discoveries in E1:
- 24 STRUCTURAL_DD (final gcd);
- 32 ACCIDENTAL (mostly gcd(r, N) on N with 14–17-bit primes, which is why all twelve 3-prime first splits
  exited accidentally);
- 8 by trial division, 5 by the perfect-power step.

E1b (diagnostic continue mode) supplied the ω = 3 DD statistics.

The toy validates the algebra and the hidden-root symmetry end to end. It cannot test relation supply (FM/LD*),
which is where the theorem's content lies.

## 4. What I checked myself
- The exact constants (`src/thmA/params.py`; tests).
- Re-derivation of the FM and LD* power counts and of every margin.
- The Theorem S′ binomial comparison.
- Kern's index via the unimodular reduction.
- Lemma W and Corollary W.
- The L2_fam class bounds (A)/(B)/(C).
- The PZ inequality with L^{o(1)} first moment.
- The full row algebra including signs (implemented and verified, P2/P3).
- (T) step by step.
- The cost minimisation.

I found items 1, 3, 10 and 21 independently; reviewers found the rest. I could not independently verify the
interior of LV's §8 proof or Harper's Proposition 2. Item 7 therefore stays open.

## 5. Literature
Sources: my WebSearch, plus a dedicated Opus literature reviewer (89 tool calls) with primary texts where
accessible.

1. **No rigorous c < 1 found**, unconditional or GRH/ERH-conditional, 1992–2026.
   - Lenstra–Pomerance (JAMS 5, 1992; c = 1) is still called the best rigorous probabilistic bound by Harvey
     (arXiv 2010.05450, Math. Comp. 2021) and by Pomerance ("Smooth numbers and the quadratic sieve", MSRI 44,
     2008).
   - Lee–Venkatesan (JNT 187, 2018; arXiv 1805.08873) and van Leeuwen (arXiv 2007.02689) prove rigorous L[1/3]
     *congruences of squares*. Factoring is conditional on LV's Conjecture 7.1, and only for p ≡ q ≡ 3 (mod 4).
   - de Boer–Pellet-Mary–Wesolowski (arXiv 2512.01588) is ERH-rigorous for class groups and descent. They say
     the factoring application is "not investigated".
   - Kleinjung–Wesolowski (JAMS 2022) is a rigorous descent, but for discrete logarithms in fixed
     characteristic.
2. **No theorem that c = 1 is a barrier.** What exists:
   - LLMP: "some people suspected that L^{1+o(1)} would be the best we could ever achieve";
   - LP 1992: their own algorithm is exactly L^{1+o(1)} and cannot be improved by re-tuning;
   - the folklore model: combining random-like auxiliary integers of size ≥ N^{1/2} into a square costs
     L_N[1/2, 1]. This is sharp for random integers: Croot–Granville–Pemantle–Tetali, Annals 2012;
     Balister–Bollobás–Morris, Annals 2018.

   The model does not bind here. The √N-size halves u₀, v₀ are completely factored (by LP, at L[1/2, 1/√2]); they
   are not relations combined in the linear algebra. The relations combined are descent norms of total size
   N^{θ}, θ ≈ 0.3–0.41 < 1/2. This is how the design pays c = (1 + θ)/√2 rather than 1. This reading is ours and
   the reviewer's, not a published statement.
3. **Prior art for the devices; the combination was not found.**
   - Fixed-degree rigorous NFS at L[1/2]: not found. LV use growing degree, and they randomise f with
     f(m) = n exactly, with no multiplier.
   - Multipliers kN in NFS polynomial selection (Murphy, Coxon, Montgomery) are heuristic. LP's discriminant
     multipliers serve provability in a different way.
   - Dixon/Rabin hidden-root fruitfulness is classical. We found no source combining it with an NFS descent.
   - Rational-reconstruction split R ≡ u/v: standard heuristic start of NFS-DL descent.
   - LP as a shape-test subroutine: not found.

   We make no novelty claim.
4. **Not accessed:** Crandall–Pomerance ch. 6, Pomerance 1987 and ICM 1994, Seysen 1987, A. K. Lenstra 1988,
   Vallée 1991 (abstract only), most of Pomerance 1982 (scan without text).

## 6. What would make this a theorem one can state publicly
1. A self-contained proof of Lemma 8.5′ Part 1 (restricted BV with absolute-character majorant over moduli ℓ·s),
   not leaning on the interior of LV §8's proof.
2. Merge the fixes above into the source notes (owners: lv6, ldbound, hypD, pcest).
3. Optionally, a growing-k diagonalisation to state 0.919 + o(1) at d = 9.
4. External human refereeing of the whole chain.
