# Paper review round r17: adjudication

Agent `paperwriter2`, 2026-10-18.

The object of review is the restructured `paper/main.tex`. It leads with Theorem AV, followed
by D′ with the fourth-moment bound, Theorem R with DCH_pair, and T1, plus a placeholder on
unbalanced ratios.

There were two independent read-only reviews:
- the **Opus critic** (subagent `critic17`), verdict **ACCEPT WITH FIXES**;
- the **GPT-6.1-sol adversarial referee** (`codex/paper_r17_review.out`, xhigh), verdict
  **MAJOR REVISION**. It found no fatal error. It marked as CORRECT:
  - Theorems 1.1 and 1.2;
  - Proposition 4.13 and every lemma of §4;
  - Theorem 5.1, Lemmas 5.2 and 5.3, and Theorem 5.4;
  - Lemmas 6.1 and 6.2;
  - Corollary 1.4, as revised;
  - Theorem 8.2.

  Its gaps concern the specification of DCL-wc's cost, the scope of T1, the proofs in
  Appendix A, the hypotheses of the PYY remark, and prose.

Both reviewers re-derived every exponent: the E_c exponents −1/12 + 5ε/3 ≤ −ε − 1/60 at
ε ≤ 1/40; the AV exponents −1/18, −1/9, 13/90 and 7/36; and E(j, c). Both verified the Carella
defect against arXiv v3 and the Ward–Maestro quotations against the report.

The source was edited while codex reviewed. It froze a snapshot and credited the fixes already
present.

## Ward–Maestro 2026 (pre-emption check, done before the reviews)

Ward–Maestro, "Floors for deterministic integer factoring…", is an iGent AI report dated
2 September 2026, 110 pp. A subagent read it in full, and I re-read the key statements myself
from the extracted text.

**Pre-empts** our full-candidate N^{1/6} covering floor (multimod Theorem A; paper Appendix A):
- Def. 2.1 requires oblivious coverage of every real p.
- Lemma 3.2 is the covering-sum (curvature) bound.
- Thm 3.6 gives T ≥ 2(c_C/2)^{1/3}N^{1/6} with free giants.
- Thm 3.4 gives N^{1/5} with materialised giants.
- Thm 9.4 treats bounded-degree forms.
- Thm 2.8 and Lemma 2.3 treat thinning by a common modulus.

Our version adds per-target lock moduli. It is not new in exponent or method. The paper now says
so in §1.5, §2 and Appendix A, and Appendix A was reduced to statements with WM credited.

**Does not pre-empt** D′, D″, the fourth-moment bound, AV, R or T1/T2:
- WM has no almost-all or average-case algorithm below N^{1/5}.
- It has no many-prime selective covering and no shared list across moduli.
- It has no true-factor counting under residue-local keys.
- Its Remark 3.5, on randomised coverage of a δ-measure, is related to T1 in spirit only.

**Prop. 35.3** (p − 1 with base 2) is unconditional, using Friedlander and Fouvry. It succeeds
on a positive proportion of the nearly balanced family x^{1−ε} < p, q ≤ x at N^{1/(4√e)+}. It is
explicitly not an average bound for a complete algorithm, and it fails on a positive proportion
for smoothness bounds up to N^{0.33}. The paper cites it as lying outside T1's model.

## Opus critic: all 15 findings accepted

| # | Finding | Action |
|---|---|---|
| 1 | ERROR: "fourth moment improves the average polylog to n^{c+o(1)}" is unprovable | Replaced by "changes only polylog factors; a factor n^{c+5/6} remains within these bounds" |
| 2 | OVERCLAIM: abstract "a run fails on n/H², n⁴/H⁴" (only E_b obeys it) | Now reads "its coverage step misses the true factor on…" |
| 3 | OVERCLAIM: "average-case theorems avoid order questions altogether" | Now says they never need an element of large order; ord_q(2) enters only in counting E_a, E_c, C(T) and in the detector |
| 4 | GAP: T1 budgets are per input or uniform n-bit, not S_n-averaged | Caveat added to the abstract, Thm 1.5 and Remark 1.6(d) |
| 5 | GAP: Corollary 1.4(b) quantifiers | Now fixes c ∈ (3/8, 1/2] and the scale Y = N^{1/(7−2c)+o(1)}; the conclusion carries +o(1); the abstract says "≤ Y^{1/2}J^{−c}, c > 3/8"; the proof is rewritten |
| 6 | WM model has a common modulus (Lemma 2.3, Thm 2.8), not "no locks" | §1.5, §2 and Appendix A reworded |
| 7 | Prop 35.3 is about the nearly balanced family, not S_n; Prop 14.2 needs M ≤ N^{1/7} | Both fixed |
| 8 | AV novelty should be called routine | Novelty item (iv) now says "routine… standard restarts and BSGS order search; the substantive input is Prop 4.13" |
| 9 | "Why this shape" cites the wrong lemma | Now cites multimod Lemma K′ |
| 10 | Prop 8.3's lower-bound sentence follows from Thm 8.2 | Reworded |
| 11 | Lemma 7.3 wording; trial division of o; w_ℓ^{−1} for A < 0 | All added |
| 12 | §11: x ≍ Y³ should be Y^{25/8} | Fixed |
| 13 | Poor explicit constants in Thm 5.1 | Sentence added before Thm 5.4 |
| 14 | Appendix A labelled "proved" while only sketched | Appendix rewritten as statements only (see codex #45) |
| 15 | Overfull hbox | Fixed |

## GPT-6.1-sol referee: adjudication

**Accepted and applied.**
- **#36.** The averaged bound is written as O(2^{n/6}n^{c+1}log n); shell fractions use
  2^{−εn}; Theorem 1.5 uses 2^{−n/4}; a note says that E_b depends on T; generator families
  are fixed per n.
- **#37.** n₀ is stated explicitly as an eventual guard (n in the thousands), with the
  experiments far below it.
- **#41.** The whole weighted error sum of Theorem 5.4 is displayed before normalisation, with
  the source rows named and the λ < 1 case made trivial. The n² integer-to-prime loss is shown
  separately.
- **#40.** The μ²/φ identity and the large-sieve hypotheses are added.
- **#15.** Lemma 4.11 now states q > 16Y and ℓ ∤ qB.
- **#42 (major).** DCL-wc gets per-block parameters ℓ_k^−, m̄_k = ⌊3HP/(ℓ_k^−)²⌋ + 2, the
  endpoint-competitive m_k and empty-block skipping. The cost uses Σ log ℓ in place of |P| log Y.
  The proof of Theorem 7.1(iii) now itemises:
  - base selection;
  - the exponent trees;
  - the γ pre-step;
  - the equal-order split (an all-but-one tree over the prime powers of E; valuations give D;
    each enumeration runs at most once);
  - depth-by-depth batched exploration (total degree ≤ m and ≤ 2T′ points per depth);
  - rejected tests.
- **#3.** Theorem 1.3's cost now names the list cost, including target construction and rejected
  tests.
- **#43.** Proposition 7.7's proof is now complete:
  - the endpoint contribution is O(λP²/n⁴ + M₀P/n²);
  - the de la Vallée Poussin PNT holds uniformly for y ≥ x/(log x)³;
  - the constants 0.2√n and 0.14√n are derived.
- **#44 (major).** T1 now states that it bounds the prescribed coverage event, not accidental
  collisions. That is said in Theorem 1.5 and Proposition 8.3. Lemma 8.1 has gcd(x, K) = 1. In
  the A = 0 swap, R′ keeps only unit classes. Canonical normalisation is qualified: an offset may
  move into the window only if it is key-determined. Proposition 8.3 is restricted to the rule
  m = m̄_k, which suffices for Theorem 1.2, and states its cost comparison.
- **#45 (major).** Appendix A is now "statements only, proved in the project notes, not
  reproduced". It has a precise model, the interval [R/√θ, R] instead of the undefined J, and
  credits Ward–Maestro. This is the referee's own suggested alternative to completing the proofs;
  nothing in the main text depends on the appendix. The pointer to Lemma K was removed from §8.
- **#46.**
  - The PYY remark restores the hypotheses: π_ℓ ∈ (0, 1/2], max π ≤ 2 min π, λ ≥ λ_{j₀}, a common
    marginal, independent thinning, and odd j₀ via their (1.1).
  - It is labelled a limited-information barrier, not an arithmetic impossibility.
  - The dilate statement now requires gcd(t, ∏ℓ) = 1.
  - The order parameter was renamed j₀, to avoid a clash with J.
- **#47.** Prose fixes:
  - Markov tails from Theorem 1.1;
  - "average case is easy" qualified;
  - "hardest natural population" changed to "balanced population";
  - the abstract's practical comparison is now "our cost proxy … favours HH";
  - the novelty item on extraction is split into its average-case and worst-case parts.
- **#49.** Hittmeir is credited (Harvey's Alg. 4.1 is almost identical to Hittmeir's Alg. 3.2).
- **#52.** The p − 1 deduction is attributed to WM Prop. 35.3, with Friedlander and Fouvry as
  its ingredients, and the population is called "nearly balanced" everywhere.
- **#53, #55 (bibliography).**
  - Provenance prose was moved out of the printed `note` fields into a comment block at the top
    of refs.bib; the "unverified by us" note on Carella is gone.
  - DOIs were added for PYY (10.1002/rsa.20329) and AKS (10.4007/annals.2004.160.781).
  - Strassen's year is now 1976/77.
  - Brown 2010 was added ("Lower Bounds for Straight Line Factoring", ePrint 2010/047, verified
    on the ePrint page), with a short paragraph "Lower bounds in restricted models".
- **#22.** The detector reports FAIL when there is no proper gcd.
- **#58.** "Confirmatory experiments were pre-registered; exploratory analyses are labelled." The
  repository URL is a placeholder ("withheld in this draft"), because the author is blank.

**Partially applied or declined, with reasons.**
- **#56 (structure).** The review history in the acknowledgements is already short. A description
  of the AI process is required for honest disclosure, so it stays. The abstract was corrected but not shortened substantially (left for the submission version).
  The superseded D′ optimisations stay as a single remark
  (the "Theorem D′" remark of §4), since they explain the numbers 1/52 and 2/11 that the notes cite.
- **#57.** The clashes X = 2^n (removed) and J (sieve denominator renamed 𝒥; the detector uses J′;
  the PYY order j₀) are fixed. The use of M(·) in the cost table is unchanged, as standard
  notation.
- **#31.** This was resolved by restricting Proposition 8.3 to m = m̄_k, rather than proving a
  window-width comparison for the endpoint rule.

## Result

The paper compiles with tectonic: 28 pp., no undefined references or citations, no overfull
boxes.

We did not run a second review round. The fixes are local, and none changes a theorem statement
except as recorded above: Corollary 1.4(b), the scope of T1, and Appendix A.

**Status:**
- Theorems 1.1–1.3 and 1.5 are proved, with AI review only.
- Appendix A is credited context.
- §9 is a placeholder pending `notes/unbalanced_results.md`.
